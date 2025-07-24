import argparse
import pandas as pd
import numpy as np
from datetime import datetime
import os
import logging
from tqdm.contrib.concurrent import process_map
import re
import random
from mpmath import mp
from functools import partial

def parse_arguments():
    parser = argparse.ArgumentParser(description="Calculate p-values and odds ratios for mutation frequencies.")
    parser.add_argument("--rho", type=int, default=100, help="Sequencing depth threshold (default: 100).")
    parser.add_argument("--theta", type=float, default=0.02, help="Allele frequency threshold (default: 0.02).")
    parser.add_argument("--ns", type=int, default=200, help="Sample size threshold (default: 200).")
    parser.add_argument("--prevalence_threshold", type=float, default=0.01, help="Mutation prevalence threshold (default: 0.01).")
    parser.add_argument("--regions", nargs="+", default=["SA", "CA", "IND"], help="Regions to process (default: ['SA', 'CA', 'IND']).")
    parser.add_argument("--input_dir", type=str, required=True, help="Root input directory (e.g., '/Volumes/public/Data').")
    parser.add_argument("--output_dir", type=str, required=True, help="Output directory path.")
    parser.add_argument("--start_month", type=str, default="2020-03", help="Start month (default: '2020-03').")
    parser.add_argument("--end_month", type=str, default="2024-05", help="End month (default: '2024-05').")
    parser.add_argument("--precision", type=int, default=100, help="Precision for p-value calculation (default: 100).")
    parser.add_argument("--max_workers", type=int, default=4, help="Maximum number of workers for parallel processing (default: 4).")
    return vars(parser.parse_args())

def generate_month_range(start_month, end_month):
    start = datetime.strptime(start_month, "%Y-%m")
    end = datetime.strptime(end_month, "%Y-%m")
    months = []
    while start <= end:
        months.append(start.strftime("%Y-%m"))
        if start.month == 12:
            start = start.replace(year=start.year + 1, month=1)
        else:
            start = start.replace(month=start.month + 1)
    return months

def random_sample_srr_files(input_dir, month, sample_threshold, region, srr_log_file):
    """
    Randomly sample SRR files for a specific month and region.
    """
    random.seed(42)

    csv_files = [file for file in os.listdir(input_dir)
                 if file.endswith('.csv') and file.startswith("SRR") and month in file]

    if len(csv_files) <= sample_threshold and len(csv_files) > 0:
        selected_csv_files = csv_files
    elif len(csv_files) > sample_threshold:
        logging.info(f"Performing random sampling for {month}, region: {region}, sample threshold: {sample_threshold}")
        selected_csv_files = random.sample(csv_files, sample_threshold)
        logging.info(f"Number of files sampled for {month} in {region}: {len(selected_csv_files)}")
    else:
        logging.warning(f"No files found for {month}, region: {region} or sample size exceeds population")
        selected_csv_files = []

    return selected_csv_files


def get_p_zs(contig_list, precision=100):
    mp.dps = precision
    within_mutated, within_nonmutated, between_mutated, between_nonmutated = [x + 0.5 for x in contig_list]
    odds_ratio = (within_mutated * between_nonmutated) / (between_mutated * within_nonmutated)
    se_ln_or = np.sqrt(1 / within_mutated + 1 / within_nonmutated + 1 / between_mutated + 1 / between_nonmutated)
    z_ln_or = np.log(odds_ratio) / se_ln_or
    p_value_ln_or = 2 * (1 - mp.ncdf(abs(z_ln_or)))
    within_total = within_mutated + within_nonmutated
    between_total = between_mutated + between_nonmutated
    RR = (within_mutated / within_total) / (between_mutated / between_total)
    return float(p_value_ln_or), odds_ratio, RR

def split_aa_column(df, column_name):
    pattern = re.compile(r'([A-Z])(\d+)([A-Z])')
    first_part = []
    middle_part = []
    last_part = []
    for aa in df[column_name]:
        match = pattern.match(aa)
        if match:
            first_part.append(match.group(1))
            middle_part.append(match.group(2))
            last_part.append(match.group(3))
    df['Ref'] = first_part
    df['Site'] = middle_part
    df['Alt'] = last_part
    return df

def get_contig_table(merge_data, total, args, region):
    depth_threshold = args['rho']
    af_threshold = args['theta']
    filter_df = merge_data.query(f"Depth >= {depth_threshold} & AF >= {af_threshold}")
    if filter_df.empty:
        return None
    original_mut_inhost = filter_df.groupby('AA').apply(lambda x: (x['AF'] * x['Depth']).sum()).round().astype(int)
    original_total_inhost = filter_df.groupby('AA')['Depth'].sum()
    original_wildtype_inhost = original_total_inhost - original_mut_inhost

    mut_between = filter_df.groupby('AA')['SRR'].nunique()
    total_between = total
    wildtype_between = total_between - mut_between

    norm_mut_inhost = original_mut_inhost * total_between / original_total_inhost
    norm_wildtype_inhost = original_wildtype_inhost * total_between / original_total_inhost

    results = pd.DataFrame({
        'region': region,
        'AA': original_mut_inhost.index,
        'mut_inhost': norm_mut_inhost.values,
        'wildtype_inhost': norm_wildtype_inhost.values,
        'mut_between': mut_between.values,
        'wildtype_between': wildtype_between.values
    })
    results.set_index('AA', inplace=True)

    norm_metrics = results.apply(lambda row: get_p_zs([row['mut_inhost'],
                                                       row['wildtype_inhost'],
                                                       row['mut_between'],
                                                       row['wildtype_between']]), axis=1)

    results = results.query(f"mut_between >= {args['prevalence_threshold'] * args['ns']}")

    if not results.empty:
        results.loc[:, 'p_value'] = norm_metrics.apply(lambda x: x[0])
        results.loc[:, 'OR'] = norm_metrics.apply(lambda x: x[1])
    else:
        return None
    return results

def process_month(month, args, region, srr_log_file):
    input_dir = os.path.join(args['input_dir'], region, f"{region}_csv")
    sample_threshold = args['ns']
    prevalence_threshold = args['prevalence_threshold']
    logging.info(f"Processing data for the month: {month}")
    selected_csv_files = random_sample_srr_files(input_dir, month, sample_threshold, region, srr_log_file)
    results_contig_table = pd.DataFrame()
    if selected_csv_files:
        try:
            merge_data = pd.concat([pd.read_csv(os.path.join(input_dir, file)).
                                    assign(SRR=file.split('_')[0]) for file in selected_csv_files], ignore_index=True)
            total = len(selected_csv_files)
            mut_contig_table = get_contig_table(merge_data, total, args, region)
            if mut_contig_table is not None:
                mut_contig_table.insert(0, 'date', month)
                mut_contig_table = mut_contig_table[
                    mut_contig_table['mut_between'] / sample_threshold >= prevalence_threshold
                ]
                mut_contig_table.reset_index(inplace=True)
                mut_contig_table = split_aa_column(mut_contig_table, 'AA')
                mut_contig_table = mut_contig_table[mut_contig_table['Ref'] != mut_contig_table['Alt']]
                results_contig_table = results_contig_table.append(mut_contig_table)
                return results_contig_table, selected_csv_files
        except pd.errors.EmptyDataError:
            logging.error(f"Empty or invalid file format for files in: {input_dir}")
    else:
        logging.warning(f"No data available for {month} in {region}")
    return pd.DataFrame(), []

if __name__ == "__main__":
    args = parse_arguments()

    os.makedirs(args['output_dir'], exist_ok=True)

    srr_log_file = os.path.join(args['output_dir'], "srr_sampling_log.csv")
    months = generate_month_range(args['start_month'], args['end_month'])
    all_regions_results = []
    all_srr_log_entries = []

    for region in args['regions']:
        process_month_with_args = partial(process_month, args=args, region=region, srr_log_file=srr_log_file)
        monthly_results = process_map(process_month_with_args, months, max_workers=args['max_workers'])
        for month, result in zip(months, monthly_results):
            srr_files = result[1] if isinstance(result, tuple) and len(result) > 1 else []
            for srr_file in srr_files:
                all_srr_log_entries.append({'region': region, 'month': month, 'SRR_file': srr_file.split('_')[0]})
        region_results = pd.concat([result[0] if isinstance(result, tuple) else result for result in monthly_results if (isinstance(result, tuple) and not result[0].empty) or (not isinstance(result, tuple) and not result.empty)], ignore_index=False)
        if not region_results.empty:
            all_regions_results.append(region_results)

    # Combine all regions' results into a single DataFrame
    if all_regions_results:
        final_results = pd.concat(all_regions_results, ignore_index=False)
        final_results.reset_index(drop=True, inplace=True)
        output_filename = os.path.join(args['output_dir'], "hypothesis_results.csv")
        final_results.to_csv(output_filename)
        logging.info(f"Final combined results saved to {output_filename}")

    srr_log_df = pd.DataFrame(all_srr_log_entries)
    srr_log_df.to_csv(srr_log_file, index=False)


    ### command line
    # python Hypothesis.py --input_dir /Volumes/public/Data --output_dir ../output