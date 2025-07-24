import os
import pandas as pd
import numpy as np
import re
from multiprocessing import Pool
from pysam import VariantFile
from Bio.SeqUtils import seq1
from datetime import datetime
from argparse import ArgumentParser


def parse_arguments():
    parser = ArgumentParser(description="Process VCF files and export genetic mutation data to CSV format.")
    parser.add_argument("--input_dir", required=True, help="Directory containing VCF files.")
    parser.add_argument("--output_dir", required=True, help="Directory where processed CSV files will be saved.")
    parser.add_argument("--metadata_path", required=True, help="Path to the metadata file with sample information.")
    parser.add_argument("--start_pos", type=int, required=True, help="Start position for filtering genetic mutations.")
    parser.add_argument("--end_pos", type=int, required=True, help="End position for filtering genetic mutations.")
    parser.add_argument("--region", required=True, help="Region code for the data processing.")
    return parser.parse_args()


def date2month(date_str: str) -> str:
    date_pattern = re.compile(r"\d{4}[-/]\d{1,2}([-/]\d{1,2})?")
    match = date_pattern.match(date_str)
    if match:
        try:
            return datetime.strptime(date_str, '%Y-%m-%d' if '-' in date_str else '%Y/%m/%d').strftime('%Y-%m')
        except ValueError:
            return np.nan
    return np.nan


def map_aa(aa: str) -> str:
    match = re.match(r'([A-Za-z]{3})(\d+)([A-Za-z]{3})', aa)
    if match:
        return f"{seq1(match.group(1)).upper()}{match.group(2)}{seq1(match.group(3)).upper()}"
    return ""


def get_processed_srr_list(output_dir: str) -> list:
    return [file.replace(".csv", "").split("_")[0]
            for file in os.listdir(output_dir)
            if file.endswith('.csv') and not file.startswith('.')]


def process_vcf_file(args):
    input_dir, output_dir, start_pos, end_pos, filename, metadata, region = args
    filepath = os.path.join(input_dir, filename)
    srr = filename.split('.')[0]
    vcf = VariantFile(filepath)

    result = []
    for rec in vcf.fetch():
        if start_pos <= rec.pos <= end_pos:
            for ann in rec.info["ANN"]:
                fields = ann.split("|")
                if fields[1] == "missense_variant":
                    aa_change = map_aa(fields[10].replace('p.', ''))
                    result.append([
                        rec.pos, rec.ref, rec.alts[0], fields[3], rec.info.get("DP"), rec.info.get("AF"),
                        aa_change, fields[1], fields[9].replace('c.', '')
                    ])

    if result:
        df = pd.DataFrame(result,
                          columns=["Nucleo_Position", "Ref", "Alt", "Protein", "Depth", "AF", "AA", "Variant_Type",
                                   "NT_Change"])
        df['Date'] = date2month(metadata.loc[srr]['collection_date'])
        df['Country'] = metadata.loc[srr]['country']

        output_path = os.path.join(output_dir, f"{srr}_{region}_{df['Date'].iloc[0]}.csv")
        df.to_csv(output_path, index=False)


def main():
    args = parse_arguments()

    metadata = pd.read_csv(args.metadata_path, sep='\t').drop_duplicates('run_accession').set_index('run_accession')
    processed_srrs = get_processed_srr_list(args.output_dir)

    to_process = [
        (args.input_dir, args.output_dir, args.start_pos, args.end_pos, filename, metadata, args.region)
        for filename in os.listdir(args.input_dir)
        if filename.endswith('.vcf') and filename.split('.')[0] not in processed_srrs
    ]

    with Pool(4) as pool:
        pool.map(process_vcf_file, to_process)


if __name__ == "__main__":
    main()

    ## Command line
    ## python Preprocessing.py --input_dir ../vcf --output_dir ../csv --metadata_path ../Data/results_read_run_tsv.tsv --start_pos 21563 --end_pos 25384 --region IND

