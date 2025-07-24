import pandas as pd
import argparse
import os


def standardized_datetime(df):
    date_column = 'Date' if 'Date' in df.columns else 'date'
    if date_column in df.columns:
        sample_date = df[date_column].dropna().iloc[0]

        if '/' in sample_date:
            df[date_column] = pd.to_datetime(df[date_column], format='%Y/%m/%d', errors='coerce')
        elif '-' in sample_date and len(sample_date) == 7:
            df[date_column] = pd.to_datetime(df[date_column].apply(lambda x: x + '-01'), format='%Y-%m-%d',
                                             errors='coerce')
        else:
            df[date_column] = pd.to_datetime(df[date_column], errors='coerce')

    return df


def calculate_month_diff(row):
    year_diff = row['Date'].year - row['date'].year
    month_diff = row['Date'].month - row['date'].month
    return year_diff * 12 + month_diff


def filter_data_by_region(mutations_data, conditions):
    filtered_pos = mutations_data[
        (mutations_data['p_value'] < conditions['p_value']) &
        (mutations_data['OR'] > conditions['OR'])
    ]

    filtered_neg = mutations_data[
        ~mutations_data['AA'].isin(filtered_pos['AA'])
    ]

    return filtered_pos, filtered_neg


def evaluate_by_answer_set(answer_set, pos_data, neg_data):
    if not pos_data.empty:
        pos_merged = pd.merge(pos_data, answer_set, left_on='AA', right_on='Mutation',
                              how='left', suffixes=('', '_answer'))
        pos_merged['LeadTime'] = pos_merged.apply(calculate_month_diff, axis=1)
    else:
        pos_merged = pd.DataFrame()

    if not pos_data.empty:
        tp = pos_merged[pos_merged['AA'].isin(answer_set['Mutation'])].drop_duplicates(subset='AA')
        fp = pos_merged[~pos_merged['AA'].isin(answer_set['Mutation'])].drop_duplicates(subset='AA')
        fn = answer_set[~answer_set['Mutation'].isin(tp['AA'])].drop_duplicates(subset='Mutation').dropna()
    else:
        tp = pd.DataFrame()
        fp = pd.DataFrame()
        fn = answer_set

    if not neg_data.empty:
        tn = neg_data[~neg_data['AA'].isin(answer_set['Mutation'])].drop_duplicates(subset='AA')
    else:
        tn = pd.DataFrame()

    return tp, fp, tn, fn


def calculate_performance_metrics(tp, fp, tn, fn):
    median_lead_time = tp[tp['LeadTime'] > 0]['LeadTime'].median() if 'LeadTime' in tp.columns and not tp.empty else None
    tpr = len(tp) / (len(tp) + len(fn)) if len(tp) + len(fn) > 0 else 0
    fpr = len(fp) / (len(fp) + len(tn)) if len(fp) + len(tn) > 0 else 0

    tpr = float(f"{tpr:.3g}")
    fpr = float(f"{fpr:.3g}")
    median_lead_time = float(f"{median_lead_time:.1g}")

    results = {
        "TP": len(tp),
        "FP": len(fp),
        "TN": len(tn),
        "FN": len(fn),
        "Sensitivity (tpr)": tpr,
        "1 - Specificity (fpr)": fpr,
        "Median Lead Time": median_lead_time
    }
    return results


def site_specific_filter(df, region_range):
    site_specific_df = df[df['Site'].apply(lambda site: site in region_range)] if not df.empty else pd.DataFrame()
    return site_specific_df


def main(data_file, answer_file, output_dir, conditions):
    regions = {
        "RBD": range(319, 542),
        "Spike": range(1, 1274)
    }

    data = pd.read_csv(data_file)
    answer_set = pd.read_csv(answer_file)

    data = standardized_datetime(data)
    answer_set = standardized_datetime(answer_set)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    all_results = []
    all_pos_data = []

    for region_name, region_range in regions.items():
        print(f"\nProcessing region: {region_name}")

        ## filter data by region
        data_filtered = site_specific_filter(data, region_range)
        answer_set_filtered = site_specific_filter(answer_set, region_range)

        if data_filtered.empty or answer_set_filtered.empty:
            print(f"No data for region {region_name}. Skipping.")
            continue

        pos_data, neg_data = filter_data_by_region(data_filtered, conditions)

        ## ouput postivie results
        if not pos_data.empty:
            pos_data['date'] = pd.to_datetime(pos_data['date'], errors='coerce')

            pos_data = pos_data.sort_values(by=['AA', 'date', 'p_value'], ascending=[True, True, True])
            pos_data = pos_data.groupby('AA', as_index=False).first()

            all_pos_data.append(pos_data)

        else:
            continue

        ## calculate performance metrics
        tp, fp, tn, fn = evaluate_by_answer_set(answer_set_filtered, pos_data, neg_data)
        metrics = calculate_performance_metrics(tp, fp, tn, fn)
        metrics.update({"region_name": region_name})

        region_results_df = pd.DataFrame([metrics])
        all_results.append(region_results_df)



        ## output true positive results
        if not tp.empty:
            tp_file_path = os.path.join(output_dir, f"tp_results_{region_name}.csv")
            tp.to_csv(tp_file_path, index=False)
            print(f"TP results for region {region_name} saved to {tp_file_path}.")
        else:
            print(f"No TP results for region {region_name}. Skipping.")
            continue


    if all_results:
        final_results = pd.concat(all_results, ignore_index=True)
        output_file = os.path.join(output_dir, 'eval_results.csv')
        final_results.to_csv(output_file, index=False)
        print(f"Metrics and combined results saved successfully to {output_file}!")
    else:
        print("No results to save.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate mutations data and calculate performance metrics.")
    parser.add_argument("--data_file", required=True, help="Path to the input data file.")
    parser.add_argument("--answer_file", required=True, help="Path to the answer set file.")
    parser.add_argument("--output_dir", required=True, help="Path to save the output results.")
    parser.add_argument("--p_value", type=float, default=3.92e-05, help="P-value threshold")
    parser.add_argument("--OR", type=float, default=20, help="Odds Ratio threshold")

    args = parser.parse_args()
    conditions = {"p_value": args.p_value, "OR": args.OR}
    main(args.data_file, args.answer_file, args.output_dir, conditions)


    ##command line
    # python Evaluation.py --data_file ../Output/hypothesis_results.csv --answer_file ../Data/answer_set_pr50.csv --output_dir ../Output --p_value 3.92E-05 --OR 20
