# BETH-2-Predicted Positive Mutations Through February 2026

This directory contains an update of **BETH-2-predicted positive Spike mutations** based on data available **from 2024-06 to 2026-02**.

The file included here is:

- `BETH-2-predicted_positive_Spike_mutations.csv`

Each mutation is listed only once, at the **first month** in which it was identified as a BETH-2-predicted positive mutation in this updated analysis.

## Column Definitions

- `Mutation`: Spike amino-acid mutation identified as a BETH-2-predicted positive mutation in the updated analysis.
- `Site`: amino-acid position of that Spike mutation.
- `Detection Date`: first month (`YYYY-MM`) in which the mutation was classified as BETH-2-predicted positive.
- `Detection Region`: region in which the mutation first met the BETH-2 positivity criteria.
- `p value`: p value at the first month when the mutation became BETH-2-predicted positive.
- `logOR`: base-10 logarithm of the odds ratio at the first month when the mutation became BETH-2-predicted positive.

## Interpretation

This update follows the same framework used in the study to define BETH-2-predicted positive mutations, based on the regional comparison across South Africa, California, and India.

A mutation was retained in this table only when it first satisfied all of the following:

- `p value < 3.92e-05`
- `logOR > 1` (equivalent to `OR > 10`)
Here, `logOR` is defined as `log10(OR)`.

If a mutation satisfied these criteria in multiple months, only the earliest month was kept in the table so that each mutation appears once as its first BETH-2 positive signal.
