#!/usr/bin/env python3
"""
Convert featureCounts tables into one TPM matrix (genes x samples).

Usage: counts_to_tpm.py sample1.featureCounts.tsv sample2.featureCounts.tsv ... > tpm_matrix.tsv

TPM for gene i in one sample:
    rate_i = count_i / (length_i / 1000)        # reads per kilobase
    TPM_i  = rate_i / sum(all rates) * 1e6      # normalise to one million
"""
import os
import sys


def read_featurecounts(path):
    """Return dicts gene_id -> length, gene_id -> gene_name, gene_id -> count."""
    lengths, names, counts = {}, {}, {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or line.startswith("Geneid"):
                continue  # skip command line comment and header
            fields = line.rstrip("\n").split("\t")
            # columns: Geneid, Chr, Start, End, Strand, Length, gene_name, count
            gene_id = fields[0]
            lengths[gene_id] = int(fields[5])
            names[gene_id] = fields[6]
            counts[gene_id] = int(fields[-1])
    return lengths, names, counts


def counts_to_tpm(counts, lengths):
    """Convert raw counts of one sample to TPM."""
    rates = {g: counts[g] / (lengths[g] / 1000) for g in counts}
    total = sum(rates.values())
    if total == 0:
        return {g: 0.0 for g in counts}  # no reads assigned at all
    return {g: rates[g] / total * 1e6 for g in rates}


samples, tpm_per_sample = [], []
for path in sorted(sys.argv[1:]):
    sample = os.path.basename(path).replace(".featureCounts.tsv", "")
    lengths, names, counts = read_featurecounts(path)
    samples.append(sample)
    tpm_per_sample.append(counts_to_tpm(counts, lengths))

# write one table: rows = genes, columns = samples
print("gene_id\tgene_name\t" + "\t".join(samples))
for gene_id in lengths:
    values = [f"{tpm[gene_id]:.4f}" for tpm in tpm_per_sample]
    print(f"{gene_id}\t{names[gene_id]}\t" + "\t".join(values))