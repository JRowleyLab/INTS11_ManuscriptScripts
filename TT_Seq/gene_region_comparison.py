import pyBigWig
import math
import csv
import argparse

def compute_region_signal(bw, chrom, start, end):
    """Return the sum of signal over the given region."""
    if start < 0:  # avoid negative start positions
        start = 0
    try:
        values = bw.values(chrom, start, end, numpy=False)
        return sum(v for v in values if v is not None)
    except:
        return 0

def process_gene_row(row, bw):
    chrom, start, stop, name, _, strand = row
    start, stop = int(start), int(stop)

    if strand == "+":
        tss = start
        tes = stop
        early_start = tss - 250
        early_end = tss + 250
        upstream_start = tss - 5000
        upstream_end = tss - 250
        body_start = tss + 250
        body_end = tes
    else:
        tss = stop
        tes = start
        early_start = tss - 250
        early_end = tss + 250
        upstream_start = tss + 250
        upstream_end = tss + 5250
        body_start = tss + 250
        body_end = tes

    # Signal
    early_signal = compute_region_signal(bw, chrom, early_start, early_end)
    upstream_signal = compute_region_signal(bw, chrom, upstream_start, upstream_end)
    body_signal = compute_region_signal(bw, chrom, min(body_start, body_end), max(body_start, body_end))

    # Avoid division by zero
    upstream_ratio = upstream_signal / early_signal if early_signal > 0 else 0
    body_ratio = body_signal / early_signal if early_signal > 0 else 0

    return row + [round(upstream_ratio, 4), round(body_ratio, 4)]

def main(bedfile, bigwigfile, outfile):
    bw = pyBigWig.open(bigwigfile)
    output_rows = []

    with open(bedfile) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) < 6:
                continue
            new_row = process_gene_row(row, bw)
            output_rows.append(new_row)

    # Write output with new columns
    with open(outfile, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        for row in output_rows:
            writer.writerow(row)

    print(f"[✓] Output written to {outfile}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate upstream/early and body/early signal ratios")
    parser.add_argument("-b", "--bed", required=True, help="Input BED file")
    parser.add_argument("-w", "--bigwig", required=True, help="Signal bigWig file")
    parser.add_argument("-o", "--output", required=True, help="Output file path")
    args = parser.parse_args()

    main(args.bed, args.bigwig, args.output)
