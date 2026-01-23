import sys

def parse_fai_longest_six(fai_file):
    """Parse the .fai file and return the longest six chromosome lengths."""
    chrom_lengths = {}
    with open(fai_file, "r") as f:
        for line in f:
            parts = line.strip().split('\t')
            chrom = parts[0]
            length = int(parts[1])
            chrom_lengths[chrom] = length

    # Sort chromosomes by length (descending) and pick the top 6
    longest_six = dict(sorted(chrom_lengths.items(), key=lambda x: x[1], reverse=True)[:6])
    return longest_six

def generate_bed(chrom_lengths):
    """Generate BED output for the first and last 2.69% of each sequence."""
    for chrom, length in chrom_lengths.items():
        region_length = int(length * 0.0269)
        print(f"{chrom}\t0\t{region_length}\t{chrom}_left")
        print(f"{chrom}\t{length - region_length}\t{length}\t{chrom}_right")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <fai_file>")
        sys.exit(1)

    fai_file = sys.argv[1]
    chrom_lengths = parse_fai_longest_six(fai_file)
    generate_bed(chrom_lengths)

