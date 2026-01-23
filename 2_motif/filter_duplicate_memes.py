import polars as pl
import sys

def filter_bed_file(file_path):
    # Define column names
    columns = [
        "chr", "min_value", "max_value", "status", "score_sum", "strand", "distance",
        "start_ref", "stop_ref", "name_ref", "score_ref", "strand_ref",
        "start_comp", "stop_comp", "name_comp", "score_comp", "strand_comp"
    ]

    # Read TSV with Polars (lazy scan to minimize memory)
    df = pl.scan_csv(file_path, separator="\t", has_header=False)
    df = df.rename({f"column_{i+1}": col for i, col in enumerate(columns)})

    # Collect after sorting
    df = df.sort("score_sum", descending=True).collect()

    used = set()
    kept_rows = []

    # Convert to Python tuples (fast iteration)
    for row in df.iter_rows(named=True):
        ref_id = (row["start_ref"], row["stop_ref"], row["strand_ref"])
        comp_id = (row["start_comp"], row["stop_comp"], row["strand_comp"])
        if ref_id not in used and comp_id not in used:
            kept_rows.append(row)
            used.update([ref_id, comp_id])

    # Convert back to Polars DataFrame
    filtered = pl.DataFrame(kept_rows)
    return filtered

if __name__ == "__main__":
    input_file = sys.argv[1]
    result = filter_bed_file(input_file)
    result.write_csv(sys.stdout, separator="\t", include_header=False)

