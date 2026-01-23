import pandas as pd
import numpy as np
import sys

bed_file = sys.argv[1]

# Read BED file
bed_df = pd.read_csv(bed_file, sep="\t", header=None, names=["chr", "start", "stop", "dot", "score", "strand", "source", "type", "dot2", "col9"])

# get name 
bed_df["name"] = bed_df["col9"].str.extract(r'ID=([^;]+)')

# Define threshold
distance_threshold = 100 

# 1. Generate all pairwise comparisons **within the same chromosome** (avoids looping)
bed_pairs = bed_df.merge(bed_df, on="chr", suffixes=("_ref", "_comp"))

# 2. Remove self-comparisons
bed_pairs = bed_pairs[bed_pairs["name_ref"] != bed_pairs["name_comp"]]

# 3. Ensure MEME-2 is always ref and MEME-1 is always comp
bed_pairs = bed_pairs[
    (bed_pairs["name_ref"].str.contains("MEME-2")) & 
    (bed_pairs["name_comp"].str.contains("MEME-1"))
]

# 4. Compute distances **vectorized** (no loops!)
bed_pairs["start2stop_distance"] = bed_pairs["start_comp"] - bed_pairs["stop_ref"]  # Raw distance
bed_pairs["stop2start_distance"] = bed_pairs["stop_comp"] - bed_pairs["start_ref"]
bed_pairs["min_value"] = bed_pairs[["start_ref", "stop_ref", "start_comp", "stop_comp"]].min(axis=1)
bed_pairs["max_value"] = bed_pairs[["start_ref", "stop_ref", "start_comp", "stop_comp"]].max(axis=1)

# 4.5 Remove overlapping features
bed_pairs = bed_pairs[~((bed_pairs["start_comp"] <= bed_pairs["stop_ref"]) & 
                        (bed_pairs["start_ref"] <= bed_pairs["stop_comp"]))]


# 5. Apply strand filtering (upstream for + strand, downstream for - strand)
bed_pairs = bed_pairs[
    ((bed_pairs["strand_ref"] == "+") & (bed_pairs["start2stop_distance"] < 0)) |  # Upstream for positive strand
    ((bed_pairs["strand_ref"] == "-") & (bed_pairs["stop2start_distance"] > 0))    # Downstream for negative strand
]

# 6. Compute sum of scores (ref + comp) and get flipped status
bed_pairs["score_sum"] = bed_pairs["score_ref"] + bed_pairs["score_comp"]
bed_pairs["status"] = np.where(bed_pairs["strand_ref"] == bed_pairs["strand_comp"], "normal", "flipped")

# 7. Apply the distance threshold filter **efficiently**
bed_pairs["distance"] = np.minimum(abs(bed_pairs["start2stop_distance"]), abs(bed_pairs["stop2start_distance"]))
filtered_pairs = bed_pairs.query("distance < @distance_threshold")

# 9. Exclude rows where start and end coordinates are identical
filtered_pairs = filtered_pairs[~((filtered_pairs["start_ref"] == filtered_pairs["start_comp"]) & 
                                  (filtered_pairs["stop_ref"] == filtered_pairs["stop_comp"]))]

# 10. Select only the relevant columns
filtered_pairs = filtered_pairs[[
    "chr", "min_value", "max_value", "status", "score_sum", "strand_ref", "distance",
    "start_ref", "stop_ref", "name_ref", "score_ref", "strand_ref",
    "start_comp", "stop_comp", "name_comp", "score_comp", "strand_comp",
]]

# 11. Print results without row and column names in TSV format
print(filtered_pairs.to_csv(sep='\t', index=False, header=False))

