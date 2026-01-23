import sys

def parse_orthogroups(file_path):
    """Parses the orthogroup file and returns dictionaries mapping sequences to orthogroups and orthogroups to sequences."""
    seq2ortho_dict = {}
    ortho_dict = {}

    with open(file_path, 'r') as orthogroups:
        for line in orthogroups:
            parts = line.rstrip("\n").split(": ")
            if len(parts) < 2:
                continue  # Skip malformed lines
            ortho, seqlist = parts[0], parts[1].split()
            for seq in seqlist:
                seq2ortho_dict[seq] = ortho
            ortho_dict[ortho] = seqlist

    return seq2ortho_dict, ortho_dict

def parse_eliminated(file_path):
    """Parses the eliminated genes file and returns a set of eliminated sequences."""
    with open(file_path, 'r') as eliminated:
        return {line.strip() for line in eliminated}

def parse_seq2locus(file_path):
    """Parses the seq2locus file and returns a dict of seq:locus."""
    with open(file_path, 'r') as seq2locus:
         seq2locus_dict = {}
         for line in seq2locus:
             seq, locus = line.split("\t")[0], line.split("\t")[1].rstrip("\n")
             seq2locus_dict[seq] = locus
    return seq2locus_dict

def parse_interproscan(file_path):
    """Parses the inteproscan file and returns a dict of seq:domain."""
    with open(file_path, 'r') as interproscan:
         interproscan_dict = {}
         for line in interproscan:
             seq, desc = line.split("\t")[0], line.split("\t")[5]
             try:
                 interproscan_dict[seq].append(desc)
             except KeyError:
                 interproscan_dict[seq] = [desc]
    return interproscan_dict

def parse_gff3(file_path):
    """Parses GFFs and returns a dict of seq:[chr, start, stop]"""
    gff3_dict = {}
    with open(file_path, 'r') as gff3:
        for line in gff3:
            if not line.startswith("#"):
               cols = line.rstrip("\n").split("\t")
               chr, start, stop = cols[0], cols[3], cols[4]
               if cols[2] == 'mRNA':
                   id = cols[8].split(";")[0].split("=")[1]
                   gff3_dict[id] = [chr, start, stop]
    return gff3_dict

def analyze_eliminated(seq, seq2ortho_dict, ortho_dict, eliminated_list, prefix, seq2locus_dict, interproscan_dict, gff3_dict):
    """Analyzes the eliminated gene and returns its associated orthogroup statistics."""
    if seq not in seq2ortho_dict:
        return None  # Sequence not found in orthogroup data
    ortho = seq2ortho_dict[seq]
    seqlist = ortho_dict.get(ortho, [])
    seq_count = len(seqlist)
    prefix_count = sum(1 for item in seqlist if item.startswith(prefix))
    eleg_list = [item for item in seqlist if item.startswith("CELEG")]
    locus_list = []
    for item in eleg_list:
        try:
            locus = seq2locus_dict[item.replace('CELEG.', '')]
            locus_list.append(locus)
        except KeyError:
            pass
    # get elegans info
    eleg_string = ";".join(eleg_list).replace('CELEG.', '')
    locus_string = ";".join(locus_list)
    eleg_count = len(eleg_list)
    # get ortho summary
    elim_count = sum(1 for item in seqlist if item in eliminated_list)
    retained_status = '-'
    if prefix_count == 1:
        copy_status = 'single'
    else:
        copy_status = 'multi'
        if prefix_count == elim_count:
             retained_status = 'eliminated_only'
        else:
             retained_status = 'copy_in_retained'
    # get interproscan output
    try:
        interproscan_string = ";".join(interproscan_dict[seq])
    except KeyError:
        interproscan_string = ''
    # get location info
    try:
        chr, start, stop = gff3_dict[seq]
    except KeyError:
        sys.exit("ERROR: " + seq + " not found in GFF3")
    return seq, chr, start, stop, ortho, copy_status, retained_status, seq_count, prefix_count, eleg_count, elim_count, eleg_string, locus_string, interproscan_string

def main():
    if len(sys.argv) != 7:
        print("Usage: python script.py <orthogroup_file> <eliminated_file> <prefix> <seq2locus_file> <interproscan_file> <gff3_file")
        sys.exit(1)

    orthogroup_file = sys.argv[1]
    eliminated_file = sys.argv[2]
    prefix = sys.argv[3]
    seq2locus_file = sys.argv[4]
    interproscan_file = sys.argv[5]
    gff3_file = sys.argv[6]

    seq2ortho_dict, ortho_dict = parse_orthogroups(orthogroup_file)
    eliminated_list = parse_eliminated(eliminated_file)
    seq2locus_dict = parse_seq2locus(seq2locus_file)
    interproscan_dict = parse_interproscan(interproscan_file)
    gff3_dict = parse_gff3(gff3_file)

    for seq in eliminated_list:
        result = analyze_eliminated(seq, seq2ortho_dict, ortho_dict, eliminated_list, prefix, seq2locus_dict, interproscan_dict, gff3_dict)
        if result:
            print("\t".join(map(str, result)))
        else:
            print(f"Warning: {seq} not found in orthogroups", file=sys.stderr)

if __name__ == "__main__":
    main()
