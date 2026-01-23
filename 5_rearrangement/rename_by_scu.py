import sys

def parse_scu_bed(scu_bed_file):
    with open(scu_bed_file, 'r') as scu_bed:
        scu_dict = {}
        chr2scu_dict = {}
        for line in scu_bed:
            cols = line.rstrip('\n').split('\t')
            scu, start, stop = cols[0], int(cols[1]), int(cols[2])
            chr = scu.split("_")[0]
            try:
                chr2scu_dict[chr].append(scu)
            except KeyError:
                chr2scu_dict[chr] = [scu]
            scu_dict[scu] = [start, stop]
        return scu_dict, chr2scu_dict
    
def rename_by_scu(oxford_table_file, sp1_scu_dict, sp1_chr2scu_dict, sp2_scu_dict, sp2_chr2scu_dict):
    with open(oxford_table_file, 'r') as oxford_table:
        for line in oxford_table:
            cols = line.rstrip('\n').split('\t')
            busco, nigon, chr1, pos1, chr2, pos2 = cols[0], cols[1], cols[2], float(cols[3]), cols[4], float(cols[5])
            # get scu of sp1 
            scu_list = sp1_chr2scu_dict[chr1]
            sp1_assigned_scu = '-'
            for scu in scu_list:
                start_stop = sp1_scu_dict[scu]
                if pos1 > start_stop[0] and pos1 < start_stop[1]:
                    sp1_assigned_scu = scu
            # get scu of sp2
            scu_list = sp2_chr2scu_dict[chr2]
            sp2_assigned_scu = '-'
            for scu in scu_list:
                start_stop = sp2_scu_dict[scu]
                if pos2 > start_stop[0] and pos2 < start_stop[1]:
                    sp2_assigned_scu = scu
            print(line.rstrip("\n") + "\t" + sp1_assigned_scu + "\t" + sp2_assigned_scu)
    

sp1_scu_bed_file = sys.argv[1]
sp2_scu_bed_file = sys.argv[2]
oxford_table_file = sys.argv[3]
sp1_scu_dict, sp1_chr2scu_dict = parse_scu_bed(sp1_scu_bed_file)
sp2_scu_dict, sp2_chr2scu_dict = parse_scu_bed(sp2_scu_bed_file)
rename_by_scu(oxford_table_file, sp1_scu_dict, sp1_chr2scu_dict, sp2_scu_dict, sp2_chr2scu_dict)