import argparse
def get_tsv(file):
    
    for line in open(file):
        line=line.strip()
        if not line or line.startswith("#"):
            continue
        yield line.split("\t")
        
def get_blast_dict(file):
    blast_dict={}
    for line in get_tsv(file):
        if line[0] not in blast_dict:
            blast_dict[line[0]]=line
    return blast_dict

def get_taxonomy(m6_file,taxonomy,prefix):
    taxid_dict={}
    blast_dict=get_blast_dict(m6_file)
    tax_file=open("%s.taxonomy.xls" % prefix,"w")
    header="\t".join(["qseqid","qlen","qstart","qend","stitle","sstart","send","pident","length","mismatch","gapopen","evalue","bitscore","qcovhsp","qstrand","taxid","superkindom","phylum","class","order","family","genus","species"])
    tax_file.write(header+"\n")
    for line in blast_dict.values():
        taxid=line[-1].split(";")[-1]
        if taxid not in taxid_dict:
            taxid_dict[taxid] = []
        taxid_dict[taxid].append([line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[-1]])

    for line in get_tsv(taxonomy):
        if line[0] in taxid_dict:
            for i in taxid_dict[line[0]]:
                superkindom=line[18] or "Unclassified"
                kindom=line[16] or "Unclassified"
                phylum=line[14] or "Unclassified"
                class_=line[12] or "Unclassified"
                order=line[10] or "Unclassified"
                family=line[8] or "Unclassified"
                genus=line[6] or "Unclassified"
                species=line[4] or "Unclassified"
                tax_file.write("\t".join([i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7],i[8],i[9],i[10],i[11],i[12],i[13],i[14],i[15],superkindom,phylum,class_,order,family,genus,species,"\n"]))
    tax_file.close()

def add_args(parser):
    parser.add_argument('-i', '--input', metavar='m6', type=str, required=True,help='Input balst comparison file.')
    parser.add_argument('-t', '--taxonomy', metavar='FILE', type=str)
    parser.add_argument('-n', '--name', metavar='STR', type=str, default= 'out',help='Output file prefix')
    return parser

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    args = add_args(parser).parse_args()
    get_taxonomy(args.input, args.taxonomy, args.name)

if __name__ == "__main__":
    main()
