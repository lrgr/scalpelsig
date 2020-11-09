# script to extract exome coordinates from GENCODE annotation file
# ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/

import pandas as pd

infile = "/home/franzese/projects/hotspot_signature_panel/data/gencode.v33.basic.annotation.gff3"
outfile = "data/gencode_exons.tsv"


df = pd.read_csv(infile, sep="\t", names=["Chromosome", "Source", "Type", "Start", "End", "Ignore6", "Ignore7", "Ignore8", "Ignore9"], skiprows=[0,1,2,3,4,5,6])

exons = df.loc[df["Type"]=="exon"]

formatted_exons = exons[["Chromosome", "Start", "End"]]

formatted_exons.to_csv(outfile, sep="\t", index=False, index_label=False)

