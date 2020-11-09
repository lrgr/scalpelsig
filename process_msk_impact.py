# process msk impact bed file into a standardized form which can be ingested by my R scripts

import pandas as pd

infile = "/home/franzese/projects/hotspot_signature_panel/data/IMPACT_v2_expanded.bed"
outfile = "data/msk_impact_regions.tsv"

df = pd.read_csv(infile, sep="\t")

# for some reason the gencode Chromosome column had "chr" prepended to each entry and I did that first
# so my R scripts which ingest that df expects it to be there, so I have to put it here as well
df['Chromosome'] = "chr" + df['Chromosome'].astype(str)

formatted_regions = df[["Chromosome", "Start", "End"]]

formatted_regions.to_csv(outfile, sep="\t", index=False, index_label=False)
