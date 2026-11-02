#!/usr/bin/env python3
# Combine ROHan .hmmrohl files into one table, convert to PLINK .hom format.
# NOT EXECUTED - test on two or three files first.
import argparse
import os
import re
import sys
import yaml

parser = argparse.ArgumentParser(description="Combine ROHan .hmmrohl files and convert to PLINK .hom format")
parser.add_argument('--rohan', required=True, nargs='+', help='ROHan .hmmrohl files')
parser.add_argument('--out', default='rohan.hom', help='Output .hom file')
parser.add_argument('--names', default='sample_names.txt', help='List of sample names')
parser.add_argument('--version', required=True, help='Nextflow process version')
args = parser.parse_args()

# Step 1: file name without the path or the .hmmrohl suffix
def sample_name(path):
    return re.sub(r'\.hmmrohl$', '', os.path.basename(path))

names = [sample_name(f) for f in args.rohan]

# Stop early if two files reduce to the same name
dupes = sorted({n for n in names if names.count(n) > 1})
if dupes:
    sys.exit(f"ERROR: duplicate sample names: {', '.join(dupes)}")

with open(args.names, 'w') as fh:
    fh.write("\n".join(names) + "\n")

# Convert to .hom format: FID IID PHE CHR SNP1 SNP2 POS1 POS2 KB NSNP DENSITY PHOM PHET
# ROHan input: SAMPLE CHROM BEGIN END ROH_LENGTH VALIDATED_SITES
# Map: FID=SAMPLE, IID=SAMPLE, PHE=-9, CHR=CHROM, SNP1=".", SNP2=".", 
#      POS1=BEGIN, POS2=END, KB=ROH_LENGTH/1000, NSNP=".", DENSITY=".", 
#      PHOM=".", PHET="."

header_written = False
with open(args.out, 'w') as out:
    # Write .hom header once
    hom_header = ["FID", "IID", "PHE", "CHR", "SNP1", "SNP2", "POS1", "POS2", 
                  "KB", "NSNP", "DENSITY", "PHOM", "PHET"]
    out.write("\t".join(hom_header) + "\n")
    
    for name, path in zip(names, args.rohan):
        with open(path) as fh:
            for line in fh:
                fields = line.split()            # handles mixed spaces/tabs
                if not fields:
                    continue                     # skip blank lines
                if fields[0].startswith('#'):    # skip header lines
                    continue
                
                # ROHan row: SAMPLE CHROM BEGIN END ROH_LENGTH VALIDATED_SITES
                # (SAMPLE is already added, so fields are: CHROM BEGIN END ROH_LENGTH VALIDATED_SITES)
                if len(fields) < 5:
                    sys.exit(f"ERROR: expected 5 columns in {path}: {line.strip()}")
                
                sample_id = name
                chrom = fields[1]
                begin = int(fields[2])
                end = int(fields[3])
                roh_length = int(fields[4])
                
                kb = round(roh_length / 1000.0, 3)
                
                hom_row = [
                    sample_id,           # FID
                    sample_id,           # IID
                    "-9",                # PHE
                    chrom,               # CHR
                    ".",                 # SNP1
                    ".",                 # SNP2
                    begin,               # POS1
                    end,                 # POS2
                    kb,                  # KB
                    ".",                 # NSNP
                    ".",                 # DENSITY
                    ".",                 # PHOM
                    "."                  # PHET
                ]
                out.write("\t".join(map(str, hom_row)) + "\n")

data = {
    'Process version': {'rohan_combine': args.version},
    'Tool version': {'python': sys.version.split()[0]}
}
with open('versions.yml', 'w') as fh:
    yaml.dump(data, fh)