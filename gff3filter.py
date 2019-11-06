#!/usr/bin/python

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -

# Todd D. Yoder
#
# Purpose: Filter GFF3 file to a single gene of interest

# #Modules/Libraries:--------------------------------------------------------------------------------------------------------------

import sys
import argparse
import pandas as pd


parser = argparse.ArgumentParser(description="Filter GFF3 by gene of interest")

parser.add_argument('-p','--path', help='path to configuration file...', required = True)
parser.add_argument('-f','--file', help='name of GFF3 file...[TSV]', required = True)
parser.add_argument('-q','--query', help='gene of interest...', required = True)

args = vars(parser.parse_args())

path = args['path']
file = args['file']
GOI = args['query']


check = []

with open(path + file) as myfile:
    for line in myfile:
        if line.startswith('#'):
            pass
        else:
            line=line.strip('\n')
            check.append(line)


gff = pd.DataFrame([sub.split("\t") for sub in check])

query = gff[gff.iloc[:,8].str.contains(GOI)]

chrom = query[0].iloc[0]
start = query[3].iloc[0]
stop = query.groupby([0])[4].transform('last').iloc[0]

gff3 = gff[gff.iloc[:,0] == chrom]
gff3 = gff3[gff3.iloc[:,3].astype(int).between(int(start), int(stop), inclusive=True)]

print(gff3)
gff3.to_csv(GOI + ".gff3", sep='\t', index=False,header=False)
