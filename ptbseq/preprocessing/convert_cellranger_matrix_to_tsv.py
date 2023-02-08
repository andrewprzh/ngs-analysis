import csv
import gzip
import os
import scipy.io
import sys

# define MEX directory
matrix_dir = sys.argv[1]
# read in MEX format matrix as table
print("Loading matrix")
mat = scipy.io.mmread(os.path.join(matrix_dir, "matrix.mtx.gz"))
 
# list of transcript ids, e.g. 'ENSG00000243485'
print("Loading features")
features_path = os.path.join(matrix_dir, "features.tsv.gz")
feature_ids = [row[0] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]
 
# list of gene names, e.g. 'MIR1302-2HG'
gene_names = [row[1] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]
 
# list of feature_types, e.g. 'Gene Expression'
feature_types = [row[2] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]

print("Loading barcodes")
barcodes_path = os.path.join(matrix_dir, "barcodes.tsv.gz")
barcodes = [row[0] for row in csv.reader(gzip.open(barcodes_path, mode="rt"), delimiter="\t")]

import pandas as pd
 
# transform table to pandas dataframe and label rows and columns
print("Coverting to table")
matrix = pd.DataFrame.sparse.from_spmatrix(mat)
matrix.columns = barcodes
matrix.insert(loc=0, column="feature_id", value=feature_ids)
matrix.insert(loc=0, column="gene", value=gene_names)
matrix.insert(loc=0, column="feature_type", value=feature_types)
 
# display matrix
#print(matrix)
# save the table as a CSV (note the CSV will be a very large file)
print("Saving to disc")
matrix.to_csv("mex_matrix.csv", index=False, sep="\t")
