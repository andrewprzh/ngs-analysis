import pandas as pd
import argparse

def merge_tsv_files(barcoded_reads_files, barcode_to_cell_file, cell_types_file, output_file):
    # Read the TSV files into DataFrames, skipping the first line (header)
    barcoded_reads = pd.concat((pd.read_csv(barcoded_reads_file, sep='\t', skiprows=1, header=None, names=['read_id', 'barcode_id', 'UMI', 'BC_score', 'valid_UMI', 'strand', 'polyT_start', 'primer_end', 'linker_start', 'linker_end', 'TSO5']) for barcoded_reads_file in barcoded_reads_files))
    barcode_to_cell = pd.read_csv(barcode_to_cell_file, sep='\t', skiprows=1, header=None, names=['cell_origin', 'barcode_id', 'cell_id'])
    cell_types = pd.read_csv(cell_types_file, sep='\t', skiprows=1, header=None, names=['cell_id', 'spot_type', 'cell_type1', 'cell_type2', 'region', 'subregion'])

    # Merge the DataFrames
    print(barcoded_reads.head(n=10).to_string(index=False))
    print(barcode_to_cell.head(n=10).to_string(index=False))
    merged_df = barcoded_reads.merge(barcode_to_cell, on='barcode_id')
    merged_df = merged_df.merge(cell_types, on='cell_id')

    # Create the new column based on the spot_type
    def create_new_column(row):
        if row['spot_type'] == 'singlet':
            return f"{row['cell_type1']}::{row['region']}::{row['subregion']}"
#        elif row['spot_type'].startswith('doublet'):
#            return f"{row['cell_type1']}_{row['cell_type2']}::{row['region']}::{row['subregion']}"
        else:
            return None

    merged_df['group_id'] = merged_df.apply(create_new_column, axis=1)

    # Select the required columns and write to the output file
    output_df = merged_df[merged_df['group_id'].notna()][['read_id', 'group_id']]
    output_df.to_csv(output_file, sep='\t', index=False, header=False)

def main():
    parser = argparse.ArgumentParser(description='Merge TSV files and output a new TSV file.')
    parser.add_argument('--barcoded_reads', nargs="+", type=str, help='Path to the barcoded reads TSV file.', required=True)
    parser.add_argument('--barcode_to_cell', type=str, help='Path to the barcode to cell ID TSV file.', required=True)
    parser.add_argument('--cell_types', type=str, help='Path to the cell types TSV file.', required=True)
    parser.add_argument('--output', type=str, help='Path to the output TSV file.', required=True)

    args = parser.parse_args()

    merge_tsv_files(args.barcoded_reads, args.barcode_to_cell, args.cell_types, args.output)

if __name__ == '__main__':
    main()
