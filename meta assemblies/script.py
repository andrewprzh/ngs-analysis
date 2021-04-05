#!/usr/bin/python
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import pysam as ps
import gffutils
import sys
import argparse
import os
from traceback import print_exc

def parse_sam(input_align):
    align = ps.AlignmentFile(input_align, "r")
    table=[]
    for read in align.fetch():
        if read.flag!=4:
            nodes=[]
            name=read.qname.split('_')[0]+'_'+read.qname.split('_')[1]
            nodes.append(name)
            nodes.append(read.reference_name)
            nodes.append(int(read.qlen))
            nodes.append(int(read.qstart))
            nodes.append(int(read.qend))
            nodes.append(int(read.reference_start))
            nodes.append(int(read.reference_end))
            nodes.append(read.flag)
            table.append(nodes)
    align.close()
    
    df=pd.DataFrame(data=table, columns=['Name','Node', 'Len', 'Qstart', 'Qstop', 'Rstart', 'Rstop', 'Flag'])
    df.sort_values(by=['Node','Rstart', 'Rstop'], inplace = True)
    return df

def check (input_df, start, stop):

    light=False
    
    for i in range(1,len(input_df)):
        if light: 
            if input_df.Rstart.iloc[i] < input_df.Rstop.iloc[i-1] or stop > input_df.Rstart.iloc[i]:
                    loc.add(i-1)
                    loc.add(i)
                    if stop<=input_df.Rstop.iloc[i]:
                        length= (sum(input_df.Rstop.iloc[list(loc)].tolist()[:-1])+stop - sum(input_df.Rstart.iloc[list(loc)].tolist()[1:])+start)/(stop-start)*100
                        light=False
                        print(input_df.Name.iloc[list(loc)].tolist())
                        return length, input_df.Name.iloc[list(loc)].tolist(), 'part'
                    
        else: 
            loc=set()
            
        #skip genes 
        if start>input_df.Rstop.iloc[i-1] or stop<input_df.Rstart.iloc[i-1]:
            continue
        
        if start>=input_df.Rstart.iloc[i-1]:
            if stop<=input_df.Rstop.iloc[i-1]:
                return 100, input_df.Name.iloc[i-1], 'full'
            
            elif input_df.Rstart.iloc[i] < input_df.Rstop.iloc[i-1]: #or stop > input_df.Rstart.iloc[i]:
                    loc.add(i-1)
                    loc.add(i)
                    
                    if stop<=input_df.Rstop.iloc[i]:
                        return 100, input_df.Name.iloc[list(loc)].tolist(), 'full'
                    
            elif stop > input_df.Rstart.iloc[i]:
                    loc.add(i-1)
                    loc.add(i)
                    
                    if stop<=input_df.Rstop.iloc[i]:
                        length=(stop+input_df.Rstop.iloc[i-1] - start - input_df.Rstart.iloc[i])/(stop-start)*100
                        return length, input_df.Name.iloc[list(loc)].tolist(), 'part'
            
                    else: light=True
                    
        if start<input_df.Rstart.iloc[i-1]:
            if stop<=input_df.Rstop.iloc[i-1]:
                return (stop - input_df.Rstart.iloc[i-1])/(stop-start)*100, input_df.Name.iloc[i-1], 'part'
            
            elif input_df.Rstart.iloc[i] < input_df.Rstop.iloc[i-1] or stop > input_df.Rstart.iloc[i]:
                    loc.add(i-1)
                    loc.add(i)
                    
                    if stop>input_df.Rstop.iloc[i]:
                        return (max(input_df.Rstop.iloc[list(loc)].tolist()) - min(input_df.Rstart.iloc[list(loc)].tolist()))/(stop-start)*100, input_df.Name.iloc[list(loc)].tolist(), 'part'
       
        if start>input_df.Rstart.iloc[i-1] and stop>input_df.Rstop.iloc[i-1]:
            return (input_df.Rstop.iloc[i-1] - start)/(stop-start)*100, input_df.Name.iloc[i-1], 'part'
        
                
            
    return 0, '', ''

def get_result(df, db):
    result=[]
    for i in db.features_of_type('CDS'):
        gene = i
        start = gene.start
        stop = gene.stop
        c_name = gene.seqid
    
        subset = df[df.Node == c_name]
        cov, r_name, type_gene = check(subset, start, stop)
   
        if cov != 0:
            result.append([c_name, r_name, start, stop, cov, gene.id, gene.strand])
    
    return result

def minimap2(contigs, transcripts):
    print("Running minimap2")
    os.system("minimap2 -a  {path_to_contigs} {path_to_transcripts}> align.sam".format(path_to_transcripts = transcripts, path_to_contigs = contigs))

def prokka(contigs):
    print("Running prokka")
    os.system("prokka {path_to_contigs} --outdir ./out_prokka --prefix PROKKA".format(path_to_contigs = contigs))


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="Output file name")
    parser.add_argument("--sam", "-f", type=str, help="initial file with alignment")
    parser.add_argument("--gff", "-g", type=str, help="initial file with annotation")
    parser.add_argument("--genome", type=str, help="initial file with genome assembly")
    parser.add_argument("--transcriptome", type=str, help="initial file with transcriptome assembly")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    sam = args.sam
    gff = args.gff
    if args.genome is not None and args.transcriptome is not None:
        minimap2(args.genome, args.transcriptome)
        sam = "./align.sam"
    if args.genome is not None:
        prokka(args.genome)
        gff = "./out_prokka/PROKKA.gff"

    gffutils.create_db(gff, 'db')
    db = gffutils.FeatureDB('db', keep_order=True)
    all_genes = get_result(parse_sam(sam), db)
    df_results=pd.DataFrame(data=all_genes, columns=['Node', 'Name', 'Start', 'Stop', 'Cov', 'Gene', 'Strand'])
    df_results.to_csv("./result.csv")



if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)


