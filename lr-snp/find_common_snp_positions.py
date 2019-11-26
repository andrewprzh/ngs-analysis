import os
import sys
import glob
import pysam

REGION_SIZE = 100

def read_vcf(vcf_file):
    positions = {}
    for l in open("vcf_file"):
        if l.startswith("#"):
            continue

        tokens = l.split()
        if len(tokens) < 2:
            continue

        chr_id = tokens[0]
        if chr_id not in positions:
            positions[chr_id] = []

        var_pos = int(tokens[1])
        if len(positions[chr_id]) > 0 and var_pos - positions[chr_id] < REGION_SIZE:
            continue
        positions[chr_id].append(var_pos)

    return positions


def get_region_list(positions):
    region_list = []
    for chr_id in positions:
        for pos in positions[chr_id]:
            region_list.append("chr" + chr_id + ":" + str(max(0, pos - REGION_SIZE)) + ":" + str(pos + REGION_SIZE))
    return region_list


DIR_10X = "/Bmo/prjbel/RNA_10x/SNP/data/10X_ShortRead_BAM/SC/EN1/"
DIR_ISOSEQ = "/Bmo/prjbel/RNA_10x/SNP/data/IsoSeq_GR/ExitedNeuron1/EN1/"
PILEUP_CMD = "samtools mpileup -I  --reference /Bmo/prjbel/RNA_10x/references/mouse/Mus_musculus.GRCm38.75.dna.cleared.short.fa -g "
BCF_CMD = "/Bmo/prjbel/tools/bcftools-1.9/bcftools call -v -A -m --ploidy 1 -o CATGCCTTCTGCTGCT.vcf CATGCCTTCTGCTGCT.pileup.bcf "
BCF_VIEW = "/Bmo/prjbel/tools/bcftools-1.9/bcftools view "

for barcode_bam in glob.glob(os.path.join(DIR_10X, "*.bam")):
    barcode = os.path.splitext(os.path.basename(barcode_bam))[0]
    os.system(PILEUP_CMD + " " + barcode_bam + " > ./pileup/" + barcode + ".pileup.bcf")
    os.system(BCF_CMD + " ./pileup/" + barcode + ".pileup.bcf -o ./vcf/" + barcode + ".bcf")
    os.system(BCF_VIEW + "  ./vcf/" + barcode + ".bcf > tmp.vcf")

    region_list = get_region_list(read_vcf("tmp.vcf"))
    isoseq_bam_file = os.path.join(DIR_ISOSEQ, barcode + ".bam")
    if not os.path.exists(isoseq_bam_file):
        print("WARN: " + isoseq_bam_file + " does not exist")

    isoseq_bam = pysam.AlignmentFile(isoseq_bam_file, "rb")
    for region in region_list:
        reads_in_region = []
        for read in isoseq_bam.fetch(region):
            reads_in_region.append(read)
        if len(reads_in_region) == 0:
            continue

        print("Found covered variant: barcode " + barcode + ", region " + region)
        out_file = pysam.AlignmentFile("./bam_regions/" + barcode + "_" + region.replace(":", "_") + ".bam", "wb", template=isoseq_bam)
        for read in reads_in_region:
            out_file.write(read)
            out_file.close()

    isoseq_bam.close()
    break


