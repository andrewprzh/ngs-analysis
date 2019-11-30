import os
import sys
import glob
import pysam

REGION_SIZE = 100

def read_vcf(vcf_file):
    positions = {}
    for l in open(vcf_file):
        if l.startswith("#"):
            continue

        tokens = l.split()
        if len(tokens) < 2:
            continue

        chr_id = tokens[0]
        if chr_id not in positions:
            positions[chr_id] = []

        var_pos = int(tokens[1])
        if len(positions[chr_id]) > 0 and var_pos - positions[chr_id][-1] < REGION_SIZE:
            continue
        positions[chr_id].append(var_pos)

    return positions


def get_region_list(positions):
    region_list = []
    for chr_id in positions:
        for pos in positions[chr_id]:
            region_list.append(chr_id + ":" + str(max(0, pos - REGION_SIZE)) + "-" + str(pos + REGION_SIZE))
    return region_list


DIR_10X = "/Bmo/prjbel/RNA_10x/SNP/data/10X_ShortRead_BAM/SC/EN1/"
DIR_ISOSEQ = "/Bmo/prjbel/RNA_10x/SNP/data/IsoSeq_GR/ExitedNeuron1/EN1/"
PILEUP_CMD = "samtools mpileup -I  --reference /Bmo/prjbel/RNA_10x/references/mouse/Mus_musculus.GRCm38.75.dna.cleared.short.fa -g "
PILEUP_ISOSEQ_CMD = "samtools mpileup -I  --reference /Bmo/prjbel/RNA_10x/references/mouse/GRCm38.p6.genome.chr.fa -g "

BCF_CMD = "/Bmo/prjbel/tools/bcftools-1.9/bcftools call -O v -v -A -m --ploidy 1 "
BCF_VIEW = "/Bmo/prjbel/tools/bcftools-1.9/bcftools view "

SPLIT_ISO = "python /Bmo/prjbel/ngs-analysis-2/lr-snp/split_isoseq_with_table.py "
SPLIT_10X = "python /Bmo/prjbel/ngs-analysis-2/lr-snp/split_10x.py "

count = 0
for barcode_bam in glob.glob(os.path.join(DIR_10X, "*.bam"))[200:210]:
    count += 1
    print("Processing " + barcode_bam)
    barcode = os.path.splitext(os.path.basename(barcode_bam))[0]
    os.system(PILEUP_CMD + " " + barcode_bam + " > ./pileup/" + barcode + ".pileup.bcf 2> /dev/null")
    os.system(BCF_CMD + " -o ./vcf/" + barcode + ".vcf ./pileup/" + barcode + ".pileup.bcf")
#    os.system(BCF_VIEW + "  ./vcf/" + barcode + ".bcf > ./vcf/" + barcode + ".vcf")

    region_list = get_region_list(read_vcf("./vcf/" + barcode + ".vcf"))
    isoseq_bam_file = os.path.join(DIR_ISOSEQ, barcode + ".bam")
    if not os.path.exists(isoseq_bam_file):
        print("WARN: " + isoseq_bam_file + " does not exist")

    print("Checking " + isoseq_bam_file + ", total regions: " + str(len(region_list)))
    isoseq_bam = pysam.AlignmentFile(isoseq_bam_file, "rb")
    for region in region_list:
        reads_in_region = []
        if region[0] not in set(['0','1','2','3','4','5','6','7','8','9']):
            continue
        for read in isoseq_bam.fetch(region="chr"+region):
            reads_in_region.append(read)
        if len(reads_in_region) == 0:
            continue

        print("Found covered variant: barcode " + barcode + ", region " + region)
        out_region_file = barcode + "_" + region.replace(":", "_") + ".bam"
        out_file = pysam.AlignmentFile("./isoseq_bam_regions/" + out_region_file, "wb", template=isoseq_bam)
        for read in reads_in_region:
            out_file.write(read)
        out_file.close()

        os.system(PILEUP_ISOSEQ_CMD + " " + "./isoseq_bam_regions/" + out_region_file + " > ./isoseq_bam_regions/" + out_region_file + ".pileup.bcf 2> /dev/null")
        os.system(BCF_CMD + " -O v -o " + "./isoseq_bam_regions/" + out_region_file + ".vcf ./isoseq_bam_regions/" + out_region_file + ".pileup.bcf")
        umi_dir = './umis/' + barcode + '/' + region.replace(":", "_") + '/isoseq/'
        if not os.path.isdir(umi_dir):
            os.makedirs(umi_dir)
        os.system(SPLIT_ISO + " --bam " + "./isoseq_bam_regions/" + out_region_file + " -o  " + umi_dir + " --table /Bmo/prjbel/RNA_10x/SNP/data/BAMFiles_ShortReadClusters/umi/ExcitNeuron1_BC_UMI --split_by UMI > /dev/null")
        umis_isoseq = set(map(lambda x: os.path.splitext(os.path.basename(x))[0], glob.glob(os.path.join(umi_dir, "*.bam"))))

        os.system("samtools view -b " + barcode_bam + " " + region + " > ./10x_bam_regions/" + out_region_file)
        os.system(PILEUP_CMD + " " + "./10x_bam_regions/" + out_region_file + " > ./10x_bam_regions/" + out_region_file + ".pileup.bcf 2> /dev/null")
        os.system(BCF_CMD + " -O v -o " + "./10x_bam_regions/" + out_region_file + ".vcf ./10x_bam_regions/" + out_region_file + ".pileup.bcf")
        umi_dir = './umis/' + barcode + '/' + region.replace(":", "_") + '/10x/'
        if not os.path.isdir(umi_dir):
            os.makedirs(umi_dir)
        os.system(SPLIT_10X + " --bam " + "./10x_bam_regions/" + out_region_file + " -o  " + umi_dir + " --split_by UMI > /dev/null")
        umis_10x = set(map(lambda x: os.path.splitext(os.path.basename(x))[0], glob.glob(os.path.join(umi_dir, "*.bam"))))

        umis_common = umis_isoseq.intersection(umis_10x)
        print("10x UMIs: " + str(len(umis_10x)) + ",Iso-seq UMIs: " + str(len(umis_isoseq)) + ", common: " +  str(len(umis_common)))
        if len(umis_common) > 0:
            print("Common UMIs:\n" + "\n".join(list(umis_common)))

    isoseq_bam.close()
    if count > 20:
        break


