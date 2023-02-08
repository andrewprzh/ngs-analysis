#!/bin/sh
# By Anoushka Joglekar 06.2022
# File to convert from IsoQuant corrected BED file to
# a corrected all-info file
# Column 6: Corrected Intron Chain
# Column 9: Corrected Exon Chain
# Column 10: Old intron chain
# Column 11: Old exon chain

# ReadID, Gene, Celltype, Barcode, UMI, IntronChain,  TSS, PolyA, ExonChain
# 0289e424-346b-43f6-befa-d9ec79950b42    ENSMUSG00000025417.8    P56::M1::Cerebellum::Glia::Oligo::MOLs  GAGTCTAGTCGAACAG        AGTCTTCTATGC    ;%;chr10_127200298_127200839_-;%;chr10_127200879_127201038_-;%;chr10_127201186_127205674_-;%;chr10_127205819_127206613_-;%;chr10_127206711_127208869_-  NoTSS   NoPolyA ;%;chr10_127200217_127200297_-;%;chr10_127200840_127200878_-;%;chr10_127201039_127201185_-;%;chr10_127205675_127205818_-;%;chr10_127206614_127206710_-;%;chr10_127208870_127208939_-    known   5
# baeb5f94-80dd-4155-9c15-1469be720a18    ENSMUSG00000039105.6    P56::M1::Cerebellum::Glia::Oligo::MOLs  AAACGCTTCTTACCAT        CACGATCCAAAA    ;%;chr4_63545000_63548584_+;%;chr4_63548686_63549944_+  chr4_63544800_63544853_+        chr4_63550702_63550709_+        ;%;chr4_63544830_63544999_+;%;chr4_63548585_63548685_+;%;chr4_63549945_63550701_+       known   2

corrBed=$1;
allInfo=$2;
allInfoCorrected=$3;

cat $corrBed | awk 'NR>1 {exonChain=""; rn=$4; chr=$1; str=$6; nE=$10; nz=split($11,bsz,","); nt=split($12,bst,","); \
if(nE>=2) {for(i=1;i<=nE;i++) {exon=chr"_"1+$2+bst[i]"_"$2+bst[i]+bsz[i]"_"str; \
exonChain=exonChain";%;"exon} print rn"\t"exonChain} }' > corrBed_ai

nCol=$(zcat $allInfo | head -1 | awk '{print NF}')

if [ $nCol == 9 ]; then
echo "processing all-info file w/ complete reads"
awk -v allInfo=$allInfo 'BEGIN{comm="cat corrBed_ai"; while(comm|getline) {ec[$1]=$2;} comm="zcat "allInfo; \
while(comm|getline) {if($1 in ec) {$11=$9; $10=$8; $9=$7; $8=$6; $7=ec[$1]; split($9,o,/_|;%;/); \
split($7,n,/_|;%;/); if(n[5]!=o[5]) {gsub(n[5],o[5],$7);} {if($7==$9) {$6=$8;} else {intronChain=""; \
nE=split($7,exChain,";%;"); for(i=2;i<=(nE-1);i++){split(exChain[i],start,"_"); split(exChain[i+1],end,"_"); \
intron=start[1]"_"start[3]+1"_"end[2]-1"_"end[4]; intronChain=intronChain";%;"intron; $6=intronChain}} } OFS="\t"; \
print }}}' | gzip -c > $allInfoCorrected
else
echo "processing all-info file w/ complete and incomplete reads"
awk -v allInfo=$allInfo 'BEGIN{comm="cat corrBed_ai"; while(comm|getline) {ec[$1]=$2;} comm="zcat "allInfo; \
while(comm|getline) {if($1 in ec) {$13=$11; $12=$10; $11=$9; $10=$6; $9=ec[$1]; split($11,o,/_|;%;/); \
split($9,n,/_|;%;/); if(n[5]!=o[5]) {gsub(n[5],o[5],$9);} {if($9==$11) {$6=$10;} else {intronChain=""; \
nE=split($9,exChain,";%;"); for(i=2;i<=(nE-1);i++){split(exChain[i],start,"_"); split(exChain[i+1],end,"_"); \
intron=start[1]"_"start[3]+1"_"end[2]-1"_"end[4]; intronChain=intronChain";%;"intron; $6=intronChain}} } OFS="\t"; \
print }}}' | gzip -c > $allInfoCorrected

fi

rm corrBed_ai

echo "Done!"
