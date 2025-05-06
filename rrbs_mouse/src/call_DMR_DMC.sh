#!/bin/bash
# yhbatch -p e9 -n 1

colnum=$1
name=$2
header_dir=header
in_dir=infile
out_dir=outfile
mkdir -p $header_dir $in_dir $out_dir
tools=~/BIGDATA2/gzfezx_shhli_1/softwares/metilene_v0.2-8


echo "###step1:according to your pheno.csv,generate the grouped matrix file"
dos2unix pheno.csv
cat pheno.csv |grep -v '^RANK'|awk -v col="$colnum" -F',' '{print $col}'|sed '1i chrom\npos'|\
tr '\n' '\t'|sed G -  > $header_dir/$name.txt

sed '1d' InGroup2.txt |\
cat $header_dir/$name.txt - |sed '/^\s*$/d' > $in_dir/InGroup2.$name.txt

##############################   for DMR analysis  ############################################
echo "###step2:use metilene to test the difference between two groups, results are generate in outfile directory"

metilene -a g1 -b g2 -t 24 $in_dir/InGroup2.$name.txt | sort -k1,1V -k2,2n  > $out_dir/InGroup2.$name.DMRs.bed

perl $tools/metilene_output.pl -q $out_dir/InGroup2.$name.DMRs.bed \
-o $out_dir/InGroup2.$name.DMRs -p 0.05 -c 5 -d 0.1

/HOME/gzfezx_shhli/gzfezx_shhlixy_2/BIGDATA2/gzfezx_shhli_2/miniconda3/bin/python calbedmeth_v5.py \
-m $in_dir/InGroup2.$name.txt -b $out_dir/InGroup2.$name.DMRs_qval.0.05.out -o $out_dir/InGroup2.$name.DMRs_qval.0.05.sampleMet.bed

##############################   for DMC analysis  ############################################

echo "###step3:use metilene to test the difference between two groups, results are generate in outfile directory for DMC analysis"

metilene -a g1 -b g2 -f 3 -t 24 $in_dir/InGroup2.$name.txt | sort -k1,1V -k2,2n  > $out_dir/InGroup2.$name.DMCs.bed

less $out_dir/InGroup2.$name.DMCs.bed |awk '($5>0.1 || $5<-0.1) && $7<0.05'|awk '{print $1,$2,$7,$5,$6,$9,$10}' OFS='\t' > $out_dir/InGroup2.$name.DMCs_pval.0.05.out


