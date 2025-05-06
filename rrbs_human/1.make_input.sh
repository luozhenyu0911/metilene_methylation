#!/usr/bin/env bash
# yhbatch -N 1 -n 64 -p mars
# yhbatch -p e9 -n 1
data_bedf=~/BIGDATA2/gzfezx_shhli_1/project/RRBS-20240821/00.data/AG20240326R1-287hsa-CSDL/bed
data_bedf1=~/BIGDATA2/gzfezx_shhli_1/WANGCR/RRBS-20220402/data1/05.bed
data_bedf2=~/BIGDATA2/gzfezx_shhli_1/WANGCR/RRBS-20220402/data2/05.bed

for i in `ls $data_bedf1/ND_case* $data_bedf2/ND_case_* $data_bedf/B_*`
do
id=$(basename $i .gz)
id2=$(basename $id .CG_5x.bedgraph)
echo $id2
done > sampleID.txt

for i in `ls $data_bedf1/Control* $data_bedf2/control_* $data_bedf/A_*`
do
id=$(basename $i .gz)
id2=$(basename $id .CG_5x.bedgraph)
echo $id2
done >> sampleID.txt

g1=`ls $data_bedf1/ND_case* $data_bedf2/ND_case_* $data_bedf/B_* | tr '\n' ','`
g2=`ls $data_bedf1/Control* $data_bedf2/control_* $data_bedf/A_* | tr '\n' ','`

echo $g1
echo $g2
perl ~/WORK2/gzfezx_shhli_1/softwares/metilene_v0.2-8/metilene_input.pl --in1 $g1 --in2 $g2 --out InGroup1.txt -b \
~/WORK2/gzfezx_shhli_1/softwares/bedtools2/bin/bedtools &
perl ~/WORK2/gzfezx_shhli_1/softwares/metilene_v0.2-8/metilene_input.pl --in1 $g1 --in2 $g2 --out InGroup2.txt -b \
~/WORK2/gzfezx_shhli_1/softwares/bedtools2/bin/bedtools --NA 0 &
wait
