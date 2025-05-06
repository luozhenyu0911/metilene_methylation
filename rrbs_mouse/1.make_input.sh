#!/usr/bin/env bash
# yhbatch -N 1 -n 64 -p mars
# yhbatch -p e9 -n 1

for i in `ls data/D1_* data/D2_*`
do
id=$(basename $i .gz)
id2=$(basename $id .CG_5x.bedgraph)
echo $id2
done > sampleID.txt

for i in `ls data/C1_* data/C2_*`
do
id=$(basename $i .gz)
id2=$(basename $id .CG_5x.bedgraph)
echo $id2
done >> sampleID.txt

g1=`ls data/D1_* data/D2_* | tr '\n' ','`
g2=`ls data/C1_* data/C2_* | tr '\n' ','`

echo $g1
echo $g2
perl ~/WORK2/gzfezx_shhli_1/softwares/metilene_v0.2-8/metilene_input.pl --in1 $g1 --in2 $g2 --out InGroup1.txt -b \
~/WORK2/gzfezx_shhli_1/softwares/bedtools2/bin/bedtools &
perl ~/WORK2/gzfezx_shhli_1/softwares/metilene_v0.2-8/metilene_input.pl --in1 $g1 --in2 $g2 --out InGroup2.txt -b \
~/WORK2/gzfezx_shhli_1/softwares/bedtools2/bin/bedtools --NA 0 &
wait
