#!/bin/bash
# yhbatch -p e9 -n 1
# source /HOME/gzfezx_shhli/gzfezx_shhlixy_2/BIGDATA2/gzfezx_shhli_2/miniconda3/etc/profile.d/conda.sh && echo conda.sh
source /HOME/gzfezx_shhli/gzfezx_shhlixy_2/BIGDATA2/gzfezx_shhli_2/miniconda3/bin/activate manta2 && echo source_activate

# columns=9
# name=ADHD_vs_Control
columns=
name=

if [ -z "$columns" ] || [ -z "$name" ]; then
    echo "Error: Variable 'columns' or 'name' is empty. Please check your input!"
    exit 1
fi

echo "columns: $columns"
echo "name: $name"

ln -s /XYFS01/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_1/BIGDATA2/gzfezx_shhli_1/WANGCR/RRBS-20220402/1.analysis/RRBS-20240827/InGroup2.txt .
ln -s /XYFS01/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_1/BIGDATA2/gzfezx_shhli_1/WANGCR/RRBS-20220402/1.analysis/RRBS-20240827/pheno.csv pheno.csv
ln -s /XYFS01/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_1/BIGDATA2/gzfezx_shhli_1/WANGCR/RRBS-20220402/1.analysis/RRBS-20240827/src/* .

files=("InGroup2.txt" "pheno.csv" "calbedmeth_v5.py" "call_DMR_DMC.sh")
missing_files=()

for file in "${files[@]}"; do
    if [ ! -f "$file" ]; then
        missing_files+=("$file")
    fi
done

if [ ${#missing_files[@]} -gt 0 ]; then
    echo "log: The following files are missing:"
    printf '%s\n' "${missing_files[@]}"
    exit 1
fi

echo "log: All the files are available. Proceed to execute the script..."

echo "log: call DMR DMC"

files=(
    "outfile/InGroup2.$name.DMCs_pval.0.05.out"
    "outfile/InGroup2.$name.DMRs_qval.0.05.out"
    "outfile/InGroup2.$name.DMRs_qval.0.05.sampleMet.bed"
)

all_files_valid=true

for file in "${files[@]}"; do
    if [ ! -s "$file" ]; then
        echo "The file is missing or empty: $file"
        all_files_valid=false
    fi
done

if $all_files_valid; then
    echo "Skip the command if all files exist and are not empty: call_DMR_DMC.sh ..."
else
    echo "Some files are missing or empty. Execute the command. call_DMR_DMC.sh ..."
    sh call_DMR_DMC.sh $columns $name
fi


Rscript plot1.annotate_DMR.R $name
Rscript plot2.gsea_DMR.R $name
Rscript plot3.heatmap_DMR.R $name
Rscript plot4.manhattan_DMC.R $name

