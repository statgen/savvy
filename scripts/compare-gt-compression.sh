#!/bin/bash
# Usage: ./compare-gt-compression.sh <input.bcf> <output_directory>
#
# This script documents the commands used to compare SAV to other variant call formats.

set -eu

bcf_file=$1
out_dir=$2

if [[ -z "$bcf_file" ]]; then
  (>&2 echo "Error: missing input BCF")
  exit 1
fi

if [[ -z "$out_dir" ]]; then
  out_dir="."
fi

vcf_file=${out_dir}/vcf/$(basename $bcf_file .bcf).vcf.gz
bgt_file=${out_dir}/bgt/$(basename $bcf_file .bcf).bgt
gds_file=${out_dir}/gds/$(basename $bcf_file .bcf).gds
pgen_file=${out_dir}/pgen/$(basename $bcf_file .bcf).pgen
sav_file=${out_dir}/pgen/$(basename $bcf_file .bcf).sav
sav_pbwt_file=${out_dir}/pgen/$(basename $bcf_file .bcf).pbwt.sav


mkdir -p ${out_dir}/vcf/
bcftools view $input -Oz -o $vcf_file
bcftools index $vcf_file

mkdir -p ${out_dir}/bgt/
bgt import -S $bgt_file $vcf_file

mkdir -p ${out_dir}/gds/
Rscript -e "library(SeqArray)" -e "seqVCF2GDS(\"${vcf_file}\", \"${gds_file}\")"

wdir=`pwd`
mkdir -p ${out_dir}/gqt/
ln -s $bcf_file ${out_dir}/gqt/$(basename $bcf_file)
cd ${out_dir}/gqt/
gqt convert bcf -i $(basename $bcf_file) # Notice: multiple processes of GQT cannot be run simultaneously from the same directory.
rm $(basename $bcf_file)
cd $wdir

mkdir -p ${out_dir}/pgen/
plink2 --threads 1 --make-pgen vzs --vcf $vcf_file --out $pgen_file

mkdir -p ${out_dir}/bcf/
sav import --phasing full -b 8192 -19 $bcf_file $sav_file

mkdir -p ${out_dir}/bcf/
sav import --phasing full -b 8192 -19 --pbwt-fields GT --sparse-threshold 0.01 $bcf_file $sav_pbwt_file
