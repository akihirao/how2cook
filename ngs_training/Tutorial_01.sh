#!/bin/bash
#Tutorial_01.sh


set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0)  && pwd)


######################################
## 1. 公開データ取得

#user_name=hogehoge #アカウント名:hogehogeの場合
main_folder=$SCRIPT_DIR/Scer
mkdir -p $main_folder
fastq_folder=$main_folder/fastq
mkdir -p $fastq_folder
cd $fastq_folder

# 酵母のリシーケンスの生リードデータERR038793の取得
#fastq-dump --split-files ERR038793

# 日本酒酵母のリシーケンスの生リードデータERR038793の取得
#fastq-dump --split-files SRR13254428

cd $main_folder


reference_folder=$main_folder/reference
mkdir -p $reference_folder
cd $reference_folder

# wgetコマンドで酵母のリファレンスゲノムを取得します。
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz

# Macではwgetコマンドの代わりにcurlコマンドを使って下さい。
#curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz

#gzを展開
gzip -d GCF_000146045.2_R64_genomic.fna.gz

#コンパクトな名前に変更
mv GCF_000146045.2_R64_genomic.fna ScerCer3.fa

cd $main_folder


######################################
## 2. リードのクオリティーコントロール（QC)

#fastpを使ってQCを実施
no_threads=3 #計算機環境に応じてスレッド数を指定
cd $fastq_folder
fastp -i ERR038793_1.fastq -I ERR038793_2.fastq -o ERR038793_1.trimmed.fastq -O ERR038793_2.trimmed.fastq -q 30 -l 40 -w $no_threads -h ERR038793.fastp.report.html


######################################
## 3. マッピング

# リファレンスのインデックス作成
cd $reference_folder
bwa index ScerCer3.fa

# マッピングの出力用フォルダを作成
bwa_out_folder=$main_folder/bwa_out
mkdir -p $bwa_out_folder
cd $bwa_out_folder

# bwaでマッピングを実施
bwa mem -t $no_threads -R "@RG\tID:ERR038793\tSM:ERR038793\tPL:Illumina" $reference_folder/ScerCer3.fa $fastq_folder/ERR038793_1.trimmed.fastq $fastq_folder/ERR038793_2.trimmed.fastq | samtools view -@ $no_threads -Sb | samtools sort -@ $no_threads > ERR038793.sorted.bam


######################################
## 4. 変異検出

samtools faidx $reference_folder/ScerCer3.fa
gatk CreateSequenceDictionary -R $reference_folder/ScerCer3.fa -O $reference_folder/ScerCer3.dict

vcf_out_folder=$main_folder/vcf_out
mkdir -p $vcf_out_folder


# 前処理
cd $bwa_out_folder
gatk MarkDuplicates -I $bwa_out_folder/ERR038793.sorted.bam -M $bwa_out_folder/ERR038793.metrics.txt -O $bwa_out_folder/ERR038793.markdup.bam

samtools view -@ no_threads -b -q 4 $bwa_out_folder/ERR038793.markdup.bam > $bwa_out_folder/ERR038793.filtered.bam
samtools index  -@ no_threads $bwa_out_folder/ERR038793.filtered.bam

# バリアントコール
gatk HaplotypeCaller -R $reference_folder/ScerCer3.fa -I $bwa_out_folder/ERR038793.filtered.bam --bam-output $bwa_out_folder/ERR038793.hpcall.bam -O $vcf_out_folder/ERR038793.raw.vcf

# SNPのみを切り分け
gatk SelectVariants -R $reference_folder/ScerCer3.fa -V $vcf_out_folder/ERR038793.raw.vcf --select-type SNP -O $vcf_out_folder/ERR038793.snp.vcf

# VcfファイルのINFO fieldを対象としたサイトベースのフィルタリング
gatk VariantFiltration -R $reference_folder/ScerCer3.fa -V $vcf_out_folder/ERR038793.snp.vcf --filter-expression "MQ < 40.0" --filter-name "MQ40" --filter-expression "QUAL < 30.0" --filter-name "QUAL30" -O $vcf_out_folder/ERR038793.snp.filtered.vcf

# VcfファイルのFORMAT fieldを対象としたサンプルベースのフィルタリング
gatk VariantFiltration -R $reference_folder/ScerCer3.fa -V $vcf_out_folder/ERR038793.snp.filtered.vcf -G-filter "GQ < 20" -G-filter-name "GQ20" -G-filter "DP < 10" -G-filter-name "DP10" -O $vcf_out_folder/ERR038793.snp.DPfiltered.vcf

cd $SCRIPT_DIR

######################################
