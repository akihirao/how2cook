#!/bin/bash
#Tutorial.single.sample.genotyping.sh


set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0)  && pwd)

#計算機環境に応じてスレッド数を指定
no_threads=3

######################################
## 1. 公開データ取得

#アカウント(hogehoge)
user_name=hogehoge
#ターミナルを立ち上げた場所の直下にScerというフォルダを作成
main_folder=$SCRIPT_DIR/Scer
mkdir -p $main_folder
#Scerの直下にfastqというフォルダを作成し、生リードの保存場所とする
fastq_folder=$main_folder/fastq
mkdir -p $fastq_folder
cd $fastq_folder


# 日本酒酵母のリシーケンスの生リードデータERR038793の取得
fastq-dump --split-files SRR5678551

# 生リードデータの一部(2M個のペアリード）を使用。あわせて分かりやすい名前に変更
head -n 8000000 SRR5678551_1.fastq | gzip > sake001_2M_1.fastq.gz
head -n 8000000 SRR5678551_2.fastq | gzip > sake001_2M_2.fastq.gz

cd $main_folder


reference_folder=$main_folder/reference
mkdir -p $reference_folder
cd $reference_folder

# wgetコマンドで酵母のリファレンスゲノムを取得します。
wget -O - 'ftp://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/chromosomes/chr*' | gunzip -c > sacCer3.fa

cd $main_folder


######################################
## 2. リードのクオリティーコントロール（QC)

cd $fastq_folder

#fastqcレポートの作成
fastqc sake001_2M_1.fastq.gz sake001_2M_2.fastq.gz

#fastpを使ってQCを実施
fastp -i sake001_2M_1.fastq.gz -I sake001_2M_2.fastq.gz -o sake001_2M_1.trimmed.fastq.gz -O sake001_2M_2.trimmed.fastq.gz -f 5 -F 5 -q 30 -l 30 -w $no_threads -h sake001.fastp.report.html

######################################
## 3. マッピング

# リファレンスのインデックス作成
cd $reference_folder
bwa index sacCer3.fa

# マッピングの出力用フォルダを作成
bwa_out_folder=$main_folder/bwa_out
mkdir -p $bwa_out_folder
cd $bwa_out_folder

# bwaでマッピングを実施
bwa mem -t $no_threads -R "@RG\tID:sacCer\tSM:sake001\tPL:Illumina" $reference_folder/sacCer3.fa $fastq_folder/sake001_2M_1.trimmed.fastq.gz $fastq_folder/sake001_2M_2.trimmed.fastq.gz | samtools view -@ $no_threads -Sb | samtools sort -@ $no_threads > sake001.sorted.bam

######################################
## 4. 変異検出

samtools faidx $reference_folder/sacCer3.fa

#Checking for gatk reference index (*.dict)
if [ ! -e $reference_folder/sacCer3.dict ]; then
	gatk CreateSequenceDictionary -R $reference_folder/sacCer3.fa -O $reference_folder/sacCer3.dict
fi

vcf_out_folder=$main_folder/vcf_out
mkdir -p $vcf_out_folder


# 前処理
cd $bwa_out_folder
# PCR duplicatesをマーク
gatk MarkDuplicates -I $bwa_out_folder/sake001.sorted.bam -M $bwa_out_folder/sake001.metrics.txt -O $bwa_out_folder/sake001.markdup.bam

# Multiple mapped readsを除去
samtools view -@ no_threads -b -q 4 $bwa_out_folder/sake001.markdup.bam > $bwa_out_folder/sake001.filtered.bam
samtools index  -@ no_threads $bwa_out_folder/sake001.filtered.bam


# バリアントコール
gatk HaplotypeCaller -R $reference_folder/sacCer3.fa -I $bwa_out_folder/sake001.filtered.bam --bam-output $bwa_out_folder/sake001.hpcall.bam -O $vcf_out_folder/sake001.raw.vcf


# SNPのみを切り分け
gatk SelectVariants -R $reference_folder/sacCer3.fa -V $vcf_out_folder/sake001.raw.vcf --select-type SNP -O $vcf_out_folder/sake001.snp.vcf

#INFO filedを-filter-expressionの閾値でマーク
gatk VariantFiltration -R $reference_folder/sacCer3.fa -V $vcf_out_folder/sake001.snp.vcf -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O $vcf_out_folder/sake001.snp.filter.vcf

#上記のコマンドでPASSしたサイトのみを抽出
grep -E '^#|PASS' $vcf_out_folder/sake001.snp.filter.vcf  > $vcf_out_folder/sake001.snp.filterPASSED.vcf

#FORMAT filedを-G-filterの閾値でマーク
gatk VariantFiltration -R $reference_folder/sacCer3.fa -V $vcf_out_folder/sake001.snp.filterPASSED.vcf -G-filter "GQ < 20" -G-filter-name "GQ20" -G-filter "DP < 10" -G-filter-name "DP10" -O $vcf_out_folder/sake001.snp.DPfilterPASSED.vcf

#FORMAT fieldでマークされたジェノタイプを無効化
gatk SelectVariants -R $reference_folder/sacCer3.fa -V $vcf_out_folder/sake001.snp.DPfilterPASSED.vcf --set-filtered-gt-to-nocall -O $vcf_out_folder/sake001.snp.DPfilterNoCall.vcf

#複数のサンプルにおいて欠損率が高いサイトを除去
gatk SelectVariants -R $reference_folder/sacCer3.fa -V $vcf_out_folder/sake001.snp.DPfilterNoCall.vcf --set-filtered-gt-to-nocall --max-nocall-fraction 0.99 --exclude-filtered -O $vcf_out_folder/sake001.snp.DPfilterNoCall.P99.vcf

cd $SCRIPT_DIR

######################################
