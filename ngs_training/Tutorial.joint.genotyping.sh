#!/bin/bash
#Tutorial.joint.genotyping.sh


set -exuo pipefail

SCRIPT_DIR=$(cd $(dirname $0)  && pwd)

#計算機環境に応じてスレッド数を指定
no_threads=3

######################################
## 1. 公開データ取得

#ターミナルを立ち上げた場所の直下にScerというフォルダを作成し、そこで作業をおこないます。
main_folder=$SCRIPT_DIR/Scer
mkdir -p $main_folder
#Scerの直下にfastqというフォルダを作成し、生リードの保存場所とする
fastq_folder=$main_folder/fastq
mkdir -p $fastq_folder
cd $fastq_folder

### sake001:日本酒酵母のリシーケンスの生リードデータSRR5678551の取得
if [ ! -e $fastq_folder/SRR5678551_1.fastq ]; then
	fastq-dump --split-files SRR5678551
fi

if [ ! -e $fastq_folder/sake001_2M_1.fastq.gz ]; then
	head -n 8000000 SRR5678551_1.fastq | gzip > sake001_2M_1.fastq.gz
	head -n 8000000 SRR5678551_2.fastq | gzip > sake001_2M_2.fastq.gz
fi

### sake002:日本酒酵母のリシーケンスの生リードデータSRR5678548の取得
if [ ! -e $fastq_folder/SRR5678548_1.fastq ]; then
	fastq-dump --split-files SRR5678548
fi

if [ ! -e $fastq_folder/sake002_2M_1.fastq.gz ]; then
	head -n 8000000 SRR5678548_1.fastq | gzip > sake002_2M_1.fastq.gz
	head -n 8000000 SRR5678548_2.fastq | gzip > sake002_2M_2.fastq.gz
fi

### sake003:日本酒酵母のリシーケンスの生リードデータSRR5678549の取得
if [ ! -e $fastq_folder/SRR5678549_1.fastq ]; then
	fastq-dump --split-files SRR5678549
fi

if [ ! -e $fastq_folder/sake003_2M_1.fastq.gz ]; then
	head -n 8000000 SRR5678549_1.fastq | gzip > sake003_2M_1.fastq.gz
	head -n 8000000 SRR5678549_2.fastq | gzip > sake003_2M_2.fastq.gz
fi


cd $main_folder


reference_folder=$main_folder/reference
mkdir -p $reference_folder
cd $reference_folder

# 酵母のリファレンスゲノムを取得していなければ、wgetコマンドでDL
if [ ! -e $reference_folder/sacCer3.fa ]; then
	wget -O - 'ftp://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/chromosomes/chr*' | gunzip -c > sacCer3.fa
fi

cd $main_folder


######################################
## 2. リードのクオリティーコントロール（QC)

cd $fastq_folder

### fastqcレポートの作成
if [ ! -e $fastq_folder/sake001_2M_1_fastqc.html ]; then
	fastqc sake001_2M_1.fastq.gz sake001_2M_2.fastq.gz
fi

if [ ! -e $fastq_folder/sake002_2M_1_fastqc.html ]; then
	fastqc sake002_2M_1.fastq.gz sake002_2M_2.fastq.gz
fi

if [ ! -e $fastq_folder/sake003_2M_1_fastqc.html ]; then
	fastqc sake003_2M_1.fastq.gz sake003_2M_2.fastq.gz
fi


### fastpを使ってQCを実施
if [ ! -e $fastq_folder/sake001.fastp.report.html ]; then
	fastp -i sake001_2M_1.fastq.gz -I sake001_2M_2.fastq.gz -o sake001_2M_1.trimmed.fastq.gz -O sake001_2M_2.trimmed.fastq.gz -f 5 -F 5 -q 30 -l 30 -w $no_threads -h sake001.fastp.report.html
fi

if [ ! -e $fastq_folder/sake002.fastp.report.html ]; then
	fastp -i sake002_2M_1.fastq.gz -I sake002_2M_2.fastq.gz -o sake002_2M_1.trimmed.fastq.gz -O sake002_2M_2.trimmed.fastq.gz -f 5 -F 5 -q 30 -l 30 -w $no_threads -h sake002.fastp.report.html
fi

if [ ! -e $fastq_folder/sake003.fastp.report.html ]; then
	fastp -i sake003_2M_1.fastq.gz -I sake003_2M_2.fastq.gz -o sake003_2M_1.trimmed.fastq.gz -O sake003_2M_2.trimmed.fastq.gz -f 5 -F 5 -q 30 -l 30 -w $no_threads -h sake003.fastp.report.html
fi


######################################
## 3. マッピング

### リファレンスのインデックス作成
cd $reference_folder
if [ ! -e $reference_folder/sacCer3.fa.bwt ]; then
	bwa index sacCer3.fa
fi

### マッピングの出力用フォルダを作成
bwa_out_folder=$main_folder/bwa_out
mkdir -p $bwa_out_folder
cd $bwa_out_folder

### bwaでマッピングを実施
if [ ! -e $bwa_out_folder/sake001.sorted.bam ]; then
	bwa mem -t $no_threads -R "@RG\tID:sacCer\tSM:sake001\tPL:Illumina" $reference_folder/sacCer3.fa $fastq_folder/sake001_2M_1.trimmed.fastq.gz $fastq_folder/sake001_2M_2.trimmed.fastq.gz | samtools view -@ $no_threads -Sb | samtools sort -@ $no_threads > sake001.sorted.bam
fi

if [ ! -e $bwa_out_folder/sake002.sorted.bam ]; then
	bwa mem -t $no_threads -R "@RG\tID:sacCer\tSM:sake002\tPL:Illumina" $reference_folder/sacCer3.fa $fastq_folder/sake002_2M_1.trimmed.fastq.gz $fastq_folder/sake002_2M_2.trimmed.fastq.gz | samtools view -@ $no_threads -Sb | samtools sort -@ $no_threads > sake002.sorted.bam
fi

if [ ! -e $bwa_out_folder/sake003.sorted.bam ]; then
	bwa mem -t $no_threads -R "@RG\tID:sacCer\tSM:sake003\tPL:Illumina" $reference_folder/sacCer3.fa $fastq_folder/sake003_2M_1.trimmed.fastq.gz $fastq_folder/sake003_2M_2.trimmed.fastq.gz | samtools view -@ $no_threads -Sb | samtools sort -@ $no_threads > sake003.sorted.bam
fi

######################################
## 4. 変異検出

### Checking for reference index (*.fai)
if [ ! -e $reference_folder/sacCer3.fa.fai ]; then
	samtools faidx $reference_folder/sacCer3.fa
fi

#Checking for gatk reference index (*.dict)
if [ ! -e $reference_folder/sacCer3.dict ]; then
	gatk CreateSequenceDictionary -R $reference_folder/sacCer3.fa -O $reference_folder/sacCer3.dict
fi

vcf_out_folder=$main_folder/vcf_out
mkdir -p $vcf_out_folder


# 前処理
cd $bwa_out_folder
# PCR duplicatesをマーク
if [ ! -e $bwa_out_folder/sake001.markdup.bam ]; then
	gatk MarkDuplicates -I $bwa_out_folder/sake001.sorted.bam -M $bwa_out_folder/sake001.metrics.txt -O $bwa_out_folder/sake001.markdup.bam
fi

if [ ! -e $bwa_out_folder/sake002.markdup.bam ]; then
	gatk MarkDuplicates -I $bwa_out_folder/sake002.sorted.bam -M $bwa_out_folder/sake002.metrics.txt -O $bwa_out_folder/sake002.markdup.bam

fi

if [ ! -e $bwa_out_folder/sake003.markdup.bam ]; then
	gatk MarkDuplicates -I $bwa_out_folder/sake003.sorted.bam -M $bwa_out_folder/sake003.metrics.txt -O $bwa_out_folder/sake003.markdup.bam
fi

# Multiple mapped readsを除去
if [ ! -e $bwa_out_folder/sake001.filtered.bam ]; then
	samtools view -@ no_threads -b -q 4 $bwa_out_folder/sake001.markdup.bam > $bwa_out_folder/sake001.filtered.bam
	samtools index  -@ no_threads $bwa_out_folder/sake001.filtered.bam
fi

if [ ! -e $bwa_out_folder/sake002.filtered.bam ]; then
	samtools view -@ no_threads -b -q 4 $bwa_out_folder/sake002.markdup.bam > $bwa_out_folder/sake002.filtered.bam
	samtools index  -@ no_threads $bwa_out_folder/sake002.filtered.bam
fi

if [ ! -e $bwa_out_folder/sake003.filtered.bam ]; then
	samtools view -@ no_threads -b -q 4 $bwa_out_folder/sake003.markdup.bam > $bwa_out_folder/sake003.filtered.bam
	samtools index  -@ no_threads $bwa_out_folder/sake003.filtered.bam
fi

# バリアントコール
gatk HaplotypeCaller -R $reference_folder/sacCer3.fa -I $bwa_out_folder/sake001.filtered.bam -O $vcf_out_folder/sake001.g.vcf.gz --emit-ref-confidence GVCF --bam-output $bwa_out_folder/sake001.g.hpcall.bam
gatk HaplotypeCaller -R $reference_folder/sacCer3.fa -I $bwa_out_folder/sake002.filtered.bam -O $vcf_out_folder/sake002.g.vcf.gz --emit-ref-confidence GVCF --bam-output $bwa_out_folder/sake002.g.hpcall.bam
gatk HaplotypeCaller -R $reference_folder/sacCer3.fa -I $bwa_out_folder/sake003.filtered.bam -O $vcf_out_folder/sake003.g.vcf.gz --emit-ref-confidence GVCF --bam-output $bwa_out_folder/sake003.g.hpcall.bam

#中間ファイルをローカルデータベースにまとめる
DB_path=$main_folder/gDB
echo -e "chrI\nchrII\nchrIII\nchrIV\nchrIX\nchrV\nchrVI\nchrVII\nchrVIII\nchrX\nchrXI\nchrXII\nchrXIII\nchrXIV\nchrXVI\nchrM" > $main_folder/intervals.list
gatk GenomicsDBImport -R $reference_folder/sacCer3.fa -V $vcf_out_folder/sake001.g.vcf.gz  -V $vcf_out_folder/sake002.g.vcf.gz  -V $vcf_out_folder/sake003.g.vcf.gz -L $main_folder/intervals.list --genomicsdb-workspace-path $DB_path

#GenotypeGVCFsコマンドを用いて3サンプルをまとめてジェノタイピング
gatk GenotypeGVCFs -R $reference_folder/sacCer3.fa -V gendb://$DB_path -O $vcf_out_folder/sake.3samples.raw.vcf.gz

cd $vcf_out_folder

# SNPのみを切り分け
gatk SelectVariants -R $reference_folder/sacCer3.fa -V $vcf_out_folder/sake.3samples.raw.vcf.gz --select-type SNP -O $vcf_out_folder/sake.3samples.snp.vcf.gz

#INFO filedを-filter-expressionの閾値でマーク
gatk VariantFiltration -R $reference_folder/sacCer3.fa -V $vcf_out_folder/sake.3samples.snp.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O $vcf_out_folder/sake.3samples.snp.filter.vcf.gz

#上記のコマンドでPASSしたサイトのみを抽出
gzip -dc $vcf_out_folder/sake.3samples.snp.filter.vcf.gz | grep -E '^#|PASS' | bgzip > $vcf_out_folder/sake.3samples.snp.filterPASSED.vcf.gz

#vdf.gzファイルのインデックス付
tabix -f -p vcf $vcf_out_folder/sake.3samples.snp.filterPASSED.vcf.gz

#FORMAT filedをG-filterの閾値でマーク
gatk VariantFiltration -R $reference_folder/sacCer3.fa -V $vcf_out_folder/sake.3samples.snp.filterPASSED.vcf.gz -G-filter "GQ < 20" -G-filter-name "GQ20" -G-filter "DP < 10" -G-filter-name "DP10" -O $vcf_out_folder/sake.3samples.snp.DPfilterPASSED.vcf.gz

#FORMAT fieldでマークされたジェノタイプを無効化
gatk SelectVariants -R $reference_folder/sacCer3.fa -V $vcf_out_folder/sake.3samples.snp.DPfilterPASSED.vcf.gz --set-filtered-gt-to-nocall -O $vcf_out_folder/sake.3samples.snp.DPfilterNoCall.vcf.gz

#複数のサンプル間で欠損率1%未満のサイトのみを抽出
gatk SelectVariants -R $reference_folder/sacCer3.fa -V $vcf_out_folder/sake.3samples.snp.DPfilterNoCall.vcf.gz --set-filtered-gt-to-nocall --max-nocall-fraction 0.99 --exclude-filtered -O $vcf_out_folder/sake.3samples.snp.DPfilterNoCall.P99.vcf.gz

gzip -dc $vcf_out_folder/sake.3samples.snp.DPfilterNoCall.P99.vcf.gz | awk '!/^#/' | wc -l

cd $SCRIPT_DIR