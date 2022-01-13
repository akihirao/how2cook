# NGSデータ解析チュートリアル


平尾　章

水産研究・教育機構 水産資源研究所
水産資源研究センター 漁業情報解析部　資源解析グループ


----

## 概要

初心者向けのチュートリアルとして、酵母 *Saccharomyces* *cerevisiae* のリシーケンス解析を題材に公開データの取得から変異検出までの解析手順の流れを学びます。 


## 本チュートリアルの流れ

1.　[公開データ取得](#公開データ取得)

2.　[クオリティーコントロール](#リードデータのクオリティーコントロール)

3.　マッピング

4.　変異検出　


## 参考情報
- Web
  - [NGSハンズオン2015: ゲノムReseq、変異検出](https://www.iu.a.u-tokyo.ac.jp/~kadota/bioinfo_ngs_sokushu_2015/20150804_amelieff_20150902.pdf): (株)アメリエフ 山口昌男氏による講義資料
  - [NGSデータから新たな知識を導出するためのデータ解析リテラシー](https://github.com/yuifu/ajacs68): 尾崎遼さんらの講義資料@AJACS68
  - [macでインフォマティクス](https://kazumaxneo.hatenablog.com): 上坂一馬さんによるNGSツールなどの紹介
  - [(Rで)塩基配列解析](http://www.iu.a.u-tokyo.ac.jp/~kadota/r_seq.html)
  - [統合TV（NGS解析だけでなくDBなども）](http://togotv.dbcls.jp)
  - [Linux標準教科書](http://www.lpi.or.jp/linuxtext/text.shtml)
- 書籍
  - [「入門者のLinux」(奈佐原顕郎著)](https://gendai.ismedia.jp/list/books/bluebacks/9784062579896):Linux初心者の方におすすめです 


## 使用NGSツールのリスト

* bedtools: a powerful toolset for genome arithmetic https://bedtools.readthedocs.io
* BWA: Burrow-Wheeler Aligner http://bio-bwa.sourceforge.net
* Bwa-mem2: the next version of the bwa-mem https://github.com/bwa-mem2/bwa-mem2) 
* fastp: an all-in-one preprocessing tool for fastq files (https://github.com/OpenGene/fastp
* fastqc: a quality control tool https://www.bioinformatics.babraham.ac.uk/projects/fastqc
* GATK: Genome Analysis Toolkit https://gatk.broadinstitute.org
* Plink: whole-genome association analysis tool https://www.cog-genomics.org/plink
* samtools: tools for manipulating NGS data https://github.com/samtools/samtools
* seqkit: an ultrafast toolkit for FASTA/Q file manipulation https://github.com/shenwei356/seqkit
* Trimmomatic: a flexible read trimming tool https://github.com/usadellab/Trimmomatic
* vcftools: a set of tools for working with VCF files https://github.com/vcftools/vcftools


## NGSツールのインストールと設定
この例では、ubuntuマシン(18.04)の/home/hogehoge/localにツール類を入れることとします。パスの設定も適宜忘れずに！


#### SRA-toolkit
Ubuntu 64 bit版を本家サイトからダウンロード https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit
```
$ cd /home/hogehoge/local
$ wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.3/sratoolkit.2.11.3-ubuntu64.tar.gz　
$ tar zxvf sratoolkit.2.11.3-ubuntu64.tar.gz
```


#### fastqc
```
$ cd /home/hogehoge/local
$ wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
$ unzip fastqc_v0.11.9.zip
```


#### fastp
```
$ cd /home/hogehoge/local
# download the latest build
$ wget http://opengene.org/fastp/fastp
$ chmod a+x ./fastp

# or download specified version, i.e. fastp v0.23.1
$wget http://opengene.org/fastp/fastp.0.23.1
$mv fastp.0.23.1 fastp
$chmod a+x ./fastp
```


#### BWA
```
$ cd /home/hogehoge/local
$ git clone https://github.com/lh3/bwa.git
$ cd bwa; make
```

---
---

## リシーケンス解析のチュートリアル
出芽酵母 *Saccharomyces* *cerevisiae* は真核生物として初めてゲノム解読されたモデル生物です。ゲノムサイズが小さい(12.1Mb)ため、塩基配列データもコンパクトで扱いやすく、高スペックの計算機でなくても解析することができます。そこで酵母のリシーケンスデータを題材にして、公開データの取得から変異検出までの解析処理の流れを学びます。

### 1. 公開データ取得
#### 1-1. シーケンスリードの取得
酵母をリシーケンスした生リードデータ(ERR038793)をDRA/SRA/ERA公共データベースからSRA-toolkitを使ってダウンロードしてみましょう。

たとえば、こんな感じで作業フォルダおよびリードデータの保管フォルダを作っておきます。
```
$ user_name=hogehoge #アカウント名:hogehoge
$ main_folder=/home/$user_name/work/Scer
$ fastq_folder=$main_folder/fastq
$ mkdir -p $fastq_folder
$ cd $fastq_folder
```
SRA-toolkitのfastq-dumpコマンドを使って、リードデータを取得
```
$ fastq-dump --split-files ERR038793 #オプション--split-filesでペアエンドを２つのfastqに分割
```

リードデータの中身確認
```
$ head ERR038793_1.fastq　#fastqの先頭部分を閲覧
```
```
@ERR038793.1 1 length=100
GGACAAGGTTACTTCCTAGATGCTATATGTCCCTACGGCCTTGTCTAACACCATCCAGCATGCAATAAGGTGACATAGATATACCCACACACCACACCCT
+ERR038793.1 1 length=100
D/DDBD@B>DFFEEEEEEEEF@FDEEEBEDBBDDD:AEEE<>CB?FCFF@F?FBFF@?:EEE:EEBEEEB=EEE.>>?=AD=8CDFFFFFEFEF@C?;DC
@ERR038793.2 2 length=100
TGGTGGTATAAAGTGGTAGGGTAAGTATGTGTGTATTATTTACGATCATTTGTTAGCGTTTCAATATGGTGGGTAAAAACGCAGGATAGTGAGTTACCGA
...
```
```
$ seqkit stats ERR038793_1.fastq　#リードデータの概要チェック
```
```
file    format  type  num_seqs      sum_len   min_len   avg_len   max_len
ERR038793_1.fastq   FASTQ DNA 739,873 73,987,300  100 100 100
```
```
$ cd $main_folder　#メインの作業フォルダに戻る
```

#### 1-2. 酵母のリファレンスゲノムの取得

リファレンスゲノムの保存フォルダの準備
```
$ reference_folder=$main_folder/reference
$ mkdir -p $reference_folder
$ cd $reference_folder
```
酵母リファレンスゲノムを取得
```
$ wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
$ gzip -d GCF_000146045.2_R64_genomic.fna.gz　#gzを展開
$ mv GCF_000146045.2_R64_genomic.fna ScerCer3.fa #コンパクトな名前に変更
```
リファレンスゲノムの中身確認
```
$ head ERR038793_1.fastq　#fastqの先頭部分を閲覧
```
```
>NC_001133.9 Saccharomyces cerevisiae S288C chromosome I, complete sequence
ccacaccacacccacacacccacacaccacaccacacaccacaccacacccacacacacacatCCTAACACTACCCTAAC
ACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCAT
TCAACCATACCACTCCGAACCACCATCCATCCCTCTACTTACTACCACTCACCCACCGTTACCCTCCAATTACCCATATC
CAACCCACTGCCACTTACCCTACCATTACCCTACCATCCACCATGACCTACTCACCATACTGTTCTTCTACCCACCATAT
...
```
```
$ seqkit stats ScerCer3.fa　#リファレンスゲノムの概要チェック
```
```
file         format  type  num_seqs     sum_len  min_len    avg_len    max_len
ScerCer3.fa  FASTA   DNA         17  12,157,105   85,779  715,123.8  1,531,933
```
```
$ cd $main_folder
```

### 2. リードデータのクオリティーコントロール
#### 2-1. リードのクオリティーチェック
NGSから出力されるリードには cutadapt アダプター配列やポリA、ポリT、低クオリティのリードが含まれている場合があります。リードのデータにそのような配列が含まれていたり、その他おかしなことがないかを確認し、必要に応じてそういった配列をFASTQファイルからアダプターを取り除く必要があります。このような操作をリードのQCと呼び、特に後者ははリードトリミングやリードフィルタリングとも呼ばれます。

FASTQファイルのクオリティを確認する代表的ツールがFastQCです。まずバージョンを確認してみましょう。
```
$　fastqc --version
FastQC v.0.11.9
```
ついでヘルプで使い方をみてみましょう。
```
$ fastqc --help
  FastQC - A high throughput sequence QC analysis too
SYNOPSIS
fastqc seqfile1 seqfile2 .. seqfileN
...
```
FastQCを実行すると、QCの結果がHTML形式でレポートが出力されます。
```
$ fastqc ERR038793_1.fastq　
```
[上記のFastQC解析のレポート例](https://github.com/akihirao/how2cook/tree/main/ngs_training/ERR038793_1_fastqc.html)

FastQCのインストール、使い方、レポートの見方 https://bi.biopapyrus.jp/rnaseq/qc/fastqc.html


#### 2-2. リードのクオリティーフィルタリング
次に低品質のリードや塩基を除去します。



### 3. マッピング

まずリファレンスのインデックスを作成します。
```
$ cd $reference_folder
$ bwa index ScerCer3.fa
```



### 4. 変異検出
