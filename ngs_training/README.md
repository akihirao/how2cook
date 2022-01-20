# NGSデータ解析チュートリアル


平尾　章

水産研究・教育機構 水産資源研究所
水産資源研究センター 漁業情報解析部　資源解析グループ


----

## 概要

初心者向けのNGSデータ解析のチュートリアルです。公開NGSデータの取得から変異検出までの解析手順の流れについて、酵母のリシーケンス解析 (illuminaショートリード) を例に学びます。



## 参考情報
- Web
  - [NGSハンズオン2015: ゲノムReseq、変異検出](https://www.iu.a.u-tokyo.ac.jp/~kadota/bioinfo_ngs_sokushu_2015/20150804_amelieff_20150902.pdf): (株)アメリエフ 山口昌男氏による講義資料
  - [NGSデータから新たな知識を導出するためのデータ解析リテラシー](https://github.com/yuifu/ajacs68): 尾崎遼さんらの講義資料@AJACS68
  - [macでインフォマティクス](https://kazumaxneo.hatenablog.com): 上坂一馬さんによるNGSツールなどの紹介
  - [(Rで)塩基配列解析](http://www.iu.a.u-tokyo.ac.jp/~kadota/r_seq.html): 門田先生らによる充実サイト
  - [統合TV（NGS解析だけでなくDBなども）](http://togotv.dbcls.jp)
  - [Linux標準教科書](http://www.lpi.or.jp/linuxtext/text.shtml)


- 書籍
  - [「入門者のLinux」(奈佐原顕郎著)](https://gendai.ismedia.jp/list/books/bluebacks/9784062579896):Linux初心者の方におすすめです


- [よく使うシェルコマンド](#よく使うシェルコマンド)

---

# リシーケンスのドライ解析のチュートリアル
出芽酵母（ *Saccharomyces* *cerevisiae* ）は、真核生物として最初にゲノムが解読されたモデル生物です。ゲノムサイズ（12.1Mb）が小さいため、塩基配列データもコンパクトで扱いやすく、高スペックの計算機でなくてもデータ解析をおこなうことができます。酵母のリシーケンスを例題として、公開データの取得から変異検出までの解析処理の流れを学びます。なお解析環境として、Ubuntuマシンに[使用NGSツールのリスト](#使用NGSツールのリスト)がインストール済みであることを想定しています。

## 本チュートリアルの流れ

1. [公開データ取得](#公開データ取得)

2. [リードのクオリティーコントロール（QC）](#リードのクオリティーコントロール（QC）)

3. [マッピング](#マッピング)

4. [変異検出](#変異検出)

5. [ゲノムビューワーによる変異の視覚化](#ゲノムビューワーによる変異の視覚化)

![](images/ngs_training_01.png)

<h2 id="公開データ取得">1.&nbsp;公開データ取得</h2>

### 1-1. シーケンスリードの取得

酵母のリシーケンスの生リードデータ [ERR038793](https://www.ncbi.nlm.nih.gov/sra/ERR038793) を公共データベースからダウンロードします。

[ERR038793](https://www.ncbi.nlm.nih.gov/sra/ERR038793) のリンクからメタ情報をみると、illumina Hiseq 2000 で解読したペアエンドリードのデータであることが分かります。

![](images/ngs_training_02.png)

ダウンロードする前に、こんな感じでメインの作業フォルダおよびリードデータの保管フォルダを作っておきます。
```
user_name=hogehoge #アカウント名:hogehogeの場合
main_folder=/home/$user_name/work/Scer
fastq_folder=$main_folder/fastq
mkdir -p $fastq_folder
cd $fastq_folder
```
ブラウザからは直接fastq形式のリードデータをダウンロードすることはできません。そこでSRA-toolkitのfastq-dumpコマンドを使って、生リードデータ [ERR038793](https://www.ncbi.nlm.nih.gov/sra/ERR038793) を[DRA/SRA/ERAデータベース](https://www.ddbj.nig.ac.jp/dra/index.html)からダウンロードします。NGSツールを使う際には、解析の再現性担保のために、そのツールのバージョンを確認することが大切です。
```
fastq-dump --version
```
```
fastq-dump : 2.8.0
```
次いでfastq-dumpの使い方をhelpで確認してみましょう。
```
fastq-dump --help
```
```
Usage:
  fastq-dump [options] <path> [<path>...]
  fastq-dump [options] <accession>
...
```
fastq-dumpコマンドにオプション--split-filesをつけて実行することで、ペアエンドのSRAデータ [ERR038793](https://www.ncbi.nlm.nih.gov/sra/ERR038793) は２つのfastqに分割して取得されます。
```
fastq-dump --split-files ERR038793
```
処理が済むと、RR038793_1.fastq とERR038793_2.fastqという２つのファイルが作成されます。ペアエンド１のRR038793_1.fastqを対象に冒頭部分を見てみましょう。
```
head ERR038793_1.fastq　#fastqの先頭部分を閲覧
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
Fastq形式について https://qiita.com/hayatak/items/0ab4f8bc3c051dd9a0d4

次いでfastq/fasta操作ツールであるseqkitを用いて、ペアエンド１の概要を確認してみましょう。
```
seqkit stats ERR038793_1.fastq　#seqkitはfasta/fastqの操作ツール
```
```
file    format  type  num_seqs      sum_len   min_len   avg_len   max_len
ERR038793_1.fastq   FASTQ DNA 739,873 73,987,300  100 100 100
```
ペアエンド１は、リード長が100bp、リード数は739,873個で、総計73,987,300bpです。

確認を終えたら、メインの作業フォルダに戻っておきましょう。

```
cd $main_folder　
```

#### 1-2. 酵母のリファレンスゲノムの取得

リファレンスゲノムの保存フォルダの準備
```
reference_folder=$main_folder/reference
mkdir -p $reference_folder
cd $reference_folder
```
wgetコマンドで酵母のリファレンスゲノムを取得します。
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
```
Macではwgetコマンドの代わりにcurlコマンドを使って下さい。
```
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
```
ダウンロードした圧縮ファイルを解答して、コンパクトな名前に付け替えます。
```
#gzを展開
gzip -d GCF_000146045.2_R64_genomic.fna.gz
#コンパクトな名前に変更
mv GCF_000146045.2_R64_genomic.fna ScerCer3.fa
```
リファレンスゲノムの中身確認
```
#fasta形式の塩基配列ファイルの先頭部分を閲覧
head ScerCer3.fa
```
```
>NC_001133.9 Saccharomyces cerevisiae S288C chromosome I, complete sequence
ccacaccacacccacacacccacacaccacaccacacaccacaccacacccacacacacacatCCTAACACTACCCTAAC
ACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCAT
..
```
リファレンスゲノムの概要チェック
```
seqkit stats ScerCer3.fa　
```
```
file         format  type  num_seqs     sum_len  min_len    avg_len    max_len
ScerCer3.fa  FASTA   DNA         17  12,157,105   85,779  715,123.8  1,531,933
```
```
cd $main_folder
```

<h2 id="リードのクオリティーコントロール（QC）">2.&nbsp;リードのクオリティーコントロール（QC）</h2>

NGSから出力されるリードにはアダプター配列や低品質のリードが含まれている場合があります。データ前処理として、リードデータの品質を確認し、ノイズとなりそうなリードやアダプター配列、塩基を取り除いておきます。前者をリードクオリティーチェック、後者をリードフィルタリング(またはリードトリミング）と呼び、これらの一連の処理をクオリティーコントロール（Quality control: QC)と呼びます。

#### 2-1. リードのクオリティーチェック

FASTQファイルのクオリティを確認する代表的ツールがFastQCです。まずFastQCのバージョンを確認してみましょう。
```
fastqc --version
```
```
FastQC v.0.11.9
```
またヘルプで使い方を確認してみましょう。
```
fastqc --help
```
```
  FastQC - A high throughput sequence QC analysis too
SYNOPSIS
fastqc seqfile1 seqfile2 .. seqfileN
...
```
FastQCを実行すると、QCの結果がHTML形式でレポート出力されます。
```
fastqc ERR038793_1.fastq ERR038793_2.fastq
```
[上記のFastQC解析のレポート例](https://github.com/akihirao/how2cook/tree/main/ngs_training/ERR038793_1_fastqc.html)

FastQCのインストール、使い方、レポートの見方について https://bi.biopapyrus.jp/rnaseq/qc/fastqc.html


#### 2-2. リードのクオリティーフィルタリング

次に低品質のリードや塩基を除去します。クオリティーのチェックからフィルタリングまでを一括して高速に処理してくれるQCツールfastpを使うことにします。ここでは平均でQ30未満のリードおよびリード長40bp未満を除去する設定を適用しますが、個々のデータに応じてフィルタリング設定を調整することが望ましいです。

fastpの使い方などについて　https://kazumaxneo.hatenablog.com/entry/2018/05/21/111947

```
no_threads=3 #計算機環境に応じてスレッド数を指定
cd $fastq_folder
fastp -i ERR038793_1.fastq -I ERR038793_2.fastq -o ERR038793_1.trimmed.fastq -O ERR038793_2.trimmed.fastq -q 30 -l 40 -w $no_threads -h ERR038793.fastp.report.html
```
```
Read1 before filtering:
total reads: 739873
total bases: 73987300
Q20 bases: 71492383(96.6279%)
Q30 bases: 66616427(90.0377%)
...
fastp v0.22.0, time used: 7 seconds
```
fastp処理が完了すると、フィルタリング前後の結果がHTML形式でレポートされます。レポートを確認して、以降の処理に進むべきか、それとも、フィルタリング設定を調整して再度QCをおこなうべきかを検討することが大事です。


<h2 id="マッピング">3.&nbsp;マッピング</h2>

マッピングとは、シーケンサーから出てきた大量の塩基配列リードついて、参照ゲノム配列の中の該当する箇所を見つける処理です。

![](images/ngs_training_4_01.png)

定番のマッピングツールであるbwaを使ってマッピングします。Bwaにはいくつかのアルゴリズがありますが、ここではbwa-memを使います。なおbwa-memアルゴリズムのネクストバーションが独立したツール [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) として公開されており、高速化に適しています（ただしメモリ使用量やインデックスサイズも大きくなるので、導入する際には留意して下さい）。

#### 3-1. リファレンスのインデックス作成

まずリファレンス (参照配列) に対してインデックス（索引）を作成します。参照配列への高速な検索を図るために事前に索引を作るといった作業になります。
```
cd $reference_folder
bwa index ScerCer3.fa
```
```
[bwa_index] Pack FASTA... 0.07 sec
[bwa_index] Construct BWT for the packed sequence...
...
[main] CMD: bwa index ScerCer3.fa
[main] Real time: 5.287 sec; CPU: 5.239 sec
```
処理が済むと、*.amb、*.ann、*.bwt、*.pac というファイルが作成され、bwaが使用するインデックスとなります。

#### 3-2. マッピング

マッピングの出力用フォルダを作成しておきます。
```
bwa_out_folder=$main_folder/bwa_out
mkdir -p $bwa_out_folder
cd $bwa_out_folder
```
bwa memコマンドの使い方を確認しましょう。
```
bwa mem

Usage: bwa mem [options] <idxbase> <in1.fq> [in2.fq]
...
```

bwa memコマンドにてマッピングを実行します
* gatk解析のためにオプション-Rでリードグループ(@RG)を指定しておきます。IDはその名のとおりID、SMはサンプル名、PLはシーケンスのプラットフォームに対応させます。gatkによるリードグループに関する解説は[こちら](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups)

```
bwa mem -t $no_threads -R "@RG\tID:ERR038793\tSM:ERR038793\tPL:Illumina" $reference_folder/ScerCer3.fa $fastq_folder/ERR038793_1.trimmed.fastq $fastq_folder/ERR038793_2.trimmed.fastq > ERR038793.sam
```
マッピングで出力されたsam形式ファイルの中身を確認
```
less ERR038793.sam
```
sam形式ファイルはテキストで記述されており、ヘッダーのメタ情報に続いて、各リードのアライメント（どこにマッピングされて、ミスマッチはいくつあったかなど）が１行毎に記載されています。
```
@SQ     SN:NC_001133.9  LN:230218
@SQ     SN:NC_001134.8  LN:813184
...
...
@SQ     SN:NC_001224.1  LN:85779
@RG     ID:ERR038793    SM:ERR038793    PL:Illumina
@PG     ID:bwa  PN:bwa  VN:0.7.17-r1198-dirty   CL:bwa mem -R @RG\tID:ERR038793\tSM:ERR038793\tPL:Illumina /home/akihirao/temp/Scer/reference/ScerCer3.fa /home/akihirao/temp/Scer/fastq/ERR038793_1.trimmed.fastq /home/akihirao/temp/Scer/fastq/ERR038793_2.trimmed.fastq
ERR038793.3     99      NC_001133.9     90      60      100M    =       386     395     ATCTAACCCTGGCCAACCTGTTTCTCAACTTACCCTCCATTACCCTACCTCCACTCGTTACCCTGTCCCATTCAACCATACCACTCCGAACCACCATCCA    HHHHHHHHHHHGHHHHGHHBHHHHHHBFHHFDGGHHGHHBHHGHHBEHEFHFFHHFBEE=FFEFB@EFFEBEFECFBFFGHFGHFHFHH?HEHEGHEHFF    NM:i:2  MD:Z:21C24G53   MC:Z:67M7D25M8S AS:i:90 XS:i:55 RG:Z:ERR038793
ERR038793.3     147     NC_001133.9     386     60      67M7D25M8S      =       90      -395    CCATGCCATAGTCACCCTCACTTGTATACTGATTTTACGTACGCACACGGATGCTACAGTATATACCACTTACCCTACTCTCATATTCCACTCCACACCC    ,?CCCEDECD?=8/?ED<BCGFFFGHFFHHGFGGGF=FEFEEEEADBDDCFCDEDCDDFFD@DC5CC;C@E;EGEFCF@GAFFDFECFFBD?GDAD?GDD    NM:i:10 MD:Z:0A9C56^ATCTCAA16G8 MC:Z:100M       AS:i:68 XS:i:53 RG:Z:ERR038793
...
```
sam形式の解説　https://bi.biopapyrus.jp/format/sam.html

samtoolsを使って、sam形式からバイナリータイプのbam形式に変換し、データサイズを圧縮しておきます。
```
samtools view -@ $no_threads -Sb ERR038793.sam > ERR038793.bam
```

次いでbamをソートします。
```
samtools sort -@ $no_threads ERR038793.bam > ERR038793.sorted.bam
```
なお中間生成ファイルによるストレージの圧迫を避けるならば、次のようにパイプを使って bwa mem から samtools sort までの処理を一度におこなってもよいでしょう。
```
#bwa mem -t $no_threads -R "@RG\tID:ERR038793\tSM:ERR038793\tPL:Illumina" $reference_folder/ScerCer3.fa $fastq_folder/ERR038793_1.trimmed.fastq $fastq_folder/ERR038793_2.trimmed.fastq | samtools view -@ $no_threads -Sb | samtools sort -@ $no_threads > ERR038793.sorted.bam
```
出来上がった bam ファイルにもインデックスをつけておきます。
```
samtools index ERR038793.sorted.bam

```
bamtoolsを使って、マッピングの概要を確認します。
```
bamtools stats -in ERR038793.sorted.bam
```
```
**********************************************
Stats for BAM file(s):
**********************************************

Total reads:       1258550
Mapped reads:      1213020	(96.3823%)
Forward strand:    651156	(51.7386%)
Reverse strand:    607394	(48.2614%)
Failed QC:         0	(0%)
Duplicates:        0	(0%)
Paired-end reads:  1258550	(100%)
'Proper-pairs':    1167167	(92.739%)
Both pairs mapped: 1210760	(96.2028%)
Read 1:            629162
Read 2:            629388
Singletons:        2260	(0.179572%)

```


<h2 id="変異検出">4.&nbsp;変異検出</h2>

GATK (Genome Analysis toolkit)を使用して、BAMファイルから変異を検出します。

まずリファレンス配列に対してGATK用のインデックスを作成しておきましょう。
```
samtools faidx $reference_folder/ScerCer3.fa
gatk CreateSequenceDictionary -R $reference_folder/ScerCer3.fa -O $reference_folder/ScerCer3.dict
```
vcf形式ファイルの出力用フォルダを作っておきます。
```
vcf_out_folder=$main_folder/vcf_out
mkdir -p $vcf_out_folder
```

#### 4-1. 前処理

シーケンスライブラリーの作成にPCRを使っている場合、マッピングされたリードの中に PCR duplication に由来する重複リードが含まれている可能性があります。このような重複リードはバリアントコールに偽陽性をもたらす可能性があるので、目印をつけておき（マーキング）、ダウンストリームで除去できるようにしておきます。なおゲノム縮約シーケンス（RAD-SeqやMig-Seq、GRAS-Di
など）では、このような重複リードの除去は不必要です。
```
cd $bwa_out_folder
gatk MarkDuplicates -I $bwa_out_folder/ERR038793.sorted.bam -M $bwa_out_folder/ERR038793.metrics.txt -O $bwa_out_folder/ERR038793.markdup.bam
```
さらに複数箇所にマッピングされたリードも偽陽性の原因となるので、除去しておきます。
```
samtools view -@ no_threads -b -q 4 $bwa_out_folder/ERR038793.markdup.bam > $bwa_out_folder/ERR038793.filtered.bam
samtools index  -@ no_threads $bwa_out_folder/ERR038793.filtered.bam
```

前処理済みのbamファイルの概要を表示します。
```
bamtools stats -in ERR038793.filtered.bam
```
```
**********************************************
Stats for BAM file(s):
**********************************************

**********************************************
Stats for BAM file(s):
**********************************************

Total reads:       1120866
Mapped reads:      1120866	(100%)
Forward strand:    560599	(50.0148%)
Reverse strand:    560267	(49.9852%)
Failed QC:         0	(0%)
Duplicates:        15838	(1.41301%)
Paired-end reads:  1120866	(100%)
'Proper-pairs':    1082195	(96.5499%)
Both pairs mapped: 1118908	(99.8253%)
Read 1:            560339
Read 2:            560527
Singletons:        1958	(0.174686%)

```
Duplicate リードが認識され、multiple mapping readsの除去によって、総リード数が減少したことが確認できます。

これらの前処理に加えて、ヒトやマウスなどの生物では既知変異データをもとに塩基スコアを再計算してBAMファイルのクオリティーを補正すること(Base Quality Score Recalibration: BQSR)の有効性が示されていますが、本チュートリアルでは省略します。

#### 4-2. バリアントコール

gatk HaplotypeCaller コマンドを使って、バリアントコールをおこないます。
* gatk ver4.0以降の HaplotypeCaller では、アクティブ領域 (各塩基のエントロピーの計算に基づいてバリアントの存在が予想される領域) を検出し、局所アッセンブルを適用することで、SNPs/INDELの検出精度が向上するという工夫が施されています。

```
gatk HaplotypeCaller -R $reference_folder/ScerCer3.fa -I $bwa_out_folder/ERR038793.filtered.bam --bam-output $bwa_out_folder/ERR038793.hpcall.bam -O $vcf_out_folder/ERR038793.raw.vcf
```
処理が済むとvcf形式のファイルが作られます。vcfファイルの中身を見てみましょう。
```
less $vcf_out_folder/ERR038793.raw.vcf
```
```
##fileformat=VCFv4.2
##FILTER=<ID=LowQual,Description="Low quality">
...
...
##contig=<ID=NC_001224.1,length=85779>
##source=HaplotypeCaller
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  ERR038793
NC_001133.9     111     .       C       T       328.64  .       AC=1;AF=0.500;AN=2;BaseQRankSum=-0.895;DP=24;ExcessHet=
3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=51.97;MQRankSum=0.167;QD=23.47;ReadPosRankSum=2.044;SOR=1.329     GT:AD:DP:GQ:PL
  0/1:5,9:14:99:336,0,158
NC_001133.9     136     .       G       A       616.06  .       AC=2;AF=1.00;AN=2;DP=28;ExcessHet=3.0103;FS=0.000;MLEAC
=2;MLEAF=1.00;MQ=51.67;QD=25.36;SOR=5.421        GT:AD:DP:GQ:PL  1/1:0,14:14:42:630,42,0
NC_001133.9     141     .       C       T       106.64  .       AC=1;AF=0.500;AN=2;BaseQRankSum=1.111;DP=27;ExcessHet=3
.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=53.61;MQRankSum=-0.488;QD=6.66;ReadPosRankSum=2.004;SOR=0.527      GT:AD:DP:GQ:PL
  0/1:12,4:16:99:114,0,413
NC_001133.9     244     .       C       CT      185.61  .       AC=1;AF=0.500;AN=2;BaseQRankSum=-0.319;DP=11;ExcessHet=
3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=55.46;MQRankSum=-0.319;QD=23.20;ReadPosRankSum=1.534;SOR=2.833    GT:AD:DP:GQ:PL
  0/1:1,7:8:19:193,0,19
...
```

vcf形式の解説　https://bi.biopapyrus.jp/gwas/vcf.html

この段階のvcfファイルにはSNPsとINDELsの変異情報が記載されています。変異の総数を確認しましょう。
```
awk '!/^#/' $vcf_out_folder/ERR038793.raw.vcf | wc -l
```
```
69151
```
計69151個の変異が検出されましたが、これらの中には偽陽性の可能性が高いものが含まれます。そこで次のステップでフィルタリングをおこないます。

#### 4-3. フィルタリング

上記のvcfファイルからSNPsの情報だけを取り出します。一般的にSNPsよりもINDELsの方が検出精度が低くなるため、それぞれを異なる閾値でフィルタリングすることが望ましいからです。
```
gatk SelectVariants -R $reference_folder/ScerCer3.fa -V $vcf_out_folder/ERR038793.raw.vcf --select-type SNP -O $vcf_out_folder/ERR038793.snp.vcf
```

gatk VariantFiltrationコマンドでフィルタリングします。
* [gatkによるフィルタリングのパラメーターと数値の意味の解説](https://gatk.broadinstitute.org/hc/en-us/articles/360035531112?id=6925)

まずvcfファイルのINFO fieldを対象としてサイトベースの
フィルタリング (オプション -filter-expression) をおこないます。。
ここでは比較的シンプルな設定を適用します。
* [gatkによる -filter/--filter-expression の解説](https://gatk.broadinstitute.org/hc/en-us/articles/360037434691-VariantFiltration#--filter-expression)

```
gatk VariantFiltration -R $reference_folder/ScerCer3.fa -V $vcf_out_folder/ERR038793.snp.vcf --filter-expression "MQ < 40.0" --filter-name "MQ40" --filter-expression "QUAL < 30.0" --filter-name "QUAL30" -O $vcf_out_folder/ERR038793.snp.filtered.vcf
```

次いでvcfファイルのFORMAT fieldを対象としてサンプルベースのフィルタリング (オプション -G-filter) をおこないます。こちらも比較的シンプルな設定です。
* [gatkによる -G-filter/--genotype-filter-expression の解説](https://gatk.broadinstitute.org/hc/en-us/articles/360037434691-VariantFiltration#--genotype-filter-expression)
```
gatk VariantFiltration -R $reference_folder/ScerCer3.fa -V $vcf_out_folder/ERR038793.snp.filtered.vcf -G-filter "GQ < 20" -G-filter-name "GQ20" -G-filter "DP < 10" -G-filter-name "DP10" -O $vcf_out_folder/ERR038793.snp.DPfiltered.vcf
```

フィルタリング後のSNPsの数を確認します。
```
awk '!/^#/' $vcf_out_folder/ERR038793.snp.DPfiltered.vcf | wc -l
```
```
63195
```
最終的に63195個のSNPsが検出されました。メインの作業フォルダに戻っておきましょう。
```
cd $main_folder
```
以上で、本チュートリアルにおけるコマンドライン操作は全て終了です。

<h2 id="ゲノムビューワーによる変異の視覚化">5.&nbsp;ゲノムビューワーによる変異の視覚化</h2>

[Integrated Genome Viewer (IGV)](https://software.broadinstitute.org/software/igv/)を使って、マッピングやバリアントコールの結果を視覚化してみましょう。

IGVの使い方について　https://bi.biopapyrus.jp/rnaseq/mapping/igv/

表示例

![](images/ngs_training_5_01.png)

<h2 id="6. その他">6.&nbsp;その他</h2>

#### 複数サンプルからの変異検出について

本チュートリアルでは1サンプルのみを対象として変異検出をおこなったので、single sample genotypingという方法を用いています。一方で、複数サンプルを対象とする場合は、すべてのサンプルをまとめてジェノタイピングする joint genotypingという方法がおすすめです。変異検出の精度が向上することが知られていますので、ぜひ適用を試みて下さい。
* [gatkによるsingle sample genotypingと joint genotypingの解説](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-)

そして複数サンプルを対象としてジェノタイピングをおこなうと、各々のサンプルによって遺伝子型が欠損するようなサイトがたびたび生じます。そのような欠損サイトの取扱や注意点いについて、ゲノム縮約シーケンスデータを対象とした岩崎貴也さんの次の講演資料が参考になります。

* [NGSのSNPデータを集団遺伝解析に使う事の利点と欠点：非モデル生物の研究で気をつけることは？](https://drive.google.com/file/d/1UK04C1IbHGvosibPWjjbKT1t1pcaTR9-/view)

---
[ページトップに戻る](https://github.com/akihirao/how2cook/tree/main/ngs_training#NGSデータ解析チュートリアル)

---
---

<h2 id="使用NGSツールのリスト">使用NGSツールのリスト</h2>

* bedtools: a powerful toolset for genome arithmetic https://bedtools.readthedocs.io
* BWA*: Burrow-Wheeler Aligner http://bio-bwa.sourceforge.net
* Bwa-mem2: the next version of the bwa-mem https://github.com/bwa-mem2/bwa-mem2)
* fastp*: an all-in-one preprocessing tool for fastq files (https://github.com/OpenGene/fastp
* fastqc*: a quality control tool https://www.bioinformatics.babraham.ac.uk/projects/fastqc
* GATK*: Genome Analysis Toolkit https://gatk.broadinstitute.org
* Plink: whole-genome association analysis tool https://www.cog-genomics.org/plink
* samtools*: tools for manipulating NGS data https://github.com/samtools/samtools
* seqkit*: an ultrafast toolkit for FASTA/Q file manipulation https://github.com/shenwei356/seqkit
* TASSEL: trait analysis by association, evolution and linkage https://www.maizegenetics.net/tassel
* Trimmomatic: a flexible read trimming tool https://github.com/usadellab/Trimmomatic
* vcftools: a set of tools for working with VCF files https://github.com/vcftools/vcftools


## NGSツールのインストールと設定
この例では /home/hogehoge/local にツール類を入れることとします。各自の環境に合わせたパスの設定を忘れずに！

以下に、いくつかの例を記しておきますが、本家サイトの最新情報を参照するように心掛けてください。

#### SRA-toolkit
Ubuntu 64 bit版を本家サイトからダウンロード https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit
```
cd /home/hogehoge/local
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.3/sratoolkit.2.11.3-ubuntu64.tar.gz　
tar zxvf sratoolkit.2.11.3-ubuntu64.tar.gz
```


#### fastqc
```
cd /home/hogehoge/local
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
unzip fastqc_v0.11.9.zip
```


#### fastp
```
cd /home/hogehoge/local
# download the latest build
wget http://opengene.org/fastp/fastp
chmod a+x ./fastp

# or download specified version, i.e. fastp v0.23.1
wget http://opengene.org/fastp/fastp.0.23.1
mv fastp.0.23.1 fastp
chmod a+x ./fastp
```


#### BWA
```
cd /home/hogehoge/local
git clone https://github.com/lh3/bwa.git
cd bwa; make
```

#### samtools
* htslibを含むソースからコンパイルする場合 (2022/01/19時点の最新版はv1.14)

```
cd /home/hogehoge/local
wget https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2
tar -xvf samtools-1.14.tar.bz2
cd samtools-1.14
./configure
make
```
#### Trimmomatic
* ver.0.39をbinaryでインストール

```
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
cd Trimmomatic-0.39
ls
## LICENSE       adapters        trimmomatic-0.33.jar
```

---
---
<h2 id="よく使うシェルコマンド">よく使うシェルコマンド</h2>

ディレクトリ操作系
```
cd ~/ #ホームディレクトリに移動
cd hogehoge #ディレクトリhogehogeに移動
cd ../ #階層が一つ上のディレクトリに移動
pwd #現在のディレクトリを表示
```
ファイル操作系
```
cp file1 file2 #file1をfile2としてコピー
rm file1 #file1を削除
mv file1 file2 #file1をfile2に移動（file１は消える）
```
ファイル圧縮解凍系
```
gzip input.fastq #input.fastqをgzに圧縮
gzip -d input.fastq.gz #解凍
```
ファイル表示系
```
cat file1 #無圧縮ファイル file1 を画面に出力
lesss file1 #無圧縮ファイル file1 をスクロールしながら見る
gzip -dc file1.gz #gzip圧縮ファイルを画面に出力
gzip -dc file1.gz | less #gzip圧縮ファイルをスクロールしながら見る
head file1 #無圧縮ファイル file1 の冒頭を画面出力
head -n 4 file1 #無圧縮ファイル file1 の冒頭の４行を画面出力
grep -c '^>' input.fasta | wc -l #fastaの配列数を表示
```

Linuxコマンド(Bash)でバックグラウンド実行する方法のまとめメモ https://qiita.com/inosy22/items/341cfc589494b8211844


[ページトップに戻る](https://github.com/akihirao/how2cook/tree/main/ngs_training#NGSデータ解析チュートリアル)
