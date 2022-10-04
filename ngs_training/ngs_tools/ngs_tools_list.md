## 使用NGSツールのリスト

* bedtools: a powerful toolset for genome arithmetic https://bedtools.readthedocs.io
* BWA: Burrow-Wheeler Aligner http://bio-bwa.sourceforge.net
* Bwa-mem2: the next version of the bwa-mem https://github.com/bwa-mem2/bwa-mem2
* fastp: an all-in-one preprocessing tool for fastq files https://github.com/OpenGene/fastp
* fastqc: a quality control tool https://www.bioinformatics.babraham.ac.uk/projects/fastqc
* GATK: Genome Analysis Toolkit https://gatk.broadinstitute.org
* Plink: whole-genome association analysis tool https://www.cog-genomics.org/plink
* samtools: tools for manipulating NGS data https://github.com/samtools/samtools
* seqkit: an ultrafast toolkit for FASTA/Q file manipulation https://github.com/shenwei356/seqkit
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

#### bedtools (a pre-compiled binary for linux)
```
wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary
mv bedtools.static.binary bedtools
chmod a+x ./bedtools
```

#### BWA
```
cd /home/hogehoge/local
git clone https://github.com/lh3/bwa.git
cd bwa; make
```

#### samtools: 2022/09/12時点の最新版はv1.16.1
HTSlib入りのソースからコンパイルする例です。HTSlib projectの一部であるtabixとbgzipは使用頻度の高いツールなので、HTSlibもビルドしておきましょう。

```
cd /home/hogehoge/local
wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2
tar -xvf samtools-1.16.1.tar.bz2
cd samtools-1.16.1
./configure --prefix=/user/hogehoge/local/
make
./samtools --version
cd htslib-1.16.1
./configure --prefix=/user/hogehoge/local/
make
./tabix --version
./bgzip --version

```
#### Trimmomatic
Ver.0.39をbinaryでインストール

```
cd /home/hogehoge/local
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
cd Trimmomatic-0.39
ls
## LICENSE       adapters        trimmomatic-0.33.jar
```

#### Plink
Ver.1.9をbinaryでインストール

```
cd /home/hogehoge/local
wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20210606.zip
unzip plink_linux_x86_64_20210606.zip
```

#### Plink2
Ver.2.0をbinaryでインストール

```
cd /home/hogehoge/local
wget https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20211217.zip
unzip plink2_linux_x86_64_20211217.zip
```

#### vcftools
Build from release tarball (https://github.com/vcftools/vcftools/releases)
```
wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz
tar zxvf vcftools-0.1.16.tar.gz
cd vcftools-0.1.16
./configure
make
```
