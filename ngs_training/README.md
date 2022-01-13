# NGSデータ解析チュートリアル


平尾　章

水産研究・教育機構 水産資源研究所
水産資源研究センター 漁業情報解析部　資源解析グループ


----

## 概要

初心者向けのチュートリアルとして、酵母 (Saccharomyces cerevisiae) のリシーケンス解析を題材に公開データの取得から変異検出までの解析手順の流れを学びます。 


## 本チュートリアルの流れ
1.　公開データ取得
2.　クオリティーコントロール
3.　マッピング　
4.　変異検出　



### 使用するNGSツール

* bedtools: a powerful toolset for genome arithmetic https://bedtools.readthedocs.io
* BWA: Burrow-Wheeler Aligner http://bio-bwa.sourceforge.net
* Bwa-mem2: the next version of the bwa-mem https://github.com/bwa-mem2/bwa-mem2) 
* fastp: an all-in-one preprocessing tool for fastq files (https://github.com/OpenGene/fastp
* fastqc: a quality control tool https://www.bioinformatics.babraham.ac.uk/projects/fastqc
* GATK: Genome Analysis Toolkit https://gatk.broadinstitute.org
* Plink: whole-genome association analysis tool https://www.cog-genomics.org/plink
* samtools: tools for manipulating NGS data https://github.com/samtools/samtools
* Trimmomatic: a flexible read trimming tool https://github.com/usadellab/Trimmomatic
* vcftools: a set of tools for working with VCF files https://github.com/vcftools/vcftools


### NGSツールのインストールについて
...



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


---
---

## リシーケンス解析チュートリアル
出芽酵母(　*Saccharomyces* *cerevisiae* )は真核生物として初めてゲノム解読されたモデル生物であり、ゲノムサイズが小さい(12.1Mb)ため、塩基配列データもコンパクトになります。そこで酵母のリシーケンスデータを用いて、変異検出までの解析処理の流れを学びます。

### 1. 公開データ取得
NGSから出力された生シーケンスデータの多くは、科学研究の再現性担保およびデータ解析による新たな発見を支えるために公共データベース（DRA/SRA/ERA)にアーカイブされています。ここでは出芽酵母のリシーケンス(ERR038793)の生リードデータをSRA-toolkitを用いてダウンロードしてみましょう。



### 2. クオリティーチェック



### 3. マッピング



### 4. 変異検出
