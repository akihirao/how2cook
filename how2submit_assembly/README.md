# ゲノムアッセンブリーのDDBJへの登録方法

新規のゲノムアッセンブリーをMass Submission Ssytem (MSS)を介してDDBJへ登録する方法のメモ

## [MSSの概要](https://www.ddbj.nig.ac.jp/ddbj/mss.html)

## 準備するもの
* ゲノムアッセンブリーの塩基配列ファイル（fasta形式: .fasta .seq.fa .fa .fna .seq）
* アノテーションファイル（プレインテキスト: .ann .annt.tsv ann.txt
-　塩基配列ファイルとアノテーションファイルの名称は、拡張子を除く名前が同一なペアになるように準備しておく

## MSS登録の流れ
- 1. ゲノムアッセンブリーの塩基配列ファイルをfasta形式で作成
- 2. アッセンブリーギャップのポジションのリストをbed形式で作成
```
gap_seq="NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
seqkit locate hogehoge.fasta -p $gap_seq --bed --only-positive-strand > hogehoge.ann.txt
```
- 3. [アノテーションファイル作成用スクリプト](make_annotation_file_DDBJ.pl)を用いてアノテーションファイルを作成
- 4.  D-wayアカウントで https://mss.ddbj.nig.ac.jp/ にログイン
- 5.  MSS Form に記入、ファイルアップロード、査定の開始

## [アノテーションファイル作成用スクリプト: make_annotation_file_DDBJ.ol](make_annotation_file_DDBJ.pl)
* 同じフォルダにゲノムアッセンブリーのfastaファイルとアセンブリーギャップのbedファイルおよび当該スクリプトを置いて、実行
```
make_annotation_file_DDBJ.pl hogehoge.fasta hogehoge.bed > hogehoge.ann.txt
```
