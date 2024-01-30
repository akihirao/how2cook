# ゲノムアッセンブリーのDDBJへの登録方法

新規のゲノムアッセンブリーをMass Submission Ssytem (MSS)を介してDDBJへ登録する方法のメモ

## [MSSの概要](https://www.ddbj.nig.ac.jp/ddbj/mss.html)

## 準備するもの
1. ゲノムアッセンブリーのfastaファイル（拡張子は .fasta, .seq.fa, .fa, .fna, .seqに対応）
2. アノテーションファイル（[様式について](https://www.ddbj.nig.ac.jp/ddbj/file-format.html#annotation)）（拡張子は .ann, .annt.tsv, ann.txtに対応)
* 塩基配列ファイルとアノテーションファイルは、拡張子を除くそれぞれの名前が同じペアになるようにすること

## MSS登録の流れ
1. ゲノムアッセンブリーの塩基配列ファイルを準備
2. アッセンブリーギャップのポジションのリストを bedで作成
```
gap_seq=$(yes 'N' | head -n 100 | tr -d '\n') # when size of assembly gap is 100 bp
seqkit locate hogehoge.fasta -p $gap_seq --bed --only-positive-strand > assembly_gap.bed
```
3. Perlスクリプト [make_annotation_file_DDBJ.pl](make_annotation_file_DDBJ.pl) を用いてアノテーションファイルを作成
```
make_annotation_file_DDBJ.pl hogehoge.fasta assembly_gap.bed > hogehoge.ann.txt
```
* SUBMITTERやREFERENCE, DATEなどの書式については、それぞれのケースに合わせてスクリプトを編集した上で、使用してください
* このスクリプトではbiological featureの記載はsource と assembly_gap のみの項目を対象としていますので、CDSなどのアノテーションを付ける場合は別途追記ください
4. D-wayアカウントで https://mss.ddbj.nig.ac.jp/ にログイン
5. MSS Form に記入、ファイルアップロード、査定の開始
