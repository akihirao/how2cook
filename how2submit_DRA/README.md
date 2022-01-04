# リードデータのDRAへの登録方法

次世代シーケンサーの生データをDDBJ Seqeunce Read Archive (DRA)に登録するためのメモ

## DDBJの公式ハンドブック
[DDBJ Seqeunce Read Archive Handbook](https://www.ddbj.nig.ac.jp/dra/submission.html)
- Submission、Study、Experiment、Sample、Runの概念と対応関係を理解しておきましょう。

## DRAへの多検体登録
他検体サンプルのDNA試料にタグを付けてシーケンスした場合、サンプルごとのデータファイルに分割してDRAに登録する必要があります。その際に、サンプルについての情報(メタデータ)を記述したXMLファイルの作成するのですが、DDBJのD-wayアカウントシステムにて200以上の他検体のデータを登録しようとすると、アップデートに非常に時間がかかりタイムアウトエラーが頻発することがあるようです([ushioさんのブログ](https://ushio-ecology-blog.blogspot.com/search?q=DRA登録%E3%80%80メタ情報))。そこでushioさんは次の方法を提案されています。

* 1回のSubmissionで登録するオブジェクトの数を200程度になるように分けて登録する
* xmlファイルを自身で作成･編集して登録する

なおxmlファイルの作成支援用ツールとしては、DDBJから[submission-excel2xml](https://github.com/ddbj/submission-excel2xml)が提供されています。

