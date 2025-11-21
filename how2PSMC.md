# PSMCソフトウェアの使い方

## ワークフローの例

1. コンセンサス配列のfastqの生成
```bash
bcftools mpileup --thread 4 -q 30 -Q 30 -Ou -f reference.fa hogehoge.bam | bcftools call  --thread 4 -c - | vcfutils.pl vcf2fq -d $min_CV -D $max_CV | gzip -c > hogehoge.fq.gz
```
$min_CV: 最低カバレッジの閾値（推奨値：平均カバレッジの1/3程度)   
$max_CV: 最大カバレッジの閾値（推奨値：平均カバレッジの２倍程度）   

現行のbcftools（2025/11/05時点）には vcfutils.plは含まれていないが、過去バージョンのbcftoolsに同包されているものを使用できる。samtoolsのmpileupを用いる場合はsamtools v.1.2を用いること！現行バージョンのsamtoolsでは出力形式の変更のためかエラーとなるので注意。


2. Pamcfa形式ファイルの生成
Psmcfaファイルには２つのハプロタイプ間で変異があるかどうか(デフォルトでは100bp間隔毎に)の情報が記録される。
```bash
fq2psmcfa hogehoge.fq.gz > hogehoge.psmcfa
```

3. PSMCの実行
```bash
psmc -N 25 -t 25 -r 5 -p "6+23*2+6" -o hogehoge.psmc hogehoge.psmcfa
```
N: number of rounds of iterations  
t: an N₀-scaled maximum coalescent time 
r: the ratio of the scaled mutationrate and the recombination rate


4. ブートストラップ用の分割データの生成
```bash
splitfapsmc hogehoge.psmcfa > hogehoge.split.psmcfa
```
デフォルトではランダム・リサンプリングの単位として、500kb毎に分割: trunk_size=500000


5. ブートストラップの実行
```bash
for(i = 1; i < 100; i++); do
	echo $i
	psmc -N30 -t15 -r5 -b -p "6+23*2+6" -o hogehoge.round.$i.psmc $hogehoge.split.psmcfa
done
```
100個のpsmc出力ファイルを集約して、Rのggplotなどを使って視覚化する。


## Atomic time intervals(p)の設定について
区間の設定は難しい問題で、ある程度の試行錯誤が必要。本家ページの設定の例はヒト用であり、そのままの設定を他の生物種に用いても望ましい結果が得られないかもしれない。信頼のおける推定結果の目安として、各区間の組み替えの回数が10以上であることが推奨されている("The `-p' and `-t' options are manually chosen such that after 20 rounds of iterations, at least ~10 recombinations are inferred to occur in the intervals each parameter spans.") 。  

[区間毎の組み替え回数の求め方に関する本家ページのQ&A](https://github.com/lh3/psmc/issues/45)

PSMCにおけるアーティファクトとして、false population size peaks (急激な集団サイズの上昇とそれに続く崩壊的な減少)がしばしば発生することが、実データおよびシミュレーションによって示されている (Hilgers et al 2025)。その対応方法として、first time intervalの分割設定（たとえば '2+2+2+23x2+6') が有効とされるが、万能の解決策でないことに注意。

## 参考資料
* Li & Durbin (2011) Inference of human population history from individual whole-genome sequences. Nature: 475, 493–496. https://doi.org/10.1038/nature10231 
* Hilgers et (2025) Avoidable false PSMC population size peaks occur across numerous studies. Current Biology: 35, 927–930. https://doi.org/10.1016/j.cub.2024.09.028
* [PSMCの本家ページ](https://github.com/lh3/psmc)   
* [「ゲノム多様性解析」](https://www.morikita.co.jp/books/mid/026171)（長田直樹編著、森北出版、2025)]  
第8章「集団サイズの推定」において、PSMCの解説とコード例が紹介。
