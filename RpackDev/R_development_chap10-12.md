---
backgroundColor: #D0DEEA
marp: true
---


# Rパッケージ開発入門

### 第10章：コンパイル済みのコード
### 第11章：インストール済みのファイル
### 第12章：その他のコンポーネント


担当：平尾章

---

# 第10章 コンパイル済みのコード


## 本章ではC++/CをRパッケージで使う方法を学びます :sparkles:
スクリプト言語としての表現の豊かさを誇るRのコードに、CやC++を組み込むことで、プログラムの実行速度の向上を保管することができます。

![width:500pt](https://github.com/akihirao/how2cook/blob/main/RpackDev/C_and_CPP_icons.png)


- Rcpp: RからC++を使うためのパッケージ
   - [本家ページ](https://www.rcpp.org/)
   - [岩嵜さんの解説](https://heavywatal.github.io/rstats/rcpp.html)
   - [みんなのRcpp by Masaki E. Tsuda](https://teuder.github.io/rcpp4everyone_ja/)

- cpp11: 主にtidyverseでC++を使うためのパッケージ
   - [本家ページ](https://www.rcpp.org/)


# 10.1 C++
Rcppを使う際は、まず最初に次のコマンドを実行
```
usethis::use_rcpp()

# 1) src/の生成
# 2) DESCRIPTIONのLinkingとImportsへのRcppの追記
# 3) .gitignoreの生成
```

パッケージ内の適当なRファイルに次のroxygenタグを追加
```
#' @useDynLib myPackage, .registration=TRUE
#' @importFrom Rcpp sourceCpp
NULL
```

document()の実行：上記のroxygenタグを介して、NAMESPACEが修正される
```
devtools::document()
```

## 10.1.1 ワークフロー

1. 新しいC++ファイルの作成

   - 例として、数を2倍にするtimesTwo関数を作成し、src/ディレクトリに適当な名前のファイル(hogehoge.cpp)で保存。
```
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int timesTwo(int x) {
  return x * 2;
}
```

   - なおRstudioにて New File > C++ File とすると、テンプレートに沿って新規C++ファイルが作成される。


2. ビルド＆ラッパー関数の作成
```
Rcpp::sourceCpp("src/hogehoge.cpp")
```


3. ラッパー関数（例：timesTwo）の動作確認

   - 引数（下記の例は20）の2倍が戻り値ならば、OK
```
timesTwo(20)
```

## 10.1.2 ドキュメント

- C++の各関数からRラッパー関数が自動生成される際に、事前にドキュメントをつけておくと便利
- ラッパー関数は基本的にR/RcppExports.Rに配置される

以下はtimeTwo関数の例
```
#include <Rcpp.h>
using namespace Rcpp;

//' 数を２倍する          #追加箇所
//'                     #追加箇所
//' @param x 1つの整数   #追加箇所
//' @RcppExports       #追加箇所
// [[Rcpp::export]]
int timesTwo(int x) {
  return x * 2;
}
```

```
devtools::document()
```
を実行すると、Rラッパー関数にroxygenタグが自動的に付与される
```
#' 数を２倍する
#'
#' @param x 1つの整数
#' @export
timesTwo <- function(x) {
    .Call(`_mypackage_timesTwo`, x)
}
```


## 10.1.3 C++コードのエクスポート

- 作成したC++のコードを他のパッケージのC++のコードから呼び出せるようにするには次の属性を追加する
```
//　[Rcpp::interfaces(r, cpp)]]
```

inst/include/mypackage.hというヘッダーファイルが生成され、他のパッケージにインクルードされる



## 10.1.4 C++コードのインポート

1. DESCRIPTIONに”LinkingTo: otherPackage”との記述を追加する

2. C++ファイル側に以下を追加
```
#include <otherPackage.h>
```

3. 外部パッケージotherPackageのC++関数の使い方
- otherPackage::foo()として使用
- using namespace otherPackageとしてグローバルに利用


## 10.1.5 ベストプラクティス

- 出力を行うときは cout << … ではなく　Rcout << …を使う
- 長時間の実行ループでは、定期的にRcpp::checkUserInterrupt()を使う
- ヘッダーやインクルードファイルの拡張子には.hを使う
- [Portable C++ for R Packages](https://journal.r-project.org/articles/RJ-2011-020/)での推奨事項に従う
- パッケージでC++コードを使うときは常に、パッケージがアンロードされた後にクリーンアップを必ず行う
- C++のコンパイルにはgccではなくclangを使ったほうがよい



# 10.2 C
### 基本的にCよりもC++の使用がおすすめ
### それでもCを選択する理由
- CのAPIを使っている古いパッケージを利用している場合
- 自分のパッケージが既存のCのライブライブと直接連携する場合

## 10.2.1 .Call()入門
.Call()とはCのコードをRから呼び出すための関数

例として、２つの数値を加算する次のCのコードをmysum.cという名前でsrc/に保存
```
#include <R.h>
#include <Rinternals.h>

SEXP add_(SEXP x_, SEXP y_){
    double x = asReal(x_);
    double y = asReal(y_);

    double sum = x + y;

    return ScalarReal(sum);
}
```

Rラッパー関数を作成して、my_add.RとしてR/に保存
```
#’ @useDynLib mypackage add_
my_add <- function(x,y) .Call(add_, x,y)
```

devtools::load_all()してから、動作確認
```
my_add(1,2)
```


## 10.2.2 .C()入門

本節のハンズオンは省略(テキストのコードの例では動作確認できなかったため)

参考までにテキストの例を以下に再掲する。

２つの数値を加算する次のCのコードをsrc.cという名前でsrc/に保存
```
void add_(double* x, double* y, double* out){
  out[0] = x[0] + y[0];
}
```

Rラッパー関数を作成して、my_add.RとしてR/に保存
```
#' @mypackage src.c add_
add <- function(x ,y) {
   .C(add_, x, y, numeric(1))[[3]]
}
```


## 10.2.3 ワークフロー

1. 新しいCファイルの作成
   
2. ラッパー関数の作成

3. ラッパー関数（例：my_add）の動作確認

   ２つの引数の合計が戻り値ならば、OK
```
my_add(24,26)
```

## 10.2.4 Cコードのエクスポート
- Cコードと対応する再配置可能なDLL（=ディスクのどこにおいても動作するDLL)を準備すること！
- R_RegisterCCallable()関数の登録でDLLの提供が可能となる
```
#include “add.h”
#include <R_ext/Rdynload.h>

viod R_init_mypackage(DllInfo * info){
    R_RegisterCCallable(info, “add_”, (DL_FUNC) &add_
}
```


## 10.2.5 Cコードのインポート
- パッケージの実装方法に依存した準備をおこなう
- これまでの説明に沿った実装ならば、DESCRIPTIONファイルにLinkingTo: otherPackageを記載し、Cのファイルの内部で#include otherPackageAPIHを使う
- そのパッケージが関数の登録をしているが、ヘッダーファイルの提供をしていない場合は、自分でラッパー関数を書く必要がある
- そのパッケージが関数の登録をしていない場合は、それらの関数を使うことはできない


## 10.2.6 ベストプラクティス
- assert(), abort(), exit()の呼び出しを避ける
- 出力には printf()ではなく、Rprintf() を使う
- 長時間の実行ループでは、定期的にRcpp::checkUserInterrupt()を使う
- Cの乱数生成器ではなく、Rの乱数生成器が持つCのAPIなどを使う
- RのマクロISNAN(x)とR_FINTE(x)を使って、NaNと無限大値を確認する
- C++と同じように、パッケージのアンロード後にはDLLもアンロードする
- Cのコンパイルにはgccではなくclangを使ったほうがよい


# 10.3 コンパイル済みのコードのデバッグ
Rstudioは使えないので、コマンドラインから実行すること！

# 10.4 Makefile
この本の対象外

# 10.5 R以外の言語
Fortran, Java

# 10.6 ライセンス
他のライブラリーを使う際は、ライセンスの互換性に配慮すること！

# 10.7 開発ワークフロー
C/C++の開発ではdevtools::load_all()ではなく、RstudioのBuild & Reloadを使うほうがよい
Cオブジェクトはload_all()してもメモリに残るため、クラッシュを起こす可能性がある

# 10.8 CRANに関する課題
C/C++を用いるパッケージはCRANへの登録の難易度が高い


---

# 第11章 インストール済みのファイル


## 本章ではinst/フォルダの取り扱いを学びます


inst/に含まれる一般的なファイル

 - inst/AUTHORとinst/COPYRIGHT

   パッケージの著作権や著作者の情報


 - inst/CITATION

   パッケージの引用方法


 - inst/extdata

   Examplesやvignettesで用いる追加データの置き場所

 - inst/doc

   古いバージョンのRではvignetteの代わりに使用．現在は推奨されず


 - inst/java, inst/pythonなど

    R以外の言語で書かれた補助スクリプトの置き場所



## inst/内のファイルの探し方

system.file()を用いる．
例えば、 SPiCTのinst/doc/spict_handbook.pdf のパスを調べるには

```
system.file(“doc”, “spict_handbook.pdf”, packaged =“spict”)
```

## 11.1 パッケージの引用

```
# Rのbaseバッケージの引用方法を表示
>citation()

# Rのパッケージ名(例としてtidyverse)を指定して、引用方法を表示
>citation("tidyverse")
```

```
bibentry(
  "Article",
  title = "Welcome to the {tidyverse}",
  author = "Hadley Wickham, Mara Averick, Jennifer Bryan, Winston Chang, Lucy D'Agostino McGowan, Romain François, Garrett Grolemund, Alex Hayes, Lionel Henry, Jim Hester, Max Kuhn, Thomas Lin Pedersen, Evan Miller, Stephan Milton Bache, Kirill Müller, Jeroen Ooms, David Robinson, Dana Paige Seidel, Vitalie Spinu, Kohske Takahashi, Davis Vaughan, Claus Wilke, Kara Woo, Hiroaki Yutani",
  year = 2019,
  journal = "Journal of Open Source Software",
  volume = 4,
  number = 43,
  pages = 1686,
  doi = "10.21105/joss.01686",
)
```

CITATIONの雛形を作成する関数
```
usethis::use_citation()
```


## 11.2 R以外の言語

 - パッケージにはR以外の言語（perlやpythonなど）で書かれた補助スクリプトがふくまれている場合がある

 - 他言語の補助スクリプトは基本的にinst/のサブディレクトリに置く（inst/perl, inst/pythonなど）

 - 多言語を用いる場合は、DESCRIPTIONにてSystemRequirementsに対応するプログラミング言語を明記する

 - Javaの場合は特別に ソースコード(java/に置き、.Rinstignoreのリストに追加）とコンパイル済みの.jar(inst/javaに置く）の両方に準備する。またIMPORTSにrJAVAを追加しておく




---

# 第12章 その他のコンポーネント

 - demo/

   パッケージのデモ用ディレクトリ
   （現在はdemoではなくvignetteの使用が推奨） 

 - exec/ 
   
   実行可能なスクリプト用ディレクトリ

 - po/   
   
   メッセージの翻訳用ディレクトリ（本教科書の対象外）

 - tools/   
   
   設定時に必要な補助ファイルや、スクリプト生成に必要なソースなど
   の補助ファイルを置くためのディレクトリ


## 12.1 デモ

```
# 利用可能なすべてのデモのリストを表示
>demo()

# graphicsのデモを表示
>demo(graphics)
```

