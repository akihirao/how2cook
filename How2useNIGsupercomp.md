# 遺伝研スパコンへのログインと共有ディレクトリーの設定


## 遺伝研スパコンへのログインの流れ
遺伝研公式ページ　https://sc.ddbj.nig.ac.jp/guides/using_general_analysis_division/  

![遺伝研スパコンへのログインの流れ](https://sc.ddbj.nig.ac.jp/assets/images/GA_division-04fd98f77e98bcc62b31c8b43d6ecb11.png)

事前にアカウントの取得および公開鍵・秘密鍵の設定を済ませておく　　

1. ゲートウェイノードへのログイン
```bash
$ ssh hogehoge@gw.ddbj.nig.ac.jp
```

2. インタラクティブノードへのログイン
ゲートウェイノードへsshログイン後、インタラクティブノード a001, a002,a003のいずれかにSSHでログインする。
```bash
ssh a001　# a001へのログインの例
```

インタラクティブノードにログイン後に、スパコン上での計算作業をおこなうこと。

***
## 共有ディレクトリーの設定
現在のディレクトリの場所を確認する。
```bash
pwd
```
インタラクティブノードにログインした直後ならば、アカウント名 hogehogeでは次のディレクトリが表示される。
```bash
/home/hogehoge
```

hogehogeフォルダの直下に共有作業ディレクトリ（たとえばkokemomoと名付ける）を新規に作成。
```bash
mkdir ~/kokemomo # ~/はホームディレクトリを示す
```
共有フォルダkokemomoでの作業を行う前に、newgrpコマンドでグループを変更してから作業する。


```bash
newgrp lingonberry　#グループ名がlingonberryの場合
```
kokemomoフォルダをグループ lingonberryで共有。
```bash
chgrp -R lingonberry ~/kokemomo
```
共有フォルダkokemomoにグループのパーミションとして、読み込み/書き出し／実行（read/write/exe）を許可。
```bash
chmod g+x ~/kokemomo
```
フォルダの実行権限を確認。
```bash
ls -la ~/kokemomo
```

***
## データ転送
遺伝研公式ページ　https://sc.ddbj.nig.ac.jp/guides/using_general_analysis_division/ga_data_transfer/

### SSH プロトコルによるファイル転送の方法 (scp, sftp)

自分のPCにて、ターミナルを介してカレントディレクトリにある your_file.txt ファイルを遺伝研スパコンに scp するには、 以下のコマンドを実行する。
```bash
scp your_file.txt hogehoge@gw.ddbj.nig.ac.jp:/home/hogehoge
```

カレントディレクトリに A001_R1.fastq.gz,A001_R2.fastq.gz, A002_R1.fastq.gz, A002_R2.fastq.gzなどと複数のリードファイルがあり、それらを一挙にスパコン上のkokemomoフォルダに転送するには、以下のコマンドを実行する。
```bash
scp *.fastq.gz hogehoge@gw.ddbj.nig.ac.jp:/home/hogehoge/kokemomo
```


***
## Slurmによるジョブの実行
遺伝研公式ページ　https://sc.ddbj.nig.ac.jp/guides/software/JobScheduler/Slurm/

バッチスクリプトの形でジョブを投げる。
```bash
man sbatch
sbatch --help
sbatch test_run.sh
```

***
## singularityを介したソフトウェアの利用

遺伝研スパコンでは、/usr/local/biotools/ 以下に各種ソフトウェアが用意されている。
利用可能なソフトウェアとバージョンを探す:
```bash
find /usr/local/biotools/ -name 'fastp*' | sort #ソフトウェアfastpの例
```

イメージとプログラム名を指定して実行:
```bash
singularity exec /usr/local/biotools/f/fastp:1.0.1--heae3180_0 fastp --help
```

fastpを用いてfastq.gz形式の生リードをQCフィルタリングするバッチスクリプトの例
```bash
#!/bin/bash
# exe_fastp.sh

##ノード数を1に指定。ジョブを1つのノードで実行
#SBATCH -N 1

#**タスク数（MPIプロセス数）**を1に指定。並列処理をしない単一プロセスのジョブ
#SBATCH -n 1

#各タスクに割り当てるCPUコア数を8に指定。マルチスレッド処理をしない前提
#SBATCH -c 8

#最大実行時間を「0日1時間0分0秒」に設定。これを超えるとジョブは強制終了
#SBATCH -t 00-00:01:00

#各CPUコアに割り当てるメモリ量を16GBに設定。1コアなので、合計1GBのメモリの割当。
#SBATCH --mem-per-cpu=24G

#ジョブ名を user_go に設定。ジョブ管理やログファイル名に使用される
#SBATCH -J fastp


#########YOUR JOB#############
n_threads=8

RawRead_folder=${HOME}"/temp/RawReads"
QCRead_folder=${HOME}"/temp/QCReads"

mkdir -p $QCRead_folder
cd $QCRead_folder


alias fastp="exec /usr/local/biotools/f/fastp:1.0.1--heae3180_0 fastp"


while read sequence_ID sample; do

	if [ ! -e $QCRead_folder/$sample ]; then
		mkdir $QCRead_folder/$sample
	fi

	R1_tag="_1"
	R2_tag="_2"
	sequence_R1=$sequence_ID$R1_tag
	sequence_R2=$sequence_ID$R2_tag
	sample_R1=$sample$R1_tag
	sample_R2=$sample$R2_tag

	fastp -i $RawRead_folder/$sequence_R1.fastq.gz\
	 -I $RawRead_folder/$sequence_R2.fastq.gz -3\
	 -o $QCRead_folder/$sample/$sample_R1.fastp.trim.fastq.gz\
	 -O $QCRead_folder/$sample/$sample_R2.fastp.trim.fastq.gz\
	 --trim_poly_g --cut_right --cut_window_size 4 --cut_mean_quality 20\
	 --average_qual 20 -q 20 -l 50 -f 5 -F 5 -t 5 -T 5\
	 -h $QCRead_folder/$sample/$sample.trimQ20.fastp.html\
	 -w $n_threads
	# -g or --trim_poly_g: polyG tail trimming
	# --cut_right: quality pruning by sliding window (--cut_right)
	#-3:  enable per read cutting by quality in tail (3'), default is disabled
	#-f: trimming how many bases in front for read1, default is 0 (int [=0])
	#-F: trimming how many bases in front for read2, default is 0 (int [=0])
	#-t: trimming how many bases in tail for read1, default is 0 (int [=0])
	#-T: trimming how many bases in tail for read2, default is 0 (int [=0])

done < ${HOME}"/temp/sample_list.txt" #list of ID
```

***
## 他の参考ページ
* 岩嵜 航さんの解説メモ　https://heavywatal.github.io/bio/nig.html


