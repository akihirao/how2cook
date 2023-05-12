# How to extract CDS from a genbank format file with using gbseqextractor
## Set up gbseqextractor

How to install gbseqextractor (in Japanese)
https://kazumaxneo.hatenablog.com/entry/2023/02/10/111950

## Execution

```
# conda activate gbseqextractor
# Extract CDS
~/.local/bin/gbseqextractor -f hogehoge.gb -type CDS -prefix hogehoge

# Extract rRNA
~/.local/bin/gbseqextractor -f hogehoge.gb -type rRNA -prefix hogehoge

# Extract tRNA
~/.local/bin/gbseqextractor -f hogehoge.gb -type tRNA -prefix hogehoge

```