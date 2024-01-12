# kraken2-rust


## 从 NCBI 下载数据
```
$ ./ncbi -h
ncbi download resource

Usage: ncbi [OPTIONS] <COMMAND>

Commands:
  taxonomy  从 NCBI 下载 taxonomy 文件 (alias: tax)
  genomes   从 NCBI 下载 genomes 数据 (alias: gen)
  help      Print this message or the help of the given subcommand(s)

Options:
  -d, --database <DATABASE>        构建数据库的目录 [default: lib]
  -n, --num-threads <NUM_THREADS>  下载时的并行大小 [default: 8]
  -h, --help                       Print help (see more with '--help')
  -V, --version                    Print version
```

```
$ ./ncbi gen -h
从 NCBI 下载 genomes 数据 (alias: gen)

Usage: ncbi genomes [OPTIONS] --group <GROUP> [COMMAND]

Commands:
  md5   仅检查文件的 md5
  fna   解析 genomic 文件，并且生成 library fna 文件
  help  Print this message or the help of the given subcommand(s)

Options:
      --site <SITE>
          从 NCBI 哪个站点目录下载（RefSeq或GenBank）

          [default: refseq]

          Possible values:
          - genbank: 下载 genbank 资源
          - refseq:  下载 refseq 资源

  -g, --group <GROUP>
          从 NCBI 站点上下载某个种类的数据信息，可以是逗号分隔的多个, archaea,bacteria,viral,fungi,plant,human,protozoa

  -h, --help
          Print help (see a summary with '-h')
```
