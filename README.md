# kraken2-rust


## 从 NCBI 下载数据
```
$ ./ncbi -h
ncbi download resource

Usage: ncbi [OPTIONS] <SITE>

Arguments:
  <SITE>  [possible values: genbank, refseq]

Options:
  -d, --database <DATABASE>        构建数据库的目录 [default: lib]
  -t, --taxonomy                   下载 taxonomy 文件
  -g, --group <GROUP>              从 NCBI 站点上下载某个种类的数据信息，必须是列表中所列名称，archaea,bacteria,viral,fungi,plant,human,protozoa
      --md5                        仅检查文件的 md5
  -n, --num-threads <NUM_THREADS>  下载时的并行大小 [default: 8]
      --fna                        解析 genomic 文件，并且生成 library fna 文件
  -h, --help                       Print help (see more with '--help')
  -V, --version                    Print version
```

