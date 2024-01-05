# kraken2-rust


## 从 NCBI 下载数据
```
$ ./ncbi -h
ncbi download resource

Usage: ncbi [OPTIONS]

Options:
  -l, --list                 列出 NCBI 站点上的种类列表信息，实时拉取
  -d, --database <DATABASE>  构建数据库的目录 [default: lib]
  -g, --group <GROUP>        从 NCBI 站点上下载某个种类的数据信息，必须是列表中所列名称，archaea,bacteria,fungi...
  -m, --md5                  仅检查文件的 md5
  -t, --threads <THREADS>    下载时的并行大小 [default: 8]
  -h, --help                 Print help (see more with '--help')
  -V, --version              Print version

```

