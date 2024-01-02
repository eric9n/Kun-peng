# kranken2-rust


## 从 NCBI 下载数据
```
cargo run --bin ncbi -- -h
   Compiling ncbi v0.1.0 (/Users/eric/kraken2-rust/ncbi)
    Finished dev [unoptimized + debuginfo] target(s) in 0.79s
     Running `target/debug/ncbi -h`
ncbi download resource

Usage: ncbi [OPTIONS]

Options:
  -l, --list                   列出 NCBI 站点上的种类列表信息，实时拉取
  -d, --database <DATABASE>    构建数据库的目录 [default: lib]
  -g, --group <GROUP>          从 NCBI 站点上下载某个种类的数据信息，必须是列表中所列名称
  -c, --check-md5 <CHECK_MD5>  检查文件的 md5 [default: true] [possible values: true, false]
  -p, --parallel <PARALLEL>    下载时的并行大小 [default: 8]
  -h, --help                   Print help (see more with '--help')
  -V, --version                Print version

```

### md5 文件校验
```
cargo run --bin ncbi_md5 -- -h
    Finished dev [unoptimized + debuginfo] target(s) in 0.07s
     Running `target/debug/ncbi_md5 -h`
ncbi check genomics file md5sum

Usage: ncbi_md5 [OPTIONS]

Options:
  -d, --database <DATABASE>  数据库的路径 [default: lib]
  -g, --group <GROUP>        从 NCBI 站点上下载某个种类的数据信息，必须是列表中所列名称
      --delete               删除校验错误的文件
  -h, --help                 Print help (see more with '--help')
  -V, --version              Print version
```
