DIR=`dirname $(realpath $0 || echo $0)`

DOWNLOADS=""
DATABASE=""
CHUNK_DIR=""


# 1. 下载 bacteria,viral 原始fna.gz格式文件和md5文件
${DIR}/ncbi -d $DOWNLOADS gen -g bacteria,viral

# 2. 下载taxonomy文件
${DIR}/ncbi -d $DATABASE taxonomy

# 3. build
${DIR}/kun_peng build -d $DATABASE --db $DATABASE

# 4. classify
./target/release/kun_peng classify --db $DATABASE --chunk-dir $CHUNK_DIR $the_sample_files
