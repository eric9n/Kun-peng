

DIR=`dirname $(realpath $0 || echo $0)`

DATABASE=""
DATABASE_CHUNK=""

# 1. 下载 bacteria,viral 原始fna.gz格式文件和md5文件
${DIR}/ncbi --db $DATABASE gen -g bacteria,viral

# 1.1 校验md5文件
${DIR}/ncbi --db $DATABASE gen -g bacteria,viral md5

# 1.2 下载taxonomy文件
${DIR}/ncbi --db $DATABASE taxonomy

# 2. 生成library.fna和prelim_map.txt子文件
${DIR}/ncbi --db $DATABASE gen -g bacteria,viral fna

# 3. 预估数据库大小
# ${DIR}/Kun estimate_capacity --db $DATABASE -k 35 -l 31

# 4. build
${DIR}/kun_peng build --db $DATABASE --chunk-dir ${DATABASE_CHUNK}
