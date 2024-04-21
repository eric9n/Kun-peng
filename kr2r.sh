

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


# 3. 生成seqid2taxid.map, 如果不指定输出文件路径,则默认保存在数据库library下
${DIR}/squid seqid2taxid --source $DATABASE

# # 4. 预估数据库大小
${DIR}/squid estimate_capacity --source $DATABASE -k 35 -l 31

# 5. build
${DIR}/squid build --source $DATABASE -H ${DATABASE}/hash.k2d -o ${DATABASE}/opt.k2d -t ${DATABASE}/taxo.k2d -m ${DATABASE}/seqid2taxid.map -n ${DATABASE}/taxonomy --chunk-dir ${DATABASE_CHUNK}
