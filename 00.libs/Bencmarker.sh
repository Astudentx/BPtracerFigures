cd /Users/astudentx/Documents/01.Work/01.Analysis/01zyz_240514_BP-Tracer2/Visualization/01.rawdata/BP2_BK/98.BKTest
# benchMarekr测试
python3 BenchmarkerNew3.py -t TprofileTax.txt -p TaxIDAbu.S -o TaxIDAbu.S.BK --normalize --info NofilterFunction
less TaxIDAbu.S.BK



# 保留平均相对丰度 > 0.01% 的物种
python3 filter_abundance.py -i TaxIDAbu.S -o TaxIDAbu.S.f1_0.0001 -t 0.0001 -f 1 --report --exclude-cols ID,Taxonomy,Lineage
# 保留平均相对丰度 > 0.1% 的物种
python3 filter_abundance.py -i TaxIDAbu.S -o TaxIDAbu.S.f1_0.001 -t 0.001 -f 1 --report --exclude-cols ID,Taxonomy,Lineage

# 保留总绝对丰度 > 1000 的物种
python3 filter_abundance.py -i TaxIDAbu.S -o TaxIDAbu.S.f2_1000 -t 1000 -f 2 --report --exclude-cols ID,Taxonomy,Lineage
# 保留缺失样本数 ≤ 5 的物种
python3 filter_abundance.py -i TaxIDAbu.S -o TaxIDAbu.S.f3_5 -t 5 -f 3 --report --exclude-cols ID,Taxonomy,Lineage
# 保留至少在一个样本中相对丰度 > 0.01% 的物种
python3 filter_abundance.py -i TaxIDAbu.S -o TaxIDAbu.S.f4_0.0001 -t 0.0001 -f 4 --report --exclude-cols ID,Taxonomy,Lineage
# 保留至少在Reads数量 > 50 的物种
python3 filter_abundance.py -i TaxIDAbu.S -o TaxIDAbu.S.f5_50 -t 50 -f 5 --report --exclude-cols ID,Taxonomy,Lineage



# 过滤策略测试
# 我的主要过滤策略为：1. 过滤掉假阳性的序列；2. 过滤掉相对丰度比较低的物种。
# 所以应该大概是 -f5(0,25,50,75,100,125,150,175,200,225,250,275,300); f1=c(0, 1e-05,1e-05,1e-06,1e-07,1e-08)

python3 filter_abundance_multi.py \                                                     
  -i TaxIDAbu.S \
  --f5 "0,25,50,75,100,125,150,175,200,225,250,275,300" \
  --f1 "0, 1e-04,1e-05,1e-06,1e-07,1e-08" \
  --order "f5,f1" \
  --prefix "BacBKStep1" \
  --exclude-cols "ID,Taxonomy,Lineage" \
  --benchmarker \
  -t TprofileTax.txt
# 合并所有结果

# 保留第一个文件的完整内容，然后从第二个文件开始跳过第一行(表头)追加
head -n 1 BacBKStep1.f5_0.f1_0.BK > Final.BacBKStep1.txt
tail -n +2 -q BacBKStep1*.BK >> Final.BacBKStep1.txt



# Batch A
python3 filter_abundance_multi.py \
  -i TaxIDAbu.S \
  --f5 "0,10,20,25,30,40,50" \
  --f1 "6e-05,8e-05,1e-04,1.2e-04,1.4e-04" \
  --order "f5,f1" \
  --prefix "BacBKStep2A" \
  --exclude-cols "ID,Taxonomy,Lineage" \
  --benchmarker \
  -t TprofileTax.txt

# Batch B
python3 filter_abundance_multi.py \
  -i TaxIDAbu.S \
  --f5 "120,125,130,140,150,160,170,175,180" \
  --f1 "6e-05,8e-05,1e-04,1.2e-04,1.4e-04" \
  --order "f5,f1" \
  --prefix "BacBKStep2B" \
  --exclude-cols "ID,Taxonomy,Lineage" \
  --benchmarker \
  -t TprofileTax.txt
  
# 保留第一个文件的完整内容，然后从第二个文件开始跳过第一行(表头)追加
head -n 1 BacBKStep1.f5_0.f1_0.BK > Final.BacBKStep1.txt
tail -n +2 -q BacBKStep*.BK >> Final.BacBKStep1.txt


  

