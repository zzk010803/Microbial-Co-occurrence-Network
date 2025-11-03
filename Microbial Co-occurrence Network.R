rm(list = ls())
getwd()
setwd("C:\\Users\\周智康\\Desktop\\Rdata")

#对OTU进行筛选
#（1）去掉平均相对丰度低于0.01%
#（2）出现次数少于总样本量1/5的OTU
otu_table <- read.csv("otu_table_nonGC.csv",header = T,row.names = 1)  #读取文件

rel_abundance <- apply(otu_table, 2, function(x) x/sum(x))  # 计算相对丰度

rel_abundance

write.csv(rel_abundance, 'non-GC_rel_abundance.csv') # 导出OTU的相对丰度表

mean_rel_abundance <- rowMeans(rel_abundance)    # 计算各个OTU在每个样本中的相对丰度

low_rel_abundance_otu <- rownames(otu_table)[mean_rel_abundance < 0.0001]  # 找到平均相对丰度小于0.01%的OTU

otu_table_filtered <- otu_table[!(rownames(otu_table) %in% low_rel_abundance_otu), ]  # 删除平均相对丰度低的OTU

freq <- apply(otu_table_filtered, 1, function(x) sum(x > 0)/length(x))

keep <- freq >= 1/5  # 根据需要改边需要的出现频率

otu_table_filt <- otu_table_filtered[keep, ] # 仅保留出现频率大于设定阈值的OTU

write.csv(otu_table_filt, 'non-GC_otu_table_filt.csv')   # 导出最终用于网络分析的OTU表





otu<-otu_table_filt

library(WGCNA)
library(psych)
library(reshape2)
library(igraph)

cor = corAndPvalue(t(otu),y=NULL,use = "pairwise.complete.obs", 
                   alternative='two.sided',method='spearman')  #OTU之间的Spearman相关系数和p值

r = cor$cor # 获取相关系数
p = cor$p                       #获取p值
p = p.adjust(p, method = 'BH')  #对p值进行BH校正

r[p > 0.001 | abs(r) < 0.60] = 0  # 对相关性进行筛选，p值>0.001或|r|<0.60的将被去除（赋0值）

write.csv(data.frame(r, check.names = FALSE), 'corr.matrix.csv')     # 将相关系数矩阵写入csv文件中          

g = graph_from_adjacency_matrix(r,mode="undirected",weighted=TRUE,diag = FALSE) #根据相关系数矩阵创建一个加权无向图
g
g = delete.vertices(g, names(degree(g)[degree(g) == 0]))        #删除度数为0的孤立节点

E(g)$corr = E(g)$weight            #为网络的边属性赋值（权重）
E(g)$weight = abs(E(g)$weight)     #为网络的边属性赋值（权重）


tax = read.csv('otu_tax.csv', row.names=1, header=T)     #读取节点分类信息  

tax = tax[as.character(V(g)$name), ]                  #为节点加上分类信息
V(g)$Kingdom = tax$Kingdom                            #界
V(g)$Phylum = tax$Phylum                              #门
V(g)$Class = tax$Class                                #纲
V(g)$Order = tax$Order                                #目
V(g)$Family = tax$Family                              #科
V(g)$Genus = tax$Genus                                #属
V(g)$Species = tax$Species                            #种

node_list = data.frame(
  label = names(V(g)),
  kingdom = V(g)$Kingdom,
  phylum = V(g)$Phylum,
  class = V(g)$Class,
  order = V(g)$Order,
  family = V(g)$Family,
  genus=V(g)$Genus,
  species = V(g)$Species)                              #创建节点列表

head(node_list)
write.csv(node_list, 'network.node_list.csv')          #并将其写入csv文件中

edge = data.frame(as_edgelist(g))                      #创建边列表
edge_list = data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  correlation = E(g)$corr
)

head(edge_list)
write.csv(edge_list, 'network.edge_list.csv')          #并将其写入csv文件中

write.graph(g, 'network.graphml', format = 'graphml')  #后续在Gephi中可视化                         


######计算网络常用的几种拓扑系数#####
nodes_num = length(V(g))                   #节点数
nodes_num

edges_num = length(E(g))                   #边数
edges_num

positive.cor_num = sum(E(g)$corr>0)        #正相关的数量
positive.cor_num

negative.cor_num = sum(E(g)$corr<0)        #负相关的数量
negative.cor_num

average_degree = mean(degree(g))           #平均度
average_degree

average_path_length = average.path.length(g, directed = FALSE)     #平均路径长度
average_path_length

network_diameter = diameter(g, directed = FALSE)                   #网络直径
network_diameter

network_density = graph.density(g)                                 #网络密度
network_density

clustering_coefficient = transitivity(g)                           #聚类系数
clustering_coefficient

network_parameter = data.frame(nodes_num, 
                               edges_num, 
                               positive.cor_num, 
                               negative.cor_num, 
                               average_degree,
                               average_path_length,
                               network_diameter, 
                               network_density,
                               clustering_coefficient                               
)

network_parameter
write.csv(network_parameter, 'network_parameter.csv')                                  

otu1 = otu
otu1[otu1>0] = 1
