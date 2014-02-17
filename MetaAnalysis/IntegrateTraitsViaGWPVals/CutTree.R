library("dynamicTreeCut")
M <- read.table('Data/JaccardArrayAlpha0.241.csv',sep=',',header=T)
D <- as.dist(1.0-M)
X <- hclust(D, method = "average", members = NULL)
# plot(X, labels = NULL, hang = 0.1, axes = TRUE, frame.plot = FALSE, ann = TRUE, main = "Cluster Dendrogram", sub = NULL, xlab = NULL, ylab = "Height")

ctDIndex        <- cutreeDynamicTree(X, maxTreeHeight = 1, deepSplit = 0, minModuleSize = 1)
Trait           <- colnames(M)
Output          <- as.data.frame(Trait)
Output$ctDIndex <- ctDIndex
write.csv(Output,'Data/TreeCutJaccardArrayAlpha0.241DeepSplit0.csv',col.names=FALSE,row.names=FALSE)

ctDIndex        <- cutreeDynamicTree(X, maxTreeHeight = 1, deepSplit = 4, minModuleSize = 1)
Trait           <- colnames(M)
Output          <- as.data.frame(Trait)
Output$ctDIndex <- ctDIndex
write.csv(Output,'Data/TreeCutJaccardArrayAlpha0.241DeepSplit4.csv',col.names=FALSE,row.names=FALSE)