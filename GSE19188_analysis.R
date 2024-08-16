#libraries:
library(MASS)
library(ggplot2)
#Read-in data
raw_data = readLines("C:/Users/Alex/Desktop/JHU Summer 2024/Gene expression data analysis and visualization/Final project/GSE19188_series_matrix.txt/GSE19188_series_matrix.txt")
data = read.delim("C:/Users/Alex/Desktop/JHU Summer 2024/Gene expression data analysis and visualization/Final project/GSE19188_series_matrix.txt/GSE19188_series_matrix.txt", comment.char = "!", header = TRUE, row.names = 1)

#prepare the metadata file
sample.title = grep("^!Sample_characteristics", raw_data, value = TRUE)
sample.title = strsplit(sample.title, "\t")[[1]][-1]
sample.names = colnames(data)
metadata = data.frame(sample.names, sample.title)
metadata$type = ifelse(grepl("healthy", metadata$sample.title), "healthy", "tumor")
names(data) = paste(colnames(data), metadata$type, sep = "-")
label = ifelse(grepl("healthy", names(data.sig)), "h", "t")

#test for outlier samples, provide visual proof (graphs)
data.cor = cor(data)
data.t = t(data)
dim(data)
#hierarchical clustering:
#GSM475776-tumor and GSM475785-tumor appear to be outliers based on dendogram, maybe gsm475802 too
data.dist = dist(data.t, method = "euclidean")
data.clust = hclust(data.dist, method = "single")
plot(data.clust, labels = names(data.t), cex = 0.55)

#remove outliers
outliers = c("GSM475776-tumor", "GSM475785-tumor")
data = data[, !colnames(data) %in% outliers]
dim(data)

#filter out genes with low expression values
gene.mean = apply(data, 1, mean)
plot(gene.mean, main = "Mean values all genes", ylab = "Mean")
abline(h = -.1)
data = data[gene.mean>-.01,]
dim(data)

#feature selection with statistical test or machine learning method
t.test.all.genes = function(x, s1, s2) {
  x1 = x[s1]
  x2 = x[s2]
  x1 = as.numeric(x1)
  x2 = as.numeric(x2)
  t.out = t.test(x1, x2, alternative = "two.sided", var.equal = T)
  out = as.numeric(t.out$p.value)
  return(out)
}

group.labels = names(data)
healthy = grep("healthy", group.labels)
tumor = grep("tumor", group.labels)
pv = apply(data, 1, t.test.all.genes, s1 = healthy, s2 = tumor)

#adjusted p-values
pv.bonferroni = p.adjust(pv, method = "bonferroni")
sum(pv.bonferroni<.05)
sum(pv<.01)
#plot scores of retained genes in histogram
pv.sort = sort(pv)
pv.bonferroni.sort = sort(pv.bonferroni)

#p-value comparison plot
plot(1:length(pv.sort), pv.sort, col = "red", type = "l", xlab = "Index", ylab = "P-value", main = "Normal vs Tumor plot: Bonferroni Adjustment")
lines(1:length(pv.bonferroni.sort), pv.bonferroni.sort, col = "blue")
legend("topleft", legend = c("Non-adjusted", "Bonferroni-adjusted"), col = c("red", "blue"), lty = 1)

#p-value after bonferroni adjustment
hist(pv.bonferroni, col = "lightblue", xlab = "Bonferroni-adjusted p-values", main = "P-value distribution between Normal and Tumor samples: Bonferroni adjustment", cex.main = 0.9)
significant = names(which(pv.bonferroni<.001))
length(significant)
#8843 genes left with a p-value cutoff of .001

#subset data by genes and use clustering to visualize
subset.data = data[significant,]
subset.pca = prcomp(t(subset.data))
pca.data = (subset.pca$x[, 1:2])
subset.kmeans = kmeans(pca.data, centers = 2)

plot(pca.data, cex = .9, main = "PCA and K-means Clustering")
points(pca.data, col = subset.kmeans$cluster, pch = 19)
legend("topright", legend = unique(subset.kmeans$cluster), col = unique(subset.kmeans$cluster), pch = 19)
text(pca.data, labels = label, pos = 4, cex = 0.7)

#use classification method to classify samples, color appropriately
datx = as.data.frame(t(subset.data))
datx = data.frame(label, datx)
dat.lda = lda(label~., datx)
dat.pred = predict(dat.lda, datx)
plot(dat.pred$x, bg = as.numeric(factor(dat.pred$class)), pch = ifelse(label == "h", 17, 22), col = 1, ylab = "Discriminant function", xlab = "Score", main = "Discriminant function for Tumor vs Healthy")
legend("bottomright", legend = c("(Actual) Tumor", "(Actual) Healthy", "(Predicted) Tumor", "(Predicted) Healthy"), pch = c(22, 24, 19, 19), col = c("black", "black", "red", "black"))

#take 5 top discriminant genes (positive and negative), look up gene information for these 10 genes
scaling = dat.lda$scaling
scaling = scaling[order(scaling[,1], decreasing = TRUE), ]
top.positive = head(names(scaling), 5)
top.negative = tail(names(scaling), 5)
list = c(top.positive, top.negative)
list
#[1] "X244291_x_at" "X242910_x_at" "X207783_x_at" "X211694_at"   "X237420_at"   "X1565602_at"  "X221015_s_at"
#[8] "X242613_at"   "X231103_at"   "X201062_at"

#another file with information on all 54,675 genes in the original dataset, used to find information on genes that DAVID could not identify
gene = read.delim("C:/Users/Alex/Desktop/JHU Summer 2024/Gene expression data analysis and visualization/Final project/gene.txt", comment.char = "#", header = TRUE, row.names = 1)
list = gsub("^X", "", list)
extracted = gene[rownames(gene) %in% list,]
extracted$GB_ACC
#[1] "AF085861"  "M81635"    "NM_017627" "AF348076"  "NM_030911" "AI656867"  "AW304248"  "AI809536"  "T89089"    "BE348646" 