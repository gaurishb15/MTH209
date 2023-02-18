df<-USArrests

df <- na.omit(df) ## clean the data by deleting missing values
df <- scale(df) ## The scale() function computes how far each observation locates from the mean.
df
d <- dist(df, method = "euclidean") ## Then we use the “Euclidean” distance as a measurement of how 
## far apart two observations are from each other. We store the distance between each pair of 
## observations in the “USArrests” data by the 50 × 50 matrix using the dist() function:

hc <- hclust(d) ##“hierarchical clustering”
hc

###To make dendrogram
dendro<-plot(hc) 

sub_grp <- cutree(hc, k = 4)
library(factoextra) ##needed for clustering
fviz_cluster(list(data = df, cluster = sub_grp))


