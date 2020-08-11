

Clust_old1 <- load_object(file = "PCA_analysis/r5/Clusters/Clusters_all_1_old.Rdata")
Clust_old2 <- load_object(file = "PCA_analysis/r5/Clusters/Clusters_all_2_old.Rdata")


Clust_new1<- load_object(file = "PCA_analysis/r5/Clusters/Clusters_all_1.Rdata")
Clust_new2<- load_object(file = "PCA_analysis/r5/Clusters/Clusters_all_2.Rdata")

arandi(Clust_old1$ClustDP1, Clust_new1$ClustDP1)
arandi(Clust_old1$ClustDP, Clust_new1$ClustDP)
arandi(Clust_old1$ClustPY1, Clust_new1$ClustPY1)

arandi(Clust_old2$ClustDP1, Clust_new2$ClustDP1)
arandi(Clust_old2$ClustDP, Clust_new2$ClustDP)
arandi(Clust_old2$ClustPY1, Clust_new2$ClustPY1)



