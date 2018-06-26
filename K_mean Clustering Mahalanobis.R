library(tidyverse)
library(broom)

set.seed(1)

# training set characteristics
centers <- data.frame( grps  = 1:5,
                       gsize = c(1000, 500, 750, 900, 800),
                       v1    = c(  -2,  -1,   0,   1,   2),
                       v2    = c(   0,   5,   1,   2,   4),
                       v3    = c(   1,   4,   3,   5,  -1),
                       v4    = c(   2,  -3,   4,  -1,   1) )

# testing set characteristics
c2      <- data.frame( grps  = 1:6,
                       gsize = c( 100, 100, 100, 100, 100, 100),
                       v1    = c(  -2,  -1,   0,   1,   2,   4),
                       v2    = c(   0,   5,   1,   2,   4,   6),
                       v3    = c(   1,   4,   3,   5,  -1,   7),
                       v4    = c(   2,  -3,   4,  -1,   1,   8) )

# training set
kd      <- centers        %>%
  group_by(grps) %>%
  do(data.frame( v1= rnorm(.$gsize[1], .$v1[1]),
                 v2= rnorm(.$gsize[1], .$v2[1]),
                 v3= rnorm(.$gsize[1], .$v3[1]),
                 v4= rnorm(.$gsize[1], .$v4[1])) ) 

# testing set
kd2     <- c2             %>%
  group_by(grps) %>%
  do(data.frame( v1= rnorm(.$gsize[1], .$v1[1]),
                 v2= rnorm(.$gsize[1], .$v2[1]),
                 v3= rnorm(.$gsize[1], .$v3[1]),
                 v4= rnorm(.$gsize[1], .$v4[1])) ) 

# kmeans for 1 to 10 clusters (training set)
kclust  <- kd                %>%
  crossing(k= 1:10) %>%
  group_by(k)       %>%
  do(clust= kmeans(select(., v1, v2, v3, v4), .$k[1], nstart=5) )

# extract kmeans results information
clusters    <- kclust  %>% tidy(clust)

assignments <- kclust  %>% augment(clust, kd)

clusterings <- kclust  %>% glance(clust)

plot(clusterings$k, clusterings$tot.withinss)

# cluster assignments for 5 cluster case (based on plot)
t           <- kclust$clust[5]
t2          <- t[[1]]
tmp         <- kd
tmp$cluster <- t2$cluster

# compare original and discovered clusters
cTable      <- table(tmp$grps, tmp$cluster)

#pc_correct <- (934+477+707+839+793)/3950

# prep for Mahalanobis distance
variability  <- function(df) {
  df$cluster <- NULL
  n          <- nrow(df)
  avg        <- apply(df, 2, mean)         # column averages
  for (i in seq_along(df) ) {
    df[[i]] <- df[[i]] - avg[i] }        # center data
  dm         <- as.matrix(df)
  sscp       <- t(dm) %*% dm               # matrix multiply
  vcv        <- (1/(n-1)) * sscp           # variance-covariance matrix
  vcvinv     <-solve(vcv)                  # inverse
  return( list(n=n, avg=avg, vcv=vcv, vcvinv=vcvinv) )
}

# calculate and retain summary information for each existing cluster
mhWork   <- tmp                       %>%
  group_by(cluster)         %>%
  do( desc=variability(.) )

clusters <- mhWork$cluster
desc     <- mhWork$desc
clusters
desc[[1]]
mhDf     <- length(desc[[1]]$avg)             # degrees of freedom for chi-sq

d2       <- matrix( -1,
                    nrow=nrow(kd2),
                    ncol=length(clusters) )

# collect Mahalanobis distance for test data across clusters
for ( i in seq_along(clusters) ) {
  wTst <- kd2
  t    <- desc[[i]]
  for ( j in seq_along(kd2) ) {
    wTst[[j]] <- wTst[[j]] - t$avg[j] }   # center the data wrt this cluster
  for ( j in 1:nrow(wTst) ) {
    tr <- as.matrix(wTst[ j, ])
    d2[j, i] <- tr %*% t$vcvinv %*% t(tr) # Mahalanobis distance squared
  }
}

# which cluster for each testing observation?
newClust      <- apply(d2, 1, which.min)
kd2$cluster   <- newClust
(cTable2      <- table(kd2$grps, kd2$cluster) )

# chi-squared cdf value for each testing observation
minStat       <- apply(d2, 1, min)
chiSq         <- dchisq(minStat, mhDf)
summary(chiSq)

kd2$minStat   <- minStat
kd2$chiSqProb <- chiSq

boxplot(chiSqProb~cluster,data=kd2)
# diagnostics via hist()
chk <- function(df, clust, bins) {
  t <- df[df$cluster == clust, ]
  hist(t$chiSqProb, breaks= bins)
}
# capture histogram for each cluster
chk(kd2, 1, 20)
chk(kd2, 2, 20)
chk(kd2, 3, 20)
chk(kd2, 4, 20)

# summarize each cluster (chiSqProb values only)
t <- kd2[kd2$cluster == 4, ]
summary(t$chiSqProb)