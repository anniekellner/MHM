library(igraph)

# Desktop working directory
#setwd("~/MHM")

# Laptop working directory
#setwd("C:/Users/Clif/Dropbox/MHM")

# Set color palette
col=c("white","grey80","grey40","black")

# Set to output 3 panel figure as PDF
pdf("Fig6.pdf", width=11, height=8.5)
par(mfrow=c(1,3))

# Import edge and weight list for median rates
rates=as.data.frame(read.csv("MHM_BEAST_median_rates.csv"))

# Create igraph network for median rates
rates.net=graph.data.frame(rates,directed=T)
V(rates.net)$color=c(rep(col[4],1),rep(col[1],8))

# View igraph network for median rates
plot.igraph(rates.net,axes=F,
            vertex.size=42,
            edge.curved=c(F, F, F, T, F, F, F, F, T, T, T),
            edge.color=c(col[3],col[2],col[3],col[3],col[2],col[2],col[3],col[2],col[4],col[4],col[2]),
            edge.width=2,
            edge.arrow.size=0.5,
            vertex.label.cex=1.3,
            vertex.label.color=c(rep(col[1],1),rep(col[4],8)),
            frame=F,layout=layout.circle)
legend(x=-0.5,y=-1.5,legend=c("<1","1 to 1.5",">1.5"),cex=1.5,
       lty=1,lwd=2,
       col=col[2:4],
       title="Median rates")

# Import edge and weight list for BF
BF=as.data.frame(read.csv("MHM_BEAST_BF.csv"))

# Create igraph network for BF
BF.net=graph.data.frame(BF,directed=T)
V(BF.net)$color=c(rep(col[4],1),rep(col[1],8))

# View igraph network for BF
plot.igraph(BF.net,axes=F,
            vertex.size=42,
            edge.curved=c(F, F, F, T, F, F, F, F, T, T, T),
            edge.color=c(col[4],col[2],col[3],col[2],col[3],col[2],col[4],col[3],col[4],col[4],col[2]),
            edge.width=2,
            edge.arrow.size=0.5,
            vertex.label.cex=1.3,
            vertex.label.color=c(rep(col[1],1),rep(col[4],8)),
            frame=F,layout=layout.circle)
legend(x=-0.5,y=-1.5,legend=c("<6","6 to 10",">10"),cex=1.5,
       lty=1,lwd=2,
       col=col[2:4],
       title="Bayes factors (2 ln K)")

# Import edge and weight list for median counts
counts=as.data.frame(read.csv("MHM_BEAST_median_counts.csv"))

# Create igraph network for median counts
counts.net=graph.data.frame(counts,directed=T)
V(counts.net)$color=c(rep(col[4],1),rep(col[1],8))

#View igraph network for median counts
plot.igraph(counts.net,axes=F,
            vertex.size=42,
            edge.curved=c(F, F, F, T, F, F, F, F, T, T, T),
            edge.color=c(col[3],col[3],col[3],col[3],col[3],col[3],col[3],col[3],col[4],col[4],col[2]),
            edge.width=2,
            edge.arrow.size=0.5,
            vertex.label.cex=1.3,
            vertex.label.color=c(rep(col[1],1),rep(col[4],8)),
            frame=F,layout=layout.circle)
legend(x=-0.5,y=-1.5,legend=c("0","1 to 3",">3"),cex=1.5,
       lty=1,lwd=2,
       col=col[2:4],
       title="Median counts")
dev.off()
