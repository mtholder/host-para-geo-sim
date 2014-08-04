dc = read.table('grand_summary_c.txt')
dn = read.table('grand_summary_n.txt')
dg = read.table('grand_summary_g.txt')
dr = read.table('grand_summary_r.txt')
er = read.table('empirical-stats.txt', stringsAsFactors=FALSE)
library(scatterplot3d)
print(dc$V5)
print(dn$V5)
print(dg$V5)
minx = min(c(dc$V3, dn$V3, dg$V3, dr$V3))
maxx = max(c(1.0, dc$V3, dn$V3, dg$V3, dr$V3))
miny = min(c(dc$V5, dn$V5, dg$V5, dr$V5))
maxy = max(c(1.0, dc$V5, dn$V5, dg$V5, dr$V5))
minz = min(c(dc$V2, dn$V2, dg$V2, dr$V2))
maxz = max(c(1.0, dc$V2, dn$V2, dg$V2, dr$V2))
pdf('scatter-3D-with-empirical.pdf')
s3d <- scatterplot3d(dc$V3, dc$V5, dc$V2,
        type="p",
        pch=1,
        axes=FALSE,
        xlim=c(minx, maxx),
        ylim=c(miny, maxy),
        zlim=c(minz, maxz),
        xlab="Mean prop. of parasite extent at 20% host tree depth",
        ylab="Mean prop. of parasite extent at 20% host tree depth",
        zlab="prop. host tree spanned"
        #cex=2
        )
s3d$points3d(dn$V3, dn$V5, dn$V2, pch=17)
s3d$points3d(dg$V3, dg$V5, dn$V2, pch=19)
s3d$points3d(dr$V3, dr$V5, dr$V2, pch=0)
s3d$points3d(er$V3[1], er$V5[1], er$V2[1], pch=2)
s3d$points3d(er$V3[2], er$V5[2], er$V2[2], pch=3)
s3d$points3d(er$V3[3], er$V5[3], er$V2[3], pch=4)
#legend3d("topright",
#        inset=0.01,
#        title="Evolutionary Scenario",
#        c(  "Co-phylogeny",
#            "Non-phylogenetic host jumping",
#            "Limited by parasite geographic dispersal",
#            "No co-speciation. High parasite extinction",
#            er$V2[[1]],
#            er$V2[[2]],
#            er$V2[[3]]),
#        pch=c(1, 17, 19, 0, 2, 3, 4),
#        pt.cex=c(2, 2, 2, 2, 2, 2, 2)
#        )
dev.off()
