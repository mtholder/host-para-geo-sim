dc = read.table('grand_summary_c.txt')
dn = read.table('grand_summary_n.txt')
dg = read.table('grand_summary_g.txt')
dr = read.table('grand_summary_r.txt')
er = read.table('empirical-stats.txt', stringsAsFactors=FALSE)

print(dc$V5)
print(dn$V5)
print(dg$V5)
r1 = er$V1[1]
r2 = er$V1[2]
r3 = er$V1[3]
minx = min(c(dc$V4, dn$V4, dg$V4, dr$V4))
maxx = max(c(1.0, dc$V4, dn$V4, dg$V4, dr$V4))
miny = min(c(dc$V5, dn$V5, dg$V5, dr$V5))
maxy = max(c(1.0, dc$V5, dn$V5, dg$V5, dr$V5))
pdf('scatter-rev-with-empirical.pdf')
plot(dc$V4, dc$V5,
        type="p",
        pch=1,
        #axes=FALSE,
        xlim=c(minx, maxx),
        ylim=c(miny, maxy),
        xlab="Proportion of parasite tree spanned",
        ylab="Mean prop. of parasite extent at 20% host tree depth",
        cex=2)
points(dn$V4, dn$V5, pch=17, cex=2)
points(dg$V4, dg$V5, pch=19, cex=2)
points(dr$V4, dr$V5, pch=0, cex=2)
points(er$V4[1], er$V5[1], pch=2, cex=2)
points(er$V4[2], er$V5[2], pch=3, cex=2)
points(er$V4[3], er$V5[3], pch=4, cex=2)
legend("topright",
        inset=0.01,
        title="Evolutionary Scenario",
        c(  "Co-phylogeny",
            "Non-phylogenetic host jumping",
            "Limited by parasite geographic dispersal",
            "No co-speciation. High parasite extinction",
            er$V1[[1]],
            er$V1[[2]],
            er$V1[[3]]),
        pch=c(1, 17, 19, 0, 2, 3, 4),
        pt.cex=c(2, 2, 2, 2, 2, 2, 2)
        )
dev.off()
