dc = read.table('grand_summary_c.txt')
dn = read.table('grand_summary_n.txt')
dg = read.table('grand_summary_g.txt')
dr = read.table('grand_summary_r.txt')
er = read.table('empirical-stats.txt', stringsAsFactors=FALSE)

print(dc$V3)
print(dn$V3)
print(dg$V3)
r1 = er$V1[1]
r2 = er$V1[2]
r3 = er$V1[3]
print(r1)
print(r2)
print(r3)
print(rownames(er$V1))
minx = min(c(dc$V2, dn$V2, dg$V2, dr$V2))
maxx = max(c(1.0, dc$V2, dn$V2, dg$V2, dr$V2))
miny = min(c(dc$V3, dn$V3, dg$V3, dr$V3))
maxy = max(c(1.0, dc$V3, dn$V3, dg$V3, dr$V3))
pdf('scatter-with-empirical.pdf')
plot(dc$V2, dc$V3,
        type="p",
        pch=1,
        #axes=FALSE,
        xlim=c(minx, maxx),
        ylim=c(miny, maxy),
        xlab="Proportion of host tree spanned",
        ylab="Mean prop. of host extent at 20% parasite tree depth",
        cex=2)
points(dn$V2, dn$V3, pch=17, cex=2)
points(dg$V2, dg$V3, pch=19, cex=2)
points(dr$V2, dr$V3, pch=0, cex=2)
points(er$V2[1], er$V3[1], pch=2, cex=2)
points(er$V2[2], er$V3[2], pch=3, cex=2)
points(er$V2[3], er$V3[3], pch=4, cex=2)
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
