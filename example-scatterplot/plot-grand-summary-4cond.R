dc = read.table('grand_summary_c.txt')
dn = read.table('grand_summary_n.txt')
dg = read.table('grand_summary_g.txt')
dr = read.table('grand_summary_r.txt')
print(dc$V3)
print(dn$V3)
print(dg$V3)
minx = min(c(dc$V2, dn$V2, dg$V2, dr$V2))
maxx = max(c(1.0, dc$V2, dn$V2, dg$V2, dr$V2))
miny = min(c(dc$V3, dn$V3, dg$V3, dr$V3))
maxy = max(c(1.0, dc$V3, dn$V3, dg$V3, dr$V3))
pdf('scatter2.pdf')
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
legend("topright",
        inset=0.01,
        title="Evolutionary Scenario",
        c(  "Co-phylogeny",
            "Non-phylogenetic host jumping",
            "Limited by parasite geographic dispersal",
            "No co-speciation. High parasite extinction"),
        pch=c(1, 17, 19, 0),
        pt.cex=c(2, 2, 2, 2)
        )
dev.off()
