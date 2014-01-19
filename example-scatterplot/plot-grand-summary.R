dc = read.table('grand_summary_c.txt')
dn = read.table('grand_summary_n.txt')
dg = read.table('grand_summary_g.txt')
print(dc$V3)
print(dn$V3)
print(dg$V3)
minx = min(c(dc$V2, dn$V2, dg$V2))
maxx = max(c(1.0, dc$V2, dn$V2, dg$V2))
miny = min(c(dc$V3, dn$V3, dg$V3))
maxy = max(c(1.0, dc$V3, dn$V3, dg$V3))
pdf('scatter.pdf')
plot(dc$V2, dc$V3,
        type="p",
        pch=1,
        #axes=FALSE,
        xlim=c(minx, maxx),
        ylim=c(miny, maxy),
        xlab="Proportion of host tree spanned",
        ylab="Mean prop. of host tree spanned at 20% parasite tree depth")
points(dn$V2, dn$V3, pch=4)
points(dg$V2, dg$V3, pch=19)
legend("topright",
        inset=0.01,
        title="Evolutionary Scenario",
        c(  "Co-phylogeny",
            "Non-phylogenetic host jumping",
            "Limited by parasite geographic dispersal"),
        pch=c(1, 4, 19)
        )
dev.off()