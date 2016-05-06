
# 4.27 fixed a major bug in the original perl script producing the 5-col file
recomb.5.col.file="galGal4_new_5-col.out"
description="Chicken linkage map (major bug fixed, Elferink et al. 2010)"

# 5.9.2015, updated data sets after manual correction
fl.bed="galGal4_relaxed_LTR-RT_chrID.bed" 
solo.bed="galGal4_relaxed_solo-ltr_chrID.bed" 

interval=1e6

recomb =read.table(recomb.5.col.file, header=FALSE, sep="\t")
colnames(recomb) = c("chrID", "rate", "position", "start.pos", "end.pos")
recomb$count = 1

chr.list= unique(recomb$chrID)
chr.list=chr.list[chr.list!="chrZ" & chr.list!="chrW" & chr.list != "chrLGE64"]
chr.list=factor(chr.list)
chr.list

all.fl = read.table(fl.bed, header=FALSE, sep="\t")
colnames(all.fl) = c("chrID", "start.bed", "end.bed")
all.fl$count = 1
all.fl$middle = (all.fl$start.bed+1 + all.fl$end.bed)/2
nrow(all.fl)

all.solo = read.table(solo.bed, header=FALSE, sep="\t")
colnames(all.solo) = c("chrID", "start.bed", "end.bed")
all.solo$count = 1
all.solo$middle = (all.solo$start.bed+ 1 + all.solo$end.bed)/2
nrow(all.solo)

all.d.recomb = data.frame()
all.c.fl = data.frame()
all.c.solo = data.frame()

for ( chrn in chr.list ){
  #  chrn = "chr1"
  chrn.recomb = subset (recomb, recomb$chrID == chrn )
  pos.cut = cut (chrn.recomb$position, seq(1, max(chrn.recomb$position)+ interval, by = interval ) )
  chrn.d.recomb = data.frame(as.data.frame( tapply (chrn.recomb$rate, list(pos.cut), sum ) ), 
                             as.data.frame( tapply (chrn.recomb$count, list(pos.cut), sum) ) )
  
  colnames(chrn.d.recomb) = c("sum", "marker.counts")
  chrn.d.recomb$cM = 1e6/interval*chrn.d.recomb$sum
  
  chrn.d.recomb$chrID = chrn
  chrn.d.recomb$lower = as.numeric( sub("[(|[](.+),.*", "\\1", levels(pos.cut) ) )  
  chrn.d.recomb$upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", levels(pos.cut) ) )
  chrn.d.recomb$middle = (chrn.d.recomb$upper + chrn.d.recomb$lower)/2
  chrn.d.recomb$key = paste(chrn.d.recomb$chrID, chrn.d.recomb$middle, sep="_" )
  
  head(chrn.d.recomb)
  
  all.d.recomb = rbind(all.d.recomb, chrn.d.recomb)
  
  ## full-length
  chrn.fl = subset (all.fl, all.fl$chrID == chrn )
  fl.pos.cut = cut(chrn.fl$middle, seq(1, max(chrn.recomb$position) + interval, by = interval) )
  chrn.c.fl = data.frame( as.data.frame( tapply( chrn.fl$count, list(fl.pos.cut), sum ) ) ) 
  colnames(chrn.c.fl) = c("fl.counts")
  chrn.c.fl$chrID = chrn
  
  chrn.c.fl[is.na(chrn.c.fl)] = 0
  chrn.c.fl$lower = as.numeric( sub("[(|[](.+),.*", "\\1", levels(fl.pos.cut) ) )
  chrn.c.fl$upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", levels(fl.pos.cut) ) )
  chrn.c.fl$middle = (chrn.c.fl$upper + chrn.c.fl$lower)/2
  chrn.c.fl$key = paste(chrn.c.fl$chrID, chrn.c.fl$middle, sep="_" )
  
  all.c.fl = rbind(all.c.fl, chrn.c.fl)
  
  ### solo
  chrn.solo = subset (all.solo, all.solo$chrID == chrn )
  solo.pos.cut = cut(chrn.solo$middle, seq(1, max(chrn.recomb$position) + interval, by = interval) )
  chrn.c.solo = as.data.frame( tapply( chrn.solo$count, list(solo.pos.cut), sum ) )
  colnames(chrn.c.solo) = "solo.counts"
  chrn.c.solo$chrID = chrn
  
  chrn.c.solo[is.na(chrn.c.solo)] = 0
  chrn.c.solo$lower = as.numeric( sub("[(|[](.+),.*", "\\1", levels(solo.pos.cut) ) )
  chrn.c.solo$upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", levels(solo.pos.cut) ) )
  chrn.c.solo$middle = (chrn.c.solo$upper + chrn.c.solo$lower)/2
  chrn.c.solo$key = paste(chrn.c.solo$chrID, chrn.c.solo$middle, sep="_" )
  
  all.c.solo = rbind(all.c.solo, chrn.c.solo)
}

fl.solo = merge(all.c.fl, all.c.solo, key = "key")
fl.solo$ratio = fl.solo$fl.counts/(fl.solo$solo.counts)
all.merge = merge(fl.solo, all.d.recomb, key="key")
all.merge$total = all.merge$solo.counts + all.merge$fl.counts  #added 2.27.2015
nrow(na.omit(all.merge))
# 345 (1-Mb); 259 (2-Mb)

### 1. windows grouped by chromosome 
### size-averaged recombination rates

all.a.grouped = data.frame() # a for average recombination rates

for ( chrn in chr.list ){
  chrn.merge = subset(all.merge, all.merge$chrID == chrn)
  chrn.grouped = data.frame(fl.counts=numeric(0), solo.counts=numeric(0), ratio=numeric(0), cM = numeric(0) )
  chrn.fl.counts = sum(chrn.merge$fl.counts)
  chrn.solo.counts = sum(chrn.merge$solo.counts)
  chrn.cM = sum(na.omit(chrn.merge$cM))
  chrn.grouped[chrn , ] = c(chrn.fl.counts, chrn.solo.counts, chrn.fl.counts/chrn.solo.counts, chrn.cM) # /chrn.size*1e6
  
  all.a.grouped = rbind(all.a.grouped, chrn.grouped)
}

all.a.grouped$size = as.numeric ( c("195276750", "148809762", "110447801", "90216835", "59580361", 
                       "34951654", "36245040", "28767244", "23441680", "19911089", 
                       "19401079", "19897011", "17760035", "15161805", "12656803", 
                       "535270", "10454150", "11219875", "9983394", "14302601", 
                       "6802778", "4081097", "5723239", "6323281", "2191139", 
                       "5329985", "5209285", "4742627", "965146" 
                       ) )

all.a.grouped$pcM = all.a.grouped$cM/all.a.grouped$size
all.a.grouped$pcM = all.a.grouped$pcM*1e6

# save data to csv
write.table( all.a.grouped, file="galGal4_chr-rr-ratio_all-chr_090415.csv", sep="," )

# check if total recomb. rates and chromosome sizes are positively correlated
plot(all.a.grouped$pcM, all.a.grouped$size) 

plot(all.a.grouped[-29,]$pcM, all.a.grouped[-29,]$ratio, xlim=c(0, 30), xlab="mean RR (cM/Mb)", ylab="f-s-ratio") # xlim=c(0,0.00003),

cor.test(all.a.grouped$pcM, all.a.grouped$ratio, method="spearman")
cor.test(all.a.grouped[-16,]$pcM, all.a.grouped[-16,]$ratio, method="spearman") # without chr. 16
cor.test(all.a.grouped[-29,]$pcM, all.a.grouped[-29,]$ratio, method="spearman") # without chrLGE22C19W28_E50C23

# plot the mean RR and f-s-ratio in chicken, with points colored/shaped by chromosome groups.
# Figure 3
macro.list = c( "chr1", "chr2", "chr3", "chr4", "chr5"  )
micro.list=c("chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", 
             "chr21", "chr22", "chr23", "chr24", "chr25", "chr26", "chr27", "chr28" )  # ,  "chrZ"
itm.list=c("chr6", "chr7", "chr8", "chr9", "chr10") # itm stands for intermediate

itm.a.grouped = all.a.grouped[itm.list,]
micro.a.grouped = all.a.grouped[micro.list,]
macro.a.grouped = all.a.grouped[macro.list,]

plot(all.a.grouped[-29,]$pcM, all.a.grouped[-29,]$ratio, xlim=c(0, 30), cex.lab=1.2,
     xlab="mean RR (cM/Mb)", ylab="f-s-ratio", cex=0.2) 

points(macro.a.grouped$pcM, macro.a.grouped$ratio, pch=18, col="royalblue3", cex=1.2)
points(itm.a.grouped$pcM, itm.a.grouped$ratio, pch=17, col="gray60")
points(micro.a.grouped$pcM, micro.a.grouped$ratio, pch=16, col="indianred1")

legend(13.5, 1.0, c("macrochromosome", "intermediate-sized", "microchromosome"), col=c("royalblue3", "gray60", "indianred1"), 
       pch=c(18, 17, 16) )

### 2.  1-Mb windows
cor.test(all.merge$cM, all.merge$ratio, use="complete", method="spearman")

### 3.  permutation test on Spearman Correlation tests, for both chromosomal level and 1 or 2-Mb windows.
## 3.1  by chromosomes; note the using of "pcM" here instead of "cM".
all.a.grouped.no.na = na.omit(all.a.grouped)

obs.rho = cor(all.a.grouped.no.na$pcM, all.a.grouped.no.na$ratio, method="spearman")
random.rhos = replicate (2000, {
  pcMstar = sample(all.a.grouped.no.na$pcM, nrow(all.a.grouped.no.na) )
  ratiostar = sample(all.a.grouped.no.na$ratio, nrow(all.a.grouped.no.na) )
  #plot(cMstar, ratiostar)
  return( cor(pcMstar, ratiostar, method="spearman" )) 
})

hist(random.rhos)
mean(random.rhos < obs.rho )
# 0.005; value flucturates each time.

##  3.2  1-Mb windows
all.merge.no.na = na.omit(all.merge)

cor.test(all.merge.no.na$cM, all.merge.no.na$ratio, method="spearman", exact=FALSE)
obs.rho = cor(all.merge.no.na$cM, all.merge.no.na$ratio, method="spearman")
random.rhos = replicate (2000, {
  cMstar = sample(all.merge.no.na$cM, nrow(all.merge.no.na) )
  ratiostar = sample(all.merge.no.na$ratio, nrow(all.merge.no.na) )
  #plot(cMstar, ratiostar)
  return( cor(cMstar, ratiostar, method="spearman" )) 
})

hist(random.rhos)
mean(random.rhos < obs.rho )
# 0.0055

## 3.3  2-Mb windows: first run the main script with interval=2e6, then re-run the above section.  
#                     P value: [1] 0.071

### 4.   telomere vs. centromere
##  4.1. centromere regions are defined by 1 Mb flanking the middle point of each centromere

cent.data = read.table("galGal4_centromere.bed")
colnames(cent.data) = c("chrID", "start", "end")

cent.data$middle = (cent.data$start + cent.data$end)/2

cent.radius = 1e6

cent.data$cstart = cent.data$middle - cent.radius
cent.data$cend = cent.data$middle + cent.radius

cent = data.frame()
for ( chrn in macro.list ){
  chrn.merge = subset(all.merge, all.merge$chrID == chrn)
  chrn.cent = subset(cent.data, cent.data$chrID == chrn)
  chrn.selected = chrn.merge[chrn.merge$middle < chrn.cent$cend & chrn.merge$middle > chrn.cent$cstart,]
  cent = rbind(chrn.selected, cent)
}

cent$loc = "cent"

##  4.2. telomeric regions are defined by x Mb from each end of each chromosome.
telo.radius = 2e6

telo1 = data.frame()
for (chrn in macro.list) {
  chrn.merge = subset(all.merge, all.merge$chrID == chrn)
  chrn.telo1 = chrn.merge[chrn.merge$middle < telo.radius, ]
  telo1 = rbind(chrn.telo1, telo1)
}
telo1$loc = "telo1"

telo2 = data.frame()
for ( chrn in macro.list ){
  chrn.merge = subset(all.merge, all.merge$chrID == chrn)
  chrn.telo2 = chrn.merge[chrn.merge$middle > ( max(chrn.merge$middle) - telo.radius ), ]
  telo2 = rbind(chrn.telo2, telo2)
}

telo2$loc = "telo2"

telo = rbind(telo1, telo2)

boxplot(na.omit(telo1$ratio), na.omit(telo2$ratio, cent), cent$ratio) #

boxplot(na.omit(telo$cM), cent$cM) 

boxplot(telo$ratio, telo$solo.counts, telo$fl.counts, cent$ratio, cent$solo.counts, cent$fl.counts,
        at=c(1,2,3, 5,6,7)) # just for fun...

# Figure 4
boxplot(na.omit(telo$ratio), na.omit(cent$ratio),
        ylab="f-s-ratio", cex.lab=1.2,
        names=c("telomere", "centromere"),
        col="gray") 

wilcox.test(telo$cM, cent$cM, alternative="greater") # p = 0.013...
wilcox.test(na.omit(telo$ratio), na.omit(cent$ratio), alternative="less") # p = 0.006... 9/17/2015

write.table( telo, file="galGal4_telo.csv", sep="," )
write.table( cent, file="galGal4_centro.csv", sep="," )

save.image("galGal4_recomb_vs_f-s-ratio.RData")
