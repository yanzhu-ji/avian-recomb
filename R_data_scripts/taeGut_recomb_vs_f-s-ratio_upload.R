# updated 4/27/15, with the bug in parsing recomb. table fixed.
recomb.5.col.file="taeGut_new_5-col_chrID.out"
description="zebra finch linkage map (bug fixed)"

# updated 5/18/2015, with the pipeline v.Final 
# re-ran 6/30/2015, 9/5/2015
fl.bed="taeGut_relaxed_LTR-RT_chrID.bed"
solo.bed="taeGut_relaxed_solo-ltr_chrID.bed"

interval=1e6

recomb =read.table(recomb.5.col.file, header=FALSE, sep="\t")
colnames(recomb) = c("chrID", "rate", "position", "start.pos", "end.pos")
recomb$count = 1

chr.list= unique(recomb$chrID)
chr.list=chr.list[chr.list!="chrZ" & chr.list!="chrW" & chr.list != "chrLGE64" & chr.list != "chr1"]
chr.list=factor(chr.list)
chr.list

all.fl = read.table(fl.bed, header=FALSE, sep="\t")
colnames(all.fl) = c("chrID", "start.bed", "end.bed")
all.fl$count = 1
all.fl$middle = (all.fl$start.bed+1 + all.fl$end.bed)/2

all.solo = read.table(solo.bed, header=FALSE, sep="\t")
colnames(all.solo) = c("chrID", "start.bed", "end.bed")
all.solo$count = 1
all.solo$middle = (all.solo$start.bed+ 1 + all.solo$end.bed)/2


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
# 127 (1-Mb); 150 (2-Mb)...
# this probably means there are more NAs in 1-Mb windows than in 2-Mb windows.

### 1. windows grouped by chromosome 

### average recomb. rates
# zebra finch chromosomes
macro.list = c( "chr1", "chr1A", "chr2", "chr3", "chr4", "chr5"  )
micro.list=c("chr4A", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", 
             "chr21", "chr22", "chr23", "chr24", "chr25", "chr26", "chr27", "chr28" )  # ,  "chrZ"
itm.list=c("chr6", "chr7", "chr8", "chr9", "chr10") # itm stands for intermediate

all.a.grouped = data.frame() # a for average recombination rates

for ( chrn in chr.list ){
  chrn.merge = subset(all.merge, all.merge$chrID == chrn)
  chrn.grouped = data.frame(fl.counts=numeric(0), solo.counts=numeric(0), ratio=numeric(0), cM = numeric(0) )
  chrn.fl.counts = sum(chrn.merge$fl.counts)
  chrn.solo.counts = sum(chrn.merge$solo.counts)
  # chrn.cM = mean(na.omit(chrn.merge$cM))
  chrn.cM = sum(na.omit(chrn.merge$cM))
  # chrn.size = max(chrn.merge$middle)
  chrn.grouped[chrn , ] = c(chrn.fl.counts, chrn.solo.counts, chrn.fl.counts/chrn.solo.counts, chrn.cM) # /chrn.size*1e6
  
  all.a.grouped = rbind(all.a.grouped, chrn.grouped)
}

all.a.grouped$size = c(20808668, 21576510, 16952381, 16419078, 14428146, # "chr10" "chr12" "chr13" "chr14" "chr15" 
                       11648728, 11201131, 11587733, 73657157, 156412533, # "chr17" "chr18" "chr19" "chr1A" "chr2"
                       15652063, 5979137, 112617285, 69780378, 62374962, # "chr20" "chr21" "chr3"  "chr4" "chr5" 
                       36305782, 39844632, 27993427, 27241186 # "chr6"  "chr7"  "chr8"  "chr9"
                       )

all.a.grouped$pcM = all.a.grouped$cM/all.a.grouped$size
all.a.grouped$pcM = all.a.grouped$pcM*1e6 # cM/bp -> cM/Mb

# save data to csv
write.table( all.a.grouped, file="taeGut_chr-rr-ratio_all-chr_090415.csv", sep="," )

plot(all.a.grouped$pcM, all.a.grouped$size) 

cor.test(all.a.grouped$pcM, all.a.grouped$ratio, method="spearman")
cor.test(all.a.grouped[-11,]$pcM, all.a.grouped[-11,]$ratio, method="spearman") # chr20 removed
# rho = -0.55, p = 0.018

# plot the mean RR and f-s-ratio in chicken, with points colored/shaped by chromosome groups.
# Figure 3
itm.a.grouped = all.a.grouped[itm.list,]
micro.a.grouped = all.a.grouped[micro.list,]
macro.a.grouped = all.a.grouped[macro.list,]

plot(all.a.grouped$pcM, all.a.grouped$ratio, xlim=c(0,5), cex=0.1, cex.lab=1.2, 
     xlab="mean RR (cM/Mb)", ylab="f-s-ratio")

points(macro.a.grouped$pcM, macro.a.grouped$ratio, pch=18, col="royalblue3", cex=1.2)
points(itm.a.grouped$pcM, itm.a.grouped$ratio, pch=17, col="gray60")
points(micro.a.grouped$pcM, micro.a.grouped$ratio, pch=16, col="indianred1")
points(all.a.grouped[11,]$pcM, all.a.grouped[11,]$ratio, cex=2) # circle the outlier of chr. 20

### 2. windows by 1-Mb or 2-Mb
# 1-Mb did not yield any correlation coefficients:
cor.test(all.merge$cM, all.merge$ratio, method="spearman", exact=FALSE)

# So...
interval=2e6

# re-run Lines 36 - 96

### 3.  permutation test on Spearman Correlation tests, for both chromosomal level and 1 or 2-Mb windows.
## 3.1  by chromosomes; note the using of "pcM" here instead of "cM".
# chromosome level
all.a.grouped.no.na = na.omit(all.a.grouped)

# 3.1.1 with the outlier (chr20, the 11th row)
cor.test(all.a.grouped$pcM, all.a.grouped$ratio, method="spearman", exact=FALSE)

cor.test(all.a.grouped.no.na$pcM, all.a.grouped.no.na$ratio, method="spearman", exact=FALSE)

obs.rho = cor(all.a.grouped.no.na$pcM, all.a.grouped.no.na$ratio, method="spearman")
random.rhos = replicate (2000, {
  pcMstar = sample(all.a.grouped.no.na$pcM, nrow(all.a.grouped.no.na) )
  ratiostar = sample(all.a.grouped.no.na$ratio, nrow(all.a.grouped.no.na) )
  #plot(cMstar, ratiostar)
  return( cor(pcMstar, ratiostar, method="spearman" )) 
})

hist(random.rhos)
mean(random.rhos < obs.rho )
# [1] 0.1165

# 3.1.2 with the outlier removed (still the 11th row after removing NAs)
cor.test(all.a.grouped[-11,]$pcM, all.a.grouped[-11,]$ratio, method="spearman", exact=FALSE)

obs.rho = cor(all.a.grouped.no.na[-11,]$pcM, all.a.grouped.no.na[-11,]$ratio, method="spearman")
random.rhos = replicate (2000, {
  pcMstar = sample(all.a.grouped.no.na[-11,]$pcM, nrow(all.a.grouped.no.na)-1 )
  ratiostar = sample(all.a.grouped.no.na[-11,]$ratio, nrow(all.a.grouped.no.na)-1 )
  #plot(cMstar, ratiostar)
  return( cor(pcMstar, ratiostar, method="spearman" ))
})

mean(random.rhos < obs.rho )
# 0.008

# 3.2 2-Mb windows in zebra finch
all.merge.no.na = na.omit(all.merge)

cor.test(all.merge.no.na$cM, all.merge.no.na$ratio, method="spearman", exact=FALSE)
obs.rho = cor(all.merge.no.na$cM, all.merge.no.na$ratio, method="spearman")
random.rhos = replicate (2000, {
  cMstar = sample(all.merge.no.na$cM, nrow(all.merge.no.na) )
  ratiostar = sample(all.merge.no.na$ratio, nrow(all.merge.no.na) )
  #plot(cMstar, ratiostar)
  return( cor(cMstar, ratiostar, method="spearman" )) 
})

mean(random.rhos < obs.rho )
# 0.52

### 4.   telomere vs. centromere: NA in taeGut as centromeres are not available.

save.image("taeGut_recomb_vs_f-s-ratio.RData") 
