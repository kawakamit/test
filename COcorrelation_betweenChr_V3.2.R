# TODO: Add comment
# 
# compare the CO density per chromosome bertween chromosomes to see if high rec individuals have high rec rate in all chromosomes
# version 2: calculate colony mean CO rate
# version 2.1: remove chromosomes where no CO was observed. By assuming at least one CO per chromosome, zero CO is like an artifact. Alternatively, one can add 1 CO in such cases
# version 2.2: cM/Mb = 0 entry was changed as cM/Mb = 1. By assuming at least one CO per chromosome, zero CO is like an artifact. 
# version 3: This version is fundamentally different from the previous version. 
#				This version is to simulate null distribution of correlation coefficient (r) between chromosomes of different individuals
#				By permuting individuals, you break-up between-chromosome correlation within individuals.
#				If between-chromosome correlation within individuals is strong, individuals with higher recombination rate should have all chromosomes with higher recombination rate. 
#				This indicates that trans-factors affect genome-wide recombination rate, resulting in higher rec rate in all chromosomes.
# version 3.1: version 3 uses pairwise comparison berween chromosomes within species for estimating observed value.
#				This has very limited power because of the small number of pairwise comparisons within colony
#				Therefore, to estimate the observed value, the genome was randomly subdivided into 2 groups (with permutation) while limiting the pair within species
#				This will be compared with the same permuted pairwise comparison BETWEEN colony
# version 3.2: family mean version of version 3.1. In addition, rather than subdividing the whole genome into two, just pick up one chromosome
# 
# Author: takikawakami
###############################################################################

#Clear the workspace
rm(list=ls())

#library(stringr)
#options(scipen=999)

#comarg    <- commandArgs()   
#co.file <- comarg[6]              # co.file <- "European.COcount.perChr.perInd.PB_v2.out"  
#size.file <- comarg[7]              # size.file <- "chrom.size.PB_v2.txt"
#fam.file <- comarg[8]				# fam.file <- "family.ind.ID.RG.list.euro"

setwd("/Users/takikawakami/Documents/Bee_Project/PBnew_Assembly/version2/CO_count_PB_v2/CO_position_euro")

shuffle.within.func <- function (co.file, size.file, nrep) {
	df.co <- read.table(co.file, header=FALSE)
	colnames(df.co) <- c("chrom", "CO_ct", "fam", "ind", "id")
	
	df.size <- read.table(size.file, header=FALSE)
	chrom.list <- as.character(df.size$V1)
	
	# ind.list <- as.character(unique(df.co$ind))
	
	# read fam file... This is needed to independently count the number of drones/colony
	#df.fam <- read.table(fam.file, header=FALSE)
	#colnames(df.fam) <- c("fam", "ind", "id", "rg")
	
	fam.list <- as.character(unique(df.co$fam))
	
	# density of CO (the number of COs/100 Mb)
	CO_dens <- vector(length=nrow(df.co))
	
	# calculate cM/Mb
	for (i in 1:nrow(df.co)) {
		chr.temp <- as.character(df.co$chrom[i])
		chr.size.temp <- subset(df.size, V1==chr.temp)$V2
		
		# If at least one CO is observed, calculate CO density
		if(df.co$CO_ct[i]>0) {
			CO_dens[i] <- df.co$CO_ct[i]/chr.size.temp*100000000
			
		} else {
			# if no CO is observed, assume that there is at least one CO per chromosome.
			CO_dens[i] <- 1/chr.size.temp*100000000			
		}

	}
	
	# add to df
	df.co$CO_dens <- CO_dens
	
	
	### mean CO rate per chromosome per family (by averagin all drones within family)
	chr.fam.out <- c()
	fam.fam.out <- c()
	co.fam.out <- c()
	
	for (fam in 1:length(fam.list)) {
		fam.temp <- fam.list[fam]	
		df.co.fam <- subset(df.co, fam==fam.temp)
		
		for (c in 1:length(chrom.list)) {
			chr.temp <- chrom.list[c]
			df.co.fam.chr <- subset(df.co.fam, chrom==chr.temp)
			co.fam.chr.mean <- mean(df.co.fam.chr$CO_dens)
			
			# store results
			chr.fam.out <- c(chr.fam.out, chr.temp)
			fam.fam.out <- c(fam.fam.out, fam.temp)
			co.fam.out <- c(co.fam.out, co.fam.chr.mean)
		}
	}
	
	# put results in df
	df.co.fam.out <- data.frame(chrom=chr.fam.out, fam=fam.fam.out, CO_dens=co.fam.out)
	

	##################
	# list files
	nchrom <- length(chrom.list)
	nchrom.top <- round(nchrom/2)
	nchrom.bot <- nchrom - nchrom.top
	##################
	
	
	#### random pick of one pair of chromosomes	
	ct <- 1
	
	# store permutation results
	cor.out <- vector(length=nrep)
	df.shuf.perm <- c()
	
	while (ct<=nrep) {
		# store shuffled result per permutation
		fam.out <- vector(length=length(fam.list))
		chr_a <- vector(length=length(fam.list))
		chr_b <- vector(length=length(fam.list))
		CO_a <- vector(length=length(fam.list))
		CO_b <- vector(length=length(fam.list))
		
		for (f in 1:length(fam.list)) {
			fam.temp <- fam.list[f]
			df.co.subFam <- subset(df.co.fam.out, fam==fam.temp)
			
			# shuffle
			df.co.subFam.shuf <- df.co.subFam[sample(nrow(df.co.subFam), replace = FALSE),]
			
			# pick a random pair of chromosomes
			fam.out[f] <- fam.temp
			chr_a[f] <- as.character(df.co.subFam.shuf$chrom[1])
			chr_b[f] <- as.character(df.co.subFam.shuf$chrom[2])
			CO_a[f] <- df.co.subFam.shuf$CO_dens[1]
			CO_b[f] <- df.co.subFam.shuf$CO_dens[2]			
		}
		
		# merge result as df
		df.shuf <- data.frame(fam=fam.out, chr_a=chr_a, chr_b=chr_b, CO_a=CO_a, CO_b=CO_b)
		df.shuf.perm <- rbind(df.shuf.perm, df.shuf)
		
		# correlation coefficient
		cor.out[ct] <- cor(df.shuf$CO_a, df.shuf$CO_b)
		
		# counter
		ct <- ct+1
	}
	
	# output
	out.list <- list("df.shuf.perm" = df.shuf.perm, "cor.out" = cor.out)
	
	return(out.list)
}


# generate observed pairwise comparison between 2 groups of chromosomes within individuals
# output is a list containing permuted combinations of chromosomes, and correlation coefficient for each permutation

df.out.mel <- shuffle.within.func("European.COcount.perChr.perInd.PB_v2.out", "chrom.size.PB_v2.txt", 1000)
df.out.scu <- shuffle.within.func("scu.COcount.perChr.perInd.PB_v2.out", "chrom.size.PB_v2.txt", 1000)
df.out.cap <- shuffle.within.func("cap.COcount.perChr.perInd.PB_v2.out", "chrom.size.PB_v2.txt", 1000)



#####################################################
# shuffle between individuals of the different family
#####################################################



shuffle.between.func <- function (co.file, size.file, nrep) {
	df.co <- read.table(co.file, header=FALSE)
	colnames(df.co) <- c("chrom", "CO_ct", "fam", "ind", "id")
	
	df.size <- read.table(size.file, header=FALSE)
	chrom.list <- as.character(df.size$V1)
	
	# ind.list <- as.character(unique(df.co$ind))
	
	# read fam file... This is needed to independently count the number of drones/colony
	#df.fam <- read.table(fam.file, header=FALSE)
	#colnames(df.fam) <- c("fam", "ind", "id", "rg")
	
	fam.list <- as.character(unique(df.co$fam))
	
	# density of CO (the number of COs/100 Mb)
	CO_dens <- vector(length=nrow(df.co))
	
	# calculate cM/Mb
	for (i in 1:nrow(df.co)) {
		chr.temp <- as.character(df.co$chrom[i])
		chr.size.temp <- subset(df.size, V1==chr.temp)$V2
		
		# If at least one CO is observed, calculate CO density
		if(df.co$CO_ct[i]>0) {
			CO_dens[i] <- df.co$CO_ct[i]/chr.size.temp*100000000
			
		} else {
			# if no CO is observed, assume that there is at least one CO per chromosome.
			CO_dens[i] <- 1/chr.size.temp*100000000			
		}
		
	}
	
	# add to df
	df.co$CO_dens <- CO_dens
	
	
	### mean CO rate per chromosome per family (by averagin all drones within family)
	chr.fam.out <- c()
	fam.fam.out <- c()
	co.fam.out <- c()
	
	for (fam in 1:length(fam.list)) {
		fam.temp <- fam.list[fam]	
		df.co.fam <- subset(df.co, fam==fam.temp)
		
		for (c in 1:length(chrom.list)) {
			chr.temp <- chrom.list[c]
			df.co.fam.chr <- subset(df.co.fam, chrom==chr.temp)
			co.fam.chr.mean <- mean(df.co.fam.chr$CO_dens)
			
			# store results
			chr.fam.out <- c(chr.fam.out, chr.temp)
			fam.fam.out <- c(fam.fam.out, fam.temp)
			co.fam.out <- c(co.fam.out, co.fam.chr.mean)
		}
	}
	
	# put results in df
	df.co.fam.out <- data.frame(chrom=chr.fam.out, fam=fam.fam.out, CO_dens=co.fam.out)
	
	
	# add a chromosome index
	chrom.list.temp <- as.character(df.co.fam.out$chrom)
	df.co.fam.out$index <- as.integer(gsub("Group", "" , chrom.list.temp))
	
	# set up a chromosome index reference
	chrom.idx.ref <- seq(1:16)
	
	##################
	# list files
	nchrom <- length(chrom.list)
	nchrom.top <- round(nchrom/2)
	nchrom.bot <- nchrom - nchrom.top
	##################
	
	
	#### random pick of one pair of chromosomes	
	ct <- 1
	
	
	# store shuffled result per permutation
	fam_a.out <- vector(length=length(nrep*length(fam.list)))
	chr_a <- vector(length=length(nrep*length(fam.list)))
	CO_a <- vector(length=length(nrep*length(fam.list)))
	fam_b.out <- vector(length=length(nrep*length(fam.list)))
	chr_b <- vector(length=length(nrep*length(fam.list)))
	CO_b <- vector(length=length(nrep*length(fam.list)))
	
	while (ct<=nrep*length(fam.list)) {
		
		# shuffle
		df.co.fam.out.shuf <- df.co.fam.out[sample(nrow(df.co.fam.out), replace = FALSE),]
		
		# select two entries whose chromosome and family are different
		entry_a <- head(df.co.fam.out.shuf, n=1)
		chrom_a <- as.character(entry_a$chrom)
		family_a <- as.character(entry_a$fam)
		
		fam_a.out[ct] <- family_a
		chr_a[ct] <- chrom_a
		CO_a[ct] <- entry_a$CO_dens
		
		for (r in 2:nrow(df.co.fam.out.shuf)) {
			entry_b <- df.co.fam.out.shuf[r,]
			chrom_b <- as.character(entry_b$chrom)
			family_b <- as.character(entry_b$fam)
			
			if (chrom_a!=chrom_b & family_a!=family_b) {
				fam_b.out[ct] <- family_b
				chr_b[ct] <- chrom_b
				CO_b[ct] <- entry_b$CO_dens
				break	
			}
		}
				
		# counter
		ct <- ct+1
	}
	
	# put results in df
	df.other.out <- data.frame(fam_a=fam_a.out, chr_a=chr_a, CO_a=CO_a, fam_b=fam_b.out, chr_b=chr_b, CO_b=CO_b)
	
	# calculate correlation coefficient for n number of entries (n = the number of families)
	nfam <- length(fam.list)
	
	# store permutation results
	cor.out <- c()
	kt <- 3
	
	
	while(kt<=nrow(df.other.out)) {
		start <- kt
		df.sub <- tail(head(df.other.out, n=kt), n=3)
		cor.out <- c(cor.out, cor(df.sub$CO_a, df.sub$CO_b))
		
		kt <- kt+3
	}
	
	
	# output
	out.list <- list("df.shuf.other" = df.other.out, "cor.out" = cor.out)
	
	return(out.list)
}


# generate observed pairwise comparison between 2 groups of chromosomes within individuals
# output is a list containing permuted combinations of chromosomes, and correlation coefficient for each permutation
# NOTE that There may some NA results in correlation analysis, which must be removed

df.out.mel.bet <- shuffle.between.func("European.COcount.perChr.perInd.PB_v2.out", "chrom.size.PB_v2.txt", 1000)
df.out.scu.bet <- shuffle.between.func("scu.COcount.perChr.perInd.PB_v2.out", "chrom.size.PB_v2.txt", 1000)
df.out.cap.bet <- shuffle.between.func("cap.COcount.perChr.perInd.PB_v2.out", "chrom.size.PB_v2.txt", 1000)





# plot.. mellifera
xlim <- c(min(df.out.mel.bet$df.shuf.other$CO_a, df.out.mel.bet$df.shuf.other$CO_b,df.out.mel$df.shuf.perm$CO_a, df.out.mel$df.shuf.perm$CO_b), max(df.out.mel.bet$df.shuf.other$CO_a, df.out.mel.bet$df.shuf.other$CO_b,df.out.mel$df.shuf.perm$CO_a, df.out.mel$df.shuf.perm$CO_b))
ylim <- xlim

plot(df.out.mel.bet$df.shuf.other$CO_a, df.out.mel.bet$df.shuf.other$CO_b, col="lightgrey", cex=0.2, pch=20, axes=FALSE, xlab="", ylab="", xlim=xlim, ylim=xlim)
lm.null.mel <- lm(df.out.mel.bet$df.shuf.other$CO_b~df.out.mel.bet$df.shuf.other$CO_a)

par(new=T)

plot(df.out.mel$df.shuf.perm$CO_a, df.out.mel$df.shuf.perm$CO_b, col="blue", cex=0.2, pch=20, xlab="Recombination rate on chromosome i (cM/Mb)", ylab="Recombination rate on chromosome j (cM/Mb)", xlim=xlim, ylim=xlim)
lm.obs.mel <- lm(df.out.mel$df.shuf.perm$CO_b~df.out.mel$df.shuf.perm$CO_a)

abline(lm.null.mel, col='darkgrey', lwd=2)
abline(lm.obs.mel, col='darkblue', lwd=2)

# test significance ... mel
co.obs <- df.out.mel$cor.out
co.null <- df.out.mel.bet$cor.out







