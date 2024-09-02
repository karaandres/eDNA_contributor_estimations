### The goal of this code is to simulate mixtures of DNA to determine
### how many individuals can be estimated from the contributor estimation
### using real goby genotypes
### Last updated 1.8.2023 by Kara Andres (akara@wustl.edu)

### Clear work environment and load packages
rm(list = ls())

library(adegenet)
library(pegas)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(splitstackshape)
library(vegan)
library(data.table)
library(devtools)
library(stringr)
library(reshape)
library(RcppAlgos)
library(gtools)
library(scales)

#################################################################################
################################## Functions ####################################
#################################################################################
# Diploid contributor estimation: calculation following Weir
# Suresh's code -- speed up computation time when # of alleles is large
MixtureLikelihood_Diploid <- function(x,p.v){
  require(RcppAlgos)
  if(length(p.v)>(2*x)) {sum.p <- 0} 
  else {
    counter <- 0
    if (length(p.v)>0){
      sum.p <- sum(p.v)^(2*x)
      if(length(p.v)>1){ # only conduct if > 1 allele is observed
        for(i in 1:(length(p.v)-1)){
          counter <- counter+1
          temp.combo <- comboGeneral(v=length(p.v), m=length(p.v)-i, repetition=FALSE, Parallel=TRUE, nThreads=10) #; i; dim(temp.combo);i<-i+1
          if(nrow(temp.combo)>15e6){
            ix <- 999999		
            row.sums.pow2x <- 1:(floor(nrow(temp.combo)/ix)+1)
            for(m in 1:length(row.sums.pow2x)){
              temp.m <- temp.combo[ ((m-1)*ix+1) : min(m*ix,nrow(temp.combo)), ]
              temp.combo.v <- c(temp.m)
              p.v.m <- matrix(p.v[temp.combo.v],nr=dim(temp.m)[1],nc=dim(temp.m)[2],byrow=F)
              row.sums.pow2x[m] <- sum(rowSums(p.v.m)^(2*x)) # reference as follows
            } # end m loop
            sum.p <- sum.p+((-1)^counter)*sum(row.sums.pow2x) } else 
            {
              # non-parsed
              temp.combo.v <- c(temp.combo)
              p.v.m <- matrix(p.v[temp.combo.v],nr=dim(temp.combo)[1],nc=dim(temp.combo)[2],byrow=F)
              row.sums.pow2x <- rowSums(p.v.m)^(2*x) # reference as follows
              sum.p <- sum.p+((-1)^counter)*sum(row.sums.pow2x)
            } # end parsing option
        }# end i loop
      } # end if
    } else sum.p <- NA # if no alleles observed, assign NA
  }
  return(sum.p)
} # end function

# loci=loci.list
# y=150

# Multi-locus diploid function across 1:y number of putative contributors
MixtureLikelihood2 <- function(y, loci){
  Likelihood_df <- data.frame()
  for (i in 1:length(loci)){ # for each locus
    for (j in 1:y) { # for each of 1:y putative contributors
      Likelihood_df[j,i] <- MixtureLikelihood_Diploid(x=j,p.v=unlist(loci[[i]][!is.na(loci[[i]])])) # function
    }  # end j loop
  } # end i loop
  # Likelihood_df$Product <- apply(Likelihood_df, 1, function(x) prod(x*10000, na.rm = TRUE)) # product across all loci, exclude NAs
  Likelihood_df$Product <- apply(Likelihood_df, 1, function(x) prod(x, na.rm = TRUE)) # product across all loci, exclude NAs
  Likelihood_df$Contributors <- c(1:y) # putative numbers of contributors
  return(tail(Likelihood_df[Likelihood_df$Product==max(Likelihood_df$Product),],n=1)) # find row with max likelihood across all loci
} # end function

# df <- melt(Likelihood_df)
# df$contributors <- 1:150
# ggplot(df, aes(contributors,value)) + geom_point() + facet_wrap(variable ~ .)

# Multi-locus diploid contributor estimation function: remove alleles that monotonically increase
MixtureLikelihood2_remove <- function(y, loci){
  Likelihood_df <- data.frame()
  for (i in 1:length(loci)){ # for each locus
    for (j in 1:y) { # for each of 1:y putative contributors
      Likelihood_df[j,i] <- MixtureLikelihood_Diploid(x=j,p.v=unlist(loci[[i]][!is.na(loci[[i]])])) # function
    }  # end j loop
    if (is.na(Likelihood_df[j,i])==FALSE){ # if locus is present in mixture
      if ((tail(Likelihood_df[,i],n=1)==tail(max(Likelihood_df[,i]),n=1))==TRUE){ # if locus exhibits monotonic increase
        Likelihood_df[,i] <- "fail"
      } # end if
    } # end if
  } # end i loop
  Likelihood_df_sub <- Likelihood_df[,colSums(Likelihood_df=="fail", na.rm=TRUE)==0]
  if (class(Likelihood_df_sub)=="numeric"){ # if exactly 1 locus remains
    Likelihood_df$Product <- Likelihood_df_sub 
  } else if (length(Likelihood_df_sub)==0){ # if no loci remain
    Likelihood_df$Product <- NA
  #} else Likelihood_df$Product <- apply(Likelihood_df_sub, 1, function(x) prod(x*10000, na.rm = TRUE)) # product across all loci, exclude NAs
  } else Likelihood_df$Product <- apply(Likelihood_df_sub, 1, function(x) prod(x, na.rm = TRUE)) # product across all loci, exclude NAs
  Likelihood_df$Contributors <- c(1:y) # putative numbers of contributors
  Likelihood_df$n_failed <- length(grep("fail", Likelihood_df[j,])) # count number of failed loci
  if (class(Likelihood_df$Product)=="logical"){ # if all loci fail
    return(data.frame(Likelihood_df[y,1:length(loci)], Product=NA, Contributors=NA, n_failed=length(loci)))
  } else return(tail(Likelihood_df[Likelihood_df$Product==max(Likelihood_df$Product),],n=1)) # find row with max likelihood across all loci
} # end function

# Haploid contributor estimation: for mitochondrial markers 
MixtureLikelihood_Haploid <- function(y,p.v){
  Likelihood_df <- data.frame()
  combo <- list(combinations(length(p.v),length(p.v)))
  if(length(p.v)>1){ # only conduct if > 1 allele if observed
    for(i in 1:(length(p.v)-1)){
      temp.combo <- list(combinations(length(p.v),length(p.v)-i))
      combo <- c(combo, temp.combo)
    } # end i loop
    for (x in 1:y){ # for x putative contributors
      if(length(p.v)>(1*x)) {
        Likelihood_df_temp <- data.frame(Contributors=x, Likelihood=0)
        Likelihood_df <- rbind(Likelihood_df, Likelihood_df_temp)
      } else {
        counter <- 0
        sum.p <- sum(p.v)^(1*x)
        for(z in 2:(length(p.v))){
          counter <- counter+1
          temp.sum <- apply(combo[[z]],MARGIN=1,function(z){sum(p.v[z])^(1*x)})
          sum.p <- sum.p+((-1)^counter)*sum(temp.sum)
        } # end z loop
        Likelihood_df_temp <- data.frame(Contributors=x, Likelihood=sum.p)
        Likelihood_df <- rbind(Likelihood_df, Likelihood_df_temp)
        if (Likelihood_df$Likelihood[x] < Likelihood_df$Likelihood[x-1]) break
      } # end if else
    } # end x loop
  } else Likelihood_df <- data.frame(Contributors=1, Likelihood=0)
  return(tail(Likelihood_df[Likelihood_df$Likelihood==max(Likelihood_df$Likelihood, na.rm=TRUE),],n=1)) # find row with max likelihood across all loci
} # end function

# Microsatellite allele frequency calculation
CalculateMsatFreqs <- function(x){
  genind.obj <- df2genind(x, sep="\\|", ploidy=2, type="codom", pop=rep("pop1",nrow(x)),
                          loc.names=colnames(x), NA.char="NA|NA")
  genpop.obj <- genind2genpop(genind.obj) # turn into genpop object
  msat_freqs <- as.data.frame(t(makefreq(genpop.obj)))
  msat_freqs$locus <- data.frame(do.call('rbind', strsplit(as.character(rownames(msat_freqs)),'.',fixed=TRUE)))[,1] # create column of loci names
  msat_freqs$allele <- data.frame(do.call('rbind', strsplit(as.character(rownames(msat_freqs)),'.',fixed=TRUE)))[,2] # create column of alleles
  msat_freqs <- msat_freqs %>% dplyr::select(locus, allele, everything()) # move locus and allele columns to the front
  colnames(msat_freqs) <- c("locus","allele","freq")
  return(msat_freqs)
}

# SNP allele frequency calculation
CalculateSNPFreqs <- function(x) {
  Freq1 <- data.frame(locus=names(glMean(x)),
                      allele=1,
                      freq=glMean(x))
  Freq2 <- data.frame(locus=names(glMean(x)),
                      allele=2,
                      freq=1-glMean(x))
  SNP_freqs <- rbind(Freq1,Freq2)
  rownames(SNP_freqs) <- NULL
  SNP_freqs <- SNP_freqs[SNP_freqs$freq!=0,] # remove fixed loci
  return(SNP_freqs)
}

# Mitochondrial haplotype frequency calculation
CalculateMitoFreqs <- function(x){
  df <- data.frame(Elor_mt=x)
  genind.obj <- df2genind(df, ploidy=1, type="codom", pop=rep("pop1",nrow(df)),
                          loc.names=colnames(df), NA.char="NA")
  genpop.obj <- genind2genpop(genind.obj) # turn into genpop object
  mito_freqs <- as.data.frame(t(makefreq(genpop.obj)))
  mito_freqs$locus <- data.frame(do.call('rbind', strsplit(as.character(rownames(mito_freqs)),'.',fixed=TRUE)))[,1] # create column of loci names
  mito_freqs$allele <- data.frame(do.call('rbind', strsplit(as.character(rownames(mito_freqs)),'.',fixed=TRUE)))[,2] # create column of alleles
  mito_freqs <- mito_freqs %>% dplyr::select(locus, allele, everything()) # move locus and allele columns to the front
  colnames(mito_freqs) <- c("locus","allele","freq")
  return(mito_freqs)
}

#################################################################################
############################# Microstatellites ##################################
#################################################################################

# Read in coral goby msat dataset, subset to Lighthouse Reef gobies, remove loci not in HWE (assumption of mixture model)
msat_genotypes <- read.delim("datasets/Goby_allmsats_Rcode_nosouth_1Aug18.txt", row.names=1)
msat_genotypes <- msat_genotypes[grepl("L", row.names(msat_genotypes)),] # subset to Lighthouse Reef gobies
genind.obj <- df2genind(msat_genotypes, sep="\\|", ploidy=2, type="codom",
                        loc.names=colnames(msat_genotypes), NA.char="NA|NA")
# remove loci with deviation from HWE and heterzygote excess
hw.output <- data.frame(hw.test(genind.obj, 1000))
hw.output$observed_heterozygosity <- summary(genind.obj)[[6]]
hw.output$expected_heterozygosity <- summary(genind.obj)[[7]]
patterns <- row.names(hw.output[hw.output$Pr.exact< 0.05/89 & hw.output$observed>hw.output$expected,]) # bonferroni corrected pvalue
patterns <- c(patterns, names(genind.obj$loc.n.all[genind.obj$loc.n.all>20]))
msat_genotypes <- msat_genotypes[ , -which(names(msat_genotypes) %in% patterns)] # remove 49 loci

# Subset to different numbers of individuals in full and small msat panels
set.seed(984)
msats_25_inds_10_loci <- CalculateMsatFreqs(msat_genotypes[sample(nrow(msat_genotypes), 25),1:10]) # 25 individuals; 10 loci, 93 alleles
msats_25_inds_all_loci <- CalculateMsatFreqs(msat_genotypes[sample(nrow(msat_genotypes), 25),]) # 25 individuals; 44 loci, 382 alleles
msats_100_inds_10_loci <- CalculateMsatFreqs(msat_genotypes[sample(nrow(msat_genotypes), 100),1:10]) # 100 individuals; 10 loci, 130 alleles
msats_100_inds_all_loci <- CalculateMsatFreqs(msat_genotypes[sample(nrow(msat_genotypes), 100),]) # 100 individuals; 44 loci, 559 alleles

msat_freqs <- list(msats_25_inds_10_loci, msats_100_inds_10_loci,
                   msats_25_inds_all_loci, msats_100_inds_all_loci)
names(msat_freqs) <- c("msats_25_inds_10_loci", "msats_100_inds_10_loci",
                       "msats_25_inds_all_loci", "msats_100_inds_all_loci")

# Estimated number of individuals in mixtures of 1:100 individuals sampled from full dataset 
t=1 # index for each of 8 datasets
v=c(10, 10, 44, 44) # number of loci to subset per dataset
for (k in msat_freqs){ # for each dataset
  simulated_mixtures <- data.frame() 
  failed_loci <- data.frame() 
  for (i in rep(seq(from=2, to=100, by=2), each=100)){ 
    sample_genotypes <- msat_genotypes[sample(nrow(msat_genotypes), i), 1:v[t]] # sample i individuals and either 10 or 44 loci from total population
    sample_freqs <- CalculateMsatFreqs(sample_genotypes) # calculate sample allele freqs
    combined_allele_freqs <- merge(k, sample_freqs, by=c("locus","allele")) # merge sampled alleles with population allele freqs
    for (j in unique(k$locus)){ # create NA row for missing loci
      if (j %in% combined_allele_freqs$locus==FALSE) 
        combined_allele_freqs <- rbind(combined_allele_freqs, data.frame(locus=j, allele=NA, freq.x=NA, freq.y=NA))
    }
    loci.list <- split(combined_allele_freqs, combined_allele_freqs$locus)
    loci.list <- lapply(loci.list, function(x) {x$freq.x}) # get population frequencies for alleles detected in mixture 
    contrib_est <- MixtureLikelihood2(y=150, loci=loci.list) # run function, retain all loci
    contrib_est$true_N <- i
    contrib_est_remove <- MixtureLikelihood2_remove(y=150, loci=loci.list) # run function, remove failed loci
    contrib_est_remove$true_N <- i
    output <- data.frame(true_N=contrib_est$true_N, 
                         est_N=contrib_est$Contributors,
                         est_N_remove=contrib_est_remove$Contributors,
                         n_failed=contrib_est_remove$n_failed,
                         est_N_allelic_rich=ceiling(max(lengths(loci.list))/2)) # number of alleles/2
    failed_loci <- rbind(failed_loci, contrib_est_remove)
    simulated_mixtures <- rbind(simulated_mixtures, output)
    print(tail(simulated_mixtures))
  }
  write.csv(simulated_mixtures, paste0("datasets/simulated_mixtures_",names(msat_freqs)[t],".csv", sep=""), row.names=FALSE)
  write.csv(failed_loci, paste0("datasets/failed_loci_",names(msat_freqs)[t],".csv", sep=""), row.names=FALSE)
  write.csv(msat_freqs[t], paste0("datasets/freqs_",names(msat_freqs)[t],".csv", sep=""), row.names=FALSE)
  t <- t+1
}

#################################################################################
################################### SNPs ########################################
#################################################################################

# Read in coral goby SNP dataset, subset to Lighthouse Reef gobies, remove loci not in HWE (assumption of mixture model)
SNP_genotypes <- read.PLINK("datasets/goby_1040_allSNPs.raw")
locNames(SNP_genotypes) <- gsub("SGOBYREF_052516_","",locNames(SNP_genotypes)) # simplify locus names
SNP_genotypes <- SNP_genotypes[grepl("L", indNames(SNP_genotypes)),] # subset to Lighthouse Reef gobies
SNP_genotypes <- SNP_genotypes[,-grep("_0",locNames(SNP_genotypes))] # remove non-SNPs

# Subset to different numbers of individuals in full and small SNP panels
set.seed(621)
SNPs_25_inds_64_loci <- CalculateSNPFreqs(SNP_genotypes[sample(nrow(SNP_genotypes), 25),1:64]) # 25 individuals; 64 loci, 128 alleles
SNPs_100_inds_64_loci <-CalculateSNPFreqs(SNP_genotypes[sample(nrow(SNP_genotypes), 100),1:64]) # 100 individuals; 64 loci; 128 alleles
SNPs_25_inds_256_loci <- CalculateSNPFreqs(SNP_genotypes[sample(nrow(SNP_genotypes), 25),1:256]) # 25 individuals, 256 loci, 512 alleles
SNPs_100_inds_256_loci <- CalculateSNPFreqs(SNP_genotypes[sample(nrow(SNP_genotypes), 100),1:256]) # 100 individuals; 256 loci, 512 alleles
SNP_freqs <- list(SNPs_25_inds_64_loci, SNPs_100_inds_64_loci, 
                  SNPs_25_inds_256_loci, SNPs_100_inds_256_loci)
names(SNP_freqs) <- c("SNPs_25_inds_64_loci", "SNPs_100_inds_64_loci", 
                      "SNPs_25_inds_256_loci", "SNPs_100_inds_256_loci")

# Estimated number of individuals in mixtures of 1:100 individuals sampled from full dataset 
t=1
v=c(64,64,256,256)
for (k in SNP_freqs){
  simulated_mixtures <- data.frame()
  for (i in rep(seq(from=2, to=100, by=2), each=100)){
    sample_genotypes <- SNP_genotypes[sample(nrow(SNP_genotypes), i),1:v[t]] # sample i individuals and either 64 or 256 loci
    sample_freqs <- CalculateSNPFreqs(sample_genotypes) # calculate sample allele freqs
    k <- k[k$freq!=1,] # remove fixed loci
    combined_allele_freqs <- merge(k, sample_freqs, by = c("locus","allele")) # merge sampled alleles with population allele freqs
    for (j in unique(k$locus)){ # create NA row for missing loci
      if (j %in% combined_allele_freqs$locus==FALSE) 
        combined_allele_freqs <- rbind(combined_allele_freqs, data.frame(locus=j, allele=NA, freq.x=NA, freq.y=NA))
    }
    loci.list <- split(combined_allele_freqs, combined_allele_freqs$locus) # split by locus
    loci.list <- lapply(loci.list, function(x) {x$freq.x}) # population frequencies for alleles detected in mixture 
    contrib_est <- MixtureLikelihood2(y=150, loci=loci.list) # run function
    contrib_est$true_N <- i
    output <- data.frame(true_N=contrib_est$true_N, 
                         est_N=contrib_est$Contributors,
                         est_N_allelic_rich=ceiling(max(lengths(loci.list))/2))
    simulated_mixtures <- rbind(simulated_mixtures, output)
    print(tail(simulated_mixtures))
  }
  write.csv(simulated_mixtures, paste0("datasets/simulated_mixtures_",names(SNP_freqs)[t],".csv", sep=""), row.names=FALSE)
  write.csv(SNP_freqs[t], paste0("datasets/freqs_",names(SNP_freqs)[t],".csv", sep=""), row.names=FALSE)
  t <- t+1
}

#################################################################################
################################ Mitogenome #####################################
#################################################################################

# Read in coral goby msat dataset, subset to Lighthouse Reef gobies, remove loci not in HWE (assumption of mixture model)
mito_genotypes <- read.delim("datasets/goby_hap_geno_mtDNA_11May22.txt", row.names=1)
mito_genotypes <- mito_genotypes[grepl("L", row.names(mito_genotypes)),] # subset to Lighthouse Reef gobies
mito_genotypes <- mito_genotypes[,grep("mt", colnames(mito_genotypes))] # subset to all seqs (10Kb) and subset (4Kb)

# Subset to different numbers of individuals
set.seed(593)
mito_25_inds_4k <- CalculateMitoFreqs(mito_genotypes[sample(nrow(mito_genotypes), 25),2]) # 25 individuals; 20 alleles
mito_100_inds_4k <- CalculateMitoFreqs(mito_genotypes[sample(nrow(mito_genotypes), 100), 2]) # 100 individuals; 39 alleles
mito_freqs <- list(mito_25_inds_4k, mito_100_inds_4k)
names(mito_freqs) <- c("mito_25_inds_4k", "mito_100_inds_4k")

# Read in coral goby msat dataset, subset to Lighthouse Reef gobies, remove loci not in HWE (assumption of mixture model)
mito_genotypes <- read.delim("datasets/goby_hap_geno_mtDNA_11May22.txt", row.names=1)
mito_genotypes <- mito_genotypes[grepl("L", row.names(mito_genotypes)),] # subset to Lighthouse Reef gobies
z <- lapply(mito_genotypes, unique)
lengths(z)
mito_genotypes <- mito_genotypes$Elor_16k_78

# Try with just a single segment: Elor_16k_78
set.seed(691)
mito_25_inds_Elor_16k_78 <- CalculateMitoFreqs(sample(mito_genotypes, 25)) # 25 individuals; 9 alleles
mito_100_inds_Elor_16k_78 <- CalculateMitoFreqs(sample(mito_genotypes, 100)) # 100 individuals; 24 allele
mito_freqs <- list(mito_25_inds_Elor_16k_78, mito_100_inds_Elor_16k_78)
names(mito_freqs) <- c("mito_25_inds_Elor_16k_78", "mito_100_inds_Elor_16k_78")

# Estimated number of individuals in mixtures of 1:100 individuals sampled from each dataset 
t=1 # index for each of 2 datasets
for (k in mito_freqs){ # for each dataset
  simulated_mixtures <- data.frame() 
  for (i in rep(seq(from=2, to=100, by=2), each=100)){ 
    sample_genotypes <- mito_genotypes[sample(nrow(mito_genotypes), i),2] # sample i individuals from total population
    #sample_genotypes <- sample(mito_genotypes, i) # sample i individuals from total population
    sample_freqs <- CalculateMitoFreqs(sample_genotypes) # calculate sample allele freqs
    combined_allele_freqs <- merge(k, sample_freqs, by=c("locus","allele")) # merge sampled alleles with population allele freqs
    p.v <- combined_allele_freqs$freq.x # get population frequencies for alleles detected in mixture 
    if (length(p.v)>0){
    contrib_est <- MixtureLikelihood_Haploid(y=150, p.v=p.v) # run function
    contrib_est$true_N <- i
    output <- data.frame(true_N=contrib_est$true_N, 
                         est_N=contrib_est$Contributors,
                         est_N_allelic_rich=max(length(p.v))) # number of haplotypes
    } else {
      output <- data.frame(true_N=i, 
                           est_N=0,
                           est_N_allelic_rich=0)
    }
    simulated_mixtures <- rbind(simulated_mixtures, output)
    print(tail(simulated_mixtures))
  }
  write.csv(simulated_mixtures, paste0("datasets/simulated_mixtures_",names(mito_freqs)[t],".csv", sep=""), row.names=FALSE)
  write.csv(mito_freqs[t], paste0("datasets/freqs_",names(mito_freqs)[t],".csv", sep=""), row.names=FALSE)
  t <- t+1
}

#################################################################################
################################## Plots ########################################
#################################################################################

### Contributor estimator plotting function
plot_cont_est <- function(x, y, nloci, N_ind_x, N_ind_y, marker, ylabel) {
  simulated_mixtures_x <- read.csv(paste0("datasets/simulated_mixtures_",x,".csv", sep=""))
  simulated_mixtures_x$N_ind <- as.factor(rep(N_ind_x))
  simulated_mixtures_y <- read.csv(paste0("datasets/simulated_mixtures_",y,".csv", sep=""))
  simulated_mixtures_y$N_ind <- as.factor(rep(N_ind_y))
  simulated_mixtures <- rbind(simulated_mixtures_x, simulated_mixtures_y)
  if (marker=="msat") {
    simulated_mixtures$est_N_remove[simulated_mixtures$est_N_remove==150] <- NA # if reaches max estimate, assume mixture failed
    simulated_mixtures$bias <- simulated_mixtures$est_N_remove-simulated_mixtures$true_N
  } else {
    simulated_mixtures$est_N[simulated_mixtures$est_N==150] <- NA # if reaches max estimate, assume mixture failed
    simulated_mixtures$bias <- simulated_mixtures$est_N-simulated_mixtures$true_N
  }
  p1 <- ggplot(simulated_mixtures) + 
    geom_point(aes(x=true_N, y=bias, color=N_ind)) +
    scale_color_manual(labels=c("25 individuals","100 individuals"), 
                       values=c(alpha("#9aab89",0.05), alpha("#58508d",0.05))) + 
    stat_summary(data=simulated_mixtures[simulated_mixtures$N_ind==N_ind_x,], 
                 aes(x=true_N, y=bias), geom="point", fun="mean", col="#9aab89", size=2, shape=19) +
    stat_summary(data=simulated_mixtures[simulated_mixtures$N_ind==N_ind_x,],
                 aes(x=true_N, y=bias), fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", color="#9aab89", width=0.2) +
    stat_summary(data=simulated_mixtures[simulated_mixtures$N_ind==N_ind_y,], 
                 aes(x=true_N, y=bias), geom="point", fun="mean", col="#58508d", size=2, shape=19) +
    stat_summary(data=simulated_mixtures[simulated_mixtures$N_ind==N_ind_y,],
                 aes(x=true_N, y=bias), fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", color="#58508d", width=0.2) +
    ylab("Bias (Estimate - True N)") + xlab("True N") +
    ylim(-100,100) + xlim(0,100) +
    geom_abline(slope=0, intercept=0) +
    theme_bw() +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=20)) +
    theme(legend.position=c(0.25, 0.88),
            legend.text=element_text(size=16),
            legend.title=element_text(size=16)) +
    labs(color = "Allele frequencies") +
    guides(colour=guide_legend(override.aes = list(alpha=1, size=3))) +
    theme(legend.background = element_rect(colour = "black", size=0.3)) +
    if (ylabel==FALSE) theme(axis.title.y=element_blank())
  density_x <- read.csv(paste0("datasets/freqs_",x,".csv", sep=""), col.names=c("locus","allele","freq"))
  density_x$inds <- as.factor(rep(N_ind_x))
  density_y <- read.csv(paste0("datasets/freqs_",y,".csv", sep=""), col.names=c("locus","allele","freq"))
  density_y$inds <- as.factor(rep(N_ind_y))
  density_comb <- rbind(density_x, density_y)
  mu <- ddply(density_comb, "inds", summarise, grp.mean=mean(freq))
  print(mu)
  p3 <- ggplot(density_comb, aes(x=freq, y = ..count../sum(..count..)*100, color=inds))+ 
    geom_density(stat='bin', cex=1, bins=20, aes(fill=factor(inds)), alpha=0.3) +
    geom_vline(data=mu, aes(xintercept=grp.mean, color=inds), linetype="dashed", cex=0.8) +
    scale_x_continuous(n.breaks=4) +
    scale_color_manual(values=c("#9aab89","#58508d")) +
    scale_fill_manual(values=c("#9aab89","#58508d")) +
    theme_bw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16)) +
    xlab("Allele frequency") + ylab("%") +
    theme(legend.position="none")
  p <- p1 + annotation_custom(ggplotGrob(p3), xmin = -5, xmax = 60, 
                              ymin = -110, ymax = -45)
  if (marker=="msat") {
    p2 <- ggplot(simulated_mixtures) +
      geom_line(aes(x=true_N, y=n_failed, color=N_ind),
                stat="summary", fun="mean") +
      scale_color_manual(values=c("#9aab89","#58508d")) +
      ylab("Failed loci") + xlab("") + ylim(0, nloci) + xlim(0,100) +
      theme_bw() +
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=16)) +
      theme(legend.position = "none")
    p <- p + annotation_custom(ggplotGrob(p2), xmin = 40, xmax = 100, 
                               ymin = 40, ymax = 105)
  }
  return(p)
}

# Save all of the plots
msat_all_loci <- plot_cont_est("msats_25_inds_all_loci","msats_100_inds_all_loci", 44, 25, 100, "msat", ylabel=TRUE)
snp_256_loci <- plot_cont_est("SNPs_25_inds_256_loci", "SNPs_100_inds_256_loci", 256, 25, 100, "snps", ylabel=FALSE)
mito_4k <- plot_cont_est("mito_25_inds_4k", "mito_100_inds_4k", 1, 25, 100, "mito", ylabel=FALSE)
msat_10_loci <- plot_cont_est("msats_25_inds_10_loci", "msats_100_inds_10_loci", 10, 25, 100, "msat", ylabel=TRUE)
snp_64_loci <- plot_cont_est("SNPs_25_inds_64_loci", "SNPs_100_inds_64_loci", 64, 25, 100, "snps", ylabel=FALSE)
mito_Elor_16k_78 <- plot_cont_est("mito_25_inds_Elor_16k_78", "mito_100_inds_Elor_16k_78", 1, 25, 100, "mito", ylabel=FALSE)
arrange_plots <- ggarrange(msat_all_loci, snp_256_loci, mito_4k,
                           msat_10_loci, snp_64_loci, mito_Elor_16k_78,
                           labels = c("Microsatellites",
                                      "SNPs",
                                      "Mitochondrial haplotypes", "","",""), 
                           ncol=3, nrow=2, hjust=-0.7, vjust=0,
                           font.label = list(size = 20))
# ggsave(filename=paste("Figure_2_simulated_mixtures_1.24.23.pdf"), plot=arrange_plots, dpi=300, width=16, height=12, units="in")

### Mean (sd) bias for table in manuscript
mean_bias <- function(x, marker, values) {
  simulated_mixtures <- read.csv(paste0("datasets/simulated_mixtures_",x,".csv", sep=""))
  if (marker=="msat") {
    simulated_mixtures$est_N_remove[simulated_mixtures$est_N_remove==150] <- NA # if reaches max estimate, assume mixture failed
    simulated_mixtures$bias <- simulated_mixtures$est_N_remove-simulated_mixtures$true_N
  } else {
    simulated_mixtures$est_N[simulated_mixtures$est_N==150] <- NA # if reaches max estimate, assume mixture failed
    simulated_mixtures$bias <- simulated_mixtures$est_N-simulated_mixtures$true_N
  }
  dataset_summary <- as.data.frame(simulated_mixtures %>% 
    group_by(true_N) %>% 
    summarise(across(bias, list(mean=mean, sd=sd), na.rm = TRUE)))
  return(dataset_summary[dataset_summary$true_N %in% values,])
}

mean_bias("msats_100_inds_all_loci", "msat", c(10,40,100))
mean_bias("msats_25_inds_all_loci", "msat", c(10,40,100))
mean_bias("msats_100_inds_10_loci", "msat", c(10,40,100))
mean_bias("msats_25_inds_10_loci", "msat", c(10,40,100))
mean_bias("SNPs_100_inds_256_loci", "SNP", c(10,40,100))
mean_bias("SNPs_25_inds_256_loci", "SNP", c(10,40,100))
mean_bias("SNPs_100_inds_64_loci", "SNP", c(10,40,100))
mean_bias("SNPs_25_inds_64_loci", "SNP", c(10,40,100))
mean_bias("mito_100_inds_4k", "mito", c(10,40,100))
mean_bias("mito_25_inds_4k", "mito", c(10,40,100))
mean_bias("mito_100_inds_Elor_16k_78", "mito", c(10,40,100))
mean_bias("mito_25_inds_Elor_16k_78", "mito", c(10,40,100))
