### The goal of this code is to simulate mixtures of DNA to determine
### how many individuals can be estimated from the contributor estimation
### using real goby genotypes
### Last updated 6.1.2022 by Kara Andres (kja68@cornell.edu)

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

#################################################################################
################################## Functions ####################################
#################################################################################
# Contributor estimator: following Weir, Haned: sum.p = 0 if # alleles exceeds 2x
# Suresh's code -- speed up computation time when # of alleles is large
MixtureLikelihood_Parsed <- function(x,p.v){
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
              row.sums.pow2x[m] <- sum(rowSums(p.v.m)^(2*x)) # reference as folows
            } # end m loop
            sum.p <- sum.p+((-1)^counter)*sum(row.sums.pow2x) } else 
            {
              # non-parsed
              temp.combo.v <- c(temp.combo)
              p.v.m <- matrix(p.v[temp.combo.v],nr=dim(temp.combo)[1],nc=dim(temp.combo)[2],byrow=F)
              row.sums.pow2x <- rowSums(p.v.m)^(2*x) # reference as folows
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

# Multi-locus function across 1:x number of putative contributors
MixtureLikelihood2 <- function(y, loci){
  Likelihood_df <- data.frame()
  for (i in 1:length(loci)){ # for each locus
    for (j in 1:y) { # for each of 1:y putative contributors
      Likelihood_df[j,i] <- MixtureLikelihood_Parsed(x=j,p.v=unlist(loci[[i]][!is.na(loci[[i]])])) # function
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

# Multi-locus contributor estimator function: remove alleles that monotonically increase
MixtureLikelihood2_remove <- function(y, loci){
  Likelihood_df <- data.frame()
  for (i in 1:length(loci)){ # for each locus
    for (j in 1:y) { # for each of 1:y putative contributors
      Likelihood_df[j,i] <- MixtureLikelihood_Parsed(x=j,p.v=unlist(loci[[i]][!is.na(loci[[i]])])) # function
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
msat_genotypes <- read.delim("/Users/kbja10/Documents/Cornell/Research/Round goby/Intraspecific_variation_with_eDNA/datasets/Goby_allmsats_Rcode_nosouth_1Aug18.txt", row.names=1)
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
msats_25_inds_10_loci <- CalculateMsatFreqs(msat_genotypes[sample(nrow(msat_genotypes), 25),1:10]) # 25 individuals; 10 loci, 91 alleles
msats_all_inds_10_loci <- CalculateMsatFreqs(msat_genotypes[,1:10]) # 127 individuals; 10 loci, 133 alleles
msats_25_inds_all_loci <- CalculateMsatFreqs(msat_genotypes[sample(nrow(msat_genotypes), 25),]) # 25 individuals; 44 loci, 402 alleles
msats_all_inds_all_loci <- CalculateMsatFreqs(msat_genotypes) # 127 individuals; 44 loci, 568 alleles
msats_50_inds_10_loci <- CalculateMsatFreqs(msat_genotypes[sample(nrow(msat_genotypes), 50),1:10]) # 50 individuals; 10 loci, 121 alleles
msats_50_inds_all_loci <- CalculateMsatFreqs(msat_genotypes[sample(nrow(msat_genotypes), 50),]) # 50 individuals; 44 loci, 460 alleles
msats_100_inds_10_loci <- CalculateMsatFreqs(msat_genotypes[sample(nrow(msat_genotypes), 100),1:10]) # 100 individuals; 10 loci, 121 alleles
msats_100_inds_all_loci <- CalculateMsatFreqs(msat_genotypes[sample(nrow(msat_genotypes), 100),]) # 100 individuals; 44 loci, 460 alleles

msat_freqs <- list(msats_25_inds_10_loci, msats_all_inds_10_loci, msats_50_inds_10_loci,msats_100_inds_10_loci,
                   msats_25_inds_all_loci, msats_all_inds_all_loci, msats_50_inds_all_loci, msats_100_inds_all_loci)
names(msat_freqs) <- c("msats_25_inds_10_loci", "msats_all_inds_10_loci", "msats_50_inds_10_loci", "msats_100_inds_10_loci",
                       "msats_25_inds_all_loci", "msats_all_inds_all_loci", "msats_50_inds_all_loci", "msats_100_inds_all_loci")

# Estimated number of individuals in mixtures of 1:100 individuals sampled from full dataset 
t=1 # index for each of 8 datasets
v=c(10, 10, 10, 10, 44, 44, 44, 44) # number of loci to subset per dataset
for (k in msat_freqs){ # for each dataset
  simulated_mixtures <- data.frame() 
  failed_loci <- data.frame() 
  for (i in rep(seq(from=2, to=100, by=2), each=10)){ 
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
  write.csv(simulated_mixtures, paste0("/Users/kbja10/Documents/Cornell/Research/Round goby/Intraspecific_variation_with_eDNA/simulated_mixtures_",names(msat_freqs)[t],".csv", sep=""), row.names=FALSE)
  write.csv(failed_loci, paste0("/Users/kbja10/Documents/Cornell/Research/Round goby/Intraspecific_variation_with_eDNA/failed_loci_",names(msat_freqs)[t],".csv", sep=""), row.names=FALSE)
  write.csv(msat_freqs[t], paste0("/Users/kbja10/Documents/Cornell/Research/Round goby/Intraspecific_variation_with_eDNA/freqs_",names(msat_freqs)[t],".csv", sep=""), row.names=FALSE)
  t <- t+1
}

#################################################################################
################################### SNPs ########################################
#################################################################################

# Read in coral goby SNP dataset, subset to Lighthouse Reef gobies, remove loci not in HWE (assumption of mixture model)
SNP_genotypes <- read.PLINK("/Users/kbja10/Documents/Cornell/Research/Round goby/Intraspecific_variation_with_eDNA/datasets/goby_1040_allSNPs.raw")
locNames(SNP_genotypes) <- gsub("SGOBYREF_052516_","",locNames(SNP_genotypes)) # simplify locus names
SNP_genotypes <- SNP_genotypes[grepl("L", indNames(SNP_genotypes)),] # subset to Lighthouse Reef gobies
SNP_genotypes <- SNP_genotypes[,-grep("_0",locNames(SNP_genotypes))] # remove non-SNPs

# Subset to different numbers of individuals in full and small SNP panels
set.seed(621)
SNPs_25_inds_64_loci <- CalculateSNPFreqs(SNP_genotypes[sample(nrow(SNP_genotypes), 25),1:64]) # 25 individuals; 64 loci, 128 alleles
SNPs_all_inds_64_loci <- CalculateSNPFreqs(SNP_genotypes[,1:64]) # 127 individuals; 64 loci; 128 alleles
SNPs_25_inds_256_loci <- CalculateSNPFreqs(SNP_genotypes[sample(nrow(SNP_genotypes), 25),1:256]) # 25 individuals, 256 loci, 512 alleles
SNPs_all_inds_256_loci <- CalculateSNPFreqs(SNP_genotypes[,1:256]) # 127 individuals; 256 loci, 512 alleles
SNP_freqs <- list(SNPs_25_inds_64_loci, SNPs_all_inds_64_loci, SNPs_25_inds_256_loci, SNPs_all_inds_256_loci)
names(SNP_freqs) <- c("SNPs_25_inds_64_loci", "SNPs_all_inds_64_loci", "SNPs_25_inds_256_loci", "SNPs_all_inds_256_loci")

SNPs_all_inds_64_loci[27,]

# Estimated number of individuals in mixtures of 1:100 individuals sampled from full dataset 
t=1
v=c(64,64,256,256)
for (k in SNP_freqs){
  simulated_mixtures <- data.frame()
  for (i in rep(seq(from=2, to=100, by=2), each=10)){
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
  write.csv(simulated_mixtures, paste0("/Users/kbja10/Documents/Cornell/Research/Round goby/Intraspecific_variation_with_eDNA/simulated_mixtures_",names(SNP_freqs)[t],".csv", sep=""), row.names=FALSE)
  t <- t+1
}

#################################################################################
################################## Plots ########################################
#################################################################################

# Microsatellite plots
plot_msats <- function(x, y, nloci, N_ind_x, N_ind_y) {
  simulated_mixtures_x <- read.csv(paste0("/Users/kbja10/Documents/Cornell/Research/Round goby/Intraspecific_variation_with_eDNA/simulated_mixtures_",x,".csv", sep=""))
  simulated_mixtures_x$N_ind <- as.factor(rep(N_ind_x))
  simulated_mixtures_y <- read.csv(paste0("/Users/kbja10/Documents/Cornell/Research/Round goby/Intraspecific_variation_with_eDNA/simulated_mixtures_",y,".csv", sep=""))
  simulated_mixtures_y$N_ind <- as.factor(rep(N_ind_y))
  simulated_mixtures <- rbind(simulated_mixtures_x, simulated_mixtures_y)
  simulated_mixtures$est_N_remove[simulated_mixtures$est_N_remove==150] <- NA # if reaches max estimate, assume mixture failed
  simulated_mixtures$est_N[simulated_mixtures$est_N==150] <- NA # if reaches max estimate, assume mixture failed
  simulated_mixtures$bias <- simulated_mixtures$est_N_remove-simulated_mixtures$true_N
  p1 <- ggplot(simulated_mixtures) + 
    geom_point(aes(x=true_N, y=bias, color=N_ind)) +
    scale_color_manual(values=c(alpha("#58508d",0.2), alpha("#9aab89",0.2))) + 
    stat_summary(data=simulated_mixtures[simulated_mixtures$N_ind==N_ind_x,], 
                 aes(x=true_N, y=bias), geom="point", fun="mean", col="#58508d", size=2, shape=19) +
    stat_summary(data=simulated_mixtures[simulated_mixtures$N_ind==N_ind_x,],
                 aes(x=true_N, y=bias), fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", color="#58508d", width=0.2) +
    stat_summary(data=simulated_mixtures[simulated_mixtures$N_ind==N_ind_y,], 
                  aes(x=true_N, y=bias), geom="point", fun="mean", col="#9aab89", size=2, shape=19) +
    stat_summary(data=simulated_mixtures[simulated_mixtures$N_ind==N_ind_y,],
                  aes(x=true_N, y=bias), fun.data=mean_sdl, fun.args = list(mult=1), 
                  geom="errorbar", color="#9aab89", width=0.2) +
    ylab("Bias (Estimate - True number of individuals)") + xlab("True number of individuals") +
    ylim(-100,100) + xlim(0,100) +
    geom_abline(slope=0, intercept=0) +
    theme_bw() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=18)) +
    theme(legend.position = c(0.13, 0.88)) +
    labs(color = "Number of \n individuals") +
    guides(colour=guide_legend(override.aes = list(alpha = 1))) +
    theme(legend.background = element_rect(colour = "black", size=0.3))
  p2 <- ggplot(simulated_mixtures) +
    geom_line(aes(x=true_N, y=n_failed, color=N_ind),
              stat="summary", fun="mean") +
    scale_color_manual(values=c("#58508d","#9aab89")) +
    ylab("N loci") + xlab("") + ylim(0, nloci) + xlim(0,100) +
    theme_bw() +
    theme(legend.position = "none")
  hist_x <- read.csv(paste0("/Users/kbja10/Documents/Cornell/Research/Round goby/Intraspecific_variation_with_eDNA/freqs_",x,".csv", sep=""), col.names=c("locus","allele","freq"))
  hist_y <- read.csv(paste0("/Users/kbja10/Documents/Cornell/Research/Round goby/Intraspecific_variation_with_eDNA/freqs_",y,".csv", sep=""), col.names=c("locus","allele","freq"))
  hist_info_x <- hist(hist_x$freq, plot = FALSE)
  hist_info_x$density <- hist_info_x$counts/sum(hist_info_x$counts)*100  
  hist_info_y <- hist(hist_y$freq, plot = FALSE)
  hist_info_y$density <- hist_info_y$counts/sum(hist_info_y$counts)*100
  p3 <- ggplot() + 
    geom_bar(aes(x=hist_info_x$mids, y=hist_info_x$density), 
             color="#58508d", fill="#58508d", stat="identity") + 
    ylim(0,100) + ylab("") + xlab("Allele frequency") + ylab("Percent") +
    theme_bw() +
    theme(axis.text=element_text(size=5),
          axis.title=element_text(size=8))
  p4 <- ggplot() + 
    geom_bar(aes(x=hist_info_y$mids, y=hist_info_y$density), 
             color="#9aab89", fill="#9aab89", stat="identity") + 
    ylim(0,100) + ylab("") + xlab("Allele frequency") + ylab("Percent")  + 
    theme_bw() +
    theme(axis.text=element_text(size=5),
          axis.title=element_text(size=8))
  p5 <- p1 + annotation_custom(ggplotGrob(p2), xmin = 50, xmax = 100, 
                  ymin = 50, ymax = 100)
  p6 <- p5 + annotation_custom(ggplotGrob(p3), xmin = 0, xmax = 25, 
                               ymin = -100, ymax = -50)
  p <- p6 + annotation_custom(ggplotGrob(p4), xmin = 25, xmax = 50, 
                              ymin =-100, ymax = -50)
  return(p1)
}

loci_10 <- plot_msats("msats_25_inds_10_loci", "msats_100_inds_10_loci", 10, 25, 100)
loci_10
loci_all <- plot_msats("msats_25_inds_all_loci","msats_100_inds_all_loci", 44, 25, 100)
loci_all

# Save all of the plots
arrange_plots <- ggarrange(loci_10, loci_all, 
                           labels = c("(A) 10 loci",
                                      "(B) 44 loci"), 
                           ncol=2, nrow=1, hjust=-0.7, vjust=0,
                           font.label = list(size = 20))
arrange_plots
ggsave(filename = paste("/Users/kbja10/Documents/Cornell/Research/Round goby/Intraspecific_variation_with_eDNA/Figure_3_simulated_mixtures_msats_5.12.22.pdf"), plot=arrange_plots, dpi=300, width=12, height=6, units="in")
# ggsave(filename = paste("/Users/kbja10/Documents/Cornell/Research/Presentations/JASM/simulated_mixtures_msats.pdf"), plot=arrange_plots, dpi=300, width=12, height=6, units="in")


# SNP plots
plot_SNPs <- function(x) {
  simulated_mixtures <- read.csv(paste0("/Users/kbja10/Documents/Cornell/Research/Round goby/Intraspecific_variation_with_eDNA/simulated_mixtures_",x,".csv", sep=""))
  simulated_mixtures$est_N[simulated_mixtures$est_N==150] <- NA # if reaches max estimate, assume mixture failed
  p <- ggplot(simulated_mixtures, aes(x=true_N, y=est_N)) + 
    geom_point(color=alpha("#58508d",0.2)) +
    stat_summary(geom="point", fun="mean", col="#58508d", size=2, shape=19) +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", color="#58508d", width=0.2) +
    geom_point(data=simulated_mixtures[is.na(simulated_mixtures$est_N),], aes(x=true_N, y=0),
               shape=4, color="#58508d") +
    ylab("Estimated number of individuals") + xlab("True number of individuals") +
    ylim(0,140) + xlim(0,100) +
    geom_abline(slope=1, intercept=0) +
    theme_bw()
  return(p)
}

A <- plot_SNPs("SNPs_25_inds_64_loci")
B <- plot_SNPs("SNPs_all_inds_64_loci")
C <- plot_SNPs("SNPs_25_inds_256_loci")
D <- plot_SNPs("SNPs_all_inds_256_loci")

# Save all of the plots
arrange_plots <- ggarrange(A, B, C, D, 
                           labels = c("(A) 64 loci, 25 individuals",
                                      "(B) 64 loci, all individuals",
                                      "(C) 256 loci, 25 individuals",
                                      "(D) 256 loci, all individuals"), 
                           ncol=2, nrow=2, hjust=-0.25)
arrange_plots
ggsave(filename = paste("/Users/kbja10/Documents/Cornell/Research/Round goby/Intraspecific_variation_with_eDNA/Figure_4_simulated_mixtures_SNPs.pdf"), plot=arrange_plots, dpi=300, width=8, height=8, units="in")


#################################################################################
################################ Mitogenome ########################################
#################################################################################

# Read in coral goby msat dataset, subset to Lighthouse Reef gobies, remove loci not in HWE (assumption of mixture model)
mito_genotypes <- read.delim("/Users/kbja10/Documents/Cornell/Research/Round goby/Intraspecific_variation_with_eDNA/datasets/goby_hap_geno_mtDNA_11May22.txt", row.names=1)
mito_genotypes <- mito_genotypes[grepl("L", row.names(mito_genotypes)),] # subset to Lighthouse Reef gobies
mito_genotypes <- mito_genotypes[,grep("mt", colnames(mito_genotypes))] # subset to all seqs (10Kb) and subset (4Kb)

# Subset to different numbers of individuals in full and small msat panels
set.seed(593)
mito_25_inds_4k <- CalculateMitoFreqs(mito_genotypes[sample(nrow(mito_genotypes), 25),2]) # 25 individuals; 16 alleles
mito_50_inds_4k <- CalculateMitoFreqs(mito_genotypes[sample(nrow(mito_genotypes), 50), 2]) # 50 individuals; 27 alleles
mito_100_inds_4k <- CalculateMitoFreqs(mito_genotypes[sample(nrow(mito_genotypes), 100), 2]) # 100 individuals; 37 alleles
mito_all_inds_4k <- CalculateMitoFreqs(mito_genotypes[,2]) # 127 individuals; 45 alleles
mito_freqs <- list(mito_25_inds_4k, mito_50_inds_4k, mito_100_inds_4k, mito_all_inds_4k)
names(mito_freqs) <- c("mito_25_inds_4k", "mito_50_inds_4k", "mito_100_inds_4k", "mito_all_inds_4k")

# Estimated number of individuals in mixtures of 1:100 individuals sampled from each dataset 
t=1 # index for each of 4 datasets
for (k in mito_freqs){ # for each dataset
  simulated_mixtures <- data.frame() 
  failed_loci <- data.frame() 
  for (i in rep(seq(from=2, to=100, by=2), each=10)){ 
    sample_genotypes <- mito_genotypes[sample(nrow(mito_genotypes), i),2] # sample i individuals and either 10 or 44 loci from total population
    sample_freqs <- CalculateMitoFreqs(sample_genotypes) # calculate sample allele freqs
    combined_allele_freqs <- merge(k, sample_freqs, by=c("locus","allele")) # merge sampled alleles with population allele freqs
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
  write.csv(simulated_mixtures, paste0("/Users/kbja10/Documents/Cornell/Research/Round goby/Intraspecific_variation_with_eDNA/simulated_mixtures_",names(mito_freqs)[t],".csv", sep=""), row.names=FALSE)
  write.csv(failed_loci, paste0("/Users/kbja10/Documents/Cornell/Research/Round goby/Intraspecific_variation_with_eDNA/failed_loci_",names(mito_freqs)[t],".csv", sep=""), row.names=FALSE)
  write.csv(mito_freqs[t], paste0("/Users/kbja10/Documents/Cornell/Research/Round goby/Intraspecific_variation_with_eDNA/freqs_",names(mito_freqs)[t],".csv", sep=""), row.names=FALSE)
  t <- t+1
}

#################################################################################
################################ Simulated data ########################################
#################################################################################

# Simulate data for 50 loci with alleles at equal frequency: 0.05, 0.10, 0.20 and w/ realistic allele frequency spectrum
set.seed(593)
p05 <- data.frame(locus=rep(1:50, each=20), allele=rep(1:20,50), freq=0.05, 
                  row.names=paste0(rep(1:50, each=20), "_", rep(1:20,50), sep=""))
p10 <- data.frame(locus=rep(1:50, each=10), allele=rep(1:10,50), freq=0.10,
                  row.names=paste0(rep(1:50, each=10), "_", rep(1:10,50), sep=""))
p20 <- data.frame(locus=rep(1:50, each=5), allele=rep(1:5,50), freq=0.20,
                  row.names=paste0(rep(1:50, each=5), "_", rep(1:5,50), sep=""))
p001 <- data.frame(locus=rep(1:50, each=20), allele=rep(1:20,50), freq=c(rep(c(0.001,0.009),each=8),0.03,0.09,0.3,0.5),
                   row.names=paste0(rep(1:50, each=20), "_", rep(1:20,50), sep=""))
hist(p001_equal$freq)

sim_freqs <- list(p05, p10, p20, p001)
names(sim_freqs) <- c("p05", "p10", "p20","p001")

# Estimated number of individuals in mixtures of 1:10 individuals sampled from each dataset 
t=1 # index for each dataset
for (k in sim_freqs){ # for each dataset
  simulated_mixtures <- data.frame() 
  # failed_loci <- data.frame()
  for (m in c(10,20,30,50)){
  for (i in rep(c(1,3,5,7,10), each=50)){ 
    sample_genotypes <- data.frame(ID=1:i)
    for (n in unique(k$locus)[1:m]){
      sample_genotypes <- cbind(sample_genotypes, paste0(sample(k[k$locus==n,"allele"], i, replace=TRUE, prob=k[k$locus==n,"freq"]),"|",sample(k[k$locus==n,"allele"], i, replace=TRUE, prob=k[k$locus==n,"freq"])))
    }
    colnames(sample_genotypes) <- c("ID", 1:m)
    sample_freqs <- CalculateMsatFreqs(sample_genotypes[,-1]) # calculate sample allele freqs
    combined_allele_freqs <- merge(k, sample_freqs, by=c("locus","allele")) # merge sampled alleles with population allele freqs
    for (j in unique(k$locus)[1:m]){ # create NA row for missing loci
      if (j %in% combined_allele_freqs$locus==FALSE) 
        combined_allele_freqs <- rbind(combined_allele_freqs, data.frame(locus=j, allele=NA, freq.x=NA, freq.y=NA))
    }
    loci.list <- split(combined_allele_freqs, combined_allele_freqs$locus)
    loci.list <- lapply(loci.list, function(x) {x$freq.x}) # get population frequencies for alleles detected in mixture 
    contrib_est <- MixtureLikelihood2(y=50, loci=loci.list) # run function, retain all loci
    contrib_est$true_N <- i
    contrib_est_remove <- MixtureLikelihood2_remove(y=50, loci=loci.list) # run function, remove failed loci
    contrib_est_remove$true_N <- i
    output <- data.frame(n_loci=m,
                         true_N=contrib_est$true_N, 
                         est_N=contrib_est$Contributors,
                         est_N_remove=contrib_est_remove$Contributors,
                         n_failed=contrib_est_remove$n_failed,
                         est_N_allelic_rich=ceiling(max(lengths(loci.list))/2)) # number of alleles/2
    # failed_loci <- rbind(failed_loci, contrib_est_remove)
    simulated_mixtures <- rbind(simulated_mixtures, output)
    print(tail(simulated_mixtures))
  }
  }
  write.csv(simulated_mixtures, paste0("/Users/kbja10/Documents/Cornell/Research/Round goby/Intraspecific_variation_with_eDNA/datasets/simulated_mixtures_",names(sim_freqs)[t],".csv", sep=""), row.names=FALSE)
  write.csv(failed_loci, paste0("/Users/kbja10/Documents/Cornell/Research/Round goby/Intraspecific_variation_with_eDNA/datasets/failed_loci_",names(sim_freqs)[t],".csv", sep=""), row.names=FALSE)
  write.csv(msat_freqs[t], paste0("/Users/kbja10/Documents/Cornell/Research/Round goby/Intraspecific_variation_with_eDNA/datasets/freqs_",names(sim_freqs)[t],".csv", sep=""), row.names=FALSE)
  t <- t+1
}

# Microsatellite plots
simulated_mixtures_1 <- read.csv(paste0("/Users/kbja10/Documents/Cornell/Research/Round goby/Intraspecific_variation_with_eDNA/datasets/simulated_mixtures_p05.csv", sep=""))
simulated_mixtures_2 <- read.csv(paste0("/Users/kbja10/Documents/Cornell/Research/Round goby/Intraspecific_variation_with_eDNA/datasets/simulated_mixtures_p10.csv", sep=""))
simulated_mixtures_3 <- read.csv(paste0("/Users/kbja10/Documents/Cornell/Research/Round goby/Intraspecific_variation_with_eDNA/datasets/simulated_mixtures_p20.csv", sep=""))
simulated_mixtures_4 <- read.csv(paste0("/Users/kbja10/Documents/Cornell/Research/Round goby/Intraspecific_variation_with_eDNA/datasets/simulated_mixtures_p001.csv", sep=""))
simulated_mixtures <- rbind(simulated_mixtures_1, simulated_mixtures_2,
                            simulated_mixtures_3, simulated_mixtures_4)
simulated_mixtures$p <- factor(rep(c("p05","p10","p20","p001"), each=1000), levels=c("p05","p10","p20","p001"))
simulated_mixtures$est_N[simulated_mixtures$est_N==50] <- NA # if reaches max estimate, assume mixture failed
simulated_mixtures$bias <- simulated_mixtures$est_N-simulated_mixtures$true_N  
ggplot(simulated_mixtures) + 
    geom_point(aes(x=true_N, y=bias), colour=alpha("gray",0.3), size=2) +
    stat_summary(data=simulated_mixtures, 
                 aes(x=true_N, y=bias), geom="point", fun="mean", col="black", size=2, shape=19) +
    stat_summary(data=simulated_mixtures,
                 aes(x=true_N, y=bias), fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", color="black", width=0.2) +
    ylab("Bias (Estimate - True number of individuals)") + xlab("True number of individuals") +
    ylim(-10,10) + scale_x_discrete(limits=c(1,3,5,7,10)) +
    geom_abline(slope=0, intercept=0) +
    theme_bw() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=18)) +
    theme(legend.position = c(0.13, 0.88)) +
    labs(color = "Number of \n individuals") +
    guides(colour=guide_legend(override.aes = list(alpha = 1))) +
    theme(legend.background = element_rect(colour = "black", size=0.3)) +
    facet_grid(vars(p), vars(n_loci))

ggsave(filename = paste("/Users/kbja10/Documents/Cornell/Research/Round goby/Intraspecific_variation_with_eDNA/Figure_3_simulated_mixtures_msats_5.12.22.pdf"), plot=arrange_plots, dpi=300, width=12, height=6, units="in")



