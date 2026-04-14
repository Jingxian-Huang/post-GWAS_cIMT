library("TwoSampleMR")
library("ieugwasr")
library("data.table")
library("plyr")
library("dplyr")
library("genetics.binaRies")
library("tzdb")
library("readr")
library("withr")
library("vroom")
library("labeling")
library("farver")
library("ggplot2")
library("stringr")
library("coloc")

# cis-MR
# Read directories of the outcome data
directory <- ""
files.out <- list.files(directory, pattern = "\\.txt$", full.names = T)

# Read directories of the exposure data
directory <- ""
files.expo <- list.files(directory, pattern = "\\.txt$", full.names = T)

# Loop through each outcome
for (file.out in files.out) {
  # Read outcome
  out <- read_outcome_data(file.out)
  name.out <- gsub(".txt", "", basename(file.out))
  out$outcome <- name.out
  
  res <- list()
  # Loop through each exposure
  for (file.expo in files.expo) {
    # Read exposure
    expo <- fread(file.expo)
    name.expo <- str_split(basename(file.expo), "_")[[1]][1]
    colnames(expo)[2:10] <- paste(colnames(expo)[2:10], ".exposure", sep = "")
    expo <-
      expo %>% mutate(
        exposure = name.expo,
        mr_keep.exposure = TRUE,
        pval_origin.exposure = "reported",
        id.exposure = name.expo,
        data_source.exposure = "textfile"
      )
    
    # Select significant SNPs
    expo$pval.exposure <- as.numeric(expo$pval.exposure)
    expo <- filter(expo, pval.exposure < 5e-8)
    # Exclude rare variants
    expo <- expo[expo$eaf.exposure > 0.01 &
                   expo$eaf.exposure < 0.99, ]
    # Exclude weak IVs by calculating F-statistics
    expo$F <- ((expo$beta.exposure) ^ 2) / ((expo$se.exposure) ^ 2)
    expo <- expo[expo$F >= 10, ]
    expo <- expo[, -16]
    
    # harmonize
    dat <- harmonise_data(expo, out)
    dat <- dat[dat$mr_keep == T, ]
    
    if (nrow(dat) > 0) {
      # clumping
      dat.clump <- mutate(dat, rsid = SNP, pval = pval.exposure)
      dat.clump <- ld_clump(
        dat.clump,
        clump_kb = 10000,
        clump_r2 = 0.001,
        plink_bin = genetics.binaRies::get_plink_binary(),
        bfile = "UKB/UKB_ref_panel"
      )
      
      dat.clump <- dat.clump[, -35:-37]
      
      # Steiger filtering
      dat.fil <- steiger_filtering(dat.clump)
      dat.fil <- dat.fil[dat.fil$steiger_dir == T, 1:34]
      
      dat.clump <- dat.fill
      
      # MR
      if (nrow(dat.clump) >= 3) {
        res.temp <- mr(
          dat.clump,
          method_list = c(
            "mr_ivw",
            "mr_egger_regression",
            "mr_weighted_median"
          )
        )
        res.plt <- mr_pleiotropy_test(dat.clump)
        res.het <- mr_heterogeneity(dat.clump, method_list = c("mr_ivw"))
        
        # Combine results
        res.temp <- mutate(
          res.temp,
          egger_intercept = res.plt$egger_intercept,
          egger_intercept_pval = res.plt$pval,
          Q = res.het$Q,
          Q_df = res.het$Q_df,
          Q_pval = res.het$Q_pval
        )
        
      } else if (nrow(dat.clump) >= 2) {
        res.temp <- mr(dat.clump, method_list = c("mr_ivw"))
        
        # Combine results
        res.temp <- mutate(
          res.temp,
          egger_intercept = NA,
          egger_intercept_pval = NA,
          Q = NA,
          Q_df = NA,
          Q_pval = NA
        )
        
      } else{
        res.temp <- mr(dat.clump, method_list = c("mr_wald_ratio"))
        
        # Combine results
        res.temp <- mutate(
          res.temp,
          egger_intercept = NA,
          egger_intercept_pval = NA,
          Q = NA,
          Q_df = NA,
          Q_pval = NA
        )
        
      }
      
      # Save IVs
      write_csv(dat.clump,
                paste0("iv/", name.expo, "_", name.out, ".csv"))
      
      # Bind results
      res <- rbind(res, res.temp)
      
    }
    
  }
  
  # Save results
  write_csv(res, paste0("~/", name.out, ".csv"))
  
}


# COLOC
res.all <- list()
for (file in list.files("out/")) {
  
  dt1 <- read_delim(file)
  dupp <- dt1[duplicated(dt1$SNP), ]
  dt1 <- dt1[(dt1$SNP %in% unique(dupp$SNP)) == F, ]
  dt1 <- dt1[dt1$A1FREQ >= 0.01 & dt1$A1FREQ <= 99.99, ]
  
  file1 <- list(
    beta = dt1$BETA,
    varbeta = (dt1$SE)^2,
    pvalues = dt1$pval,
    N = dt1$N,
    snp = dt1$SNP,
    MAF = dt1$A1FREQ,
    type = "quant"
  )
  
  dt2 <- read_delim(paste0("cis_pqtl/", file))
  dupp <- dt2[duplicated(dt2$SNP), ]
  dt2 <- dt2[(dt2$SNP %in% unique(dupp$SNP)) == F, ]
  dt2 <- dt2[dt2$eaf >= 0.01 & dt2$eaf <= 99.99, ]
  
  file2 <- list(
    beta = dt2$beta,
    varbeta = (dt2$se)^2,
    pvalues = dt2$pval,
    N = dt2$samplesize,
    snp = dt2$SNP,
    MAF = dt2$eaf,
    type = "quant"
  )
  
  res = coloc.abf(file1, file2)
  res <- data.frame(t(res$summary))
  
  res$dat1 <- "cimt"
  res$dat2 <- str_split(file, "_")[[1]][1]
  
  res.all <- rbind(res.all, res)

}

write_csv(res.all, "/coloc/res.cimt.csv")