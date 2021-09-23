## original author: D. R. Mani
## further modified by KK

library (pacman)
p_load (cmapR, dplyr, tidyr, gPCA, glue)

source('https://raw.githubusercontent.com/karstenkrug/R-code/main/my_plots.r')

check.batch.effect <- function (f, ## path to GCT v1.3 file
                                batch.var='Experiment', ## batch variable in cdesc object
                                pval.thresh=0.05
                                ) {
  d <- parse_gctx (f)
  dx <- d@mat %>% data.frame %>% drop_na
  # use gPCA to check for batch effect
  # see Reese et. al., 2013. Bioinformatics 29(22):2877-2883
  b <- gPCA.batchdetect (t (dx), batch=d@cdesc[, batch.var])
  print ( glue ("\n\nBatch Effect {ifelse (b$p.val < pval.thresh, 'PRESENT', 'ABSENT')} for {f}:\n  delta={round (b$delta, digits=3)}, p-value={b$p.val}\n\n"))
}

####################################
## perform test on first principal
## component using the batch as class
## variable
myPCA_batchdetect <- function(
                  m, 
                  b,
                  label='test',
                  scale=T
                  ){
  m <- m %>% data.frame %>% drop_na
  if(scale) m <- scale(m)

  pca <- prcomp(t(m))
  
  pc1 <- pca$x[,1]
  
  dat_plot <- data.frame(PC1=pc1, batch=b)
  
  b_levels <- unique(b)
  comps <- list()
  cc <- 1
  for(i in 1:(length(b_levels)-1))
    for(j in (i+1):length(b_levels)){
      comps[[cc]] <- c(i, j)
      cc <- cc + 1
    }
  comps=NULL
  p <- ggboxplot(dat_plot, x = 'batch', y='PC1',add = 'jitter') +
   ggtitle(label) +  
   stat_compare_means(comparisons=comps)
  
  plot(p)
  
  
  res <- kruskal.test(dat_plot$PC1, g=dat_plot$batch)
  return(res)
}

##################################################################
##
## use gPCA to check for batch effect
## see Reese et. al., 2013. Bioinformatics 29(22):2877-2883 
##
check_batch_effect <- function (f, ## path to GCT v1.3 file
                                batch.var=c('Experiment'), ## batch variable in cdesc object
                                method=c('gPCA', 'myPCA'),
                                plot=T,
                                label='', ## label to include in the plot
                                filt=NULL ## (optional) the number of features to retain after applying a variance filter. 
) {
  method <- match.arg(method)
  
  d <- parse_gctx (f)
  dx <- d@mat %>% data.frame %>% drop_na
  
  b_res <- lapply(batch.var, function(v){
    label2 <- paste(v,label)
    cat(label2, ' ----- \n')
    
    # use gPCA to check for batch effect
    # see Reese et. al., 2013. Bioinformatics 29(22):2877-2883
    if(method == 'gPCA'){
      
      ## batch needs to be numeric
      batch_num <- d@cdesc[, v] %>% as.factor %>% as.numeric
      res <- gPCA.batchdetect (t (dx), batch=batch_num, filt=filt)
      if(plot){
        par(mfrow=c(1,2))
        PCplot(res, ug='guided', main=paste0(label2, ' (guided PCA)'))
        PCplot(res, ug='unguided', main=paste0(label2, ' (regular PCA)'))
      }
    }
    if(method == 'myPCA'){
      res <- myPCA_batchdetect(dx, d@cdesc[, v], label=label2)
    }
    res
  })
  names(b_res) <- batch.var
  
  ########################
  ## summary table
  if(method == 'gPCA'){
    summary_tab <- lapply( b_res, function(x) x[c( "p.val", "delta", "n", "p", "b", "varPCu1", "varPCg1", "nperm")])
  }
  if(method == 'myPCA'){
    summary_tab <- lapply( b_res, function(x) x[c( "p.value", "statistic", "method")])
  }
  summary_tab <- Reduce(rbind, summary_tab)
  rownames(summary_tab) <- batch.var
  
  
  return(list(summary=summary_tab, all=b_res))
  }
