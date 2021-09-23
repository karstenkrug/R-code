
## ################################################################
## - create site centric table from Spectrum Mill VM reports
##
## - tab    : table in Spectrum Mill VM export format
## - id.col : header id of column containing SM site ids, e.g.
## - mod    : string appended to the site ids, e.g. 
##                -p for phospho
##                -ub for ubiquitylation
##                -ac for lysine-acetylation
##
smVM2sites <- function(tab, id.col='id', acc.type='uniprot'){
  
  mod='-p'
  
  ## just in case
  rownames(tab) <- make.names(rownames(tab))
  
  ## #################################################
  ## SM site ids
  site.ids <- gsub(' ','', tab[, id.col])
  
  # In case of RefSeq database, remoove all non-RefSeq accession numbers,
  # e.g. UniProt accessions of lab contaminants
  if(acc.type == 'refseq'){
    idx <- grep('^(NP_|XP_|YP_|ZP_)', site.ids)
    site.ids <- site.ids[idx]
    tab <- tab[idx, ]
  }

  ## index of fully localized sites 
  loc.idx <- which( sapply (strsplit(sub('^.*_([1-9]_[0-9])_.*', '\\1', site.ids), '_'), function(x)length(unique(x)) ) == 1)
  
  ## update table and ids
  tab <- tab[loc.idx, ]
  site.ids <- site.ids[ loc.idx ]
  
  ## the code below works for UniProt accession numnbers. Leave as is for the time being
  ## and merge at a later point
  if(acc.type == 'uniprot'){
    
    ## variable sites
    site.var <- sub('^(.*?_.*?)_.*', '\\1', site.ids)
  
    ## accession number plus modified residue
    all.sites <- lapply( strsplit(site.var, '(s|t|y)'), function(x){ 
      prot.acc=sub('_.*', '', x[1])
      x=sub('.*?_', '', x)
      paste(prot.acc, grep( '(S|T|Y)', x, value=T), sep=';' )
    
    })
    names(all.sites) <- rownames(tab)
  }
  ## RefSeq accesion numbers
  if(acc.type == 'refseq'){
    
    ## variable sites
    site.var <- sub('^(.*?_.*?_.*?)_.*', '\\1', site.ids)
    
    ## accession number plus modified residue
    all.sites <- lapply( strsplit(site.var, '(s|t|y)'), function(x){ 
      prot.acc=sub('_[S|T|Y].*', '', x[1])
      x=sub('^.*_([S|T|Y].*)$', '\\1', x)
      paste(prot.acc, grep( '(S|T|Y)', x, value=T), sep=';' )
      
    })
    names(all.sites) <- rownames(tab)
  }  
  
  ## Exclude sites on multiply phosphorylated peptides for which we
  ## also detect singly phosphorylated version
  n.sites <- sapply(all.sites, length)
  ## separate singly and multiply phosphorylated peptides
  all.sites.mult <- all.sites[which(n.sites > 1)]
  all.sites.sing <- all.sites[which(n.sites == 1)]
  all.sites.mult <- lapply(all.sites.mult, function(x){
    rm.idx=which(x %in% unlist(all.sites.sing))
    if(length(rm.idx)>0)
      x=x[-rm.idx]
    x
    })
  ## all sites as list
  all.sites <- append(all.sites.mult, all.sites.sing)
  
  ## all sites as string, for faster searches  
  all.sites.str <- sapply(all.sites, paste, collapse='|')
  names(all.sites.str) <- names(all.sites)
  
  ## unique sites
  all.sites.unique <- unique( unlist(all.sites) )
 
  ## site centric table
  tab.site.centric <- matrix(NA, length(all.sites.unique), ncol=ncol(tab), dimnames=list( all.sites.unique, colnames(tab)))
  
  count <- 0
  prev <- 0 ## to monitpr progress
  N <- length(all.sites.unique)
  
  cat('Assembling', N, 'unique sites...\n')
  
  ## loop over each site
  for(i in all.sites.unique){
      
      ## row indices in the input table that contains the site 
      ##row.idx <- names(all.sites)[ which( all.sites == i) ]
      #row.idx <- names(all.sites)[ grep(paste('^',i, '$', sep=''), all.sites) ]
      
     # row.idx <- names(all.sites)[ which(sapply(all.sites, function(x) ifelse(i %in% x, T, F))) ]
      
      row.idx <- names(all.sites)[ grep(paste('\\||^', i, '\\||$'), all.sites.str)]
      
      ## check whether there are multiple rows, i.e. the same site has been quantified on different
      ## multiply phosphorylated peptide versions
      #cat(i, length(row.idx), '\n')
      if(length(row.idx) > 1){
         
        #cat(row.idx, '\n')
        ## select the peptide version with lowest number of modifications
        numb.sites.tmp <- sapply(all.sites[row.idx], length)
        
        row.idx <- row.idx[ which.min(numb.sites.tmp) ] 
        
      }
     tab.site.centric[i, colnames(tab)] <- unlist( tab[row.idx, ] )
  
      # progress
     count <- count + 1
     curr <- round(count/N*100)
     if(curr %% 5 == 0 ){
       if(curr != prev){
        prev <- curr
        cat(curr, '% ')
       }
       
     }
     }
  ## rownames of orig. table
  #site.index=unlist( lapply(names(all.sites), function(x) rep(x, length(all.sites[[x]]))) )
  
  ## single site
  #single.sites <- unlist(all.sites)
  
  ## skeleton for site-centric table
 # tab.site.centric <- matrix(NA, length(site.index), ncol=ncol(tab), dimnames=list(make.unique(site.index), colnames(tab)))
  
  ## fill new table
  #for(i in 1:length(site.index))
  #  tab.site.centric[i, colnames(tab)] <- unlist( tab[site.index[i], ] )
  
  ## add single site description
  tab.site.centric <- data.frame(tab.site.centric, stringsAsFactors=F)
  
  ## ###############################################
  ## accession number
  #site.acc <- sub('^(.*?)_.*', '\\1', site.ids)
  
  #varSiteUP <-sub('^(.*?)_.*', '\\1', tab.site.centric[, id.col])
  
  ##########################################################
  ## site id: UniProt + position
  #site.id <- paste(paste(varSiteUP, single.sites, sep=';'), mod, sep='')
  #tab.site.centric <- data.frame(site_id=site.id, tab.site.centric, stringsAsFactors = F)
  tab.site.centric <- data.frame(site_id=paste( all.sites.unique, mod, sep=''), tab.site.centric, stringsAsFactors = F)
  cat('\n')
  return(tab.site.centric) 
}

## ##############################################################
## map PTM sites between organisms using PhosphositePlus
## site-group IDs
## ###############################################################
mapSites <- function(gct.str, # character, path to input GCT 1.3 file
                     psp.str='c:/Users/Karsten/Dropbox/Manuscripts/20170630_PTM-GSEA_manuscript/PSP_Oct2017/Phosphorylation_site_dataset', # character, path to input PSP table 
                     src.org='mouse',
                     dest.org='human',
                     ofile='site-centric-HUMANIZED'){
  require(pacman)
  p_load(cmapR)
  
  ## import tables
  gct <- parse.gctx(gct.str) 
  #tab <- read.delim(tab.str, stringsAsFactors = F, skip=2)
  psp <- read.delim(psp.str, stringsAsFactors = F, skip=2)
  
  ## separate psp
  psp.src <- psp[which(psp$ORGANISM == src.org), ]
  psp.dest <- psp[which(psp$ORGANISM == dest.org), ]
  
  ## extract ids
  #site.id <- tab[, id.col]
  site.id <- gct@rid
  site.acc <- sub('^(.*?);.*', '\\1', site.id)
  site.res <- sub('^.*?;(.*?)$', '\\1', site.id)
  
  ## vector to store the mapping results
  site.map <- data.frame( matrix('', nrow=length(site.res), ncol=3), stringsAsFactors=F)
  colnames(site.map) <-c(src.org, 'SITE_GRP_ID', dest.org)
  site.map[, src.org] <- paste(site.acc,site.res, sep=';')
  
  
  ## loop over site ids
  for(i in 1:length(site.id)){
    
      
    ## check protein accession
    if(site.acc[i] %in% psp$ACC_ID){
      
      ## extract all sites in PSP of that protein 
      psp.map <- psp.src[which(psp.src$ACC_ID == site.acc[i]), ]
      
      ## check modified residue
      if(site.res[i] %in% psp.map$MOD_RSD){
        
        psp.map <- psp.map[which(psp.map$MOD_RSD == site.res[i]), ]
        
        site.grp <- psp.map$SITE_GRP_ID
        
        ## map to destination organism
        if(site.grp %in% psp.dest$SITE_GRP_ID){
          psp.map.dest <- psp.dest[which(psp.dest$SITE_GRP_ID == site.grp), ]
          
          ## if there are multiple sites (isoforms)
          ## pick the canonical sequence
          if(nrow(psp.map.dest) > 1){
            idx <- grep(' iso', psp.map.dest$PROTEIN)
            if(length(idx) > 0 & length(idx) < nrow(psp.map.dest))
              psp.map.dest <- psp.map.dest[-idx,]
            psp.map.dest <- psp.map.dest[sample( nrow(psp.map.dest), 1),]
          }  
          site.map[i, dest.org] <- paste( psp.map.dest$ACC_ID,  psp.map.dest$MOD_RSD, sep=';')
          site.map[i, 'SITE_GRP_ID'] <- site.grp
          
        }
      }
    }
    
  }
  
  ## export mapping table
  #write.table(site.map, sep='\t', row.names = F, quote=F, file=paste(src.org, '-', dest.org, '_PSPmapped_sites.txt', sep=''))
  
  
  ## extract mapped sites
  site.map.mapped <- site.map[ which(nchar(site.map[, dest.org]) > 0), ]
  site.map.mapped <- site.map.mapped[ which(!duplicated(site.map.mapped[ , dest.org])), ]
  
  #map mapped sites to original table  
  idx <- match( site.map.mapped[, src.org], gct@rid)
  
  #update 
  mat <- matrix(gct@mat[idx, ])
  rid <- gct@rid[idx ]
  rdesc <- gct@rdesc[idx, ]
  
  #ids mapped to destination organisms
  site.id.dest <- site.map.mapped[ , dest.org]
  
  #update gct
  rdesc <- data.frame(site.map.mapped, rdesc)
  rownames(mat) <- rownames(rdesc) <- rid <- site.id.dest 
  gct@rid <- rid
  gct@mat <- mat
  gct@rdesc <- rdesc
  
  #export
  write.gct(gct, ofile=ofile)
  
  return(site.map)
}



## ############################################################
## - Convert between different scales.
## - Useful to make plots with two data axes
##  
## https://stackoverflow.com/questions/929103/convert-a-number-range-to-another-range-maintaining-ratio
convert.scales <- function(x, NewMin, NewMax, OldMin=NULL, OldMax=NULL){
  
  if(is.null(OldMin))
    OldMin <- min(x, na.rm = T)
  if(is.null(NewMin))
    OldMax <- max(x, na.rm = T)
  
  NewValue <- sapply(x, function(OldValue) (((OldValue - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)) + NewMin )
  
  #NewValue = (((OldValue - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)) + NewMin
  return(NewValue)
}

#################################################
##   Given a string and a number of characters
##   the function chops the string to the
##   specified number of characters and adds
##   '...' to the end.
## parameter
##   string     - character
##   nChar      - numeric
## value
##   string of 'nChar' characters followed
##     by '...'
##################################################
chopString <- function(string, nChar=10, add.dots=T)
{

    string.trim <- strtrim(string, nChar)

    if(add.dots)
        string.trim[ which(nchar(string) > nChar) ] <-  paste(string.trim[which(nchar(string) > nChar) ], '...')
    if(!add.dots)
        string.trim[ which(nchar(string) > nChar) ] <-  paste(string.trim[which(nchar(string) > nChar) ])

    return(string.trim)

}



################################################################################################
#                                significance A
#
#
#
# changelog:   20100420 implementation
#              20110531 finally changed to the formula without 0.5 ...
################################################################################################
significanceA <- function( r, log =T, log.base=2 ){

    ###################################################
    #  tranform to log scale
    ###################################################
    if(log)
        r <- log(r, log.base)

    quant.idx <- which(!is.na(r))
    ##################################################
    # robust and asymetrical estimate of the sd
    #
    ##################################################
    estimate.sd <- quantile(r[quant.idx], c(0.1587, 0.5, 0.8413))

    # the notation used in the MQ supplement
    r.minus1 <- estimate.sd[1]
    r.zero <- estimate.sd[2]
    r.plus1 <- estimate.sd[3]

    ##################################################
    # z tranformation
    ##################################################
    r.z <- ifelse(r >= r.zero, (r-r.zero)/(r.plus1-r.zero),  (r.zero-r)/(r.zero-r.minus1) )


    ##################################################
    # complementary error function
    ##################################################
    #sigA <- 0.5*erfc( r.z/sqrt(2) )
    sigA <- erfc( r.z/sqrt(2) )

return(sigA)
}


#############################################################
#                 Significance B
# r        - ratios
# i        - corresponding intensities
# nbin     - numeric, number of proteins in an intensity bin
#          - each bin contains at least that number of proteins
#
# changelog:  20100423 implementation
#             20110616 if length(r) < nbin sigA is applied
#############################################################
significanceB <- function(  r, i, nbin = 300, log=T, log.base=2  ){

    names(r) <- names(i) <- 0:(length(r)-1)

    # original order
    org.order <- names(r)

    ###################################################
    #  tranform to log scale
    ###################################################
    if(log)
        r <- log(r, log.base)

    ###################################################
    # order according to intensities
    ###################################################
    order.idx <- order(i, decreasing=T)
    i <- i[ order.idx ]
    r <- r[ order.idx ]

    ###################################################
    # index of all proteins having a ratio
    ###################################################
    quant.idx <- which(!is.na(r))

    ###################################################
    # get the ratios
    ###################################################
    r.quant <- r[quant.idx]
    i.quant <- i[quant.idx]

   ###################################################
    # bin the intensities
    ###################################################
    #bin <- round(seq(0, length(i[quant.idx]), length.out=ceiling(length(i[quant.idx])/nbin)  ))
    bin <- floor(seq(0, length(i.quant), length.out=ceiling(length(i.quant)/nbin)  ))
    #bin <- round(seq(0, length(i[quant.idx]), by=nbin  ))
    #bin <- round(seq(0, length(i), length.out=ceiling(length(i)/nbin)  ))


    ###################################################
    # loop over the bins and apply significance A
    ###################################################
    # vector to store the values
    sigB <- rep(NA, length(r))
    names(sigB) <- names(r)

    ###################################################
    #
    ###################################################
    if(length(r.quant)  <= nbin){

      sigB[ names(r.quant) ] <-  significanceA(r.quant, log=F)

      sigB <- sigB[ org.order ]

      return(list(sigB))
    }



    # store some informations on the bins
    info <-vector("list", length(bin)-1)
    names(info) <- paste("bin", 1:(length(bin)-1))

    # loop over the bins
    for(b in 1:(length(bin)-1)){

        # get the ratios in the current bin
        r.quant.bin <-  r.quant[ (bin[b]+1) : bin[b+1] ]

        # apply sigA
        sigB.tmp <- significanceA ( r.quant.bin, log=F)
        sigB[names(r.quant.bin)] <- sigB.tmp

        # some additional information
        info.tmp <- vector("list", 4)
        names(info.tmp ) <- c("proteins", "bin", "ratio", "intensity")
        info.tmp[[1]] <- length(r.quant.bin)
        info.tmp[[2]] <- names(r.quant.bin)
        info.tmp[[3]] <- r[names(r.quant.bin)]
        info.tmp[[4]] <- i[names(r.quant.bin)]

        info[[b]] <- info.tmp
    }

    ####################################################
    #   restore the original ordering
    ####################################################
    sigB <- sigB[org.order]

    out <- list()
    out[[1]] <- sigB
    out[[2]] <- info

    return(out)

}

#########################################################
#
#        complementary error function assuming
#        a normal distribution
#
#########################################################
erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)
