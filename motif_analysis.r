
#################################################################################################################
##
##                             aa composition of motifs
##
## - ph       - phosphosite table
## - column   - character, the column name of the sequence window
## - items    - all possible items that can be present in 'column'
##            - usually all amino acids, but also for example classes of amino acids
##
##
## changelog:   20100415 parameter column and items
##              20110307 estimation of relative frequencies:
##                          - absolute numbers diveded by number of sites
##                            (and not by the sum of all aa at that position)
##              20160216 adapted for usage at the Broad Institute
##                       parameter'use.all': - boolean, if multiple sequence windows are present
##                                           - SM doesn't report sites but (multiply phosphoryalted) peptides
#################################################################################################################
motif.comp <- function( ph, use.all=F, column="Sequence.Window", items="ACDEFGHIKLMNPQRSTVWY" ){

    ## sequence window of the site
    seq.win <- as.character(ph[, column])

    ## deteremin length of sequence window
    seqwinL <- unique(nchar(as.character( unlist( strsplit(seq.win[1], ';|\\|'))  ) ))

    ## make the list unique, i.e. use the first one
    if(!use.all)
        seq.win <- sub("(;|\\|).*", "", seq.win)
    if(use.all)
       seq.win <- unlist(strsplit(seq.win, ';|\\|' ))

    Nsites=length(seq.win)

    ## split
    seq.win.split <- strsplit(seq.win, "")

    ## matrix
    ## seq.win.split.mat <- matrix(unlist(seq.win.split), ncol=length( seq.win.split[[1]] ), byrow=T, dimnames=list( rownames(ph), c("-6", "-5", "-4", "-3", "-2", "-1", "0", "+1", "+2", "+3", "+4", "+5", "+6")  )  )
     seq.win.split.mat <- matrix(unlist(seq.win.split), ncol=seqwinL, byrow=T, dimnames=list( 1:length(seq.win.split), -((seqwinL-1)/2):((seqwinL-1)/2)  )  )

    ##############################
    ## items to count
    ###############################
    items <- unlist(strsplit(items, ""))

    ##############################
    ## motif aa composition matrix
    ##############################
    motif.comp <- apply(seq.win.split.mat, 2, table)
    motif.comp.mat <- motif.comp.mat.rel <- matrix( 0, ncol=length(items), nrow=length( seq.win.split[[1]] ), dimnames=list( names(motif.comp), items  ) )

      ################################
      # now fill the matrix
       ###############################
      for(i in names(motif.comp)){

         pos <- motif.comp[[i]]

         for(j in setdiff(names(pos), c("_", "X"))){
            motif.comp.mat[i,j] <- pos[j]
         }
       }
       # relative occurence
       for(i in 1:dim(motif.comp.mat)[1])
           #motif.comp.mat.rel[i,] <- as.numeric(motif.comp.mat[i,])/sum(as.numeric(motif.comp.mat[i,]), na.rm=T)
           motif.comp.mat.rel[i,] <- as.numeric(motif.comp.mat[i,])/Nsites

       # output, absolute and relative counts
       out <- vector("list", 2)
       names(out) <- c("count", "percent")

       out[[1]] <- motif.comp.mat
       out[[2]] <- motif.comp.mat.rel

       return( out )
}

##########################################################################
#
#     hyper.test - hypergeometrical test (Fisher's exact test)
#
#                        drawn
#                       y     n
#                     ----------
#            white    | x  |    |  m
#            black    |    |    |  n
#                     ----------
#                       k
#
#   x   - the number of white balls drawn without replacement from the urn
#         (the number of successes)
#   m   - the number of white balls in the urn
#   n   - the number of black balls in the urn (N - m)
#   k   - the number of balls drawn from the urn
#
# changelog:   20100415 imlementation
##########################################################################
hyper.test <- function(x, m, n, k, alternative=c("two.sided", "less", "greater")){

    alternative <- match.arg(alternative)

    ################################
    # one sided: depletion
    ################################
    if(alternative == "less")
        p <- sum( dhyper(0:x, m, n, k)  )

    ################################
    # one sided: enrichment
    ################################
    if(alternative == "greater")
        p <- sum( dhyper(x:m, m, n, k)  )

    ################################
    # two sided
    ################################
    if(alternative == "two.sided"){

         p.x <- dhyper(x, m, n, k)
         p.all <- dhyper(0:m, m, n, k)

         p <- sum(p.all[ p.all <= p.x  ])
    }


    return(p)
}

######################################################################################################################################
##
##
##
##  test.sites    - part of phospho table containing the test sites
##  bg.sites      - part of phospho table containing the background dataset
##                - or a character string specifying the path to a protein database in fasta format
##  p.motif       - numeric, p-value defining 'significant'  enrichment/depletion
##  rm.zero       - logical, if TRUE the modification position will be removed from the matrix
##  dist          - character, specifies the type of distribution used for testing
##                - can either be binomial or hypergeometric
##  main          - character, title for the heatmap
##  alternative   - character, specifies the direction of the test
##                - any other value that 'greater' has only an effect on the binomial test
##  col           - color scheme to use, e.g. "redgreen", "greenred", "heatcolors"
##  p.min         - numeric, minimal p-value, i.e. p-values that are below this value will be set to this value before logarithmize the
##                  value. this is only done to increase the 'dynamic range' of the heat map and has no effect on the real p-values
##                  that are returned by this function!
##  p.max         - numeric, this value determines the color scale of the heat map, i.e. if the most significant p-value in the data is 1e-2
##                  the colors are still scaled from -p.max to p.max
##                - this shall avoid any overinterpretation of not very significant p-values only because they jump out of the matrix...
##  column        - character, the column name holding the sequence motif
##  items         - character, all possible items that can occur in the column
##  cn            - logical, if TRUE the p values are shown as cell notes in the heatmap
##
##  ...           - further arguments passed to 'binomial.test'
##
## changelog:    20100323 - reimplementation based on the original version
##               20100324 - some doumentation
##               20100325 - added a legend that summarizes the test parameters
##               20100401 - the color scale is now always symmetric
##                        - parameter 'pmax'
##                        - swapped the matrix
##               20100415 - added the calculation of 'two.sided' and 'less' for the hypergeometric test
##                        - parameter column and items
##               20100422 - changed parameter 'col'
##               20100504 - added parameter 'cellnote'
##               20110307 - frequencies can now be compared against overall amino acid frequencies in the protein
##                          database using the binomial test
##               20110919 - anti-motif is spelled out in the figure legend
##               20160216 - changed visuale style of the heatmap
##               20160217 - parameter 'legend' to change info given in legend
##
#####################################################################################################################################

motif.enrich.new <- function(test.sites, bg.sites, p.motif = 1e-6, rm.zero=F, dist=c("binomial", "hypergeometric"), main="", alternative=c("greater", "less", "two.sided"), col=c("darkgreen", "white", "darkred"), p.min=1e-10, p.max=1e-5, column="Sequence.Window", items="AILVFWYNCQMSTDERHKGP", cellnote=T, keysize=.8, legend=c('all', 'compact', 'none'), ...){

    ## legend style
    legend <- match.arg(legend)

    ##########################################################
    # Test set:
    # count the amino acids occuring at each offset position
    test <- motif.comp(test.sites, column=column, items=items)


    ################################################
    # if 'bg.sites' is a character string, assume
    # it is the filename of a fasta file
    ################################################
    if(length(bg.sites) == 1){
        require(seqinr)

        # import fasta
        db.bg <- read.fasta(bg.sites, seqtype="AA")

        # determine the overall freuencies of all amino acids
        db.bg <- unlist(db.bg)
        bg.aa <- table(db.bg)

        # background matrix
        bg <- test
        for(aa in names(bg.aa)){
            # absolute counts
            bg[[1]][, aa] <- rep(bg.aa[aa], dim(bg[[1]])[1] )

            # relative frquencies
            bg[[2]][, aa] <- rep(bg.aa[aa]/sum(bg.aa), dim(bg[[1]])[1] )
        }

        # set the test to 'binomial'
        dist="binomial"
    } else{
        bg <- motif.comp(bg.sites, column=column, items=items)
    }

    # remove the the offset zero position, i.e. the site itself
    if(rm.zero){
        z.idx <- grep("^0$", rownames(test[[1]]))
        if(length(z.idx) > 0){
          test[[1]] <- test[[1]][-z.idx, ]
          test[[2]] <- test[[2]][-z.idx, ]

          bg[[1]] <- bg[[1]][-z.idx, ]
          bg[[2]] <- bg[[2]][-z.idx, ]
        }
    }

    # matrix storing the p-values
    p.mat <- matrix(NA, nrow=dim(test[[1]])[1], ncol=dim(test[[1]])[2], dimnames=list( rownames(test[[1]]), colnames(test[[1]]) ) )

    #################################
    # determine the test distribution
    #################################
    dist <- match.arg(dist)

    alternative <- match.arg(alternative)

    if(dist == "binomial" ){
       # list storing the whole output from binom.test
       binom.list <- vector("list", dim(test[[1]])[1])
       names(binom.list) <- rownames(test[[1]])
    } else {
        binom.list = NULL
    }

    ############################################################
    # binomial/hypergeometric test of each aa at each position
    ############################################################
    # loop over the offset positions
    for(i in rownames(test[[1]])){

        if(dist == "binomial"){
           binom.list[[i]] <- vector("list", dim(test[[1]])[2])
           names(binom.list[[i]]) <- colnames(test[[1]])
        }

        # loop over the amino acids
        for(j in colnames(test[[1]])){

            # binomial test
            if(dist == "binomial"){

                 # parameter
                 n.success <- test[[1]][i,j]                # number of successe
                 n.draws <- dim(test.sites)[1]              # number of draws
                 #p.bg <- bg[[1]][i,j]/dim(bg.sites)[1]     # hypothesized probality, i.e. relative occurence of the aa at this position in the bg dataset
                 p.bg <- bg[[2]][i,j]

                 # test
                 binom.list[[i]][[j]] <- binom.test( n.success, n.draws, p.bg, alternative=alternative, ...  )
                 # p-value
                 p.mat[i, j] <- binom.list[[i]][[j]]$p.value
             }
                # hypergeometrical test
            if(dist == "hypergeometric"){

                 # parameter
                 n.success <- test[[1]][i,j]   # number of with balls drawn
                 n.wb      <- bg[[1]][i, j]    # number of white balls ( number of occurences in the bg dataset  )
                 n.bb      <- dim(bg.sites)[1] - n.wb
                 #n.bb      <- sum(bg[[1]]["0", ]) - n.wb
                 n.draws   <- dim(test.sites)[1]


                 # test
                 if(alternative == "greater"){
                     p.mat[i, j] <- sum( dhyper( n.success:n.wb, n.wb, n.bb, n.draws  )  )
                 } else if(alternative == "less" ){
                     p.mat[i, j] <- sum( dhyper( 0:n.success, n.wb, n.bb, n.draws  )  )
                 } else if(alternative == "two.sided"){

                     # probability of success
                     p.success <- dhyper( n.success, n.wb, n.bb, n.draws  )

                     # all probabilities
                     p.hyper <- dhyper(0: n.wb, n.wb, n.bb, n.draws  )

                     # two sided p-value: any configuration that is equal or more extreme than observed
                     p.mat[i,j] <- sum( p.hyper[p.hyper <= p.success  ]  )
                 }
            }

        }
    }
    ##########################
    # cell notes
    ##########################
    cn <- matrix(paste(test[[1]], "/",  bg[[1]], "\n(", round(test[[2]], 2), "/",round(bg[[2]],2),")\np=", formatC(p.mat, digits=3),sep=""), nrow=dim(test[[1]])[1], byrow=F)


    #####################################
    # replace p-values of 0 by 1e-300
    # 1e-400 is already zero
    #####################################
    p.mat.plot <- p.mat
    p.mat.plot[p.mat.plot < p.min] <- p.min


    ##########################
    # the values to plot
    ##########################
    p.mat.log <- -log(p.mat.plot,10)
    if(alternative == "two.sided"){
        for(i in rownames(p.mat.log))
            for(j in colnames(p.mat.log))
                if( (test[[2]][i,j])  <  (bg[[2]][i,j])  ) p.mat.log[i,j] <- -1*p.mat.log[i,j]
                #if( (test[[1]][i,j]/dim(test.sites)[1])  <  (bg[[1]][i,j]/dim(bg.sites)[1])  ) p.mat.log[i,j] <- -1*p.mat.log[i,j]
    }

    ###################################################################################
    ##
    ## extract the possible motif / anti motif
    ##
    ##################################################################################
    motif.string <- "enriched: "
    motif.string.depl <- "depleted: "

       for(i in rownames(p.mat)){

           # get aa with p-value lower than 'p.motif'
           aa <- colnames(p.mat)[ which(p.mat[i, ] < p.motif & p.mat.log[i,] > 0) ]
           aa.depl <- colnames(p.mat)[ which(p.mat[i, ] < p.motif & p.mat.log[i,] < 0) ]

           #########################
           # enrichment
           #########################
           if(length(aa) == 0)
               aa.add <- "x-"
           if(length(aa) == 1)
               aa.add <- paste(aa, '-', sep='')
           if(length(aa) > 1 )
               aa.add <- paste("[",paste(aa, collapse="")  ,"]-", sep="")
           if(i != "0")
               motif.string <- paste( motif.string, aa.add, sep=""  )

           # add the site itself [STY]
           if(rm.zero){
              if(i == "-1")
                 motif.string <- paste( motif.string, "[ST]-", sep=""  )
           } else {
              if(i == "0")
                 motif.string <- paste( motif.string, "[ST]-", sep=""  )

           }
           rm(aa.add)
           #########################
           # depletion
           #########################
           if(length(aa.depl) == 0)
               aa.add <- "x-"
           if(length(aa.depl) == 1)
               aa.add <-paste(aa.depl,'-',sep='')
           if(length(aa.depl) > 1 )
               aa.add <- paste("[",paste(aa.depl, collapse="")  ,"]-", sep="")

           if(i != "0")
               motif.string.depl <- paste( motif.string.depl, aa.add, sep=""  )

           # add the site itself [STY]
           if(rm.zero){
              if(i == "-1")
                 motif.string.depl <- paste( motif.string.depl, "[ST]-", sep=""  )
           } else {
              if(i == "0")
                 motif.string.depl <- paste( motif.string.depl, "[ST]-", sep=""  )

           }

       }
    motif.string <- sub('-$','', motif.string)
    motif.string.depl <- sub('-$','', motif.string.depl)
    motif.string <- paste( paste(motif.string, sep=""), "\n",paste(motif.string.depl, sep=""), collapse="")


    #####################################
    ## 'breaks': mapping data to colors
    #####################################

    ## maximal absolute data value
    ##p.log.scale <- max( abs(p.mat.log) )

    p.log.scale = -log(p.min, 10)

    if(-log(p.max,10) > p.log.scale) p.log.scale = -log(p.max,10)

    p.log.scale <- min( p.log.scale, -log(p.min, 10)  )

    breaks.param <- seq( -p.log.scale, p.log.scale, 0.5  )

    #####################################
    # colors
    #####################################
    require(gplots)

    if(length(col) == 3)
        col <- colorpanel( length(breaks.param)-1, col[1], col[2], col[3]  )


    #####################################
    ## RowSideColors: AA properties
    if(items=="AILVFWYNCQMSTDERHKGP")
        RowSideColors=c(rep('grey40', 4), rep('chocolate4', 3), rep('dodgerblue4', 6), rep('green3', 2), rep('indianred3', 3), rep('black', 2) )

    #######################
    # plot
    #######################
    if(items=="AILVFWYNCQMSTDERHKGP") {

        ## with cellnotes
        if(cellnote){
            heatmap.2(  t(p.mat.log) , Rowv=F,Colv=F, trace="n", dendrogram="none", symkey=ifelse(alternative=="two.sided",T,F), cellnote=t(cn), notecex=.5, notecol="black", sub=paste(motif.string, "\n(p < ", p.motif,")"), main=paste( main, "(", dim(test.sites)[1],"sites )" ), col=col, breaks=breaks.param, colsep=1:ncol(t(p.mat.log)), rowsep=1:nrow(t(p.mat.log)), sepcolor='grey',  density.info="none", key.title='p-value', key.xlab='log10(P-value)', keysize=keysize, RowSideColors=RowSideColors)

        } else { # without cellnotes
            heatmap.2(  t(p.mat.log) , Rowv=F,Colv=F, trace="n", dendrogram="none", symkey=ifelse(alternative=="two.sided",T,F), sub=paste(motif.string, "\n(p < ", p.motif,")"), main=paste( main, "(", dim(test.sites)[1],"sites )" ), col=col, breaks=breaks.param , colsep=1:ncol(t(p.mat.log)), rowsep=1:nrow(t(p.mat.log)), sepcolor='grey',  density.info="none", key.title='p-value', key.xlab='log10(P-value)', keysize=keysize, RowSideColors=RowSideColors)
        }

        ###########################################
        ## add legend
        if(legend=='all')
            legend("top", legend=c(paste("background:", dim(bg.sites)[1], "sites"), paste("minimal p-value used for plotting: ", p.min), paste("test:", dist), paste("alternative:", alternative)), bty="n")
        if(legend=='compact')
            legend("top", legend=c(paste("background:", dim(bg.sites)[1], "sites"), paste("test:", dist, '/', alternative)), bty="n")

        ###########################################
        ## add AA anno
        mtext('Hydrophobic\naliphatic', side=2, line=3, at=.82, adj=0, las=2, cex=.8, col='grey40')
        mtext('Hydrophobic\naromatic', side=2, line=3, at=.65, adj=0, las=2, cex=.8, col='brown')
        mtext('Polar\nneutral', side=2, line=3, at=.44, adj=0, las=2, cex=.8, col='darkblue')
        mtext('Charged\nacidic', side=2, line=3, at=.24, adj=0, las=2, cex=.8, col='darkgreen')
        mtext('Charged\nbasic', side=2, line=3, at=.11, adj=0, las=2, cex=.8, col='darkred')
        mtext('Unique', side=2, line=3, at=0, adj=0, las=2, cex=.8, col='black')


    } else { ## no AA annotations
        if(cellnote){
            ## cellnotes
            heatmap.2(  t(p.mat.log) , Rowv=F,Colv=F, trace="n", dendrogram="none", symkey=ifelse(alternative=="two.sided",T,F), cellnote=t(cn), notecex=.5, notecol="black", sub=paste(motif.string, "\n(p < ", p.motif,")"), main=paste( main, "(", dim(test.sites)[1],"sites )" ), col=col, breaks=breaks.param, colsep=1:ncol(t(p.mat.log)), rowsep=1:nrow(t(p.mat.log)), sepcolor='grey',  density.info="none", key.title='p-value', key.xlab='log10(P-value)', keysize=keysize, RowSideColors=RowSideColors)
        } else { # without cellnotes

            heatmap.2(  t(p.mat.log) , Rowv=F,Colv=F, trace="n", dendrogram="none", symkey=ifelse(alternative=="two.sided",T,F), sub=paste(motif.string, "\n(p < ", p.motif,")"), main=paste( main, "(", dim(test.sites)[1],"sites )" ), col=col, breaks=breaks.param , colsep=1:ncol(t(p.mat.log)), rowsep=1:nrow(t(p.mat.log)), sepcolor='grey',  density.info="none", key.title='p-value', key.xlab='log10(P-value)', keysize=keysize)
        }
        if(legend=='all')
            legend("top", legend=c(paste("background:", dim(bg.sites)[1], "sites"), paste("minimal p-value used for plotting: ", p.min), paste("test:", dist), paste("alternative:", alternative)), bty="n")
        if(legend=='compact')
            legend("top", legend=c(paste("background:", dim(bg.sites)[1], "sites"), paste("test:", dist, '/', alternative)), bty="n")


    }


    ########################
    # output
    ########################
    out <- list()
    out[[1]] <- p.mat
    out[[2]] <- p.mat.log
    out[[3]] <- binom.list
    out[[4]] <- motif.string

    names(out) <- c("p-values", "p-plot", "binom.test", "motif")

    return(out)
}
