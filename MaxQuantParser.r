###################################################################################################
# Karsten Krug, 20090506
# 20111116   - outsourced all functions necessary to prodeuce the summary pdf file
#              to a speparate R file
#
#                                      MaxQuantParser.r
#
# This file contains all functions that will be included in the MaxQuant parser R package.
# A major advantage over the previous implementations in the introduction of C code to speed
# up some computations.
#
# heatmap.go
# getProteinState
# motif.comp
# hyper.test
# annotation.enrich
# motif.enrich.new
# sitesOnProteinPlot
# getNRSites
# siteOverlap
# singlePeptideFilter
# patternMatch2
# patternMatch
# compareExperiments
# splitMultipleProteinIDsInEvidence      - help function
# checkSeparation
# summary
# removeConRev
# extractRatios
# mergeDataSet
# timeProfileClusterPlot
# TimeProfilePlotTwoGroups
# getTopMembers
# topGOAnalysis
# ORA
# ORA.C
# Z.score
# getSILAClabel
# colorGradient
# chopString
# myScan
# DP.Mass.Difference
# DPmasshist
# vennDiag
# getMassAcc
# getAverageRT
#
# termPieChart
# fancyPlot
# fancyDensPlot
# getMW
#
# sigRatioPlot
# sigRatioPlot.ph
# sigRatioPlotLF
#
# GOresultPlot
# GOplot                                  - help function
# GOslim.analyis
# PAI
# observablePeptides
# spectralCount
#
# importBLAST
# collapseBLASTresults
# chromosomeMap
#
###################################################################################################
#source("c:/Users/karsten/Dropbox/Devel/MQparse/pdfscript.r")


##############################################################################
#
# parse result psipred result files
#
#
# changelog: 20120614 implementaion
##############################################################################
parse.psipred <- function(file){

   ss.tab <- read.table(file, stringsAsFactors=F, row.names=1)

   colnames(ss.tab) <- c( "aa", "ss", "prob.C", "prob.H", "prob.E")

   return(ss.tab)
}


################################################################################
#
# list - list of length N containing tables returned by function
#        annotation.enrich
#
# changelog: 20121017 implementation
#            20121218 changed default colors
################################################################################
heatmap.go <- function(list, p=0.05, pcol="p.adj", margins=c(10,25), col=c(rgb(0,178,43, maxColorValue=255) , rgb(255,252,171, maxColorValue=255), rgb(178,27,37, maxColorValue=255))
, pmin=1e-10, keysize=0.8, ...){

    dircol="direction"

    # store original tables, i.e. unfiltered
    list.org <- list

    # filter the tables
    for(l in 1:length(list))
        list[[l]] <- list[[l]][ which(list[[l]][, pcol] < p),  ]

    # get significant terms
    term.sig <- unique(unlist(lapply( list, function(x)rownames(x))))


    # minimal p-value used for plotting
    #p.min.plot <- abs(log(p.min, 10))

    # build a matrix
    gomat <- gomat.p <- matrix(NA, nrow=length(term.sig), ncol=length(list), dimnames=list(term.sig, names(list)))

    # fill the matrix
    for(tt in term.sig){
        for(l in 1:length(list)){

            if(tt %in% rownames(list.org[[l]])){

                # check whether the p-value exceeds the minimal p-value used for plotting
                p.tmp <- list.org[[l]][ tt, pcol ]
                if(p.tmp < pmin)
                    p.tmp = pmin

                # put to log
                gomat[tt, l] <- log(p.tmp , 10)*ifelse(list.org[[l]][tt, dircol] == "enrichment", -1, 1)

                # orginal p-value
                gomat.p[tt, l] <- list.org[[l]][tt, pcol]
            }
        }

    }


    require(gplots)
    ################################
    # colors
    p.scale = max(abs(gomat) )

    p.scale.breaks = seq(-p.scale, p.scale, 1)
    p.col <- colorpanel( length(p.scale.breaks)-1, col[1], col[2], col[3] )

    ################################
    # heatmap
    ################################
    heatmap.2( gomat , Colv=F, Rowv=T, dendrogram="row", trace="n", margins=margins, symkey=T, col=p.col, breaks=p.scale.breaks, keysize=keysize, density.info="none",...)

    return(gomat.p)
}


##########################################################################################################
#
#                                peak time index calculation
# - see Olsen et al. 2010
#
# arguments:
#  - data, numeric vector containing expression ratios on raw scale, i.e. no log
#
# changelog:  20120321 implementation
#
##########################################################################################################
peakTimeIndex <- function(data){

    r <- data/max(data)
    n <- length(data)
    imax <- which.max(r)

    if(imax >= 2 & imax <= (n-1))
        tpeak <- ( (imax-1)* r[imax-1] + imax*r[imax] + (imax+1)*r[imax+1]  )/(r[imax-1] + r[imax] + r[imax+1])

    if(imax == 1)
        tpeak <- (imax*r[imax] + (imax+1)*r[imax+1]  )/(r[imax] + r[imax+1] + r[n])

    if(imax == n)
        tpeak <- ( (imax-1)*r[imax-1] + imax*r[imax] + (n+1)*r[1]  )/( r[imax-1] + r[imax] + r[1]  )

    return(tpeak)

}

##########################################################################################################
#
#   given a protein sequence (imported with 'read.fasta' and as.string=F) the function
#   returns all possible sequence windows around the specified amino acids
#
# - prot.seq   - character vector, amino acid sequence
# - size       - numeric, size of the window, i.e. -size to +size
# - AA         - character, amino acids
#
# changelog: 20111208 implementation
##########################################################################################################
getSequenceWindows <- function(prot.seq, size=6, AA=c("S", "T", "Y")){

   # prot.seq <- toupper(prot.seq)
    #prot.seq <- split(prot.seq,"")

    seq.win.all <- c()

    # loop over specified aa
    for(aa in AA){

        # position
        idx <- grep(aa, prot.seq)

        # loop over the positions and get the sequence window
        if(length(idx) > 0){

            windows.tmp <- c()

            for(p in idx){

                seq.win.tmp <- paste( prot.seq[max((p-size),1):min((p+size), length(prot.seq))], collapse="")

                if( (p-size) <= 0 ){
                    seq.win.tmp <- paste( paste(rep("_", abs(p-size)+1), collapse=""), seq.win.tmp, sep="" )
                }
                if( (p+size > length(prot.seq)) ){
                        seq.win.tmp <- paste(  seq.win.tmp, paste(rep("_", (p+size)-length(prot.seq)), collapse=""), sep="" )
                }

                windows.tmp <- c(windows.tmp, seq.win.tmp)
            }

            seq.win.all <- c(seq.win.all, windows.tmp)

        }

    }

    return(seq.win.all)
}

############################################################################################################
#
# changelog: 20111209 implementation
#            20120613 added follwing motifs which I got from SIlke Hauf (S phase specific for S. pombe)
#                     ATR/ATM like kinases
#                     Cdc7/DDK
#
############################################################################################################
siteMotifs <- function(){

    kinases <- c("PKA", "CK1", "CK2", "GSK3", "CDK2", "CAMK2", "ERK/MAPK", "PKA/AKT", "PKC", "PKD", "LCK", "ABL", "SRC", "ALK", "EGFR", "CDK1", "Proline-directed", "AURORA", "AURORA-A", "PLK", "PLK1", "NEK6", "CHK1/2", "CHK1", "PDK1", "NIMA", "PIM1/2", "F box bTrCP", "14-3-3 binding", "WW GroupIV", "FHA1 Rad53p", "FHA2 Rad53p", "FHA KAPP", "ATR ATM", "CDC7 DDK")

    siteMotifs <- vector("list", length(kinases))
    names(siteMotifs) <- kinases

    ########################
    # PKA
    ########################
    siteMotifs[["PKA"]] <- c( "^....R.(S|T)......$",
                              "^....(R|K).(S|T)......$",
                              "^...R(R|K).(S|T)......$"
                             )
    ########################
    # CK1
    ########################
    siteMotifs[["CK1"]] <- c( "^...S..(S|T)......$",
                              "^..(S|T)...S......$"
                             )

    ########################
    # CK2
    ########################
    siteMotifs[["CK2"]] <- c( "^......(S|T)..(D|E)...$"
                             )

    ########################
    # GSK3
    ########################
    siteMotifs[["GSK3"]] <- c( "^......S...S..$"
                             )

    ########################
    # CDK2
    ########################
    siteMotifs[["CDK2"]] <- c( "^......(S|T)P.(K|R)...$"
                             )

    ########################
    # CAMPK2
    ########################
    siteMotifs[["CAMK2"]] <- c( "^...R..(S|T)......$",
                                "^...R..(S|T)V.....$"
                             )
    ########################
    # ERK/MAPK
    ########################
    siteMotifs[["ERK/MAPK"]] <- c( "^....P.(S|T)P.....$",
                                   "^....V.(S|T)P.....$",
                                   "^....PE(S|T)P.....$"
                             )

    ########################
    # PKA/AKT
    ########################
    siteMotifs[["PKA/AKT"]] <- c( "^...R(R|S|T).(S|T).(S|T)....$",
                                   "^.R.R..(S|T)......$"
                             )

    ########################
    # PKC
    ########################
    siteMotifs[["PKC"]] <- c( "^...R..(S|T).R....$"
                             )


    ########################
    # PKD
    ########################
    siteMotifs[["PKD"]] <- c( "^.(L|V|I).(R|K)..(S|T)......$"
                             )


    ########################
    # LCK
    ########################
    siteMotifs[["LCK"]] <- c( "^.....(I|E|V)Y(E|G)(E|D|P|N)(I|V|L)...$"
                             )

    ########################
    # ABL
    ########################
    siteMotifs[["ABL"]] <- c( "^.....(I|V|L)Y..(P|F)...$"
                             )

    ########################
    # SRC
    ########################
    siteMotifs[["SRC"]] <- c( "^...(E|D)..Y..(D|E|A|G|S|T)...$"
                             )

    ########################
    # ALK
    ########################
    siteMotifs[["ALK"]] <- c( "^......Y..(I|L|V|M)...$"
                             )

     ########################
    # EGFR
    ########################
    siteMotifs[["EGFR"]] <- c( "^....(D|P|S|A|E|N).Y(V|L|D|E|I|N|P).....$"
                             )

    ########################
    # CDK1
    ########################
    siteMotifs[["CDK1"]] <- c( "^......(S|T)P.(K|R)...$",
                               "^......(S|T)P(K|R)....$"
                             )

    ########################
    # Proline-directed
    ########################
    siteMotifs[["Proline-directed"]] <- c( "^......(S|T)P.....$"
                                          )

    ########################
    # AURORA
    ########################
    siteMotifs[["AURORA"]] <- c( "^....(R|K).(S|T)(I|L|V).....$"
                                )

    ########################
    # AURORA-A
    ########################
    siteMotifs[["AURORA-A"]] <- c( "^...(R|K|N)R.(S|T)(M|I|L|V).....$"
                                )

    ########################
    # PLK
    ########################
    siteMotifs[["PLK"]] <- c( "^....(D|E).(S|T)(V|I|L|M).(D|E)...$"
                                )

    ########################
    # PLK1
    ########################
    siteMotifs[["PLK1"]] <- c( "^....(D|E).(S|T)(F|L|I|Y|W|V|M).....$"
                                )

    ########################
    # NEK6
    ########################
    siteMotifs[["NEK6"]] <- c( "^...L..(S|T)......$"
                                )


    ########################
    # CHK1/2
    ########################
    siteMotifs[["CHK1/2"]] <- c( "^.L.R..(S|T)......$"
                                )

    ########################
    # CHK1
    ########################
    siteMotifs[["CHK1"]] <- c( "^.(M|I|L|V).(R|K)..(S|T)......$"
                                )

    ########################
    # PDK1
    ########################
    siteMotifs[["PDK1"]] <- c( "^..F..F(S|T)(F|Y).....$"
                                )

    ########################
    # NIMA
    ########################
    siteMotifs[["NIMA"]] <- c( "^...(F|L|M)(R|K)(R|K)(S|T)......$"
                                )

    ########################
    # PIM1/2
    ########################
    siteMotifs[["PIM1/2"]] <- c( "^.(R|K).RH.(S|T)......$"
    )

    ########################
    # F box bTrCP
    ########################
    siteMotifs[["F box bTrCP"]] <- c( "^.....D(S|T)G..(S|T)..$",
                                      "^.D(S|T)G..(S|T)......$"
    )

    ########################
    # 14-3-3 binding
    ########################
    siteMotifs[["14-3-3 binding"]] <- c( "^...RS.(S|T).P....$",
                                      "^..R.(Y|F).(S|T).P....$"
    )

    ########################
    # Polo box
    ########################
    siteMotifs[["Polo box"]] <- c( "^.....S(S|T)P.....$"
    )

    ########################
    # WW GroupIV
    ########################
    siteMotifs[["WW GroupIV"]] <- c( "^....P(P|A|V|L|I|M|C|F|Y|W|H|K|R|Q|N|E|D|S|T)(S|T)P.....$"
    )

    ########################
    # FHA1 Rad53p
    ########################
    siteMotifs[["FHA1 Rad53p"]] <- c( "^......T..D...$"
    )

    ########################
    # FHA2 Rad53p
    ########################
    siteMotifs[["FHA2 Rad53p"]] <- c( "^......T..(I|L)...$"
    )

    ########################
    # FHA KAPP
    ########################
    siteMotifs[["FHA KAPP"]] <- c( "^......T..(S|A)...$"
    )

    ########################
    # ATR ATM
    ########################
    siteMotifs[["ATR ATM"]] <- c( "^......[S|T]Q.....$"
    )

    ########################
    # CDC7 DDK
    ########################
    siteMotifs[["CDC7 DDK"]] <- c( "^...[E|S]S[D|E]S......$"
    )





    return(siteMotifs)
}

##########################################################################################################
#
#      given a MSMS or evidence table the function returns a non redundant version of
#      peptide sequences filtered according to 'col'
#
# tab
# col
# lower.better
#
# changelog: 20120329 implementation
#            20120410 'sequence' -> 'Sequence'
##########################################################################################################
table.nr <- function(tab, col="PEP", lower.better=T){

    if( "sequence" %in% colnames(tab) )
        colnames(tab)[ grep( "^sequence$", colnames(tab))  ] <- "Sequence"

    if( !(col %in% colnames(tab)))
        stop(paste("cannot find column name", col))

    index.tmp <- paste(rownames(tab), tab[, col], sep=";")


    index <- tapply( index.tmp, tab$Sequence, function(x){
                      xxx<-matrix( as.numeric(unlist( lapply(x, function(xx) unlist(strsplit(xx, ";") )) )), ncol=2, byrow=T);
	              return(ifelse( lower.better,  xxx[which.min(xxx[,2]), 1],xxx[which.max(xxx[,2]), 1]  ))
                      }
                    )
    index <- sub( " *$", "", sub("^ *","", as.character(format(index, scientific=F))))

    return( tab[index, ] )
}


##########################################################################################################
#
#
#
# p      - character, protein group id
# groups - character vector, maps the colnames of 'data' to groups
# data   - part of the protein groups table containing the data
#         e.g. intensities
#
# changelog: 2011      implementation
#            20120223  parameter 'plot.title' -> title of plots if 'plot=T'
#
##########################################################################################################

ANOVA <- function(p, groups, data, plot=F, plot.title="", conf.level=0.99){

    dat.frame  <- data.frame( g = factor(groups), y = as.numeric(data[p, ]) )
    #dat.frame  <- data.frame( g = factor(groups), y = data[p, ] )
    # anova
    AOV <- aov( dat.frame$y ~ dat.frame$g , dat.frame)

    # get the p-value of AOV
    AOV.p <- unlist(summary(AOV))["Pr(>F)1"]

    # Tukey HSD test
    AOV.hsd <- TukeyHSD(AOV, order=T, conf.level=conf.level)

    if(plot){

        par(mfrow=c(2,2))
        require(gplots)
        # plot means with ci
        plotmeans( dat.frame$y ~ dat.frame$g, xlab="group", ylab = "response", ci.label=F, digits=2, p=conf.level, use.t=T)
        # boxplots
        plot(dat.frame$g, dat.frame$y, xlab="group", ylab="response")
        # plot results of Zukey HSD
        plot(AOV.hsd)

        #
        textplot(plot.title)

    }

    ###########################
    # output
    ###########################
    out <- list()
    out[[1]] <- AOV
    out[[2]] <- AOV.p
    out[[3]] <- AOV.hsd
    names(out) <- c("aov", "p-value", "tukeyhsd")

    return(out)
}
################################################################################################
#
#                             aa composition of motifs
#
# - ph       - phosphosite table
# - column   - character, the column name of the sequence window
# - items    - all possible items that can be present in 'column'
#            - usually all amino acids, but also for example classes of amino acids
#
#
# changelog:   20100415 parameter column and items
#              20110307 estimation of relative frequencies:
#                          - absolute numbers diveded by number of sites
#                            (and not by the sum of all aa at that position)
#
################################################################################################
motif.comp <- function( ph, column="Sequence.Window", items="ACDEFGHIKLMNPQRSTVWY" ){

       Nsites <- dim(ph)[1]

       # sequence window of the site
       seq.win <- ph[, column]

       # make the list unique
       seq.win <- sub(";.*", "", seq.win)

       # split
       seq.win.split <- strsplit(seq.win, "")

       # matrix
       seq.win.split.mat <- matrix(unlist(seq.win.split), ncol=length( seq.win.split[[1]] ), byrow=T, dimnames=list( rownames(ph), c("-6", "-5", "-4", "-3", "-2", "-1", "0", "+1", "+2", "+3", "+4", "+5", "+6")  )  )

       #############################
       # items to count
       #############################
       items <- unlist(strsplit(items, ""))

       ##############################
       # motif aa composition matrix
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

#################################################################################################
#
#
#  prot.test   - numeric vector of the same length as rows in 'tab'
#              - the entries have to be either 0 (protein not in the test set) or 1 (protein
#                is in the test set)
#  tab         - the table containing the annotation, e.g. protein groups
#
#
# changelog:     20100419 implementation
#                20100714 - parameters of hypergeometric model are returned
#                         - associated proteins are returned
#                         - ONLY TERMS THAT ARE AMONG THE TEST SET ARE NOW TESTED!!! is this only valid when NOT testing for depletion??
#                20110209 - several columns columns can be submitted, e.g. KEGG, GO, Pfam
#                20110418 - bugfix: if the the terms in the table are separated by '; ' and multiple columns
#                                   were used, some terms could not be found when checking the column the terms
#                                   was derived from
#                20110523 - added parameter 'replace.ids'
#                         - if TRUE the names of input vector 'prot.test' are replaced by the first protein
#                           of column 'Protein.IDs' in the final output table
#                20110830 - 'min.x' minimal number of proteins in the test set associated with a term
#                         - if it's below that number, the term is disregarded
#                20121009 - in case of 'depletion' or 'two-sided' test, the function loops over all terms
#                           in the background table
#                         - in case of 'greater', the loop only concerns terms present in the test dataset
#                20121018 - changed parameter 'min.x' to 'min.m', i.e. the number of proteins in the
#                           background dataset is controlled
#
##################################################################################################
annotation.enrich <- function( prot.test ,tab, anno.col="KEGG.Pathways", desc.col=NULL, replace.ids=F, alternative=c( "greater", "two.sided", "less"), adjust=c("BH", "none"), min.m=5  ){

    ############################
    # get the arguments
    ############################
    alternative = match.arg(alternative)
    adjust = match.arg(adjust)

    ##############################
    # check whether the specified columns are present
    ##############################

    if( sum( anno.col %in% colnames(tab)) != length(anno.col) )
        stop(paste("There is no column of name ",anno.col[which(!(anno.col %in% colnames(tab)))], " in the table!\n", collapse=" "))

    ############################
    # get annotation terms
    ############################
    anno <- tab[, anno.col]             # vector, for each protein all terms separated by ';'

    # collapse several annotation columns
    if( length(anno.col) > 1 ){
        anno <- apply(anno, 1, function(x) paste(x[which(nchar(x)>0)], collapse=";"))
    }
    # make it character
    if(mode(anno) != "character") anno <- as.character(anno)


    anno <- strsplit(anno, ";")                        # list, for each protein a vector of all terms
    anno <- lapply(anno, function(x) sub("^ ", "", x)) # remove leading blanks


    anno <- lapply(anno, unique)                       # make the annotation terms non redundant

    names(anno) <- rownames(tab)


    #####################################################
    # parameters for Fisher's test
    # k   - number of proteins in the test set, i.e. the
    #       number of draws
    # N   - number of all proteins, i.e. the sum of
    #       white and black balls in the urn
    #####################################################
    k = sum(prot.test)
    N = dim(tab)[1]

    #################################################
    # divide the list of proteins and their terms
    # into the test set and the background
    #################################################
    anno.test <- anno[ which(prot.test == 1) ]
    #anno.bg <- anno[ which(prot.test == 0)  ]

    ###################################
    # non redundant vector of all terms
    # to test
    if(alternative == "greater")
        anno.nr <- unique(unlist(anno.test))        # only terms present in the test set
    else
        anno.nr <- unique(unlist(anno))             # all terms in the background; necesary to test for depletion

    # is it correct to test only terms that are present in the test set?
    # anno.nr <- unique(unlist(anno.test))

    # use only annotated terms, i.e. remove ""
    zero.char.idx <- which(unlist(lapply(anno.nr, nchar)) == 0)
    if(length(zero.char.idx) > 0)
        anno.nr <- anno.nr[ -zero.char.idx ]
    #return(anno.nr)

    #########################################
    # output vectors
    #########################################
    p.vec <- m.vec <- x.vec <- vector("numeric", length(anno.nr))         # raw p-value; m, x parameters of hypergeometrical model
    direction.vec <- desc.vec <- prots.vec <- cat.vec <- vector("character", length(anno.nr))     # direction & description & annotation categories
    names(p.vec) <- names(direction.vec) <- names(desc.vec) <- names(m.vec) <- names(x.vec) <- names(prots.vec) <- names(cat.vec) <- anno.nr


    ##############################################
    #    loop over all annotation terms
    ##############################################
    cat("\ntesting", length(anno.nr), "terms:\n" )
    count = 1


    for(tt in anno.nr){

        cat(".")
        if((count %% 100) == 0) cat(count, "\n")
        count <- count  + 1


        # number of all proteins that are annotated with the current term
        m <- sum(unlist(lapply(anno, function(x) as.numeric(tt %in% x) )))
        m.vec[[tt]] <- m

        # get  proteins in test set associiated with the current term
        term.prots.idx <- unlist(lapply(anno.test, function(x) as.numeric(tt %in% x) ))

        # number of proteins in the test set that are annotated with the current term
        x <- sum( term.prots.idx )

        # now get the protein ids
        term.prots <- ""
        if( x > 0 )
            term.prots <- paste( names(anno.test)[which(term.prots.idx == 1)], collapse=";")
        prots.vec[[tt]] <- term.prots
        x.vec[[tt]] <- x


        ############################
        # now test
        ############################
        pval <- hyper.test( x, m, N-m, k, alternative=alternative)

        # check whether the term is enriched or depleted
        enrich.direction <- ifelse( (x/k) > (m/N), "enrichment", "depletion" )

        # some additional information
        if( !is.null(desc.col) ){
            term.desc <- desc.nr[tt]
        }

        # store the results
        p.vec[tt] <- pval
        direction.vec[tt] <- enrich.direction

        ##########################################################
        # determine the annotation column of the current term
        ##########################################################
        if(length(anno.col) > 1){

            # loop over the annotation terms
            for(cc in anno.col){

                # flag all brackets in the term description by \\
                tt.exp <- gsub("\\(", "\\\\(", tt)
                tt.exp <- gsub("\\)", "\\\\)", tt.exp)
                tt.exp <- gsub("\\[", "\\\\[", tt.exp)
                tt.exp <- gsub("\\]", "\\\\]", tt.exp)

                if( length( grep(paste("(^|;|; )", tt.exp,"($|;)", sep=""), tab[, cc]  )  ) > 0 ) {
                    cat.vec[tt] <- cc
                    break;
                }
                cat.vec[tt] <- "TEST"
            }

        } else {
            cat.vec[tt] <- anno.col
        }

        if(!is.null(desc.col))
            desc.vec[tt] <- term.desc

    }
    cat("\n")



    ########################
    # adjust the p-values
    ########################
    p.vec.adj <- p.adjust(p.vec, adjust)

    ##############################################################
    # replace protein group ids by protein names (leading protein)
    ##############################################################
    if(replace.ids){
        prots.vec <-  unlist(lapply( prots.vec, function(x) paste(sub(";.*", "", tab[ as.character(unlist(strsplit(x, ";"))), "Protein.IDs"]), collapse=";") ))

    }


    ##############################################################
    # put it all together
    ##############################################################
    tab.dat <- data.frame(  p=p.vec, p.adj=p.vec.adj, direction=direction.vec, m=m.vec, x=x.vec, N=rep(N, length(p.vec)), k=rep(k, length(p.vec)), Category=cat.vec, proteins=prots.vec )

    # order according to adjusted p-value
    tab.dat <- tab.dat[order(tab.dat$p), ]



    ##############################################################
    # - remove out terms that occur less than 'min.m' in test dataset
    # - adjust p-values again
    # - add to output
    ##############################################################
    if(min.m > 1){

        # remove rows
        rm.idx <- which(tab.dat$m < min.m)
        if(length(rm.idx) > 0)
            tab.filt <- tab.dat[-rm.idx, ]

        # update adjusted p-values
        tab.filt$p.adj<- p.adjust( tab.filt$p, adjust )

    } else{
        tab.filt <- tab.dat
    }


    ######################################
    #        output
    ######################################

    # list with single vectors
    out <- vector("list", 10)
    names(out) <- c("p", "p adj", "direction", "m", "x", "N", "k", "proteins", "table", "table.filt")
    out[["p"]] <- p.vec                         # raw p-values
    out[["p adj"]] <- p.vec.adj                 # adjusted p-values
    out[["direction"]] <- direction.vec         # direction
    out[["m"]] <- m.vec                         # number of all proteins annotated with the terms
    out[["x"]] <- x.vec                         # number of proteins in the test set with the term
    out[["N"]] <- N                             # number of proteins in the universe
    out[["k"]] <- k                             # number of proteins in test set
    out[["proteins"]] <- prots.vec              # proteins associated with the terms


    out[["table"]] <- tab.dat                   # as table
    out[["table.filt"]] <- tab.filt             # table, min. proteins  per term in test set


    return(out)
}

#####################################################################################################################################################

annotation.enrich.par <- function( prot.test ,tab, anno.col="KEGG.Pathways", desc.col=NULL, cores=2, alternative=c( "greater", "two.sided", "less"), adjust=c("none", "BH")  ){

    require(doSMP)

    ############################
    # get the arguments
    ############################
    alternative = match.arg(alternative)
    adjust = match.arg(adjust)

    ##############################
    # check whether the specified columns are present
    ##############################
    if( sum( anno.col %in% colnames(tab)) == 0 ) stop(paste("There is no column of name ",anno.col," in the table!\n"))
    if(!is.null(desc.col))
        if( sum( anno.col %in% colnames(tab)) == 0 ) stop(paste("There is no column of name ",desc.col," in the table!\n"))


    ############################
    # get annotation terms
    ############################
    anno <- tab[, anno.col]             # vector, for each protein all terms separated by ';'
    if(mode(anno) != "character") anno <- as.character(anno)
    anno <- strsplit(anno, ";")         # list, for each protein a vector of all terms
    anno <- lapply(anno, unique)        # make the annotation terms non redundant
    #names(anno) <- rownames(anno)
    names(anno) <- rownames(tab)


    #####################################################
    # parameters for Fisher's test
    # k   - number of proteins in the test set, i.e. the
    #       number of draws
    # N   - number of all proteins, i.e. the sum of
    #       white and black balls in the urn
    #####################################################
    k = sum(prot.test)
    N = dim(tab)[1]

    #################################################
    # divide the list of proteins and their terms
    # into the test set and the background
    #################################################
    anno.test <- anno[ which(prot.test == 1) ]
    anno.bg <- anno[ which(prot.test == 0)  ]

    # non redundant vector of all terms
    anno.nr <- unique(unlist(anno))

    # is it correct to test only terms that are present in the test set?
    anno.nr <- unique(unlist(anno.test))


    #########################################
    # output vectors
    #########################################
    p.vec <- m.vec <- x.vec <- vector("numeric", length(anno.nr))         # raw p-value; m, x parameters of hypergeometrical model
    direction.vec <- desc.vec <- prots.vec <- vector("character", length(anno.nr))     # direction & description
    names(p.vec) <- names(direction.vec) <- names(desc.vec) <- names(m.vec) <- names(x.vec) <- names(prots.vec) <- anno.nr


    ##############################################
    #    loop over all annotation terms
    ##############################################
    workers = startWorkers(cores)
    registerDoSMP(workers)

    #for(tt in anno.nr){
    ttt <- foreach(tt=anno.nr ) %dopar% {

        # number of all proteins that are annotated with the current term
        m <- sum(unlist(lapply(anno, function(x) as.numeric(tt %in% x) )))
        m.vec[[tt]] <- m

        # get  proteins in test set associiated with the current term
        term.prots.idx <- unlist(lapply(anno.test, function(x) as.numeric(tt %in% x) ))

        # number of proteins in the test set that are annotated with the current term
        x <- sum( term.prots.idx )
        term.prots <- ""
        if( x > 0 )
            term.prots <- paste( names(anno.test)[which(term.prots.idx == 1)], collapse=";")
        prots.vec[[tt]] <- term.prots
        x.vec[[tt]] <- x


        ############################
        # now test
        ############################
        pval <- hyper.test( x, m, N-m, k, alternative=alternative)

        # check whether the term is enriched or depleted
        enrich.direction <- ifelse( (x/k) > (m/N), "enrichment", "depletion" )

        # some additional information
        if( !is.null(desc.col) ){
            term.desc <- desc.nr[tt]
        }

        # store the results
        p.vec[tt] <- pval
        direction.vec[tt] <- enrich.direction
        if(!is.null(desc.col))
            desc.vec[tt] <- term.desc

    }
    stopWorkers(workers)
    ########################
    # adjust the p-values
    ########################
    p.vec.adj <- p.adjust(p.vec, adjust)

    ######################################
    #        output
    ######################################

    # list with single vectors
    out <- vector("list", 9)
    names(out) <- c("p", "p adj", "direction", "m", "x", "N", "k", "proteins", "table")
    out[["p"]] <- p.vec                         # raw p-values
    out[["p adj"]] <- p.vec.adj                 # adjusted p-values
    out[["direction"]] <- direction.vec         # direction
    out[["m"]] <- m.vec                         # number of all proteins annotated with the terms
    out[["x"]] <- x.vec                         # number of proteins in the test set with the term
    out[["N"]] <- N                             # number of proteins in the universe
    out[["k"]] <- k                             # number of proteins in test set
    out[["proteins"]] <- prots.vec              # proteins associated with the terms

    # the same as dataframe
    tab.dat <- data.frame(p=p.vec, p.adj=p.vec.adj, direction=direction.vec, m=m.vec, x=x.vec, N=rep(N, length(p.vec)), k=rep(k, length(p.vec)), proteins=prots.vec )
    # order according to p-value
    tab.dat <- tab.dat[order(tab.dat$p.adj), ]

    out[["table"]] <- tab.dat

    return(out)
}
######################################################################################################################################
#
#
#
# changelog: 20111027 implementation
######################################################################################################################################
fold.enrichment.plot <- function(tab, p=0.01, ...){

    # get all terms with p < p
    tab.filt <- tab[which(tab$p.adj < 0.01), ]

    # calculate fold enrichment
    fold.enrich <-  apply(tab.filt, 1, function(x) (as.numeric(x["x"])/as.numeric(x["m"])) / (as.numeric(x["k"])/as.numeric(x["N"])))
    ord.idx <- order(fold.enrich)

    # p-values color gradient
    p.val <- tab.filt$p.adj

    # color gradient
    cgr <- colorpanel(length(p.val), "red", "orange", "yellow")

    #######################
    # plot
    #######################
    par(mar=c(5, 20, 5, 5))
    barplot(fold.enrich[ord.idx], horiz=T, las=1, xlab="fold enrichment", col=cgr[ord.idx], ...)

    #fancyPlot(fold.enrich, p.val)

return(fold.enrich[ord.idx])
}

####################################################################################################################################
#                               barcharts of enriched terms test vs. background
#
#
# table.filt  - matrix, results table from function 'annotation.enrich'
# color       - vector of length two
# xlab        - character
# main        - character
#
# changelog: - 20120229 implementation
#              20120305 bugfix: if there is only one significant term, the figure will be drawn correctly
#
####################################################################################################################################
GO.chart <- function(tab.filt, xlim=NULL, xlab="Proportion of proteins in percent", color=c(colors()[430], colors()[135]), legend.txt=c("Down regulated", "Detected phospho proteome"), main=""){

    barplot.mat <- matrix(0, nrow=dim(tab.filt)[1], ncol=2, dimnames=list(rownames(tab.filt), c("test", "background")))

    for(term in rownames(tab.filt)){

        ratio.test <- (tab.filt[term, "x"]/tab.filt[term, "k"]) * 100
        ratio.bg <- (tab.filt[term, "m"]/tab.filt[term, "N"]) * 100

        barplot.mat[term, "test"] <- ratio.test
        barplot.mat[term, "background"] <- ratio.bg
    }

    par(mar=c(5, 20, 5, 5))
    if(dim(barplot.mat)[1] == 1)
        barplot(t(barplot.mat), horiz=T, beside=T, las=1, xlab=xlab, xlim=xlim, col=color, main=main, ylim=c(0.5, 3))
    else
        barplot(t(barplot.mat[order(barplot.mat[, "test"], decreasing=F), ]), horiz=T, beside=T, las=1, xlab=xlab, xlim=xlim, col=color, main=main)
    legend("bottomright", legend=legend.txt, fill=color, bty="n" )

}


######################################################################################################################################
#
#
#
#  test.sites    - part of phospho table containing the test sites
#  bg.sites      - part of phospho table containing the background dataset
#                - or a character string specifying the path to a protein database in fasta format
#  p.motif       - numeric, p-value defining 'significant'  enrichment/depletion
#  rm.zero       - logical, if TRUE the modification position will be removed from the matrix
#  dist          - character, specifies the type of distribution used for testing
#                - can either be binomial or hypergeometric
#  main          - character, title for the heatmap
#  alternative   - character, specifies the direction of the test
#                - any other value that 'greater' has only an effect on the binomial test
#  col           - color scheme to use, e.g. "redgreen", "greenred", "heatcolors"
#  p.min         - numeric, minimal p-value, i.e. p-values that are below this value will be set to this value before logarithmize the
#                  value. this is only done to increase the 'dynamic range' of the heat map and has no effect on the real p-values
#                  that are returned by this function!
#  p.max         - numeric, this value determines the color scale of the heat map, i.e. if the most significant p-value in the data is 1e-2
#                  the colors are still scaled from -p.max to p.max
#                - this shall avoid any overinterpretation of not very significant p-values only because they jump out of the matrix...
#  column        - character, the column name holding the sequence motif
#  items         - character, all possible items that can occur in the column
#  cn            - logical, if TRUE the p values are shown as cell notes in the heatmap
#
#  ...           - further arguments passed to 'binomial.test'
#
# changelog:    20100323 - reimplementation based on the original version
#               20100324 - some doumentation
#               20100325 - added a legend that summarizes the test parameters
#               20100401 - the color scale is now always symmetric
#                        - parameter 'pmax'
#                        - swapped the matrix
#               20100415 - added the calculation of 'two.sided' and 'less' for the hypergeometric test
#                        - parameter column and items
#               20100422 - changed parameter 'col'
#               20100504 - added parameter 'cellnote'
#               20110307 - frequencies can now be compared against overall amino acid frequencies in the protein
#                          database using the binomial test
#               20110919 - anti-motif is spelled out in the figure legend
#####################################################################################################################################
motif.enrich.new <- function(test.sites, bg.sites, p.motif = 1e-6, rm.zero=F, dist=c("binomial", "hypergeometric"), main="", alternative=c("greater", "less", "two.sided"), col=c("green","black", "red"), p.min=1e-10, p.max=1e-5, column="Sequence.Window", items="AGILMVCNPQSTFYWHKRDE", cellnote=T, ...){

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
    #cn <- matrix(paste(test[[1]], "/",  bg[[1]], "\n(", round(test[[1]]/n.draws, 2), "/",round(bg[[1]]/dim(bg.sites)[1],2),")\np=", formatC(p.mat, digits=3),sep=""), nrow=13, byrow=F)

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

    ##############################################
    # extract the possible motif / anti motif
    ##############################################
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
               aa.add <- "X"
           if(length(aa) == 1)
               aa.add <- aa
           if(length(aa) > 1 )
               aa.add <- paste("[",paste(aa, collapse="")  ,"]", sep="")

           if(i != "0")
               motif.string <- paste( motif.string, aa.add, sep=""  )

           # add the site itself [STY]
           if(rm.zero){
              if(i == "-1")
                 motif.string <- paste( motif.string, "[ST]", sep=""  )
           } else {
              if(i == "0")
                 motif.string <- paste( motif.string, "[ST]", sep=""  )

           }
           rm(aa.add)
           #########################
           # depletion
           #########################
           if(length(aa.depl) == 0)
               aa.add <- "X"
           if(length(aa.depl) == 1)
               aa.add <- aa.depl
           if(length(aa.depl) > 1 )
               aa.add <- paste("[",paste(aa.depl, collapse="")  ,"]", sep="")

           if(i != "0")
               motif.string.depl <- paste( motif.string.depl, aa.add, sep=""  )

           # add the site itself [STY]
           if(rm.zero){
              if(i == "-1")
                 motif.string.depl <- paste( motif.string.depl, "[ST]", sep=""  )
           } else {
              if(i == "0")
                 motif.string.depl <- paste( motif.string.depl, "[ST]", sep=""  )

           }

       }
     motif.string <- paste( paste(motif.string, sep=""), "\n",paste(motif.string.depl, sep=""), collapse="")
    #####################################
    # 'breaks': mapping data to colors
    #####################################

    # maximal absolute data value
    p.log.scale <- max( abs(p.mat.log) )

    if(-log(p.max,10) > p.log.scale) p.log.scale = -log(p.max,10)

    p.log.scale <- min( p.log.scale, -log(p.min, 10)  )
    breaks.param <- seq( -p.log.scale, p.log.scale, 0.5  )

    #####################################
    # colors
    #####################################
    require(gplots)

    if(length(col) == 3)
        col <- colorpanel( length(breaks.param)-1, col[1], col[2], col[3]  )


    #######################
    # plot
    #######################

    # with cellnotes
    if(cellnote){
        heatmap.2(  t(p.mat.log) , Rowv=F,Colv=F, trace="n", dendrogram="none", symkey=ifelse(alternative=="two.sided",T,F), cellnote=t(cn), notecex=.5, notecol="black", sub=paste(motif.string, "\n(p < ", p.motif,")"), main=paste( main, "(", dim(test.sites)[1],"sites )" ), col=col, breaks=breaks.param)
    } else { # without cellnotes
        heatmap.2(  t(p.mat.log) , Rowv=F,Colv=F, trace="n", dendrogram="none", symkey=ifelse(alternative=="two.sided",T,F), sub=paste(motif.string, "\n(p < ", p.motif,")"), main=paste( main, "(", dim(test.sites)[1],"sites )" ), col=col, breaks=breaks.param)
    }
    legend("top", legend=c(paste("background:", dim(bg.sites)[1], "sites"), paste("minimal p-value used for plotting: ", p.min), paste("test:", dist), paste("alternative:", alternative)), bty="n")



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
################################################################################################
#
#                               motif.enrich.iter
# do enrichment analysis in several iterations
#
# 1) do 'normal' motif analysis
# 2) get significantly enriched AA at specific positions
# 3) get subset of test-sites having the AA
# 4) do analysis on them
#
# changelog: 20110919 implementation
################################################################################################
motif.enrich.iter <-  function(test.sites, bg.sites, p.motif = 1e-6, rm.zero=F, dist=c("binomial", "hypergeometric"), main="", alternative=c("greater", "less", "two.sided"), col=c("green","black", "red"), p.min=1e-10, p.max=1e-5, column="Sequence.Window", items="AGILMVCNPQSTFYWHKRDE", cellnote=T, mult.win=F, ...){

    ##########################################################
    # do motif analysis on the whole test set
    ##########################################################
    motif.ana <- motif.enrich.new(test.sites=test.sites, bg.sites=bg.sites, p.motif=p.motif, rm.zero=rm.zero, dist=dist, main=main, alternative=alternative, col=col, p.min=p.min, p.max=p.max, column=column, items=items, cellnote=cellnote, ...  )

    # get matrix of p-values
    p.mat <- motif.ana[['p-values']]
    p.plot <- motif.ana[['p-plot']]

    # set everything that is not significant to 1
    p.mat[p.mat >= p.motif] <- 1
    p.mat[p.plot < 0] <- 1      # set all significant p-avlues for depletion to insignificant

    # check if there are some significant AAs at certein position
    if(sum(p.mat) < (dim(p.mat)[1]*dim(p.mat)[2])){
         #####################################################################
         # loop over matrix entries and test every AA separately
         #####################################################################
         for(pos in rownames(p.mat)){

             ####################################################
             # check combinatoions of  AAs at a specific position
             ####################################################
             if( sum( p.mat[pos, ] < 1  ) > 1  ){

                    grep.vec <- rep(".",dim(p.mat)[1])
                    names(grep.vec) <- rownames(p.mat)

                    grep.vec[pos] <- paste( "[", paste(colnames(p.mat)[ which(p.mat[pos, ] < 1 ) ], collapse="|"), "]",sep="")
                    grep.vec <- paste( "^", paste(grep.vec, collapse=""), "$", sep="")

                    test.tmp <- test.sites[ grep(grep.vec, test.sites[, column]) ,]

                    if(mult.win) X11()
                    motif.tmp <- motif.enrich.new(test.sites=test.tmp, bg.sites=bg.sites, p.motif=p.motif, rm.zero=rm.zero, dist=dist, main=paste(paste(colnames(p.mat)[ which(p.mat[pos, ] < 1)  ], collapse="/"), "at", pos), alternative=alternative, col=col, p.min=p.min, p.max=p.max, column=column, items=items, cellnote=cellnote, ...  )


             }
             ###########################################
             # check AAs separately
             ###########################################
             for(aa in colnames(p.mat)){

                # check whether entry is significant
                if(p.mat[pos,aa] < 1){

                    grep.vec <- rep(".",dim(p.mat)[1])
                    names(grep.vec) <- rownames(p.mat)

                    grep.vec[pos] <- aa
                    grep.vec <- paste( "^", paste(grep.vec, collapse=""), "$", sep="")

                    ############################
                    # get all sites with that aa at that position
                    ############################
                    test.tmp <- test.sites[ grep(grep.vec, test.sites[, column]) ,]

                    # do analysis again
                    if(mult.win) X11()
                    motif.tmp <- motif.enrich.new(test.sites=test.tmp, bg.sites=bg.sites, p.motif=p.motif, rm.zero=rm.zero, dist=dist, main=paste(aa, "at", pos), alternative=alternative, col=col, p.min=p.min, p.max=p.max, column=column, items=items, cellnote=cellnote, ...  )

                } # end if

            } # end for aa
        } #end or pos


    } # end if
}



################################################################################################
#                         pairwise sequence alignment
#
# changelog:   20100421  implementation
#              20100906  switched to Biostrings package
#
################################################################################################
myPairAlign <- function( seq1, seq2, info1="Sequence1", info2="Sequence2", gapOpening=-10, gapExtension=-4, alignment=c("global", "local", "overlap")){
    #require(pairseqsim)
    require(Biostrings)


    #options(warn=-1)
    data(BLOSUM62)
    #options(warn=0)

    alignment=match.arg(alignment)

    # create sequence objects
    #seq1.obj <- new("AASequence", seq1, info=info1)
    #seq2.obj <- new("AASequence", seq2, info=info2)

    # peform the alignment
    #align <- salign( seq1.obj, seq2.obj, EPAM110, delta=delta, gapext=gapext, alignment=alignment, scoring=scoring  )
    align <- try(pairwiseAlignment( seq1, seq2, substitutionMatrix=BLOSUM62, gapOpening=gapOpening, gapExtension=gapExtension, type=alignment ), silent=T)


    return(align)
}

################################################################################################
#
#                  mapp KO (KEGG orthologs) ids to kegg pathway ids
#
# ko         - vector of KO ids
#
# changelog: 20100924 implementation
################################################################################################
KO2KEGG <- function( ko ){

  ko.map <- KO.MAP()
  title.map <- MAP.TITLE()


  # get the KEGG ids of the KO ids
  KEGG.ids <- ko.map[ rownames(ko.map) %in%  ko,]
  names(KEGG.ids) <- rownames(ko.map)[rownames(ko.map) %in%  ko]

  # output matrix
  out <- matrix("", ncol=3, nrow=length(ko), dimnames=list(ko, c("KEGG.ko", "KEGG", "KEGG.Names")))

  # based on them get the KEGG pathway names
  for(id in names(KEGG.ids)){

      kegg <- unlist(strsplit(KEGG.ids[id], " "))

      kegg.path <- paste(title.map[ rownames(title.map) %in% kegg,   ], collapse=";")

      out[id, "KEGG.ko"] <- id
      out[id, "KEGG"] <- paste(kegg, collapse=";")
      out[id, "KEGG.Names"] <- kegg.path
  }

 return(out)
}
################################################################################################
#
#                               get the kegg ortholog map
#
# changelog: 20100924 implementation
#################################################################################################
KO.MAP <- function(ftp="ftp://ftp.genome.jp/pub/kegg/pathway/ko/ko_map.tab"){

    map <- read.delim(ftp, row.names=1, stringsAsFactors=F, header=F)

    return(map)
}

#################################################################################################
#
#                         get the id -> name map
# changelog: 20100924 implementation
#################################################################################################
MAP.TITLE <- function(ftp="ftp://ftp.genome.jp/pub/kegg/pathway/map_title.tab"){

    map <- read.delim(ftp, row.names=NULL,header=F, stringsAsFactors=F)

    # fix the KEGG ids, i.e "20" -> "00020"
    ids <- map[,1]
    ids <- unlist(lapply(ids, function(x) ifelse( (nc <- nchar(x)) < 5, paste(c(rep("0", 5-nc), x), collapse=""), x ) ))
    map <- map[, -c(1)]
    dim(map) <- c( length(ids), 1  )
    rownames(map) <- ids

    return(map)

}


################################################################################################
#
# pg.id     - protein group id
# ph        - phospho site table
# pg        - protein groups table
# pfam      - character, filename
# fraction  - numeric between 1 to 50
#           - specifies a binning of the sequence; afterwards the number of sites per bin
#             is reported
#
# changelog:  20100304 implementation
#             20100330 pfam
#################################################################################################
sitesOnProteinPlot <- function(pg.id, ph, pg, pfam=NULL, label.column="Motifs", label="Aurora" ,fraction=10, rel=F ){


    ##################################
    # - am annotation downloaded from
    # Sanger S.pombe website
    # - most likely this won't work
    #   with other files
    ##################################
    if(!is.null(pfam) & (mode(pfam) == "character") ){

        pfam <- read.delim(pfam,  stringsAsFactors=F, comment.char="#", header=F)
        # make it tab-separated
        pfam.tab <- apply(pfam,1, function(x) gsub(" {1,20}", "\t", x))
        pfam.tab <- lapply(pfam.tab, strsplit,"\t")
        pfam <- matrix( unlist(pfam.tab), byrow=T, ncol=15  )
        colnames(pfam) <- c("seq id","alignment start", "alignment end", "envelope start", "envelope end", "hmm acc", "hmm name", "type", "hmm start", "hmm end", "hmm length", "bit score", "E-value", "significance", "clan")

    }

    # relative positions of sites
    #rel=T

    #################################
    # get all sites on that protein
    #################################
    #ph.pg <- ph[ph[, "Protein.Group.IDs"] ==pg.id,  ]
    ph.pg <- ph[ grep(paste("(^|;)", pg.id, "($|;)", sep="")  , ph$Protein.Group.IDs)  ,]

    ##############################################################
    # if pg is  not a table assume that it is a fasta database
    ##############################################################
    if(!is.null(dim(pg))){

       ##################################
       # get the row of the protein group
       ##################################
       pg.pg <- pg[pg.id, ]

       prot.id <- sub( ";.*" ,"" ,pg.pg["Protein.IDs"])

       #################################
       # protein length
       #################################
       protL <- as.numeric(pg.pg["Sequence.Length"])

       #################################
       # protein description
       #################################
       protDesc <- pg.pg["Protein.Descriptions"]

    } else {

        prot.id <- ph.pg$Protein[1]
        protL = length(pg[[prot.id]])
        protDesc <- ph.pg$Protein.Descriptions[1]

    }
    #################################
    # if there is at least 1 sites....
    #################################
    if(dim(ph.pg)[1] > 0){

        # get the positions
        ph.pos <- as.numeric(ph.pg[, "Position"])
        names(ph.pos) <- rownames(ph.pg)

        # relative position
        if(rel) ph.pos <- round(ph.pos/protL,3)

        #########################
        # plot
        #########################
        if(rel){ xlim= c(0,1)
        } else {
            xlim=c(1, protL )
        }
        ylim=c(0,1)

        plot.new()
        plot.window(xlim=xlim, ylim=ylim)

        # the protein
        lines(xlim, c(.5,.5), lwd=2 )
        text(xlim[1], 0.45, xlim[1], srt=45)
        text(xlim[2], 0.45, xlim[2], srt=45)
        text(xlim[1], 0.8, paste( protDesc, "  ( sequence length: ", as.numeric(protL), " aa )"), pos=4)

        # pfam annotation
        if(!is.null(pfam)){

            # check if the protein has any pfam annotation
            if( prot.id %in% pfam[,1]){

                pfam.pg <- pfam[which(pfam[, 1] ==  prot.id),   ]

                # if there is only one pfam anno
                if(is.null(dim(pfam.pg))){
                    pf.coords <- as.numeric(pfam.pg[ c("envelope start", "envelope end")])
                    if(rel) pf.coords <- pf.coords/protL

                    lines( pf.coords, c(0.52, 0.52), lwd=2, col="red"  )
                    text( pf.coords[1] + ((pf.coords[2]-pf.coords[1])/2), 0.53, paste(pfam.pg["hmm name"], " (",pfam.pg[ "type"], ")", sep=""), col="red", pos=3 )
                } else {

                    for(pf in 1:dim(pfam.pg)[1]){
                        pf.coords <- as.numeric(pfam.pg[pf,  c("envelope start", "envelope end")])
                        if(rel) pf.coords <- pf.coords/protL

                        lines( pf.coords, c(0.52, 0.52), lwd=2, col="red"  )
                        text( pf.coords[1] + ((pf.coords[2]-pf.coords[1])/2), 0.53, paste(pfam.pg[pf, "hmm name"], " (",pfam.pg[pf, "type"], ")", sep=""), col="red", pos=3 )

                    }

                }

            }

        } # end if is.null pfam

        # the sites
        for(s in names(ph.pos)){

           # colors
           #col <- ifelse( ph.pg[s, "class.down.BH"] == "class1", "green", "black"  )
           col <- ifelse( ph.pg[s, "Localization.Prob"] >= .75, "green", "black"  )


           lines(rep(ph.pos[s],2), c(0.5,0.55) , col=col  )
           text(ph.pos[s], 0.45, ph.pos[s], srt=45)

           # this works only with the Aurora table
           #text(ph.pos[s], 0.61, ph.pg[s, "class.down.BH"], srt=90, pos=3)
           #text(ph.pos[s], 0.3, paste(ph.pg[s, c( "Exp1.BH","Exp2.BH", "Exp3.BH" )], collapse=" "), srt=90)

        }
        ##################################
        # count the number of sites per
        # segements of 10%
        ###################################
        site.count <- vector("numeric", round(100/fraction))
        tmp <- seq(0,100, fraction)
        names.site.count <- c()
        for(f in 1:(length(tmp)-1))
            names.site.count <- c(names.site.count, paste(tmp[f], "-", tmp[f+1], sep=""))
        names(site.count) <- names.site.count

        frac.rel <- fraction/100
        count = 1
        for(i in seq(0, 1-frac.rel, frac.rel ) ){
            site.count[ count ] <- sum( ph.pos > i & ph.pos <= (i+frac.rel)  )
            count <- count + 1
        }

        site.rel <- site.count/length(ph.pos)
        ##################################
        # some output
        ##################################
        out <- vector("list", 2)
        names(out) <- c(paste(pg.id, "sites absolute"), paste(pg.id,"sites relative"))
        out[[1]] <- site.count
        out[[2]] <- site.rel

        return(out)
    }

}
#################################################################################################
#
#
#
#################################################################################################


##################################################################################################
#                                  singlePeptideFilter
#
# determines the minimal PEP value of a reverse peptide hit and removes each
# single peptide pg with a PEP value greater or equal than this number.
#
# arguments:
#   proteinGroups      - 'proteinGroups.txt'
#   peptides           - 'peptides.txt'
#
# value:
#   proteinGroups      - filtered 'proteinGroups'
#   PEPthr             - minimal PEP value of a reverse peptide hit
#   filtered PG        - number of pg that were filtered out
#
# changed:   20090513 - implementation
#
##################################################################################################
singlePeptideFilter <- function( proteinGroups, peptides )
{
    # get the best PEP value of a reverse peptide hit
    PEPthr <- min(as.numeric(peptides[ which(peptides[, "Reverse"]== "+" ), "PEP"]   ))

    # determine all single peptide protein groups AND that have a PEP greater or equal than 'PEPthr'
    spPGidx <- which( (as.character(proteinGroups[, "Peptides..seq."]) == "1") & (proteinGroups[, "PEP"] >= PEPthr) )

    # ... and filter them out
    proteinGroups <- proteinGroups[-spPGidx, ]


    ############################
    #     output
    ############################
    output <- vector("list", 3)
    names(output) <- c("proteinGroups", "PEPthr", "filtered PG")

    output[[1]] <- proteinGroups
    output[[2]] <- PEPthr
    output[[3]] <- length(spPGidx)


    return(output)

}
##################################################################################################
#
#
#
#
patternMatch2 <- function( s, regexpr)
{
    # check if there are some flagged brackets
    brackets <- (grep( "\\(", regexpr) == 1) & (grep("\\)", regexpr) == 1)

    if( brackets )
    {
        bL <- regexpr("\\(", regexpr, fixed=T)

        bracketsLeft <- substr( regexpr, 1, bL[1]+unlist(attributes(bL))-3)
        cat(bL, "   ", bracketsLeft,"\n")

        bR <- regexpr("\\)", regexpr, fixed=T)
        bracketsRight <- substr( regexpr, bR[1]+3, nchar(regexpr))

        cat(bR, "   ", bracketsRight,"\n")

        if(nchar(bracketsLeft) > 0)
            s <- unlist(strsplit(s, bracketsLeft, perl=T))
        if(nchar(bracketsRight) > 0)
            s <- unlist(strsplit(s, bracketsRight, perl=T))

        #return(s)

        #bracketExpr <- unlist(strsplit( regexpr, "\\(", fixed=T  ))
        #bracketExpr <- unlist(lapply( bracketExpr, function(x) unlist( strsplit( x, "\\)", fixed=T))))

        #regexpr <- paste(bracketExpr, collapse="")
        #bracketExpr<- bracketExpr[2]
        return(s[which.max(nchar(s))])
    }

    # check if 'regexpr' matches 's'
    m <- regexpr( regexpr, s, extended=F)

    if(m[1] != -1)
    {
        subString <- substr( s, m[1], m[1]+unlist(attributes(m)))

        if(brackets)
        {    tmp <- regexpr(bracketExpr, subString, extended=F)
             return( substr(subString, tmp[1], tmp[1]+unlist(attributes(tmp)))  )
        }
        else
            return(subString)

    }
    else
        return(-1)
}
###########################################################
patternMatch <- function(s, regexpr)
{
    tmp <- regexpr(regexpr, s)

    return( substr( s, tmp[1], tmp[1] + unlist(attributes(tmp))-1  ) )
}




############################################################################################################
#                                       summaryPDFGUI
# template: 'guiDlgFunction'
#
#
# changelog:  20110921 implementaion
#############################################################################################################
summaryPDFGUI <- function(){

    require(svDialogs)


    fun="summaryPDF"
    execfun = getOption("guiExecFun")
    width=40
    labelwidth=15

    if (!exists(fun, where = 1, mode = "function"))
        stop(fun, "does not exist or is not a function!")

    # get all arguments of 'summaryPDF'
    arg.fun <- formals( get(fun, pos = 1, mode = "function"))


    ######################################
    # initialize gui object
    ######################################
    Tpl <- list( list( fun=fun, title="MaxQuantParser", sep=NULL, width=width, labelwidth=labelwidth, message=NULL, help=NULL  )   )

    # parameters to change
    Tpl[[2]] <- list( type="entry", argname = "title", default=deparse(arg.fun[[1]]), message="title for pdf file")                                                      # title
    Tpl[[3]] <- list( type="entry", argname = "separation", default=deparse(arg.fun[[grep("separation", names(arg.fun))]]), message="should separation be checked?")     # separation
    Tpl[[4]] <- list( type="entry", argname = "norm.ratio", default=deparse(arg.fun[[grep("norm.ratio", names(arg.fun))]]), message="plot normalized SILAC ratios?")     # normalized SILAC ratio?
    Tpl[[5]] <- list( type="entry", argname = "ic", default=deparse(arg.fun[[grep("ic", names(arg.fun))]]), message="should SILAC incorporation check be performed?")    # incorporation check


    #
    class(Tpl) <- c("guiDlg", "gui")

    pdf=display(Tpl)
    if (is.null(execfun))
        execfun <- "guiEval"

    if (exists(execfun, where = -1, mode = "function")) {
            get(execfun, pos = -1, mode = "function")(pdf)
        }
        else warning(execfun, " not found!")


    return(invisible(pdf))
}





#####################################################################################
#
# generate latex code: - description of MQ tables
#
#
#
#####################################################################################
summaryPDF.table.desc <- function(){


    SweaveFile <- paste("\\section{Description of result tables}\n",  sep="")

    ###############################################
    #  Summary table
    ###############################################
    SweaveFile <- paste(SweaveFile, "\n\\begin{frame}
         \\frametitle{Summary table}

          The summary file contains summary information for all the raw files processed with a single MaxQuant run.
          The summary information consists of some MaxQuant parameters, information of the raw file contents, and
          statistics on the peak detection. Based on this file a quick overview can be gathered on the quality of the data
          in the raw file.
          The last row in this file contains the summary information for each column on each of the processed files.

         \\end{frame}\n", sep="")



    SweaveFile <- paste(SweaveFile, "\n\\begin{frame}
          \\frametitle{Summary table}

          \\begin{table}[ht]
            \\begin{center}
            \\tiny
            \\begin{tabular}{p{3.2cm}|p{7.5cm}}
            \\hline
            Name & Description \\\\
            \\hline
           Raw File & The raw file processed. \\\\
           \\hline
           Protease & The protease used to digest the protein sample. \\\\
           \\hline
           Protease first search & The protease used for the first search. \\\\
           \\hline
           Use protease first search &  Marked with '+' when a different protease setup was used for the first search. \\\\
           \\hline
           Fixed modifications & The fixed modification(s) used during the identification of peptides. \\\\
           \\hline
           Variable modifications & The variable modification(s) used during the identification of peptides. \\\\
           \\hline
           Variable modifications first search & The variable modification(s) used during the first search. \\\\
           \\hline
           Use variable modifications first search & Marked with '+' when different variable modifications were used for the first search. \\\\
           \\hline
           Multiplicity & The number of labels used. \\\\
           \\hline
           Max. missed cleavages & The maximum allowed number of missed cleavages. \\\\
           \\hline
           Labels0 & The labels used in the SILAC experiment. Allowed values for X: 0=light; 1=medium; 2=heavy SILAC partner. \\\\
           \\hline
           LC-MS run type & The type of LC-MS run. Usually it will be 'Standard' which refers to a conventional shotgun proteomics run with datadependent MS/MS. \\\\
           \\hline
           Time-dependent recalibration & When marked with ‘+’, time-dependent recalibration was applied to improve the data quality. \\\\
           \\hline
           MS & The number of MS spectra recorded in this raw file. \\\\
           \\hline
           MS/MS & The number of tandem MS spectra recorded in this raw file. \\\\
           \\hline
           MS/MS Submitted & The number of tandem MS spectra submitted for analysis. \\\\
           \\hline
           MS/MS Submitted (SIL) & The number of tandem MS spectra submitted for analysis, where the precursor ion was detected as part of a SILAC cluster. \\\\
           \\hline
           MS/MS Submitted (ISO) & The number of tandem MS spectra submitted for analysis, where the precursor ion was detected as an isotopic pattern. \\\\
           \\hline
           MS/MS Submitted (PEAK) & The number of tandem MS spectra submitted for analysis, where the precursor ion was detected as a single peak. \\\\


          \\end{tabular}
          \\end{center}
          \\end{table}\n
         \\end{frame}\n", sep="")


    return(SweaveFile)
}



#####################################################################################
#                                  extractRatios
# The function extracts the ratios from the MaxQuant output.
#
# Assuming you have performed a SILAC experiment, i.e. you have protein ratios respective
# to the labels, and you have
#
# arguments:
#    proteinGroups   - matrix representing the output file
#                         'proteinGroups.txt'
#    groupLabel      - either 'NULL' or a character vector specifying the group labels used
#                      depending on the experiment type
#                    - group labels are specified in 'experimentalDesign.txt'
#    silacLabel      - character that specifies the SILAC labels used
#                    - either 'none', 'doublets', 'triplets'
#    sigB            - significance threshold (significance B)
#                    - only proteins with at least one sinificant ratio
#                      are considered
#    adjust          - logical, if TRUE the sigB values are adjusted for multiple comparisons
#                    - NOTE: currently (20090305) I am not sure if one has to adjust the p-values of the unfiltered table
#                        or if it is sufficient to adjust the p-values of the remaining protine groups (i.e. proteins that
#                        are quantified in all timepoints)
#
#    method          - character, method used to adjust the p-values
#                    - default: BH - Benjamini & Hochberg (1995)
#
#####################################################################################
extractRatios <- function(proteinGroups, groupLabel=NULL, silacLabel = c("none", "doublets", "triplets") , sigB=.05, adjust=T, method="BH", dot=T ){

      # vector for storing the colnames of of SILAC ratios
      # this depends on parameters'groupLabel' and 'silacLabel'
      ratioIdx <- c()
      ratioSigBIdx <- c()

      # matrix that contains some informations about the filtering
      summary <- matrix(0, nrow=6, ncol=1)
      rownames(summary) <- c("protein groups", "group label", "silac label", "completely quantified ", paste("sigB >", sigB) , "# remaining")

      # check if  there is a group label
      if(is.null(groupLabel)) groupLabel=""

      summary[1, 1] <- dim(proteinGroups)[1]
      summary[2, 1] <- paste(groupLabel, collapse=",")
      summary[3, 1] <- silacLabel


      # determine the experiment type
      silacLabel = match.arg(silacLabel)
      silacLabel = getSILAClabel(silacLabel, dot=dot)

      ##############################################################
      # produce the respective colnames for indexing the columns
      # the parameter 'dot' specifies if the dots within the
      # column names are present or were removed
      ##############################################################
      count=1
      for(grLa in groupLabel)       # loop over the group labels
          for(siLa in silacLabel)   # loop over the SILAC labels
          {
              if(dot)
               {   ratioIdx[count] <- paste( "Ratio.",siLa,".Normalized.", grLa, sep="")
                   ratioSigBIdx[count] <- paste( "Ratio.", siLa, ".Significance.B..", grLa, sep="" )
               }
              if(!dot)
              {
                   ratioIdx[count] <- paste( "Ratio",siLa,"Normalized", grLa, sep="")
                   ratioSigBIdx[count] <- paste( "Ratio", siLa, "SignificanceB", grLa, sep="" )
              }
              count=count+1
          }

      #################################################
      # extract the columns containing all ratios     #
      #################################################
      ratios <- proteinGroups[, ratioIdx]

      ##################################################################################################
      #      now extract only protein groups that have a ratio in ALL group/SILAC label combinations   #
      ##################################################################################################

      # determine the number of 'NA' value per protein group. if you observe an 'NA' in one of the ratios,
      # the protein group was not quantified in the respective group. here we consider only proteins
      # that were quantified in all groups
      NAnumb <- apply(ratios, 1, function(x) sum(is.na(x)))
      NAidx <- which(NAnumb == 0)

      # extract the respective protein groups
      proteinGroups <- proteinGroups[NAidx, ]
      ratios <- ratios[NAidx, ]

      # number of proteins that remain
      summary[4,1] <- sum(NAnumb == 0)

      ##################################################################################################
      #     extract protein groups having at least one significant ( significance B ) ratio over
      #     the time points
      ##################################################################################################

      #################################################
      # extract the columns containing significance B #
      #################################################
      significanceB <- proteinGroups[, ratioSigBIdx]

      #################################################
      # adust the p-values (Significance B)
      #################################################
      if(adjust){
          significanceB <- apply( significanceB, 2, p.adjust, method=method)  #p.adjust( significanceB, method=method )
      }


      # determine the number of 'significant' ratios per protein group
      sigBNumb <- apply( significanceB, 1, function(x) sum( x <= sigB) )

      # consider only protein groups with at least one significant ratio
      proteinGroups <- proteinGroups[ which(sigBNumb >= 1) , ]
      ratios <- ratios[which(sigBNumb >= 1), ]

      summary[5,1] <- sum(sigBNumb == 0)
      summary[6,1] <- dim(proteinGroups)[1]

      ########################################
      #     generate the output
      ########################################
      output <- vector("list", 2)
      names(output) <- c("ratios", "summary")

      output[["ratios"]] <- ratios
      output[["summary"]] <- summary

      return(output)
}


########################################################################################################
#
#
# arguments
#    data          - matrix returned by function 'extractRatios'
#    silacLabel    - doublets, triplets
#    groupLabel    - character vector specifying different groups
#                  - have to be consistent with the column names of 'data'
#                  - note that the group labels have to be in the RIGHT ORDER
#                    regarding to the time points, i.e. 024, 0824 (NOT: 0824, 024 )!!
#    Tp            - numeric vector specifying the timepoints, e.g. c(0,2,4,8,24)
#    TPLabel       - character vector of same length as 'Tp' specifiying the labels
#                    corresponding to the time points
#    commomTp      - common time point, currently NOT USED!
#
#
########################################################################################################
mergeDataSet <- function(data, silacLabel="triplets", groupLabel=c("A024", "A0824"), Tp=c(0, 2, 4, 8, 24), TpLabel=c("L", "M", "H", "M", "H"), commonTP=0, commonTPlabel="L"  , dot=T)
{
    ####################################################
    #
    ####################################################
    TpRelative = TpLabel[1]

    ####################################################
    # get the column names of respective SILAC  ratios
    ####################################################
    silacLabel = getSILAClabel( silacLabel  , dot=dot)

    colIdx <- c()
    count=1

    ##################################################
    #   matrix to store the time profiles
    ###################################################
    timeProfile <- matrix(0, nrow=dim(data)[1], ncol=length(Tp), dimnames=list( rownames(data), Tp ))


    ####################################################
    #  loop over the time points
    ####################################################
    count=1
    groupCount=1
    for(tp in Tp)
    {
        ##################################################
        # determine the respective group
        # NOTE: the group labels have to be in the correct order!!
        ##################################################
        if( count > ceiling(length(Tp)/2)  )
            groupCount = 2

        ##################################################
        # ratio in time point zero is always one ???!!!
        # maybe there is a better solution
        if(tp == 0)
            timeProfile[, as.character(tp)] <- rep(1, dim(timeProfile)[1])
        else{

            # get the column label for the current time point
            if(!dot) label <- paste("Ratio", TpLabel[count], TpRelative, "Normalized", groupLabel[groupCount], sep=""  )
            else stop("IMPLEMENT THE DOT!!!\n\n")

            timeProfile[, as.character(tp)] <- data[, label]
        }

        count <- count+1
    }
    return(timeProfile)

}
####################################################################################################################
#
#
# changelog: 20111215 implementation
####################################################################################################################
my.cmeans <- function(x, centers=NULL, iter.max=1000, maxK=NULL, verbose=FALSE, dist="euclidean", method="cmeans", m=NULL, rate.par=NULL, weights=1, control=list(reltol=(.Machine$double.eps))){


    # x <- rbind( matrix(rnorm(100,sd=0.3),ncol=2), matrix(rnorm(100,sd=0.3, mean=3),ncol=2), matrix(rnorm(100,sd=0.3, mean=10), ncol=2), matrix(rnorm(100,sd=0.3, mean=5),ncol=2) )


    require(e1071)

    #########################################################
    # estimate optimal value for fuzzyfication parameter 'm'
    # Schwämmle et.al 2010
    #########################################################
    if(is.null(m)){

        D <- dim(x)[2]
        N <- dim(x)[1]

        mt <- 1 + ( (1418/N) + 22.05) * D^-2 + ( (12.33/N) + 0.243 )*D^(-0.0406*log(N, exp(1))-0.1134)

        cat("optimal fuzzyfication parameter 'm' estimated as:", mt, "\n")
    }
    #########################################################
    # determine the number of clusters
    #########################################################
    if(is.null(centers)){

        # maximal cluster number to test
        if(is.null(maxK))
            K = floor( sqrt(dim(x)[1]) )
        else
            K = maxK

        # loop over different number of clusters
        min.centroid.dist <- vector("numeric", K-1)
        names(min.centroid.dist) <- 2:K

        clustering <- vector("list", K-1)
        names(clustering) <- 2:K

        for(cc in 2:K){

            # cluster the data
            cl.tmp <- cmeans(x, centers=cc, iter.max=iter.max, verbose=verbose, dist=dist, method=method, m=mt, rate.par=rate.par, weights=weights, control=control)
            clustering[[as.character(cc)]] <- cl.tmp

            # get centroids
            centroids.tmp <- cl.tmp$centers
            rownames(centroids.tmp) <- paste( "center_", rownames(centroids.tmp), sep=""  )

            # estimate distances to the centroids
            #dist.centroids <- as.matrix( dist(rbind(x, centroids.tmp), method=dist ))
            dist.centroids <- as.matrix( dist(centroids.tmp, method=dist ))
            diag(dist.centroids) <- NA

            # minimal centroid distance
            min.centroid.dist[as.character(cc)] <- min( dist.centroids^2, na.rm=T )
        }
        ####################################
        # calculate K_opt
        ####################################
        dist.decay <- vector("numeric", length(min.centroid.dist)-1)

        for(i in 1:(length(min.centroid.dist)-1))
            dist.decay[i] <- abs(min.centroid.dist[i] - min.centroid.dist[i+1])


        X11()
        # par(mfrow=c(1,2))
        plot(2:K, min.centroid.dist, main="Minimum centroid distance", xlab="cluster number", type="b", pch=20)
        #barplot(dist.decay)

        #
        dist.decay.norm <- dist.decay/max(dist.decay)
        #copt <- max(which(dist.decay.norm > 0.05)) + 1

        copt <- (dist.decay.norm > 0.05)
        copt[min(which(copt == FALSE)):length(copt) ] <- FALSE
        copt <- max(which(copt == TRUE)) + 1


        cat("optimal number of clusters:", copt, "\n")


    }
    return(clustering[[as.character(copt)]])
    #return(copt)
    #return(dist.centroids)
    #return(dist.decay.norm)
    #return(min.centroid.dist)

}

####################################################################################################################
#
#
# parameter
#    data - feature matrix without NAs, raw scale
#    norm - either 'NULL' or a character specifying the normalization method to use
#           transformation:
#               Zscore      - transform to Z-score row-wise
#           distance measures:
#               euclidean   - usual square distance bewteen two rows ( 2 norm )
#               maximum     - maximum distance between two components of the rows ( supremum norm )
#               canberra
#           correlation:
#               spearman    - Spearman's rank correlation coefficient
#               pearson     - pearson correlation coefficient
#    log       - logical indicating whether the data should be transformed to log-scale
#    logBase   - numerical, the base of logarithm
#    plotOnRawScale   - logical, if TRUE the y-axis of the time profile plot is the non-normalized scale of 'data'
#    ylim
#    algorithm
#    nClust
#    file
#    pdf
#    nrows     - numeric, number of rows in the plot
#    ...
# value
#   the object returned by
#   the cluster algorithm
#
#
# changelog:
#             20100721 added parameter 'hc'; can be an 'hclust' object
#                      Zscore transformation is done by function 'scale'
#             20120320 parameter 'alpha': for non-colorgradients, e.g. hclust, kmeans, ...
#             20120221 parameter 'legend': if TRUE a legend will be plotted, at least in case
#                      of 'cmeans' and 'conensusclust'
#############################################################################################################
timeProfileClusterPlot <- function(data, hc=NULL, norm=c("none", "Zscore"), log=T, logBase=2, plotOnRawScale=F, addAverageProfile=T, y.lim=NULL, algorithm=c("cmeans", "mycmeans", "kmeans", "consensus", "hclust"), nClust=3, file="clusterTimeProfiles.pdf", pdf=F, ylab="Change in expression", col=c("green", "magenta", "darkred"), nrows=2, minMembShip=0.5, cmeansm=2, alpha=40, legend=T,...)
{
    alpha=80

    if(!is.null(hc) && class(hc) == "hclust")
        algorithm="hclust"
    if(!is.null(hc) && class(hc) == "fclust")
        algorithm="mycmeans"

    ####################################
    #  numeric matrix
    ####################################
    if(mode(data) != "numeric")
        data <- data.matrix(data)

    ####################################
    #      log transformation
    ####################################
    if(log)
        data <- log(data, logBase)
    # data before normalization
    dataRaw <- data

    #####################################
    #        normalize the data
    #####################################
    norm=match.arg(norm)

    if(norm != "none")
    {
       ################################
       #    Z-score
       ################################
       if(norm == "Zscore")
           data <- t( scale(t(data), center=T, scale=T))
       else if(norm == "pearson" | norm == "spearman")
           data <- cor(t(data), method=norm)
       else{
           data <- dist(data, method=norm)
       }
    }
    ######################################################
    #
    #  different colors for each cluster
    #
    ######################################################
    if(nClust == length(col)){
           COLORS=col
    } else {
           COLORS=rep("black", nClust)
    }
    for(ii in 1:length(COLORS))
           COLORS[ii] <- my.col2rgb( as.character(COLORS[ii], alpha=alpha)  )

    #######################################################################################
    #
    #                            perform the clustering
    #
    #######################################################################################
    if(is.null(hc)){

       # determine cluster algorithm
       algorithm=match.arg(algorithm)

       ################################
       # kmeans consensus clustering
       ################################
       if(algorithm == "consensus")
       {
           require(compdiagTools)
           clustering <- consensusCluster( t(data),  nclass=nClust, ...)
       }
       ################################
       # kmeans consensus clustering
       ################################
       if(algorithm == "consensus+") {
           require(ConsensusClusterPlus)
       }
       #################################
       #
       #################################
       if(algorithm == "mycmeans"){

           clustering <- my.cmeans(data, ...)
       }
       #################################
       #         cmeans
       #################################
       if(algorithm == "cmeans"){
           require(e1071)
           #clustering <- eval(parse( text=paste(algorithm, "(data, centers=", nClust,", ...)", sep=""  )  ))

           N <- dim(data)[1]
           D <- dim(data)[2]

           cmeansm <-  1 + ( (1418/N) + 22.05) * D^-2 + ( (12.33/N) + 0.243 )*D^(-0.0406*log(N, exp(1))-0.1134)

           cat("optimal fuzzyfication parameter 'm'=", cmeansm, "\n")


           clustering = cmeans(data, centers=nClust, m=cmeansm, ...)

       }
       #################################
       #        k-means
       #################################
       if(algorithm == "kmeans"){
           clustering = kmeans(data, centers=nClust, ...)
       }

       ########################################
       # cluster assigment
       ########################################
       cl <- clustering$cluster

    } else{ # end if is.null(hc)

       ##################################################
       #         hierarchical clustering
       ##################################################
       if(class(hc) == "hclust"){
           cl <- cutree(hc, nClust)
       }
       if(class(hc) == "fclust"){
           cl <- hc$cluster
           clustering=hc
       }

    }

    ##########################################################
    #
    #                   plot the stuff
    #
    ##########################################################

    # determine the number of clusters
    nC <- length(levels(as.factor(cl)))

    # determine whether the normalized data shall be plotted
    if( plotOnRawScale )
        dataPlot <- dataRaw
    else
        dataPlot <- data

    #################################
    # ylim, equal for all clusters
    #################################
    if(is.null(y.lim)){
            y.lim = rep(max(abs(dataPlot), na.rm=T), 2)
            y.lim[1] <- -1*y.lim[1]
        }



    #################################
    # open pdf file
    #################################
    if(pdf) pdf(file=file, width=ceiling(nC/nrows)*6, height= 5*nrows)

    par(mfrow=c(nrows, ceiling(nC/nrows)))

    #################################
    #        loop over the cluster
    #################################
    ccount=1
    for( cluster in levels(as.factor(cl)) ){

        # determine the members of the current cluster
        clusterMembers = rownames(dataPlot)[ which( cl == cluster ) ]


        ################################################
        #
        ################################################
        #COLOR <- my.col2rgb( "black", alpha=alpha )
        COLOR = COLORS[ccount]

        ################################################
        # color gradient respective to membership
        ################################################
        if(algorithm=="cmeans" | algorithm=="mycmeans"){

            item.score.range <- c(0,1)
            COLOR <- colorGradient( max(clustering$membership[ clusterMembers[1], ]) , mi=item.score.range[1], ma=item.score.range[2], scheme=col )

            COLOR <- my.col2rgb( COLOR, alpha=alpha )

            # calculate the mean cluster membership score
            maxScore <- apply( clustering$membership[clusterMembers, ],1, max  )

        }

        ################################################
        # color gradient respective to membership
        ################################################
        if(algorithm=="consensus")
        {

            # min and max item consensus
            #item.score.range <- range(clustering$itemConsensus[ clusterMembers ])
            item.score.range <- range(clustering$itemConsensus)
            #item.score.range <- c(0,1)
            COLOR <- colorGradient( clustering$itemConsensus[ clusterMembers[1] ] , mi=item.score.range[1], ma=item.score.range[2], scheme=col )

            COLOR <- my.col2rgb( COLOR, alpha=clustering$itemConsensus[ clusterMembers[1] ] * alpha   )


        }


        # plot the first profile....
        plot(dataPlot[clusterMembers[1], ], type="l", ylim=y.lim, main=paste("cluster",cluster, "(",length(which(cl==cluster)),")") , col=COLOR, xaxt="n", ylab=ylab, lwd=2, xlab="" )

        if((algorithm=="cmeans" | algorithm=="mycmeans") & legend)
            legend("topright", legend=c( paste( "mean: ", round(mean(maxScore),2) ) , paste("median: ", round(median(maxScore),2)), paste("N <", minMembShip, ":", sum(maxScore < minMembShip)) ), title="Membership" )

        if(algorithm=="consensus" & legend)
            legend("topright", legend=c( paste( "cluster consensus: ", round(clustering$clusterConsensus[as.numeric( cluster) ], 3) ) ) )

        axis(1, at=seq(1,dim(data)[2]), labels=colnames(data), las=2)


        # if there are other cluster memebrs....
        if( length(clusterMembers) > 1 ) {
        # ... and loop over the remaining members
           for( i in 2:length(clusterMembers)){

                if(algorithm=="cmeans" | algorithm=="mycmeans"){
                   COLOR <- colorGradient( max(clustering$membership[clusterMembers[i],]), mi=item.score.range[1], ma=item.score.range[2], scheme=col )
                      COLOR <- my.col2rgb( COLOR, alpha=alpha )
                }
                if(algorithm=="consensus"){
                   COLOR <- colorGradient( clustering$itemConsensus[clusterMembers[i]], mi=item.score.range[1], ma=item.score.range[2], scheme=col )
                   COLOR <- my.col2rgb( COLOR, alpha= clustering$itemConsensus[clusterMembers[i]] * alpha )
                 }
                lines(dataPlot[clusterMembers[i], ], type="l", main=paste("cluster",cluster, "(",length(which(cl==cluster)),")") , col=COLOR, lwd=2)

              }

           # add the average profile
           if(addAverageProfile)
                  lines( apply(dataPlot[clusterMembers, ], 2, mean, na.rm=T), lwd=2, col="red" )
        }
        ccount = ccount + 1
    }
    par(mfrow=c(1,1))

    # close pdf file
    if(pdf)dev.off()



    if(is.null(hc)){
        return(clustering)
    } else {
        return(cl)
    }
}
###############################################################################################################
#
#
# dataGroup1     - time profile group 1, e.g. normal
# dataGroup2     - time profile group 2, e.g. tumor
# clustering     - result of cmeans clustering of group 1
# col            - color scheme to use
# file           - pdf filenames
# pdf            - logical
#
###############################################################################################################
TimeProfilePlotTwoGroups <- function( dataGroup1, dataGroup2, clustering, col="topo", file="timeProfilesNormalvsHeLa.pdf", pdf=T, addAverageProfile=T )
{
    # Problem!! Ratios are only extracted from proteins that were quantified in all timepoints, and thus have a
    # ratio in all time point. However, this can vary between the groups
    # get all common proteins -> the number of considered proteins decreases

    commonProteins <- intersect( rownames(dataGroup1), rownames(dataGroup2) )

    dataGroup1 <- dataGroup1[commonProteins, ]
    dataGroup2 <- dataGroup2[commonProteins, ]



    # get the cluster membership
    cl <- clustering$cluster[commonProteins]

    ######################################
    #               plot
    ######################################
    # open pdf file
    if(pdf) pdf(file=file, width=8, height=4*length(levels(as.factor(cl))))
    par(mfrow=c( length(levels(as.factor(cl))), 2  ), mar=c(2.5,1.5,1.5,2.5))

    # loop over the cluster
    for(cluster in levels(as.factor(cl)))
    {
         # determine the members of the current cluster
        clusterMembers = rownames(dataGroup1)[ which( cl == cluster ) ]

        # color gradient respective to membership
        COLOR <- colorGradient( max( clustering$membership[clusterMembers[1], ] ), scheme=col )

        # calculate the mean cluster membership score
        maxScore <- apply( clustering$membership[clusterMembers, ],1, max  )

        # determine ylim
        y.lim = range(dataGroup1[clusterMembers, ])

         # plot the first profile....
        plot(dataGroup1[clusterMembers[1], ], type="l", ylim=y.lim, main=paste("cluster",cluster, "(",length(which(cl==cluster)),")") , col=COLOR, xaxt="n" , ylab="change in regulation")
        #legend("topright", legend=c( paste( "mean: ", round(mean(maxScore),2) ) , paste("median: ", round(median(maxScore),2)) ) )
        axis(1, at=seq(1,dim(dataGroup1)[2]), labels=colnames(dataGroup1))

        # ... and loop over the remaining members
        for( i in 2:length(clusterMembers)){
             COLOR <- colorGradient( max( clustering$membership[clusterMembers[i], ] ), scheme=col )

            lines(dataGroup1[clusterMembers[i], ], type="l", main=paste("cluster",cluster, "(",length(which(cl==cluster)),")") , col=COLOR)

        }

         # add the average profile
        if(addAverageProfile)
            lines( apply(dataGroup1[clusterMembers, ], 2, mean), lwd=2, col="red" )

        #####################################################################
        #  plot the time profiles of the same proteins in the second group
        #####################################################################
        plot(dataGroup2[clusterMembers[1], ], type="l", ylim=y.lim, main=paste("cluster",cluster, "(",length(which(cl==cluster)),")") , xaxt="n", col=COLOR , ylab="change in regulation")
         axis(1, at=seq(1,dim(dataGroup2)[2]), labels=colnames(dataGroup2))

        for( i in 2:length(clusterMembers))
        {
             COLOR <- colorGradient( max( clustering$membership[clusterMembers[i], ] ), scheme=col )
             lines(dataGroup2[clusterMembers[i], ], type="l", main=paste("cluster",cluster, "(",length(which(cl==cluster)),")"), col=COLOR)

         }
         # add the average profile
        if(addAverageProfile)
            lines( apply(dataGroup2[clusterMembers, ], 2, mean), lwd=2, col="red" )
    }

    par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)

    if(pdf)
        dev.off()


    output <- vector( "list",3 )
    names(output) <- c("group1", "group2", "cluster")

    output[["group1"]] <- dataGroup1
    output[["group2"]] <- dataGroup2
    output[["cluster"]] <- cl

    return(output)
}

#########################################################################
#   get a table containing the cluster members and some annotation
#   obtained from the MaxQuant output
#
# parameter
#   - clustering       - object returned by the cmeans cluster algorithm
#   - proteinGroups    - MaxQuant output
#   - scoreThr         - numeric, indicates the minimum membership score
#                        a member must have in order to be returned
#   - n                - numeric, specifies how many proteins per cluster
#                        will be returned
#                      - if  'NULL', all remaining members will be returned
#
# value
#   - a table
#########################################################################
getTopMembers <- function( clustering, proteinGroups, n=NULL,  scoreThr = 0.8, digits=NULL, chop=NULL, dot=F, anno=c("ProteinNames", "GOBPNames"))
{
    ############################
    #  annotoation
    #############################
    #if(dot)
    #    anno <- c( "Protein.Names" )
    #else
    #    anno <- c(  "ProteinNames"  )

    # get the number of cluster
    K <- length( levels( as.factor( clustering$cluster  )  )  )

    # list to store the particular tables
    members <- vector("list", K)
    names(members) <- paste( "cluster", K )


    # loop over the cluster
    #  - get the mebers
    #  - sort and filter according to membership score
    #  - build up the table
    for( cl in 1:K )
    {
        ############################
        # get the members
        ############################

        # determine the cluster members
        MembIdx <- which(clustering$cluster == cl  )

        # get the membership score of each member in the
        # cluster
        MembScore <- clustering$membership[MembIdx, cl]

        # now filter the members according to a defined membership score threshold
        members[[cl]] <- MembScore[ which(MembScore >= scoreThr)  ]

        # sort them according to membership score
        members[[cl]] <- members[[cl]][ order(members[[cl]], decreasing=T)  ]

        # return only the 'n' top hits
        if( !is.null(n) )
        {
            if(length( members[[cl]] ) > 0)
                members[[cl]] <- members[[cl]][1:min(n, length(members[[cl]]))]
        }

        ##############################
        # build the table
        ##############################
        if( cl == 1  ){
            table <- cbind(  rep(cl, length(members[[cl]])) ,members[[cl]],  proteinGroups[ names(members[[cl]]),  anno ]  )

        }
        else
            table <- rbind(table , cbind(  rep(cl, length(members[[cl]])) ,members[[cl]], proteinGroups[ names(members[[cl]]),  anno ]  ) )
    }


    # add column names
    colnames(table)=c("Cluster", "Score", anno)


    ########################
    #    format the table
    ########################
    if( !is.null(digits) )  # round the score
         table[, "Score"] <- round( as.numeric( table[, "Score"]), digits  )
    if( !is.null(chop) )
    {
        if(length(anno) == 1)
            table[, anno] <- unlist( lapply( table[, anno], chopString, nChar=chop )  )

        else
            for(a in anno)
                table[, a] <- unlist( lapply( table[, a], chopString, nChar=chop )  )

    }

 return(table)

}



#############################################################################################
#
#                  Gene Ontology Enrichment Analysis based on the 'topGO' package
#
#
# arguments:
#  proteins       - character vector, the protein ids of 'interesting' proteins
#  GOmap          - dataframe or matrix,
#                 - rows are a proteins; columns are GO terms.
#                 - If there are multiple GO terms for a single protein, the
#                   terms have to be separated by ';'.
#                 - rownames are the protein ids!
#                 - there have to be at least three columns, namly 'BP.term', 'MF.term', 'CC.term'
#                 - e.g.
#                                      BP.term                 BP.description                  CC.term  ...
#                   Contig0.snap1   GO:0040010;GO:0018991    "positive regulation..."        GO:0005634
#  test           - character, specifies the direction of the test, i.e.
#                      overrepresentation, underrepresentation, two.sided
#  algorithm      - character, specifies the algorithm to use
#                 - 'classicCount' - determine the contingency table on all GO terms
#                 - 'elimCount'    - eliminates redundancy to get the most specific term
#                 - 'weightCount'
#  pvalThr        - numeric, specifies the maximal p-value of the top nodes from which
#                   the resulting graphs are build
#  label          - character, the label is used to create a diretory where the results are stored
#  adjust         - character, specifies the method used for multiple testing correction
#
#  ...            - further arguments passed to classes for the test statistics,
#                   e.g. 'cutOff' -> 'elimCount': nodes are significant below this p-value cut off
#                        'sigRatio' -> KS-test
#
#
# value:           none
#
#
# changelog:  20090810   - implementation based on a script
#             20090812   - some documentation
#             20090820   - fdr
#             20090825   - two sided fisher test
#                        - added a column indicating the direction of regulation, i.e. over/under
#             20090911   - the function now returns the tables for each category
#             20091102   - annotation summary, i.e. # proteins in test set, # annotated proteins in
#                          test set, # proteins in background
#
##############################################################################################
topGOAnalysis <- function( proteins, GOmap, test=c("over", "under", "two.sided"), algorithm=c( "classicCount", "elimCount", "weightCount" ), adjust=c("BH", "BY", "fdr", "holm", "hochberg", "hommel", "none"), pvalThr=.01, label="test", ...  )
{
    require(topGO)
    require(Rgraphviz)

    #######################################################
    #  test statistics
    #######################################################
    test <- match.arg(test)

    #######################################################
    # algorithm
    #######################################################
    algorithm <- match.arg(algorithm)

    #######################################################
    # method for adjusting p-values
    #######################################################
    adjust <- match.arg(adjust)

    #######################################################
    # specify a directory where the results will be stored
    #######################################################
    resDir = paste( label, algorithm, test, sep="_")
    dir.create(resDir)

    #######################################################
    # define the funtions for testing depletion
    # and two sided test
    #######################################################
    if(test != "over" )
    {
        ###############################################
        # define a 2-sided test statistic
        ###############################################

        ## define the test statistic which will detect over- and underrepresentation
        if(!isGeneric("GOFisherTestTwoSided"))
            setGeneric("GOFisherTestTwoSided", function(object) standardGeneric("GOFisherTestTwoSided"))


        ## set the method
        setMethod("GOFisherTestTwoSided", algorithm,
          function(object) {

           contMat <- contTable(object)
                    if(all(contMat == 0))
             p.value <- 1
           else
             p.value <- fisher.test(contMat, alternative = "two.sided")$p.value    ## "greater" is for over-, "less" for under-, and
                                                                                   ## "two.sided" is for both alternatives

           return(p.value)
         })

        ###############################################
        # underepresentation
        ###############################################
        ## define the test statistic which will detect underrepresentation
        if(!isGeneric("GOFisherTestUnder"))
            setGeneric("GOFisherTestUnder", function(object)
        standardGeneric("GOFisherTestUnder"))

        # set the method
        setMethod("GOFisherTestUnder", algorithm,
         function(object) {

           contMat <- contTable(object)
                    if(all(contMat == 0))
             p.value <- 1
           else
             p.value <- fisher.test(contMat, alternative = "less")$p.value


           return(p.value)
         })



    }
    #################################################
    #             test statisitics
    #################################################
    if(test == "over") test.statistic <- new( algorithm, testStatistic = GOFisherTest, name = "Fisher test Enrichment", ...)
    if(test == "under") test.statistic <- new( algorithm, testStatistic = GOFisherTestUnder, name = "Fisher test Depletion", ...)
    if(test == "two.sided") test.statistic <- new( algorithm, testStatistic = GOFisherTestTwoSided, name = "Fisher test Two-Sided", ...)


    #####################################################
    #                  get GO terms
    #####################################################
    BP <-strsplit( GOmap[ , "BP.term"], ";")
    MF <-strsplit( GOmap[ , "MF.term"], ";")
    CC <-strsplit( GOmap[ , "CC.term"], ";")

    ######################################################
    # vector indicating wether a protein belongs to the
    # test set or to background
    ######################################################
    proteins.In <- factor(as.numeric( rownames(GOmap) %in% proteins  ))
    names(proteins.In) <- rownames( GOmap )

    ######################################################
    # determine how many proteins are annotated in both
    # sets (test-set and background)
    ######################################################
    nProtAnnotated.test <- sum( rownames(GOmap) %in% proteins )

    ######################################################
    # determine how many proteins are annotated at all,
    # i.e. the background
    ######################################################
    nProtAnnotated.bg <- dim( GOmap  )[1]

    annot.summary <- c( length(proteins), nProtAnnotated.test, nProtAnnotated.bg  )
    names(annot.summary) <- c("test.set", "test.set.annotated", "background")

    #######################################################
    #         build the topGOdata objects
    #######################################################
    GOdata.BP <- new( "topGOdata", ontology="BP", allGenes=proteins.In, annot=annFUN.gene2GO, gene2GO=BP )
    GOdata.CC <- new( "topGOdata", ontology="CC", allGenes=proteins.In, annot=annFUN.gene2GO, gene2GO=CC )
    GOdata.MF <- new( "topGOdata", ontology="MF", allGenes=proteins.In, annot=annFUN.gene2GO, gene2GO=MF )

    #######################################################
    #                run the analysis
    #######################################################
    resultBP <-  getSigGroups( GOdata.BP, test.statistic  )
    resultCC <-  getSigGroups( GOdata.CC, test.statistic  )
    resultMF <-  getSigGroups( GOdata.MF, test.statistic  )


    ######################################################
    #                 adjust p-values
    #######################################################
    BP.pval.adjust <- p.adjust(  score(resultBP), adjust  )
    resultBP.fdr <- resultBP
    score(resultBP.fdr) <- BP.pval.adjust

    CC.pval.adjust <- p.adjust(  score(resultCC), adjust  )
    resultCC.fdr <- resultCC
    score(resultCC.fdr) <- CC.pval.adjust

    MF.pval.adjust <- p.adjust(  score(resultMF), adjust  )
    resultMF.fdr <- resultMF
    score(resultMF.fdr) <- MF.pval.adjust


    #######################################################
    #              export the results
    #######################################################

    ####################
    #     tables
    ####################

    ## BP
    # get the table of all terms
    resultBP.mat <- GenTable(GOdata.BP, resultBP, topNodes=length(score(resultBP)) )
    # add adjusted p-values
    resultBP.mat <- cbind(resultBP.mat, BP.pval.adjust[resultBP.mat[, "GO.ID"]]  )
    # add the direction of regulation
    resultBP.mat <- cbind(resultBP.mat, ifelse(resultBP.mat[, "Significant"] > resultBP.mat[, "Expected"], "over", "under"  ) )
    colnames(resultBP.mat)[(dim(resultBP.mat)[2]-2):dim(resultBP.mat)[2]  ] <- c("p-value", paste("adjusted p-value (", adjust,")", sep=""), "over/under")
    # export
    write.table(resultBP.mat, file=paste(resDir, "/signifGO_BP_",algorithm,"_",test,".txt", sep=""),  sep="\t", quote=F, row.names=F  )

    # CC
    resultCC.mat <- GenTable(GOdata.CC, resultCC, topNodes=length(score(resultCC)) )
    resultCC.mat <- cbind(resultCC.mat, CC.pval.adjust[resultCC.mat[, "GO.ID"]]  )
    resultCC.mat <- cbind(resultCC.mat, ifelse(resultCC.mat[, "Significant"] > resultCC.mat[, "Expected"], "over", "under"  ) )
    colnames(resultCC.mat)[(dim(resultCC.mat)[2]-2):dim(resultCC.mat)[2]  ] <- c("p-value", paste("adjusted p-value (", adjust,")", sep=""), "over/under")
    write.table(resultCC.mat, file=paste(resDir, "/signifGO_CC_",algorithm,"_", test,".txt", sep=""),  sep="\t", quote=F, row.names=F  )

    # MF
    resultMF.mat <- GenTable(GOdata.MF, resultMF, topNodes=length(score(resultMF)) )
    resultMF.mat <- cbind(resultMF.mat, MF.pval.adjust[resultMF.mat[, "GO.ID"]]  )
    resultMF.mat <- cbind(resultMF.mat, ifelse(resultMF.mat[, "Significant"] > resultMF.mat[, "Expected"], "over", "under"  ) )
    colnames(resultMF.mat)[(dim(resultMF.mat)[2]-2):dim(resultMF.mat)[2]  ] <- c("p-value", paste("adjusted p-value (", adjust,")", sep=""), "over/under")
    write.table(resultMF.mat, file=paste(resDir, "/signifGO_MF_",algorithm,"_", test,".txt", sep=""),  sep="\t", quote=F, row.names=F  )

    #####################
    #    GO graph
    #####################

    # determine the number of 'top N' nodes according
    # to a p-value threshold, 'pvalThr'
    if(pvalThr < 1){
        topBP <- length(which( score(resultBP.fdr) < pvalThr  ))
        topCC <- length(which( score(resultCC.fdr) < pvalThr  ))
        topMF <- length(which( score(resultMF.fdr) < pvalThr  ))
    }
    else {
        topBP <- topCC <- topMF <- pvalThr
    }
    printGraph(GOdata.BP, resultBP.fdr, firstSigNodes=topBP, pdfSW=T, useInfo="all", fn.prefix=paste(resDir, "/BP_adj_pVal_", pvalThr, sep=""))
    printGraph(GOdata.CC, resultCC.fdr, firstSigNodes=topCC, pdfSW=T, useInfo="all", fn.prefix=paste(resDir, "/CC_adj_pVal_", pvalThr ,sep=""))
    printGraph(GOdata.MF, resultMF.fdr, firstSigNodes=topMF, pdfSW=T, useInfo="all", fn.prefix=paste(resDir, "/MF_adj_pVal_", pvalThr ,sep=""))


    #########################################
    # some output
    #########################################
    output <- list()
    output[[1]] <- resultBP.mat
    output[[2]] <- resultMF.mat
    output[[3]] <- resultCC.mat

    output[[4]] <- annot.summary

    names(output) <- c("BPtable", "MFtable", "CCtable", "Annotation.Summary")

    return(output)
}



######################################################
ORA.C <- function(proteins, clustering, test="exact", verbose=T)
{

    # get the number of clusters
    #K <- length( levels( as.factor(clustering) )  )
    K <- length(unique( as.character(clustering)  ))

    # get all annotation terms of all proteins
    annoTerms <- unique( unlist( lapply( proteins, function(x) unlist(strsplit( as.character(x), ";"  ))  )))

    # get all proteins that were annoted by at least one term
    protIdx <- which( unlist( lapply(proteins, function(x) nchar(as.character(x)))  ) > 0 )
    proteins <- proteins[ protIdx  ]
    clustering <- clustering[ protIdx ]

    # number of all proteins that were annotated
    N <- length(proteins)

    if(verbose)
        cat("considering", N,"annotated protein groups\n")


    # output
    pValMat <- matrix(0.0, nrow=length(annoTerms), ncol=K, dimnames=list( annoTerms, paste("cluster",1:K )))

    tmp <- .C("ora", as.character(annoTerms), as.character(proteins), as.character(clustering), as.integer(length(annoTerms)), as.integer(length(proteins)), as.integer(K), as.character(test), matrix(0.0, nrow=length(annoTerms), ncol=K), as.double( rep(0, length(annoTerms)*K) )  )


    #####################################
    # matrix containing the raw p-values
    #####################################
    pValMat <- matrix( tmp[[9]], nrow=length(annoTerms), ncol=K, byrow=T )
    dimnames(pValMat) <- list( annoTerms, paste("cluster",1:K ))



    return(pValMat)
}

##################################################
#                  Z-score-transformation
#
# x   - vector of raw scores
#
# the values of x are tranformed such that the mean
# of the resulting distribution is 0 and sd equals
# one
##################################################
Z.score <- function(x, na.rm=F){

    return( (x - mean(x, na.rm=na.rm))/sd(x, na.rm=na.rm))

}

####################################################################
#                     getSILAClabel
# 'translates' a string characterizing a SILAC experiment to the
# column names of interest of MaxQuant output
#
####################################################################
getSILAClabel <- function( silacLabel = c("none", "doublets", "triplets") , dot=T)
{
      silacLabel = match.arg(silacLabel)
      if(silacLabel == "none" ) silacLabel = ""
      if(silacLabel == "doublets") silacLabel = c("H.L")
      if(silacLabel == "triplets") silacLabel = c("M.L", "H.L", "H.M")

      if(!dot) silacLabel=gsub( ".","", silacLabel, fixed=T )

      return(silacLabel)
}



#######################################################################
#                        colorGradient
#
# arguments
#   x          numeric, value between 'mi' and 'ma'
#   scheme     color scheme to use
#
# value
#   RGB color corresponding to x
#
#
#######################################################################
colorGradient <- function(x, mi=0, ma=1, scheme=c("heat", "topo", "terrain", "grey", "cm"), eps=0.02 )
{
    # intervall
    int <- seq(mi, ma, length.out=1000  )

    # color scheme
    if(length(scheme) == 1){

      scheme=match.arg(scheme)

      if(scheme=="heat")
         cols <- heat.colors(length(int))
      if(scheme=="topo")
        cols <- topo.colors(length(int))
      if(scheme=="terrain")
        cols <- terrain.colors(length(int))
      if(scheme=="grey")
        cols <- grey.colors(length(int))
      if(scheme=="cm")
        cols <- cm.colors(length(int))
    }
    if(length(scheme) == 3){
        require(gplots)
        cols <- colorpanel(length(int), scheme[1], scheme[2], scheme[3])
    }
    if(length(scheme) > 3){

        int <- seq(mi, ma, length.out = length(scheme))
        cols <- scheme

    }

    return( cols[ max( which( x >= (int-eps) ) )  ] )
}

#################################################
#   Given a string and a number of characters
#   the function chops the string to the
#   specified number of characters and adds
#   '...' to the end.
#
# parameter
#   string     - character
#   nChar      - numeric
#
# value
#   string of 'nChar'´characters followed
#     by '...'
#
##################################################
chopString <- function(string, nChar)
{
    string <- as.character(string)

    if(nchar(string) <= nChar)
        return(string)

    return( paste( paste( unlist( strsplit( string, "" )  )[1:nChar], collapse=""), "..." ))

}

########################################################################################################################
#
#
#
#######################################################################################################################
myScan <- function(file, sep="\t", what="character", row.names=T, col.names=T,...)
{
    # get the number of columns
    #nCol <- length(scan(file=file, sep=sep, what=what, nlines=1, quiet=T))
    nCol <- dim(read.delim(file=file, sep=sep, row.names=1, nrows=1))[2]


    # now import all
    data <-scan(file=file, sep=sep, what=what, quiet=T,...)

    data <- matrix(data, ncol=nCol, byrow=T)

    data <- as.data.frame(data, stringsAsFactors=F, row.names=data[, 1])

    if(col.names){
        colnames(data) <- make.names(data[1, ])
        data <- data[-c(1),]
    }

    return(data)
}
########################################################################################################################
#
#
#
########################################################################################################################
importMQtable <- function(file, table="msms", what=NULL){

    if(table == "msms"){

    }

}


####################################################################################################################
#
#
# changelog: ??        implementation
#            20120423  decoy hist are removed
#                      second legend showing the number of dp and decoy hits
####################################################################################################################
DPmassHist <- function(allPeptides, massRange=c(0,100), sigNumb=100, main="Dependent Peptide Mass Histogram", ...)
{

    ###############################################
    # get all dependent peptides
    ###############################################
    dp <- allPeptides[!is.na(allPeptides[, "DP.Mass.Difference"]),  ]

    # number of dp decoy hits
    dp.decoy <- grep("\\+", dp$DP.Decoy)

    dp.all <- dp

    if(length(dp.decoy)>0)
        dp <- dp[-dp.decoy, ]

    # sort accoring to mass difference
    dp <- dp[order(dp[, "DP.Mass.Difference"]), ]

    # filter accordig to mass range
    dp <- dp[  dp[, "DP.Mass.Difference"] >= massRange[1] & dp[, "DP.Mass.Difference"] <= massRange[2],  ]


    ##############################################
    # get all dependent peptide clusters
    ##############################################
    DpClustIDs <- unique(as.character(dp[, "DP.Cluster.Index"]))


    ##############################################
    #  create a matrix containing all necessary
    #  information for the histogram
    ##############################################
    dpMassHistMatrix <- matrix("", nrow=3, ncol=length(DpClustIDs), dimnames=list( c("DP.Mass.Difference", "DP.Modification", "counts"), DpClustIDs ))

    # loop over the dp cluster ids
    for(id in DpClustIDs)
    {
        # get the index of all peptides in the curretn dp cluster
        dpClustidx <- which( dp[, "DP.Cluster.Index"] == id  )

        # number of peptides in cluster
        idCount <- length( dpClustidx )

        # average mass difference
        idMass <- mean( dp[ dpClustidx,  "DP.Mass.Difference"]  )

        # dp modification
        idMod <- paste( unique(dp[dpClustidx, "DP.Modification"   ]), collapse=";"  )

        # store all this stuff in the matrix
        dpMassHistMatrix[, id] <- c( idMass, idMod, idCount  )
    }


    ##########################################################
    #    make the plot
    ##########################################################
    plot( as.numeric( dpMassHistMatrix["DP.Mass.Difference", ] ), as.numeric(dpMassHistMatrix["counts", ]), type="h", axes=F, xlab="Mass Difference (Da)", ylab="Number of Spectra", main=main, ylim=c(0, max(as.numeric(dpMassHistMatrix["counts",]))+max(as.numeric(dpMassHistMatrix["counts",]))*0.2  ),...  )

    # add the modifications
    dpMassHistMatTmp <- dpMassHistMatrix[ , which(as.numeric(dpMassHistMatrix["counts", ]) >= sigNumb) ]
    # remove unknowns
    dpMassHistMatTmp <- dpMassHistMatTmp[ , which(nchar(dpMassHistMatTmp[2,]) > 0)  ]
    # if there are several modification (or aa substitution/truncation/...) choose only one for plotting
    mods <- strsplit(dpMassHistMatTmp[2, ],";")
    #mods <- unlist(lapply(mods, function(x){ idx <- grep("Peptide",x);if(length(idx)>0)return(x[-idx]);if(length(idx)==0)return(x) } ))

    mods <- unlist(lapply(mods, function(x){ idx <- grep("Peptide",x); if(length(idx)==length(x))return(paste(x, collapse=";")); if(length(idx)==0)return(x); return(paste(x[-idx], collapse=";"))} ) )

    dpMassHistMatTmp[2,] <- mods

    # if there are any mass differences meeting the criteria 'sigNumb'
    if(!is.null(dim(dpMassHistMatTmp)))
    {
        for(i in 1:dim(dpMassHistMatTmp)[2] )
            text( as.numeric(dpMassHistMatTmp["DP.Mass.Difference" , i]), as.numeric(dpMassHistMatTmp["counts" , i]), ifelse(nchar(dpMassHistMatTmp["DP.Modification", i]) > 30, chopString(dpMassHistMatTmp["DP.Modification" , i], 30), dpMassHistMatTmp["DP.Modification" , i]), cex=1.2, pos=3, srt=90, offset=3  )

    ############
    # axes
    ############
    # x-axis, tickmarks at 'massRange' and at each mass whose frequency is >= sigNumb
        axis(side=1, at=round( as.numeric( c(massRange[1], dpMassHistMatTmp["DP.Mass.Difference", ], massRange[2])),2  ))
        # y-axis
        maxCount <- max(as.numeric(dpMassHistMatTmp[ "counts", ]))
        scale <- eval(parse(text=paste("1e-", (nchar(as.character(maxCount))-1), sep="")))
        ymax <- ceiling( maxCount*scale ) / scale
        #cat(ymax,"\n")
        axis(side=2, at= round(seq(0, ymax,  length.out=10)), las=2 )
    }
    ##########################################################
    # legend
    ##########################################################
    legendTXT <- apply(dpMassHistMatTmp,2, function(x) paste(x[2],"  (", x[3], "spectra )")  )
    legend("topright", legend=legendTXT, bty="n")

    legend("topleft", legend=c(paste("# total dependent peptides:", nrow(dp.all)), paste("# decoy hits:", length(dp.decoy), "(", round( length(dp.decoy)/nrow(dp.all)*100, 1 ),"%)"), paste("# dp between", massRange[1],"and", massRange[2], "Da:", nrow(dp)) ), bty="n")

    return(dpMassHistMatrix)
}


##########################################################################################
#
#                                  Venn Diagram
#
#
# 20090604   - implementation
# 20090618   - added code for three sets
##########################################################################################
vennDiag <- function(x, y, z=NULL, names=c("set 1", "set 2"), ...)
{
    require(limma)

    # check if there are 2 or 3 sets
    n = ifelse(is.null(z),  2, 3)

    if(n == 2)
    {
        allItems <- union(x,y)

        allItemsMat <- matrix(0, ncol=n, nrow=length(allItems), dimnames=list(allItems, names)  )

        allItemsMat[ , 1] <- as.numeric( allItems %in%  x )
        allItemsMat[ , 2] <- as.numeric( allItems %in%  y )

        vennDiagram(vennCounts( allItemsMat  ), ...)

        text(-1, -2, paste(length(unique(x))), cex=1.5)
        text(1, -2, paste(length(unique(y))), cex=1.5)
        text(0, -2.3, length(union(x,y)), cex=1.5)

        # overlay the anoying null
        points(2.3,-2.1, col="white", cex=3, bg="white", pch=21)

    }
    if(n == 3)
    {
        allItems <- union(x, y)
        allItems <- union(allItems, z)

        allItemsMat <- matrix(0, ncol=n, nrow=length(allItems), dimnames=list(allItems, names)  )

        allItemsMat[ , 1] <- as.numeric( allItems %in%  x )
        allItemsMat[ , 2] <- as.numeric( allItems %in%  y )
        allItemsMat[ , 3] <- as.numeric( allItems %in%  z )

        vennDiagram(vennCounts( allItemsMat  ), ...)

        text(0, 3, paste( "total:",length(allItems) ) , cex=1.5)
        points(2.5,-3, col="white", cex=4, bg="white", pch=21)


    }



}


###################################################
#  given a peptide IDs the function determines all
#  associated evidences, their PEP values and returns
#  the mass error of the evidence having the lowest
#  PEP values
#
#
# pepID   - character, peptide ID
# ev      - evidence.txt
# pep     - peptides.txt
###################################################
getMassAcc <- function(pepID, peptides, evidence){

    # get the evidence IDs of the current peptide ID
    evIDs <- unlist(strsplit( as.character(peptides[pepID, "Evidence.IDs"]) , ";"))

    # get the mass accuracy of the evidence having the lowest PEP
    ma <- evidence[ evIDs[which.min( evidence[evIDs, "PEP"] ) ] , "Mass.Error..ppm." ]

    return(ma)

}
####################################################
#   given a peptide ID the function calculates
#   the average retention time PER EXPERIMENT if
#   specified. Otherwise the RT is averaged over
#   all evidences
#
# arguments:
#  pepID     - character, peptide id
#  peptides  - dataframe, character
#  evidence  - dataframe, character
#
# value:
#  list containing average retention time
#
# changelog: 20091015  - implementation
####################################################
getAverageRT <- function(pepID, peptides, evidence){

    # get all evidence ids for the given peptide id
    evIDs <- unlist(strsplit( as.character(peptides[pepID, "Evidence.IDs"]) , ";"))


    # check if there was an experimental design file specified
    if( length(grep("^Experiment$", colnames(evidence))) > 0  ){

        # get the experiments
        ex <- unique(evidence[, "Experiment"])

        rt.list <- vector("list", length(ex))
        names(rt.list) <- ex

        # remove all other evidences
        evidence <- evidence[evIDs, ]

        # now calculate the average RTs for each experiment
        for(e in ex){
            rt.list[[e]] <- mean( evidence[grep(paste("^",e,"$", sep=""), evidence[, "Experiment"]  ), "Retention.Time"]  )
        }
    } else {
        # remove all other evidences
        evidence <- evidence[evIDs, ]

        rt.list <- vector("list", 1)
        names(rt.list) = "mean.rt"

        rt.list[[1]] <- mean(evidence[, "Retention.Time"])

    }
    return(rt.list)
}



################################################
#
#                 termPieChart
#
#
# term          - vector containing the terms, e.g.
#               GO, KEGG, ...
# first.only    - logical, if there are more terms seperated by 'sep'
#               only the first will be used
# rm.empty      - logical, discard all empty elements of 'terms'
# sep           - character, symbol used as seperator if there
#                 are more terms per element
#
#
# changelog:
# 20090703     - implementation
#
################################################

termPieChart <- function(terms, first.only=T, rm.empty=F, sep=";")
{
    # convert to character
    terms <- as.character(terms)

    # take only the first term
    if(first.only)
        terms <- sub(paste(sep, ".*",sep=""), "", terms)

    # remove empty elements
    if(rm.empty)
        terms <- terms[nchar(terms) > 0  ]

    # count the occurences
    termTab <- table(terms)


    ########################
    #       pie chart
    ########################
    pie(termTab, clockwise=TRUE)
    legend("topleft", legend=names(termTab))


    return(terms)
}
###################################################################
#                  Bland-Altman-Plot
#
#
###################################################################
ba.plot <- function( x, y, ... ){

    diff.xy <- x - y
    mean.xy <- (x+y)/2
    #mean.xy <- y

    fancyPlot(mean.xy, diff.xy, ...)


}

#########################################################################
#
#                         fancy 3D pie charts
#
# changelog:  20101110 implementation
##########################################################################
fancyPie3D <- function(x, labels=x, col=rainbow(length(x)), explode=.05, ... ){

   require(plotrix)

   par(mfrow=c(1,2))
   pie3D(x,  labels=labels, col=col, explode=explode, ...  )
   plot.window(xlim=c(-1,1), ylim=c(-1,1))
   plot.new()
   legend("top", legend=paste( names(x), " (",  x,")", sep=""  ), bty="n", ncol=1, fill=rainbow(length(x)) )


}

##########################################################################
#
#                                 fancyPlot
#
#  produce fancy scatterplot
#
#
# x
# y
# rug
# grid
# loess
# cor
#
#
# changelog:
# 20090703     - implementation
# 20100216     - checking missing values also for NAN and infinity
# 20100311     - added transparancy
# 20100325     - number of datapoints is now reported
# 20100510     - added parameter 'groups'
# 20101011     - changed parameter 'loess' to 'reg' and
#                added linear regression
#              - parameter 'cex.text': controls 'cex' for lagend texts like
#                correlation or regression
# 20110615     - changed color vector: in previous R version one had to
#        ´       supply a color vector starting with an arbitrary color followed
#                by the actual colors (function: scatterpot)
#              - now the color vector starts with the actual color
# 20121004     - included 'reset.par' ain the argument list
#              - before it was set to FALSE as default
############################################################################
fancyPlot <- function(x, y, rug=F, grid=T, reg=c("none", "linear", "loess"), cor=T, cor.method = c( "pearson", "spearman", "kendall"), pch=19, col="black", alpha=60, boxplots=c("xy","x", "y", "n"), groups=1, cex.text=1, reset.par=F, ...)
{
    library(car)

    # regression
    reg = match.arg(reg)

    ############################
    # color vector?
    ############################
    col.asis = F
    if(length(col) == length(x))
        col.asis = T

    ############################
    # color
    ############################
    if(!col.asis){

        color <- list()

        for(i in 1:length(col)){
          color[[i]] <- col2rgb(col[i])
          color[i] <- rgb(color[[i]][1], color[[i]][2], color[[i]][3], alpha, maxColorValue=255)
        }
        #col=c("black", unlist(color))
        col=c(unlist(color)) # 20110615
    }

    # plotting characters: the same number as groups
    if(length(pch) == 1 & length(unique(groups)) > 1) pch <- rep(pch,length(unique(groups)))

    #############################
    # correlation coefficient
    #############################
    cor.method=match.arg(cor.method)

    x <- as.numeric(x)
    y <- as.numeric(y)
    ##############################
    # remove missing values
    ##############################
    #na.idx <- union(union( which(is.na(x)), which(is.na(y)) ), union( which(is.nan(x)), which(is.nan(y)) ))
    na.idx <- union( which(is.na(x)), which(is.na(y)) )
    if(length(na.idx) > 0){
        x <- x[-na.idx]
        y <- y[-na.idx]

        if(length(groups)>1) groups <- groups[-na.idx]
        if(col.asis) col <- col[-na.idx]
    }
    inf.idx <- union( which(is.infinite(x)), which(is.infinite(y))  )
    if(length(inf.idx)>0){
        x <- x[-inf.idx]
        y <- y[-inf.idx]
        if(length(groups) > 1) groups <- groups[-inf.idx]
        if(col.asis) col <- col[-inf.idx]
    }

    #if(col.asis) col <-  c("black", col)
    if(col.asis) col <-  c("black", col) # 20110615

    ###############################
    # boxplot on the axes
    ###############################
    boxplots=match.arg(boxplots)


    palette.org <- palette()
    ########################################
    #            plot
    #########################################
    scatterplot(x,y, col=col, pch=pch, boxplots=boxplots, reg.line=F, smooth=F, sub=paste("N=", length(x),sep=""), groups=groups, reset.par=reset.par ,legend.plot=ifelse(length(groups) == 1, F, T), ...)

    #
    if(grid) grid()
    if(rug) rug(x, col="indianred")

    # regression
    if(reg == "loess"){
        LOESS <- loess( y~x  )
        lines(x, LOESS$fitted, col="red", type="p", pch=20)
        legend("topright", legend="loess", fill="red", bg="white", cex=cex.text)
    } else if(reg == "linear"){

        fit = lm( formula=y~x)

        #abline( fit$coefficients[1], fit$coefficients[2], col="red", lwd=1.5 )
        abline(fit, col="red", lwd=1.5)

        legend("bottomright", legend=paste("y = ", round(fit$coefficient[2], 3), "x ", ifelse(fit$coefficient[1] < 0, "- ", "+ "), round(abs(fit$coefficient[1]), 3), sep="" ) , bty="n", text.col="red", cex=cex.text)
    }


    # correlation coefficient
    if(cor){
        minY <- min(y, na.rm=T)
    legend("top", legend=paste( toupper(substr(cor.method, 1,1)), substring(cor.method, 2) ,"-correlation = ",  round(cor(x,y, method=cor.method), 2), sep=""), bty="n", cex=cex.text)

    }
}


################################################
#              getMW
#
# sequence    - character vector, each element
#               represents an aa
#             - string
#
################################################
getMW <- function(sequence){


    ###############################
    # aa annotation
    ###############################
    AA <- data.frame(
     letter1=c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "B", "Z", "J", "X"),
     letter3=c("Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His","Ile","Leu","Lys","Met","Phe","Pro","Ser","Thr","Trp","Tyr","Val","Asx","Glx","Xle","Xaa"),
     MW=c(71.03712,156.10112,114.04293,115.02695,103.00919,128.05858,129.04260,57.02147,137.05891,113.08407,113.08407,128.09497,131.04049,147.06842,97.05277,87.03203,101.04768,186.07932,163.06333,99.06842,114.53490,128.55050,113.08400,118.80568),

     row.names=1)

    ################################
    #  monoiso. mass water
    ################################
    mwH2O <- 18.01056

    ################################
    #  check if vector or string
    ###############################
    if(length(sequence)==1)
        sequence = unlist(strsplit(sequence, ""))

    ################################
    #  calculate  MW
    ################################
    mw <- 0
    for(aa in 1:length(sequence))
   {     mw <- mw + AA[ toupper(sequence[aa]), "MW"]
         #cat(AA[ toupper(sequence[aa]), "MW"],"\n")
    }
    # substract mw of peptide bonds
    mw = mw + mwH2O

    return(mw)

}

#######################################################################
#
# plots the (log) SILAC protein ratios vs intensities, marks significant
# ratios (sig A or B) and adds protein ids of significant hits
#
# pg         - character, protein groups table
# norm       - logical, if TRUE the normalized ratios are usd
# rmCon      - logical, if TRUE contaminants are filtered out
#            - by default they have no SILAC label thus they
#              are heavily downregulated, if Re-quantify was check in
#              the identify module
# silac.state- character, which which SILAC label
# p          - numeric, significance threshold
# experiment - character, name of experiment
#
# changelog: - 20090905 implementation
#              20091202 rmCon, silac.state added
#                       reverse hits are now filterd by default
#              20091203 some documentation
#
#######################################################################

sigRatioPlot <- function(pg, norm=T, rmCon=F, silac.state=c("H.L", "M.L", "H.M"), p=0.05, adjust.p=c("none", "BH"), experiment="", signif="B", sig.col="red", label.col="blue", add.label=T, label.top=5, main="Ratio vs. Intensity", pch=20, col="black", alpha=50, ...){

    # remove contaminants
    if(rmCon){
        con.idx <- grep("^\\+$", pg[, "Contaminant"])
        if(length(con.idx)>0) pg <- pg[-con.idx, ]
    }
    # remove revererse hits
    rev.idx <- grep("^\\+$", pg[, "Reverse"])
    if(length(rev.idx)>0) pg <- pg[-rev.idx, ]

    # get the silac label
    silac.state <- match.arg(silac.state)

    # if there is an experiment defined...
    if(nchar(experiment) > 0){
        ratioColumn <- ifelse(norm, paste("Ratio.", silac.state,".Normalized.", experiment, sep=""),paste("Ratio.",silac.state,".", experiment, sep=""))
        intColumn <- paste("Intensity.", experiment, sep="")
        sigColumn <-  paste("Ratio.", silac.state, ".Significance.",signif,"..", experiment, sep="")

    } else {
        ratioColumn <- ifelse(norm, paste("Ratio.", silac.state,".Normalized", sep=""),paste("Ratio.",silac.state, sep=""))
        intColumn <- paste("Intensity", sep="")
        sigColumn <-  paste("Ratio.",silac.state,".Significance.",signif,".", sep="")
    }

    ########################################
    # get ratios and intensities
    ########################################
    ratio <- log( pg[, ratioColumn], 2)   # always log 2

    intens <- log( pg[, intColumn ], 10)  # alway log 10
    intens[intens == "-Inf"] <- NA

    # significance values
    sig <- pg[, sigColumn]

    # adjust significance values for multiple hypthesis testing
    adjust.p <- match.arg(adjust.p)
    sig <- p.adjust(sig, adjust.p)


    # median ratio
    median.ratio <- median(ratio, na.rm=T)

    #############
    # regulation
    #############
    ratioRegulated <- ratio[sig < p]
    intensRegulated <- intens[sig < p]

    ########################################
    # plot
    ########################################

    ###################################
    # color
    ###################################
    col.org <- col
    sig.col.org <- sig.col
    col <- col2rgb(col)
    col <- rgb(col[1], col[2], col[3], alpha, maxColorValue=255)
    sig.col <- col2rgb(sig.col)
    sig.col <- rgb(sig.col[1], sig.col[2], sig.col[3], alpha, maxColorValue=255)

    # set the color palette: the first entry is just a dummy
    palette( c("black", col, sig.col, label.col))


    # color
    #col <- col2rgb(col)
    #col <- rgb(col[1], col[2], col[3], alpha, maxColorValue=255)

    plot( ratio[sig >= p], intens[sig >= p], ylim=range(intens, na.rm=T), xlim=c( -max(abs(ratio), na.rm=T), max(abs(ratio), na.rm=T) ), main=main, pch=pch,
         xlab=expression(log[2](Ratio)), ylab=expression(log[10](Intensity)), col=col, ... )

    points( ratioRegulated, intensRegulated, col=sig.col, pch=20 )


    if(add.label){
        # protein accessions, leading protein
        protAcc <- sub(";.*", "", pg[, "Protein.IDs"])

        # row index of the top n proteins
        topIdx <- order(sig)[1:label.top]

        top.x <- ratio[topIdx]
        top.y <- intens[topIdx]
        top.id <- protAcc[topIdx]

        text(top.x, top.y + max(intens, na.rm=T)*.05, top.id, cex=.55)

    }
    legend("topright", legend=c(paste("singificance", signif, "<", p)), fill="red")

    ###########################################
    #              output
    # - number of pg having a ratio
    # - number of significantly regulated pg
    ###########################################
    output <- vector("list", 3)
    names(output) <- c("number.ratios","number.down", "number.up" )
    output[["number.ratios"]] <- sum( !is.na(ratio) )
    output[["number.down"]] <- sum(ratioRegulated < median.ratio, na.rm=T )  #sum(ratio[sig < p] < 1, na.rm=T)
    output[["number.up"]] <- sum(ratioRegulated >= median.ratio, na.rm=T )    #sum(ratio[sig < p] > 1, na.rm=T)


    return(output)
}


#####################################################################
#
#   GOlist   - list returned by 'topGOanalysis'
#   pval     - 0.01
#
# changelog:  20090911  - implementation
#             20100217  - only the first three entries of the list are
#                         taken containing the results of the three
#                         categories
#             20100218  - a short summary is now included in the plot
#####################################################################
GOresultPlot <- function(GOlist,  pval=.01, adj.p = T, category=c("all", "BP", "MF", "CC"),  label="", mfrow=NULL, titles=NULL, ...){

       # take only the first three entries
       annot <- GOlist[[4]]
       GOlist <- GOlist[1:3]


       category <- match.arg(category)

       if(category != "all"){

           # get the GO category to plot
           cat.name <- names(GOlist)[grep(paste("^", category,sep=""), names(GOlist))]
           GOlist <- GOlist[[cat.name]]

           # title for the plot
           if(is.null(titles)){
               titles=category
           } else {
               if(length( titles) > 1 ){
                   titles <- titles[1]
                   warning("\n\nMore titles than GO categories provided! I took the first title to annotate the plot!\n\n")
               }

           }
           GOplot(GOlist, pval, titles)

       } else {

            if(is.null(mfrow))
                par(mfrow=c(2,ceiling(length(GOlist)/2  )))
            else par(mfrow=mfrow)

            # plot titles
            if(is.null(titles))
                titles <- names(GOlist)


           for(i in 1:length(GOlist))
           {
               # get the table
               tab <- GOlist[[i]]

               GOplot(tab, pval, titles[i], adj.p)

           }
            ########################
            # add some information
            ########################
            par(mar=rep(0,4))
            plot.window(xlim=c(0,1), ylim=c(0,1))
            plot.new()

            # title
            text(x=.5, y=.9, "Summary", cex=1.5)

            # label
            if(nchar(label)>0) text(x=.1, y=.8, paste("label: ", label), cex=1, pos=4)

            # number of proteins in the test and reference set
            text(x=0.1, y=.7, paste("test set: ", annot[1], sep=""), pos=4 )
            text(.1, 0.65, paste("annotated: ", annot[2], sep=""), pos=4)
            text(.1, 0.6, paste("reference set: ", annot[3], sep=""), pos=4)

            # p-values
            text(.1, 0.5, paste("p-values: ", ifelse(adj.p, paste(colnames(tab)[grep("^adjusted p-value", colnames(tab))]), "unadjusted"), sep=""), pos=4)
            text(.1, 0.45, paste("p-value cutoff: ", pval, sep=""), pos=4)

            par(mfrow=c(1,1))
   }

}

###########################################
#  helper function for 'GOresultPlot'
#
# this function does the plotting
#
#
# changelog:  20100407 p-values are now transformed
#                      to +/- log10(p) instead of
#                         +/- 10*log2(p)
###########################################
GOplot <- function(tab, pval, titles, adj.p=T ){


           # determine the column containing the adjusted p values
           if(adj.p){
               pvalCol <- grep("^adjusted p-value", colnames(tab))
           } else pvalCol <- grep("^p-value$", colnames(tab))

           # determine the column containing the direction of significance, i.e. 'over/under'
           directionCol <- grep("^over/under$", colnames(tab))


           # get the terms with p-val < 'pval'
           tab <- tab[ which(as.numeric(tab[, pvalCol]) < pval),  ]

           # if there are any significant p values ...
           if(dim(tab)[1] > 0 ){
               width <- ifelse(dim(tab)[1] == 1, 0.1, 1)

              # log p-values
              pval.log <-  ifelse( tab[, directionCol] == "over", -log(as.numeric(tab[, pvalCol]),10), log(as.numeric(tab[,pvalCol]),10)  )
              names(pval.log) <- tab[, 2]

              # reverse vector such that the most significant terms are on top
              pval.log <- pval.log[length(pval.log):1]

              if(length(pval.log) == 1) pval.log <- c(0, pval.log)

              # colors
              cols <- ifelse(pval.log < 0, "green", "red")

              # plot
              par(mar = c(5,15,4,2), cex.axis=.7)
              barplot(pval.log, xlim=c(-max( abs(pval.log) ), max(abs(pval.log))  ),

                   xlab=expression(paste( "+/-", log[10](p)  )),
                   main=titles,
                   width=width,
                   horiz=T, las=2, cex.lab=.9, cex.sub=.7, space=.5, col=cols)
              grid(ny=0)
              legend("bottomright", legend=c("Enrichment"), bty="n", text.col="red")
              legend("bottomleft", legend=c("Depletion"), bty="n", text.col="green")
              abline(v=0, lwd=1.5)
           } else {

               plot(0,0, col="white", main=titles, axes=F, xlab="",ylab="")
           }
}

#########################################################################
#             hypergeometrical test
#
#  test    - vector of protein ids
#  goslim  - go slim terms downloaded from AmiGO using GOslimmer
#          - the file has to be imported with 'readLines(file)'
#  pcut    - numeric, p-value for cutoff
#  adjust  - p-values correction
#
# changelog:  20100316 implementation
#             20100407 protein ids are now made unique before testing
##########################################################################
GOslim.analyis <- function( test, goslim, p.cut=0.01, adjust="BH", alternative=c("two.sided", "greater", "less"), plot=T, ...  ){

    alternative <- match.arg(alternative)

    # count each protein only once
    test <- unique(test)

    goslim.split <- strsplit(goslim, "\t")

    # get the terms and names
    go.slim.terms <- unlist(lapply(goslim.split, function(x)x[1]  ))
    go.slim.names <- unlist(lapply(goslim.split, function(x)x[2]  ))
    names(goslim.split) <-  names(go.slim.names) <- go.slim.terms


    # all annotated protein ids
    prot.anno <-  sub("^GeneDB_Spombe:", "", unique(unlist(lapply(goslim.split, function(x) unlist(strsplit(x[4]," "))))))


    # vector storing the p-values
    p <- dir <- vector("numeric", length(go.slim.terms))
    names(p) <- names(dir) <- go.slim.terms


    ################################################
    # Fisher's exact test
    #
    ################################################
    for(go in go.slim.terms){

        # all proteins that are annotated with current GO term
        prot.anno.GO <- sub("^GeneDB_Spombe:", "", unlist(strsplit(goslim.split[[go]][4], " ")))


        # parameters for hypergeometric distribution
        x <- sum( test %in% prot.anno.GO)
        k <- sum( test %in% prot.anno )

        N <- length(prot.anno)
        m <- length(prot.anno.GO)

        # now do the test
        p[go] <- hyper.test(x,m,N-m,k, alternative=alternative)
        dir[go] <- ifelse(x/k >= m/N  , "enrichment", "depletion"  )

    }


    ##################################
    # adjust p-values
    ##################################
    p.BH <- p.adjust(p, adjust)

    ##################################
    # get significantly enriched terms
    ##################################
    goslim.enriched <- p.BH[ p.BH <= p.cut ]
    goslim.enriched.names <- go.slim.names[names(goslim.enriched)]
    goslim.enriched.direction <- dir[names(goslim.enriched)]

    if(plot & (length(goslim.enriched) > 0)){
      if(alternative == "two.sided"){

          p.plot <- unlist( lapply( names(goslim.enriched), function(x) ifelse( goslim.enriched.direction[x] == "enrichment", -log(goslim.enriched[x], 10), log(goslim.enriched[x], 10)  )   ) )

      }
      if(alternative == "greater")
          p.plot <- -log(goslim.enriched,10)

      if(alternative == "less")
          p.plot <- -log(goslim.enriched,10)

      # order according to p-values
      goslim.enriched <- goslim.enriched[ order(goslim.enriched, decreasing=T)  ]
      goslim.enriched.names <- goslim.enriched.names[ names(goslim.enriched) ]
      goslim.enriched.direction <- goslim.enriched.direction[ names(goslim.enriched) ]
      p.plot <- p.plot[names(goslim.enriched) ]

      col <- ifelse( p.plot < 0, "green", "red"  )

      ####################################
      #   plot
      ####################################
      par(mar=c(5,25,5,5))
      barplot( p.plot, names.arg=goslim.enriched.names, las=1, xlab=expression(+- log[10](p)), horiz=T, xlim=c(-max(abs(p.plot)),max(abs(p.plot)) ), col=col
       ,   space=0,...)
      text( ifelse(p.plot < 0, -0.05, 0.05), (1:length(p.plot))-.5, paste("p=", formatC(goslim.enriched)), pos=ifelse(p.plot < 0, 2, 4), cex=0.8  )
    }

    ####################################
    # output
    ####################################
    out <- vector("list", 2)
    names(out) <- c("p", "direction")

    out[[1]] <- p.BH
    out[[2]] <- dir

    return(out)
}

###########################################################################################
#            given a list of GO term ids or GO term names the function
#            returns a higher level term at the specified level
#
# terms      - character, vector of GO ids or names
# go.level   - numeric, GO level in graph
#
# changelog:  20101119 implementation
#             20101201 parameter 'rm.higher.levels'
############################################################################################
getGOterm <- function(terms, go.level=3, rm.higher.terms=T ){

    require(GO.db)

    # all GO terms
    GOTERMS <- as.list(GOTERM)

    #######################################
    # check whether ids or terms
    #######################################
    if( length(grep("GO:", terms[1])) == 0   ){
        # get all information about the GO terms
        term.GOTERM <- unlist(lapply( terms, function(x) GOTERMS[which(x == lapply(GOTERMS, Term))]))
    } else {
        term.GOTERM <- GOTERMS[terms]
    }
    # get the Ontologies
    ontology <- names( table( unlist(lapply(term.GOTERM, Ontology)) ))

    # create objects for each Ontology
    for(o in 1:length(ontology))
        assign( paste( ontology[o], "PARENTS", sep=""),  as.list(eval(parse(text=paste("GO", ontology[o], "PARENTS", sep="") ))) )


    ####################################################
    # create some objects to store the results
    ####################################################
    gomap.mat <- matrix("", nrow=length(terms), ncol=2, dimnames=list(terms, c(paste("Term level", go.level), "Ontology")) )

    goparent.list <- vector("list", length(terms))
    names(goparent.list) <- terms

    ###########################################################################
    # now loop over these terms and get all parent terms on the way to the root
    for(term.id in names(term.GOTERM)){

        # ontology of that term
        term.ontology <- Ontology(term.GOTERM[[term.id]])
        parent.list <- list()

        count <- 1

        # go to the root
        node = term.id
        while(  !("all" %in% node)){

             # get the parent node
             node = eval(parse(text= paste(term.ontology, "PARENTS[[node]]", sep="")) )

             # consider only 'is a' relationships
             node <- node[grep("is_a", names(node))][1]

             # store the node
             parent.list[[count]] <- node
             count <- count + 1
       }

       # node levels as names
       names(parent.list) <- paste( "level",(length(parent.list)-1):0, sep="_"   )


       #############################
       # store the results
       ##############################
       term.name <-   Term(term.GOTERM[[term.id]])

       goparent.list[[term.name]] <- parent.list

       #############################################
       # if the term is on a higher level than specified
       # either no term or the original terms is returned
       # depending on paramter 'rm.higher.levels'
       if(rm.higher.terms){

           if(paste("level", go.level, sep="_" ) %in% names(parent.list)){

               gomap.mat[term.name, paste("Term level", go.level)] <- Term(GOTERMS[[ parent.list[[paste("level", go.level, sep="_")]]]] )

           } else {
               # check whether the term is already on the specified level
               if( as.numeric(sub("level_", "", names(parent.list))[1]) == (go.level-1)  ){
                   gomap.mat[term.name, paste("Term level", go.level)] <- term.name
               } else {
                   warning(  paste(names(parent.list)[1],"term", term.name, "has been removed!\n")  )
               }
           }
       } else {

           gomap.mat[term.name, paste("Term level", go.level)] <- ifelse( paste("level", go.level, sep="_" ) %in% names(parent.list),  Term(GOTERMS[[ parent.list[[paste("level", go.level, sep="_")]]]] ), term.name)
       }


       gomap.mat[term.name, "Ontology"] <- term.ontology

    }

    ###############################
    # output
    ###############################
    out <- list()
    out[[1]] <- gomap.mat
    out[[2]] <- goparent.list
    names(out) <- c("term.matrix", "parent.list")

    return(out)
}

############################################################################################
#               given a table of GO terms as returned by function
#               'annotation.enrich' the function maps the GO terms
#               to more general ones ('go.level'), if possible.
#
# go.table         - data frame, returned by function 'annotation.enrich'
#                    rownames have to be GO terms or GO ids
# go.level         - numeric, level in GO graph
# rm.higher.terms  - logical,if TRUE all GO terms found in 'go.table' that are on a higher
#                    level than 'go.level' are removed
#                  -> all terms in the pie are on the same level
# ...              - further arguments passed to 'fancyPie3D'
#
# changelog:  20101122 implementation
#             20101201 parameter 'rm.higher.terms'
############################################################################################
makeGOpie <- function(go.table, go.level=3, rm.higher.terms=T, ...){

    # rownames of the table should be ther terms
    terms <-  rownames(go.table)

    # get the GO terms at specified level
    terms.level <- getGOterm( terms, go.level=go.level, rm.higher.terms=rm.higher.terms  )
    terms.level.mat <- terms.level$term.matrix
    if(rm.higher.terms){

        # check if any terms have o be removes
        rm.idx <- which(nchar(terms.level.mat[, paste("Term level", go.level)]) == 0)
        if(length(rm.idx) > 0){
            terms.level.mat <- terms.level.mat[-rm.idx, ]
            terms <- terms[-rm.idx]
            go.table <- go.table[-rm.idx, ]
        }
    }

    # append to the original table
    go.table <- cbind(go.table, terms.level.mat )

    # count proteins associated with each non-redundant term at level 'go.level'
    terms.level.nr <- unique(as.character(go.table[, paste('Term level', go.level)]))
    terms.level.nr.prot <- unlist(lapply( terms.level.nr, function(x) length(unique(unlist(strsplit(as.character( go.table[ which(go.table[, paste('Term level', go.level)] ==x)  , "proteins"]), ";" ))))  ))
    names(terms.level.nr.prot) <- terms.level.nr

    # pie chart
    fancyPie3D( terms.level.nr.prot, ...  )

    return(go.table)
    return(terms.level.nr.prot)
}

###########################################################################################
#
#
#
#
###########################################################################################
insilico.digest <- function(db,  enzyme=c("trypsin", "lysine") ){

    enzyme = match.arg(enzyme)

    # upper case
    #seq <- toupper(seq)

    # remove the asterisk at the end
    #seq <- sub("\\*$", "", seq)

    if(enzyme == "trypsin")
        cleav = "K|R"

    if(enzyme == "lysine")
        cleav = "K"

    # number of observable peptides without missed cleavages
    peptides <- unlist(strsplit(seq.dig, " "))


    # remove 'zero length' peptides
    pepL <- unlist( lapply(peptides, nchar)  )
    if( sum( pepL == 0  ) > 0  )
        peptides <- peptides[ -which(pepL  == 0)  ]


    return(peptides)
}

############################################################################################
#
#                                          PAI
#
# calculates for a given protein group the PAI value respective to its leading protein
#
# arguments:
#   pg.row     - vector representing a single row of the 'proteinGroups.txt' table
#              - vector names are column names of that table
#   mc         - numeric, number of missed cleavages
#   enzyme     - character
#   minAA      - numeric, minimal peptide length
#   db         - fasta db, read in by function 'read.fasta'
#
# changelog: 20090925   implementation
############################################################################################
PAI.obsolete <- function(pg.row, mc=2, enzyme=c("trypsin", "lysc"), minAA=6, db, emPAI=T  ){

    enzyme = match.arg(enzyme)

    # get the protein sequence from the database
    protID <- sub(";.*", "", pg.row["Protein.IDs"])
    sequenceProtein <- paste(db[[protID]], collapse="")
    sequenceProtein <- sub("\\*$", "", sequenceProtein)

    if(nchar(sequenceProtein) == 0) stop("Protein not found in database!!\n")

    # number of observed peptides (leading protein, i.e. first protein in the groups)
    nPepObserved <- as.numeric(sub(";.*", "", pg.row["Peptide.Counts"]))

    # number of observable peptides
    nPepObservable <- observablePeptides( sequenceProtein, enzyme, mc, minAA  )

    PAI = nPepObserved/nPepObservable

    # emPAI: PAI values show a linear relationship with the logarithm
    # of protein concentration. emPAI values are exponetially modified
    # PAI values that are proportional to protein content in a mixture
    if(emPAI) PAI <- 10^PAI - 1


    return(PAI)
}


######################################################
#             observablePeptides
#
# given a protein sequence, an enzyme specificity
# and the number of missed cleavages, this function
# calculates the number of theorethically
# observable peptides.
#
# arguments:
#    seq            - character, the protein sequence
#    enzyme         - character
#    mc             - numeric, number of allowd missed
#                     cleavages
#
#
#######################################################
observablePeptides.obsolete <- function(seq, enzyme=c("trypsin", "lysine"), mc=2, minAA=6){

    enzyme = match.arg(enzyme)
    if(enzyme == "trypsin"){
        cut <- "R|K"
    }

    # upper case
    seq <- toupper(seq)

    # remove the asterisk at the end
    seq <- sub("\\*$", "", seq)

    # number of observable peptides without missed cleavages
    peptides <- unlist(strsplit(seq, cut))

    nPepZeroMC <- length( peptides[ nchar(peptides) >= (minAA-1) ]   ) # minus one, because the AA after which the sequence
                                                                       # is cutted is removed by function 'strsplit'

    nPep <- nPepZeroMC
    # missed cleaved peptides
    if(mc > 0){
      for(missedC in 1:mc){
          nPep <- nPep + max( (nPepZeroMC-missedC), 0)
      }
    }

    return(nPep)

}
############################################################################################
#
#                                          PAI
#
# calculates for a given protein group the PAI value respective to its leading protein
#
# arguments:
#   pgID       - protein group id
#   evidence           - evidence table, there are the observed precursor masses
#   enzyme     - character
#   mass.range - numeric vector of length 2 specifying the precursor mass range in Th
#   db         - fasta db, read in by function 'read.fasta'
#
# changelog: 20091110   implementation
#            20111010   added parameter 'experiment'
############################################################################################
PAI <- function(pgID, evidence, experiment=NULL, enzyme=c("trypsin", "lysc"), mass.range=c( 600, 6000  ), db, emPAI=T  ){

    enzyme = match.arg(enzyme)


    if(!is.null(experiment)){

        # check whether the experiment is defined
        if( !(experiment %in% evidence$Experiment)){
            stop(paste("cannot find experiment '", experiment, "' in the evidence table...\n", sep=""))
        } else {
            exp.idx <- grep(experiment, evidence$Experiment)
        }
    } else {
        exp.idx <- 1:dim(evidence)[1]
    }

    ######################################################
    # get all evidences associated with protein group id
    ######################################################
    evidence <- evidence[ grep( paste("(^|;)",pgID,"($|;)", sep="" ), evidence[ exp.idx, "Protein.Group.IDs"]  ),   ]
    #evidence <- evidence[ grep( paste("^",pgID,"$", sep="" ), evidence[, "Proteins"]  ),   ]

    ######################################################
    # leading protein
    ######################################################
    protID <- sub(";.*", "",evidence[1, "Proteins"])

    # get the protein sequence from the database
    sequenceProtein <- paste(db[[protID]], collapse="")
    sequenceProtein <- sub("\\*$", "", sequenceProtein)
    #if(nchar(sequenceProtein) == 0) stop("Protein not found in database!!\n")
    if(nchar(sequenceProtein) == 0) return(NA )



    # number of observed peptides (leading protein, i.e. first protein in the groups)
    nPepObserved <- length( unique( evidence[, "m.z"]  )  )

    # number of observable peptides
    nPepObservable <- observablePeptides( sequenceProtein, enzyme,  mass.range  )

    PAI = nPepObserved/nPepObservable

    # emPAI: PAI values show a linear relationship with the logarithm
    # of protein concentration. emPAI values are exponetially modified
    # PAI values that are proportional to protein content in a mixture
    if(emPAI) PAI <- 10^PAI - 1


    return(PAI)
}


######################################################
#             observablePeptides
#
# given a protein sequence, an enzyme specificity
# and the number of missed cleavages, this function
# calculates the number of theorethically
# observable peptides.
#
# arguments:
#    seq            - character, the protein sequence
#    enzyme         - character
#
# changelog:  20111007 added 'lysc'
#######################################################
observablePeptides <- function(seq, enzyme=c("trypsin", "lysc"), mass.range=c(600, 6000), charge.range=c(2,4)  ){

    enzyme = match.arg(enzyme)

    # upper case
    seq <- toupper(seq)

    # remove the asterisk at the end
    seq <- sub("\\*$", "", seq)

    if(enzyme == "trypsin"){

        seq.dig <- gsub( "K", "K ", seq   )
        seq.dig <- gsub( "R", "R ", seq.dig)
    }
    if(enzyme == "lysc"){

        seq.dig <- gsub( "K", "K ", seq   )
    }

    # number of observable peptides without missed cleavages
    peptides <- unlist(strsplit(seq.dig, " "))

    # remove 'zero length' peptides
    pepL <- unlist( lapply(peptides, nchar)  )
    if( sum( pepL == 0  ) > 0  )
        peptides <- peptides[ -which(pepL  == 0)  ]


    # calculate the theoretical masses for each peptide
    pepMW <- unlist( lapply(peptides, getMW)  )


    # remove all peptides falling outside the mass range
    pepMWrange <- which( pepMW >= mass.range[1] & pepMW <= mass.range[2]  )


    return( length(pepMWrange) )


    #return(pepMW[pepMWrange])

}
###################################################################
#                          split.fasta
#
#
# changelog:  20101126 implementation
###################################################################
split.fasta <- function(fasta, nseq=5000, seqtype=c("AA", "DNA")){

    require(seqinr)

    seqtype=match.arg(seqtype)

    # import database
    db <- read.fasta(fasta, seqtype=seqtype, as.string=T)
    db.names <- names(db)
    db.anno <- unlist(lapply(db, function(x) attributes(x)["Annot"]))

    # remove duplicates
    if(length(unique(db.names)) < length(db.names)){
        dups.idx <- which(duplicated(db.names))

        db <- db[-dups.idx]
        db.names <- db.names[-dups.idx]
        db.anno <- db.anno[-dups.idx]
    }

    int <- c(seq(0, length(db),by=nseq), length(db))

    #return(int)

    count=1;
    if(length(int) > 1){
        for(i in 2:length(int)){

            write.fasta(db[(int[i-1]+1): int[i]], file.out=paste( fasta, "_part_",count,".fasta", sep=""  ), names=db.names[(int[i-1]+1):int[i] ] )

            count=count+1

        }
    }

    return(db)
}



###################################################################
#
#                    spectralCount
#
# counts the number of spectra per peptide sequence, i.e. regardless
# of any modifications
#
# arguments:
#   evidence          - 'evidence.txt' table
#
#
# changelog:    20090928 implementation
####################################################################
spectralCount <- function( evidence  ){

    # get all protein group ids
    pgID <- unlist( strsplit( evidence[, "Protein.Group.IDs"], ";"  )  )

    pgIDunique <- unique(pgID)

    nSpec <- unlist( lapply( pgIDunique,  function(x){
        length(grep(paste("^", x, "$", sep=""), pgID )   )

         } ) )

    names(nSpec) <- pgIDunique

    return(nSpec)
}


#############################################################################
#  import BLAST result files from our BLAST server
#  currently only txt-format is supported
#
#  file  - file of Blast output in hit table format
#  fix   - logical, if T a new column will be added indicating
#          the reading direction and the start/stop positions will be fixed
#
# changelog:  20091216  implementation
#             20100202  parameter 'fix'
#############################################################################
importBLAST <- function( file, format="txt", fix=T ){

    format = match.arg( format)

    if(format == "txt"){

        # get the header, i.e. read in the first ten lines and check
        # if they equal a typical BLAST header
        header <- readLines(file, n=150)
        header <- header[ grep("^# Fields:", header) ][1]
        header <-  sub("^# Fields: ", "", header)


        header <- unlist(strsplit(header, ","))
        header <- sub("^ ", "", header)

        # now import the BLAST table
        table <- read.delim(file, comment.char="#", header=F, stringsAsFactors=F)
        colnames(table) <- header


    }
    if(fix){
        # the positions
        start <- as.numeric(apply(table, 1, function(x) ifelse(x["s. start"]< x["s. end"], x["s. start"], x["s. end"] )  ))
        end <- as.numeric(apply(table, 1, function(x) ifelse(x["s. start"] > x["s. end"], x["s. start"], x["s. end"] )  ))

        # determine the strand
        strand <- apply(table, 1, function(x) ifelse(x["s. start"]< x["s. end"], "fwd", "rev" )  )

        # now add it to the table
        table[, "s. start"] <- start
        table[, "s. end"] <- end

        table <- cbind(table, strand)
        colnames(table)[dim(table)[2]] <- "strand"


    }

    return(table)
}

#####################################################################################################
#
# the function takes as argument the table returned by
# function 'importBLAST' and determines for each query
# the BLAST hit having the highest similarity score.
# for each query a string is produced conatining the
# query id, subject id, similarity score, ...
# this string is then used for the fasta header of the
# collapsed database
#
# input:
#  - table returned by function 'importBLAST'
#
#
# changelog:  20100120 implementation
#             20200210 bugfix, if the sequence with maximal score was
#                      not reported on the first position, it was
#                      not recognized
#             20100408 'no hits found' is removed
#             20100413 'highestBitScore' - if T the match with the highest bit score is returned
#                      problem with simscore: although the similiarity is 100% this does not
#                      mean that the entire sequence was matched but only a few AAs ...
#             20100414 'fasta.query' if specified the query length will be added to the table
#             20100722 'delta.score' is now returned, i.e. the difference between best and second best
#                       hit
#             201012   peptide.match: if TRUE, the alignment length has to be the same as the query length AND
#                                     similarity has to be 100%, i.e. no mismatches are allowed
#                      - if there are several matches like tat,they are all returned, i.e. the resulting
#                        table is not neccessarily non-redundant
#
##########################################################
collapseBLASTresults <- function( table, simscore=100, peptide.match=F, highestBitScore=T, fasta.query=NULL, par=F ){

    # perfectmatch
    if(peptide.match){
        highestBitScore=F

        if(is.null(fasta.query))
            stop("\nYou have to specify a fasta file containing the query sequences in order to find perfect matches!\n")
    }
    # replace spaces by dots in the header
    colnames(table) <- gsub(" ", ".", colnames(table))

    ################################
    #  output table
    ################################
    output.mat <- matrix("", ncol=dim(table)[2], nrow=0 )
    colnames(output.mat) <- colnames( table )

    # the query ids
    query.ids <-  unique(table[, "query.id"])

    ############################################
    # if there is a fasta file specified
    # store the query length in a vector
    ############################################
    if(!is.null(fasta.query)){

        require(seqinr)
        fasta.query <- read.fasta(fasta.query )

        query.length <- vector("numeric", length(query.ids))
        names(query.length) <- query.ids
    }


    # check if there is a 'no hits found'
    no.hit.idx <- grep("No hits found", query.ids)
    if(length(no.hit.idx) > 0) query.ids <- query.ids[ -no.hit.idx ]

    ######################################
    # vector to store delta bit score
    #######################################
    delta.bit <- vector("numeric", length(query.ids))
    names(delta.bit) <- query.ids

    #######################################
    #
    #######################################
    #if(par){
    #      w <- startWorkers(cores)
    #registerDoSMP(w)
    #######################################
    # loop over all query ids
    #######################################
    #tmp.mat <- foreach( query = query.ids, .combine=rbind )  %dopar% {

    #}

    #######################################
    # loop over all query ids
    #######################################
    cat("getting 'best' blast hit\nlooping over ", length(query.ids), " queries...\n")
    count = 1
    for( query in query.ids ){
        cat(".")
        if((count %% 200) == 0) cat(" ", count, "\n")

        # all hits of current query
        table.tmp <- table[ which(table[, "query.id"]  == query) , ]

        ##############################
        # add query length and delta bit
        # vectors (only a dummy for
        # delta bit)
        ##############################
        #if(!is.null(fasta.query)){
        #    query.length.tmp <- ""
        #    if( query %in% names(fasta.query) )
        #        query.length.tmp <- length( unlist( fasta.query[query] ))
        #
        #    table.tmp <- cbind(table.tmp, rep("", dim(table.tmp)[1]), rep( query.length.tmp,  dim(table.tmp)[1] ) )
        #    colnames(table.tmp)[(dim(table.tmp)[2]-1):dim(tample.tmp)[2]  ] <-c( "delta.bit", "query.length")
        #}

        ##############################
        # query length
        ##############################
        if(!is.null(fasta.query)){
            if(query %in% names(fasta.query)){
                query.length[query] <- length( unlist( fasta.query[query] ) )
            }
        }



        #####################################
        # if there are multiple hits...
        #####################################
        if(dim(table.tmp)[1] > 1){

            ##########################################
            # order according to bit score
            ##########################################
            table.tmp.ord <- table.tmp[ order( table.tmp$bit.score, decreasing=T),  ]

            ##########################################
            # get all hits with the required simscore
            ##########################################
            if(!highestBitScore & !peptide.match){

              table.tmp <- table.tmp[  which(table.tmp[, "%.identity"] >= simscore), ]



              if(dim(table.tmp)[1] > 0)
                  output.mat <- rbind(output.mat, table.tmp)
                 #if(dim(table.tmp)[1] > 1 ){
                 #    output.mat <- rbind( output.mat, table.tmp[ which.max(table.tmp[, "%.identity"])  ,] )
                 #} else if(dim(table.tmp)[1] == 1){
                 #    output.mat <- rbind( output.mat, table.tmp )
                 #}
            }
            ###########################################
            # filter according to the highest bit score
            ###########################################
            if(highestBitScore & !peptide.match){

                # maximal bit score
                max.bit <- max(table.tmp[, "bit.score"])

                # check whether there are multiple hits having that score
                max.bit.hits <- table.tmp[which(table.tmp[, "bit.score"] == max.bit ), ]

                # if there is only one hit with that score just append it
                if(dim(max.bit.hits)[1] == 1 ){

                    output.mat <- rbind( output.mat, max.bit.hits )

                } else { # if there are multiple hits, choose the one with highest similarity score

                    # maximal similarity score
                    max.sim <- max(max.bit.hits[, "%.identity"])

                    max.sim.hits <- max.bit.hits[which(max.bit.hits[, "%.identity"] == max.sim) , ]

                    output.mat <- rbind( output.mat, max.sim.hits)

                }

                #output.mat <- rbind( output.mat, table.tmp[ which.max(table.tmp[, "bit.score"])  ,] )
            }
            ###########################################
            # perfect match
            ###########################################
            if(peptide.match & !highestBitScore){

                # get the peptide length: has to be the same as 'alignment length'
                #pepL <- nchar(sub(".*;", "", query))

                pepL <- query.length[query]

                output.mat <- rbind(output.mat, table.tmp[ which( (table.tmp[, "alignment.length"] == pepL) & (table.tmp[, "%.identity"] == 100) )  , ] )

            }

            #######################
            # delta bit score
            ########################
            delta.bit[query] <- table.tmp.ord[1, "bit.score"] - table.tmp.ord[2, "bit.score"]

        } else { # if there was only on blast hit...

            ##################################
            # ... check its similarity score
            ##################################
            if(!highestBitScore & !peptide.match){
               if( table.tmp[1, "%.identity"] >= simscore){
                 output.mat <- rbind( output.mat, table.tmp[1,] )
               }
            }
            #############################
            # bit score: just append it as it is
            #############################
            if(highestBitScore & !peptide.match){
                output.mat <- rbind( output.mat, table.tmp[1,] )
            }

            #############################
            # peptide match
            #############################
            if(peptide.match & !highestBitScore){

                # get the peptide length: has to be the same as 'alignment length'
                pepL <- query.length[query]

                if((pepL == table.tmp[1, "alignment.length"]) & (table.tmp[1, "%.identity"] == 100) )
                    output.mat <- rbind( output.mat, table.tmp[1,] )
            }

            #######################
            # delta bit score
            ########################
            delta.bit[query] <- NA
        }

        count = count + 1
    }

    ###########################
    # append query length
    ###########################
    if(!is.null(fasta.query)){

        #output.mat <-cbind(output.mat, rep("", dim(output.mat)[1]))
        #colnames(output.mat)[dim(output.mat)[2]] <- "query.length"

        query.length <- query.length[ match(output.mat[, "query.id"], names(query.length)) ]
        output.mat <- cbind(output.mat, query.length)

    }


    ############################
    # append delta bit score
    ############################
    #output.mat <-cbind(output.mat, rep("", dim(output.mat)[1]))
    #colnames(output.mat)[dim(output.mat)[2]] <- "delta.bit"

    delta.bit <- delta.bit[ match(output.mat[, "query.id"], names(delta.bit)) ]
    output.mat <- cbind(output.mat, delta.bit)

    return(output.mat)
}

#########################################################################################################
#
#                              parallelized version
#
#
#########################################################################################################
collapseBLASTresults.PAR <- function( table, simscore=100, highestBitScore=T, fasta.query=NULL, cores=2 ){

    require(doSMP)

    colnames(table) <- gsub(" ", ".", colnames(table))

    ################################
    #  output
    ################################
    output.mat <- matrix("", ncol=dim(table)[2], nrow=0 )
    colnames(output.mat) <- colnames( table )

    # the query ids
    query.ids <-  unique(table[, "query.id"])

    ############################################
    # if there is a fasta file specified
    # store the query length in a vector
    ############################################
    if(!is.null(fasta.query)){

        require(seqinr)
        fasta.query <- read.fasta(fasta.query )

        query.length <- vector("numeric", length(query.ids))
        names(query.length) <- query.ids
    }


    # check whether there is a 'no hits found'
    no.hit.idx <- grep("No hits found", query.ids)
    if(length(no.hit.idx) > 0) query.ids <- query.ids[ -no.hit.idx ]

    ######################################
    # vector to store delta bit score
    #######################################
    delta.bit <- vector("numeric", length(query.ids))
    names(delta.bit) <- query.ids


    w <- startWorkers(cores)
    registerDoSMP(w)
    #######################################
    # loop over all query ids
    #######################################
    tmp.mat <- foreach( query = query.ids, .combine=rbind )  %dopar% {


        # subset of the table
        table.tmp <- table[ which(table[, "query.id"]  == query) , ]
        #####################################
        # if there are multiple hits...
        #####################################
        if(dim(table.tmp)[1] > 1){

            ##########################################
            # order according to bit score
            ##########################################
            table.tmp.ord <- table.tmp[ order( table.tmp$bit.score, decreasing=T),  ]

            ##########################################
            # get all hits with the required simscore
            ##########################################
            if(!highestBitScore){
              table.tmp <- table.tmp[  which(table.tmp[, "%.identity"] >= simscore), ]

                 if(dim(table.tmp)[1] > 1 ){
                     out <- table.tmp[ which.max(table.tmp[, "%.identity"])  ,]

                 } else if(dim(table.tmp)[1] == 1){
                     out <- table.tmp
                 }
            }
            ###########################################
            # filter according to the highest bit score
            ###########################################
            if(highestBitScore){
                    out <- table.tmp[ which.max(table.tmp[, "bit.score"])  ,]
            }
            #######################
            # delta bit score
            ########################
            delta.bit <- table.tmp.ord[1, "bit.score"] - table.tmp.ord[2, "bit.score"]

            out <- c(out, delta.bit)

        } else { # if there was only on blast hit...

            ##################################
            # ... check its similarity score
            ##################################
            if(!highestBitScore){
               if( table.tmp[1, "%.identity"] >= simscore){
                 output.mat <- rbind( output.mat, table.tmp[1,] )
               }
            }
            #############################
            # just append it as it is
            #############################
            if(highestBitScore){
                out <- table.tmp[1,]
            }

            #######################
            # delta bit score
            ########################
            delta.bit <- NA
            out <- c(out, delta.bit)

        }
        names(out)[length(out)] <- "delta.bit"

        ##############################
        # query length
        ##############################
        if(!is.null(fasta.query)){
            if(query %in% names(fasta.query)){
                query.length <- length( unlist( fasta.query[query] ) )
                out <- c(out, query.length)
                names(out)[length(out)] <- "query.length"
            }
        }

        out
    }

    stopWorkers(w)

    return(tmp.mat)
}



#########################################################################################################
#
#  peptides        - - blast table, peptides hitting the chromosome
#  peptides.hit    - blast table, peptides hitting a protein and chromosome
#  proteins        - blast table, protein database
#
#  peptides.hit.id - character, regular expression used to annotated peptides hitting a protein and chromosome
#
# changelog:  20100201 implementation
#########################################################################################################
chromosomeMap <- function(peptides, peptides.hit=NULL, proteins=NULL, chromosome=NULL, xlim=NULL, x.at=NULL, add.names=F, label=NULL, peptides.hit.id="^chunk_", main="", cex.main=2, ... ){

    #################################################################
    # replace every blank in the column names of the tables by a dot
    #################################################################
    colnames(peptides) <- gsub(" ", ".", colnames(peptides))
    colnames(peptides) <- tolower(colnames(peptides))

    if(!is.null(peptides.hit)) colnames(peptides.hit) <- gsub(" ", ".", colnames(peptides.hit))
    if(!is.null(proteins)) colnames(proteins) <- gsub(" ", ".", colnames(proteins))

    # determine maximal genomic postion
    if(is.null(xlim)){
        x.max <-max(peptides[, "s..end"])
        xlim <- c(0, x.max)
    # user defined xlim
    } else {
        x.max <- xlim[2]
        # filter the tables such that only peptides/proteins within these coordinates are present
        peptides <- peptides[ ((peptides[, "s..end"] >= xlim[1]) & (peptides[, "s..start"] <= xlim[2])),  ]
        if(!is.null(proteins))
            proteins <- proteins[ ((proteins[, "s..end"] >= xlim[1]) & (proteins[, "s..start"] <= xlim[2])),  ]

        if(!is.null(peptides.hit))
            peptides.hit <- peptides.hit[ ((peptides.hit[, "s..end"] >= xlim[1]) & (peptides.hit[, "s..start"] <= xlim[2])),  ]
    }
    ####################################################
    # make the plot
    ####################################################
    plot.new()
    plot.window(xlim=xlim, ylim=c(0,2), ...)
    title(main, cex.main=2, ...)

    # draw the x-axis
    if(is.null(x.at) ){
        axis(1 )
    } else{
        axis(1, at=seq(0, x.max, x.at) )
    }
    # y axis
    axis(2, labels=F, at=c(0,2))


    #################################################
    #  peptides only hitting the chromosome
    #################################################
    if(dim(peptides)[1] > 0){
       for(i in 1: dim(peptides)[1]){

           # get the genomic region where the peptide matches
           pep.coord <- c(max(peptides[i, "s..start"], xlim[1]), min(peptides[i, "s..end"], xlim[2]))

           lines( y=c(0.5,0.5), x=pep.coord, lwd=4  )


           if(add.names){
               text( x=(pep.coord[1]+round((pep.coord[2]-pep.coord[1])/2)), y=1, sub(";.*", "", peptides[i, "query.id"]), cex=0.9, srt=90  )

               lines( y=c(0.48, 0.52), x=c(pep.coord[1], pep.coord[1]) )
               lines( y=c(0.48, 0.52), x=c(pep.coord[2], pep.coord[2]))

           }
       }
       #############################################
       # if there are at most 5 peptides in that
       # region add a legend containing the
       # corresponding blast table
       #############################################
       if(dim(proteins)[1]  < 6){

          # add some informatino about the quality of the blast hits
          pep.info <- apply(peptides, 1, function(x)paste( x[2:length(x)], collapse="-"))
          legend("topright", legend=pep.info, title=paste("BLAST peptides"), bty="n")

           }


    }
    ##########################################################
    # peptides hitting the chromosome and annotated proteins
    ##########################################################
    if(!is.null(peptides.hit)){

      if(dim(peptides.hit)[1] > 0){

       for(i in 1: dim(peptides.hit)[1]){

           # get the genomic region where the peptide matches
           pep.hit.coord <- c(max(peptides.hit[i, "s..start"], xlim[1]), min(peptides.hit[i, "s..end"], xlim[2]))

           lines( y=c(0.6,0.6), x=pep.hit.coord, lwd=4, col="green"  )

           if(add.names){
               name <- unlist( strsplit( peptides.hit[i, "query.id"], ";") )

               name <- unlist(name[ grep(peptides.hit.id, name)] )[1]


               text( x=(pep.hit.coord[1]+round((pep.hit.coord[2]-pep.hit.coord[1])/2)), y=1.1, labels=name, cex=0.9, srt=90  )


            lines( y=c(0.58, 0.62), x=c(pep.hit.coord[1], pep.hit.coord[1]) )
            lines( y=c(0.58, 0.62), x=c(pep.hit.coord[2], pep.hit.coord[2]) )

           }
       }


     }
    }
    ##########################################################
    # the annotated proteins
    ##########################################################
    if(!is.null(proteins)){

        if(dim(proteins)[1] > 0){
            # order the proteins according to start position
            proteins <- proteins[order(proteins[, "s..start"]),]
            for(i in 1: dim(proteins)[1]){

                if(i %% 2 == 0){
                    y.coord = 0.35
                } else{
                    y.coord = 0.4
                }
                # get the genomic coordinates of the protein
                prot.coord <- c(max(proteins[i, "s..start"], xlim[1]), min(proteins[i, "s..end"], xlim[2]))

                lines(  y= c(y.coord, y.coord), x=prot.coord, lwd=4, col="red"  )


                if(add.names){
                    text( x=(prot.coord[1]+round((prot.coord[2]-prot.coord[1])/2)), y=y.coord-0.1, proteins[i, "query.id"], cex=0.9  )


                    lines( y=c(y.coord - .02, y.coord + .02), x=c(prot.coord[1],prot.coord[1]) )
                    lines( y=c(y.coord - .02, y.coord + .02), x=c(prot.coord[2], prot.coord[2]))

                }
            }
            #############################################
            # if there are at most 5 peptides in that
            # region add a legend containing the
            # corresponding blast table
            #############################################
            if(dim(proteins)[1]  < 6){

               prot.info <- apply(proteins, 1, function(x)paste( x[2:length(x)], collapse="-"))
               legend("topleft", legend=prot.info, title=paste("BLAST proteins"), bty="n")

           }


        }
    }

    ###########################################################
    # add a label if specified
    ###########################################################
    if(!is.null(label))
        text( xlim[1]+(xlim[2]-xlim[1])/2, 1.9, label, cex=1  )

}

##########################################################################################################
chromosomeMap2 <- function(peptides, proteins=NULL, chromosome=NULL, xlim=NULL, x.at=NULL, add.names=F){

    # determine maximal genomic postion of a peptide
    if(is.null(xlim)){
        x.max <- 0

        # loop over all peptide blast tables....
        for(i in 1:length(peptides))
            x.max <- max(x.max, peptides[[i]][, "s. end"])
        xlim <- c(0, x.max)

    } else {
        x.max <- xlim[2]

     #           if(!is.null(proteins))
    #        proteins <- proteins[ ((proteins[, "s. end"] >= xlim[1]) & (proteins[, "s. start"] <= xlim[2])),  ]

    #    if(!is.null(peptides.hit))
     #       peptides.hit <- peptides.hit[ ((peptides.hit[, "s. end"] >= xlim[1]) & (peptides.hit[, "s. start"] <= xlim[2])),  ]
    }


    plot.new()
    plot.window(xlim=xlim, ylim=c(0, length(peptides) + length(proteins)))

    # draw the x-axis
    if(is.null(x.at) ){
        axis(1 )
    } else{
        axis(1, at=seq(0, x.max, x.at) )
    }
    # y axis
    axis(2, labels=F, at=c(0, length(peptides) + length(proteins)))

    #################################################
    # loop over all peptide tables
    #
    #################################################
    for(pep.tab in 1:length(peptides)){

       for(i in 1: dim(peptides[[pep.tab]])[1]){

           # and filter the tables such that only peptides/proteins within these coordinates are present
           peptides[[pep.tab]] <- peptides[[pep.tab]][ ((peptides[[pep.tab]][, "s. end"] >= xlim[1]) & (peptides[[pep.tab]][, "s. start"] <= xlim[2])),  ]


           # get the genomic region where the peptide matches
           X.coord <- c(max(peptides[[pep.tab]][i, "s. start"], xlim[1]), min(peptides[[pep.tab]][i, "s. end"], xlim[2]))

           # plot in the middle
           lines( y=c( (pep.tab + 1)+0.5, (pep.tab + 1 )+0.5), x=X.coord, lwd=4  )


           if(add.names){
               text( x=(X.coord[1]+round((X.coord[2]-X.coord[1])/2)), y=1, sub(";.*", "", peptides[[pep.tab]][i, "Query id"]), cex=0.9, srt=90  )

               lines( y=c(pep.tab-0.05, pep.tab + 0.05), x=c(X.coord[1], X.coord[1]) )
               lines( y=c(pep.tab-0.05, pep.tab + 0.05), x=c(X.coord[2], X.coord[2]) )

           }
       }
    }

    ##########################################################
    # the annotated proteins
    ##########################################################
    if(!is.null(proteins)){

        # order the proteins according to start position
        proteins <- proteins[order(proteins[, "s. start"]),]
        for(i in 1: dim(proteins)[1]){

           if(i %% 2 == 0){
               y.coord = 0.35
           } else{
               y.coord = 0.4
           }
       # get the genomic coordinates of the protein
       prot.coord <- c(max(proteins[i, "s. start"], xlim[1]), min(proteins[i, "s. end"], xlim[2]))

       lines(  y= c(y.coord, y.coord), x=prot.coord, lwd=4, col="red"  )


       if(add.names){
            text( x=(prot.coord[1]+round((prot.coord[2]-prot.coord[1])/2)), y=y.coord-0.1, proteins[i, "Query id"], cex=0.9  )


            lines( y=c(y.coord - .02, y.coord + .02), x=c(prot.coord[1],prot.coord[1]) )
            lines( y=c(y.coord - .02, y.coord + .02), x=c(prot.coord[2], prot.coord[2]))

        }
    }


    }


}

####################################################################
#
#  peptides   - matrix returned by function 'collapseBLASTresult'
#               that contains all peptides hitting only the chromosome
#
#  proteins   - the annotated/predicted coding regions
#
# the function tries to find discrepancies between the predicted coding
# regions and the observed peptides, e.g.
#    - peptides spanning N/C terminus of a protein
#    - peptides that are very close to annotated regions
#    - ...
#
####################################################################
findHits <- function( peptides, pep.prot=NULL, proteins, gap.dist=50, plot=F, file="interestingHits.pdf"  ){

    # lists to store indices
    ol.N <- ol.C <- ad.N <- ad.C <- list()

    # counting variables
    ol.N.count <- ol.C.count <- ad.N.count <- ad.C.count <- 1

    # loop over the peptides
    for( p in 1:dim(peptides)[1]  ){

        # check if the current peptide overlaps with a N-terminus of an annotated protein
        ol.N.idx <- which(( peptides[p, "s. start"] < proteins[, "s. start"]) & (peptides[p, "s. end"] > proteins[, "s. start"]) )
        # check if the current peptide overlaps with a C-terminus of an annotated protein
        ol.C.idx <- which( ( peptides[p, "s. start"] < proteins[, "s. end"]) & (peptides[p, "s. end"] > proteins[,"s. end"]) )
        # check if the current peptide is adjacent to a N-terminus
        ad.N.idx <- which( ( proteins[, "s. start"] - peptides[p, "s. end"] ) < gap.dist &  ( proteins[, "s. start"] - peptides[p, "s. end"] ) > 0  )
        # check if the current peptide is adjacent to a C-terminus
        ad.C.idx <- which( ( peptides[p, "s. start"] - proteins[, "s. end"] ) < gap.dist & ( peptides[p, "s. start"] - proteins[, "s. end"] ) > 0  )

        if( length( ol.N.idx) > 0){
            ol.N[[ol.N.count]] <-ol.N.idx
             names(ol.N)[length(ol.N)] <- peptides[p, "Query id"]
            ol.N.count <- ol.N.count + 1
        } else if( length( ol.C.idx) > 0){

            ol.C[[ol.C.count]] <-  ol.C.idx
            names(ol.C)[length(ol.C)] <- peptides[p, "Query id"]
            ol.C.count <- ol.C.count + 1

        } else if(length(ad.N.idx) > 0){

            ad.N[[ad.N.count]] <-  ad.N.idx
            names(ad.N)[length(ad.N)] <- peptides[p, "Query id"]
            ad.N.count <- ad.N.count + 1
        } else if(length(ad.C.idx) > 0 ){

            ad.C[[ad.C.count]] <- ad.C.idx
            names(ad.C)[length(ad.C)] <- peptides[p, "Query id"]
            ad.C.count <- ad.C.count
        }
    }

    #######################################
    # find clusters of peptides matching
    # only the genome
    #######################################
    pep.clust <- findPeptideCluster(peptides, gap.dist=gap.dist)

    #######################################
    # the output list
    #######################################
    out <- list()

    out[[1]] <- lapply(ol.N, function(x) proteins[x, ]  )    #proteins[ unlist(ol.N), ]
    out[[2]] <- lapply(ol.C, function(x) proteins[x, ]  ) #proteins[ unlist(ol.C), ]
    out[[3]] <- lapply(ad.N, function(x) proteins[x, ]  ) #proteins[ unlist(ad.N), ]
    out[[4]] <- lapply(ad.C, function(x) proteins[x, ]  ) #proteins[ unlist(ad.C), ]
    out[[5]] <- pep.clust

    names(out) <- c("overlap.left", "overlap.right", "adjacent.left", "adjacent.right", "peptide.cluster")

    ###########################################
    #    plot
    ###########################################
    if(plot){
        if( is.null(pep.prot) )warning("blubb\n")

        pdf( paste(file, sep="" ), height=8, width=16 )
        for(i in 1:length(out)){
            for(j in 1:length(out[[i]])){

                # margins for the figure
                if(i ==1 | i == 3) {l.tol <- 5e2; r.tol = 1e2
                                }else { l.tol=1e2; r.tol=5e2 }
                x.min <- min( out[[i]][[j]][ , "s. start"] ) - l.tol
                x.max <- max( out[[i]][[j]][ , "s. end" ] ) + r.tol


                chromosomeMap(peptides, pep.prot, proteins, xlim=c( x.min, x.max), add.names=T, label=paste( names(out)[i], ",max. gap", gap.dist ) )

            }
        }
             dev.off()
    }

    return( out  )
}

########################################################################
#
#
#
#
findPeptideCluster <- function(pep.blast, gap.dist=50 ){

    # sort according to start position
    pep.blast <- pep.blast[order(pep.blast[, "s. start"]), ]

    #
    pep.clust <- list()
    clust.count = 1
    pep.count = 1
    ok=T

    while(ok){

       # cat("pep:",pep.count, "\n")

        if(pep.blast[pep.count+1, "s. start"] - pep.blast[pep.count, "s. end"]  < gap.dist){

            out <- pep.blast[pep.count:(pep.count+1), ]

            #cat(out[1,"s. start"], "\n", out[2, "s. start"], "\n")

            # check if you can extend this cluster

            extend = T
            ext.count = pep.count+1

            while(extend){

                if(pep.blast[ext.count+1, "s. start"] - pep.blast[ext.count, "s. end"]  < gap.dist){
                    out <- rbind(out, pep.blast[ext.count:(ext.count + 1), ])

                    #cat("ext:", ext.count, "\n")

                    #cat(out[dim(out)[1]-1,"s. start"], "\n", out[dim(out)[1], "s. start"], "\n")

                    ext.count <- ext.count + 1
                } else {
                    extend=F

                    #cat("\n")
                }

            }

        pep.clust[[clust.count]] <- out
        clust.count <- clust.count + 1

        pep.count <- ext.count + 1

        } else {
         pep.count <- pep.count + 1

        }

        if(pep.count >= dim(pep.blast)[1])
            ok=F
    }

    return(pep.clust)
}


######################################################################################################################
#
# - based on 'heatmap.2'
#
# changes:
#   - hclust.method         - character, method used in function 'hclust'
#   - dist.method           - character, method used in function 'dist'
#
######################################################################################################################
heatmap.3 <- function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
    distfun = dist, hclustfun = hclust, dendrogram = c("both",
        "row", "column", "none"), symm = FALSE, scale = c("none",
        "row", "column"), na.rm = TRUE, revC = identical(Colv,
        "Rowv"), add.expr, breaks, symbreaks = min(x < 0, na.rm = TRUE) ||
        scale != "none", col = "heat.colors", colsep, rowsep,
    sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, notecex = 1,
    notecol = "cyan", na.color = par("bg"), trace = c("column",
        "row", "both", "none"), tracecol = "cyan", hline = median(breaks),
    vline = median(breaks), linecol = tracecol, margins = c(5,
        5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr),
    cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL,
    key = TRUE, keysize = 1.5, density.info = c("histogram",
        "density", "none"), denscol = tracecol, symkey = min(x <
        0, na.rm = TRUE) || symbreaks, densadj = 0.25, main = NULL,
    xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL,
    hclust.method="average", dist.method="euclidian",


    ...)
{
    require(gplots)

    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale))
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col))
        col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none"))
        warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv))
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv))
        Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv))
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote))
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv))
                dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv))
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x, dist.method), hclust.method)
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x, dist.method), hclust.method)
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc)
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm)
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm)
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow))
        labRow <- if (is.null(rownames(x)))
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
        labCol <- if (is.null(colnames(x)))
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) <
        1) {
        if (missing(col) || is.function(col))
            breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks)
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function")
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei))
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid))
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
        if (!missing(ColSideColors)) {
            if (!is.character(ColSideColors) || length(ColSideColors) !=
                nc)
                stop("'ColSideColors' must be a character vector of length ncol(x)")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] +
                1)
            lhei <- c(lhei[1], 0.2, lhei[2])
        }
        if (!missing(RowSideColors)) {
            if (!is.character(RowSideColors) || length(RowSideColors) !=
                nr)
                stop("'RowSideColors' must be a character vector of length nrow(x)")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) -
                1), 1), lmat[, 2] + 1)
            lwid <- c(lwid[1], 0.2, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
    if (length(lhei) != nrow(lmat))
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    if (!symm || scale != "none") {
        x <- t(x)
        cellnote <- t(cellnote)
    }
    if (revC) {
        iy <- nr:1
        if (exists("ddr"))
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 +
        c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col,
        breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr"))
        retval$rowDendrogram <- ddr
    if (exists("ddc"))
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) {
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
            col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
        cex.axis = cexCol)
    if (!is.null(xlab))
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
        cex.axis = cexRow)
    if (!is.null(ylab))
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
        eval(substitute(add.expr))
    if (!missing(colsep))
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0,
            length(csep)), xright = csep + 0.5 + sepwidth[1],
            ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1,
            col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) +
            1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) +
            1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1,
            col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol,
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote))
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
            col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main))
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x))
            tmpbreaks[length(tmpbreaks)] <- max(abs(x))
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }
        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row")
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column")
            mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, "Value", line = 2)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
                lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
                col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        else title("Color Key")
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
        high = retval$breaks[-1], color = retval$col)
    invisible(retval)
}


















########################################################################################################################




###############################################
## Plot a 2-Way, 3-Way or 4-Way Venn Diagram ##
###############################################
## Author: Thomas Girke
## Last update: Nov 6, 2008
## Utility: Plots a non-proportional 2-, 3- or 4-way venn diagram based on overlaps among data sets (vectors)
## Detailed instructions for running this script are available on this page:
##     http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/R_BioCondManual.html#R_graphics_venn

## Define venndiagram function
venndiagram <- function(x=x, y=y, z=z, w=w, unique=T, title="Venn Diagram", labels=c("x", "y", "z", "w"), lines=1, lcol=1, tcol=1, diacol=1, plot=T, type="3", printsub=TRUE, ...) {
	## Remove duplicates and NA fields in x, y, z and w
	if(unique==T) {
		x <- unique(x); x <- as.vector(na.omit(x))
		y <- unique(y); y <- as.vector(na.omit(y))
		if(!missing("z")) {
			z <- unique(z); z <- as.vector(na.omit(z))
		}
		if(!missing("w")) {
			w <- unique(w); w <- as.vector(na.omit(w))
		}
	}

	## Check valid type selection
	if(!type %in% c("2", "2map", "3", "3map", "4", "4map", "4el", "4elmap")) {
		return("Error: the 'type' argument can only be set to one of these values: 2, 2map, 3, 3map, 4, 4map, 4el, 4elmap.")
	}

	## Plot a 2-way venn diagram
	if(type=="2") {
		# Define ovelap queries
		q1 <- x[x %in% y]
		q2 <- x[!x %in% y]
		q3 <- y[!y %in% x]

		## Store query vectors in list
		qlist <- list(q1=q1, q2=q2, q3=q3)

		## Perfom query counts
		count <- unlist(lapply(qlist, length))
		countDF <- data.frame(query=names(count) , count=as.vector(count))
		olDF <- data.frame(x=c(5.0, 3.1, 7.0), y=c(6.1, 6.1, 6.1), count=countDF$count)
		if(printsub==TRUE) {mysub <- paste(paste("N unique: xy =", length(unique(c(x,y)))), paste("; x =", length(unique(x))), paste("; y =", length(unique(y))), sep="")} else {mysub <- ""}
		if(plot==T) {
			## Plot the 2-way venn diagram
			symbols(x=c(4, 6), y = c(6, 6), circles=c(2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main=title, sub=mysub, xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...);
			text(olDF$x, olDF$y, olDF$count, col=tcol, ...); text(c(2.0, 8.0), c(8.8, 8.8), labels[1:2], col=lcol, ...)
		}

		## Return query list
		return(qlist)
	}

	## Plot 2-way mapping venn diagram
	if(type=="2map") {
		olDFdebug <- data.frame(x=c(5.0, 3.1, 7.0), y=c(6.1, 6.1, 6.1), count=paste("q", 1:3, sep=""), ...)
		symbols(x=c(4, 6), y = c(6, 6), circles=c(2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main="Mapping Venn Diagram", xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...);
		text(olDFdebug$x, olDFdebug$y, olDFdebug$count, col=tcol, ...); text(c(2.0, 8.0), c(8.8, 8.8), paste(labels[1:2], "=", c("x","y")), col=lcol, ...)
	}

	## Plot a 3-way venn diagram
	if(type=="3") {
		## Define ovelap queries
		q1 <- x[x %in% y & x %in% z]
		q2 <- x[x %in% z]; q2 <- q2[!q2 %in% y]
		q3 <- y[y %in% z]; q3 <- q3[!q3 %in% x]
		q4 <- x[x %in% y]; q4 <- q4[!q4 %in% z]
		q5 <- x[!x %in% y]; q5 <- q5[!q5 %in% z]
		q6 <- y[!y %in% z]; q6 <- q6[!q6 %in% x]
		q7 <- z[!z %in% x]; q7 <- q7[!q7 %in% y]

		## Store query vectors in list
		qlist <- list(q1=q1, q2=q2, q3=q3, q4=q4, q5=q5, q6=q6, q7=q7)

		## Perfom query counts
		count <- unlist(lapply(qlist, length))
		countDF <- data.frame(query=names(count) , count=as.vector(count))
		olDF <- data.frame(x=c(5.0, 3.8, 6.3, 5.0, 3.0, 7.0, 5.0), y=c(5.6, 4.6, 4.6, 6.9, 6.5, 6.5, 3.0), count=countDF$count)
		if(printsub==TRUE) {mysub <- paste(paste("N unique: xyz =", length(unique(c(x,y,z)))), paste("; x =", length(unique(x))), paste("; y =", length(unique(y))), paste("; z =", length(unique(z))), sep="")} else { mysub <- "" }
		if(plot==T) {
			## Plot the 3-way venn diagram
			symbols(x=c(4, 6, 5), y = c(6, 6, 4), circles=c(2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main=title, sub=mysub, xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...);
			text(olDF$x, olDF$y, olDF$count, col=tcol, ...); text(c(2.0, 8.0, 6.0), c(8.8, 8.8, 1.1), labels[1:3], col=lcol, ...)
		}

		## Return query list
		return(qlist)
	}

	## Plot 3-way mapping venn diagram
	if(type=="3map") {
		olDFdebug <- data.frame(x=c(5.0, 3.8, 6.3, 5.0, 3.0, 7.0, 5.0), y=c(5.6, 4.6, 4.6, 6.9, 6.5, 6.5, 3.0), count=paste("q", 1:7, sep=""), ...)
		symbols(x=c(4, 6, 5), y = c(6, 6, 4), circles=c(2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main="Mapping Venn Diagram", xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...);
		text(olDFdebug$x, olDFdebug$y, olDFdebug$count, col=tcol, ...); text(c(2.0, 8.0, 6.0), c(8.8, 8.8, 1.1), paste(labels[1:3], "=", c("x","y","z")), col=lcol, ...)
	}

	## Overlap queries for 4-way venn diagram
	if(type=="4" | type=="4el" | type=="4elmap") {
		## Define ovelap queries
		xy <- x[x %in% y]; xz <-x[x %in% z]; xw <- x[x %in% w]; yz <- y[y %in% z]; yw <- y[y %in% w]; zw <- z[z %in% w]
		q1 <- xy[xy %in% zw]
		q2 <- xw[xw %in% z]; q2 <- q2[!q2 %in% y]
		q3 <- yz[yz %in% w]; q3 <- q3[!q3 %in% x]
		q4 <- yz[yz %in% x]; q4 <- q4[!q4 %in% w]
		q5 <- xw[xw %in% y]; q5 <- q5[!q5 %in% z]
		q6 <- xy[!xy %in% z]; q6 <- q6[!q6 %in% w]
		q7 <- zw[!zw %in% x]; q7 <- q7[!q7 %in% y]
		q8 <- xz[!xz %in% y]; q8 <- q8[!q8 %in% w]
		q9 <- yw[!yw %in% x]; q9 <- q9[!q9 %in% z]
		q10 <- x[!x %in% c(y,z,w)]
		q11 <- y[!y %in% c(x,z,w)]
		q12 <- z[!z %in% c(x,y,w)]
		q13 <- w[!w %in% c(x,y,z)]
		q14 <- xw[!xw %in% y]; q14 <- q14[!q14 %in% z]
		q15 <- yz[!yz %in% x]; q15 <- q15[!q15 %in% w]

		## Store query vectors in list
		qlist <- list(q1=q1, q2=q2, q3=q3, q4=q4, q5=q5, q6=q6, q7=q7, q8=q8, q9=q9, q10=q10, q11=q11, q12=q12, q13=q13, q14=q14, q15=q15)

		## Perfom query counts
		count <- unlist(lapply(qlist, length))
		countDF <- data.frame(query=names(count) , count=as.vector(count))
		olDF <- data.frame(x=c(4.8, 3.9, 5.7, 3.9, 5.7, 4.8, 4.8, 3.0, 6.5, 3.0, 6.5, 3.0, 6.5, 4.8, 4.8), y=c(5.2, 4.2, 4.2, 6.3, 6.3, 7.2, 3.2, 5.2, 5.2, 7.2, 7.2, 3.2, 3.2, 1.0, 0.4), count=countDF$count)

		if(printsub==TRUE) {mysub <- paste(paste("N unique: xyzw =", length(unique(c(x,y,z,w)))), paste("; x =", length(unique(x))), paste("; y =", length(unique(y))), paste("; z =", length(unique(z))), paste("; w =", length(unique(w))), sep="") } else { mysub <- "" }

	## Plot 4-way venn diagram as circles
		if(plot==T & type=="4") {
			symbols(x=c(4, 5.5, 4, 5.5), y = c(6, 6, 4.5, 4.5), circles=c(2, 2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main=title, sub=mysub, xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...);
			text(olDF$x[1:13], olDF$y[1:13], olDF$count[1:13], col=tcol, ...) # rows 14-15 of olDF are printed in last step
			text(c(2.0, 7.5, 2.0, 7.5), c(8.3, 8.3, 2.0, 2.0), labels, col=lcol, ...)
			text(c(3.8, 3.8), c(1.0, 0.4), c(paste("Only in ", labels[1], " & ", labels[4], ": ", olDF$count[14], sep=""), paste("Only in ", labels[2], " & ", labels[3], ": ", olDF$count[15], sep="")), col=diacol, ...)
		}

	## Plot 4-way venn diagram as ellipses
	if(plot==T & (type=="4el" | type=="4elmap")) {
		olDF <- data.frame(x=c(5.0, 4.2, 6.4, 3.6, 5.8, 2.9, 7.1, 3.1, 6.9, 1.5, 3.5, 6.5, 8.5, 5.0, 5.0), y=c(2.8, 1.4, 4.0, 4.0, 1.4, 5.9, 5.9, 2.2, 2.2, 4.8, 7.2, 7.2, 4.8, 0.7, 6.0), count=countDF$count)
		## Plot ellipse
		plotellipse <- function (center=c(1,1), radius=c(1,2), rotate=1, segments=360, xlab="", ylab="", ...) {
			angles <- (0:segments) * 2 * pi/segments
			rotate <- rotate*pi/180
			ellipse <- cbind(radius[1] * cos(angles), radius[2] * sin(angles))
			ellipse <- cbind( ellipse[,1]*cos(rotate) + ellipse[,2]*sin(rotate), ellipse[,2]*cos(rotate) - ellipse[,1]*sin(rotate) )
			ellipse <- cbind(center[1]+ellipse[,1], center[2]+ellipse[,2])
			plot(ellipse, type = "l", xlim = c(0, 10), ylim = c(0, 10), xlab = "", ylab = "", ...)
		}
		## Plot ellipse as 4-way venn diagram
		ellipseVenn <- function(lines=lines, olDF, title=title, labels=labels, sub=mysub, main, lcol=lcol, tcex=1.3, ...) {
			split.screen(c(1,1))
			plotellipse(center=c(3.5,3.6), radius=c(2,4), rotate=-35, segments=360, xlab="", ylab="", col=lines[1], axes=FALSE, main=title, sub=mysub, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(4.7,4.4), radius=c(2,4), rotate=-35, segments=360, xlab="", ylab="", col=lines[2], axes=FALSE, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(5.3,4.4), radius=c(2,4), rotate=35, segments=360, xlab="", ylab="", col=lines[3], axes=FALSE, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(6.5,3.6), radius=c(2,4), rotate=35, segments=360, xlab="", ylab="", col=lines[4], axes=FALSE, ...)
			text(olDF[1:15,1], olDF[1:15,2], olDF[1:15,3], col=tcol, ...)
			text(c(0.4, 2.8, 7.5, 9.4), c(7.3, 8.3, 8.3, 7.3), labels, col=lcol, ...)
			close.screen(all=TRUE)
		}
		## Plot 4-way ellipse venn diagram
		if(type=="4el") {
			ellipseVenn(olDF=olDF, lcol=lcol, lines=lines, labels=labels, title=title, ...)
		}

		## Plot 4-way ellipse mapping venn diagram
		if(type=="4elmap") {
			olDFdebug <- data.frame(x=c(5.0, 4.2, 6.4, 3.6, 5.8, 2.9, 7.1, 3.1, 6.9, 1.5, 3.5, 6.5, 8.5, 5.0, 5.0), y=c(2.8, 1.4, 4.0, 4.0, 1.4, 5.9, 5.9, 2.2, 2.2, 4.8, 7.2, 7.2, 4.8, 0.7, 6.0), count=paste("q", 1:15, sep=""), ...)
			ellipseVenn(olDF=olDFdebug, lcol=lcol, lines=lines, labels=paste(labels, "=", c("x","y","z","w")), title="Mapping Venn Diagram", ...)
		}
	}

	## Return query list
	return(qlist)
	}

	## Plot 4-way circle mapping venn diagram
	if(type=="4map") {
		olDFdebug <- data.frame(x=c(4.8, 3.9, 5.7, 3.9, 5.7, 4.8, 4.8, 3.0, 6.5, 3.0, 6.5, 3.0, 6.5, 4.8, 4.8), y=c(5.2, 4.2, 4.2, 6.3, 6.3, 7.2, 3.2, 5.2, 5.2, 7.2, 7.2, 3.2, 3.2, 1.0, 0.4), count=paste("q", 1:15, sep=""), ...)
		symbols(x=c(4, 5.5, 4, 5.5), y = c(6, 6, 4.5, 4.5), circles=c(2, 2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main="Mapping Venn Diagram", xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...);
		text(olDFdebug$x[1:13], olDFdebug$y[1:13], olDFdebug$count[1:13], col=tcol, ...); text(c(2.0, 7.5, 2.0, 7.5), c(8.3, 8.3, 2.0, 2.0), paste(labels, "=", c("x","y","z","w")), col=lcol, ...)
		text(c(3.8, 3.8), c(0.97, 0.36), c(paste("Only in ", labels[1], " & ", labels[4], ": ", olDFdebug$count[14], sep=""), paste("Only in ", labels[2], " & ", labels[3], ": ", olDFdebug$count[15], sep="")), col=tcol, ...)
	}

}

## Generate overlap reports
olReport <- function(qlist=qlist, missing=".", type=1) {
	## Check valid type selection
	if(!type %in% c(1, 2, 3, 4)) {
		return("Error: the 'type' argument can only be set to the values: 1, 2 or 3.")
	}

	## Output data frame with overlap keys in separate columns
	if(type==1) {
		ids <- sort(unique(as.vector(unlist(qlist))))
		qDF <- matrix(ids, nrow=length(ids), ncol=length(qlist), dimnames=list(1:length(ids), names(qlist)))
		lqDF <- as.data.frame(lapply(names(qlist), function(x) qDF[,x] %in% qlist[[x]]))
		colnames(lqDF) <- colnames(qDF)
		lqDF <- as.matrix(lqDF)
		qDF[!lqDF] <- missing
		qDF <- data.frame(IDs=ids, qDF)
		return(qDF)
	}

	## Output data frame with overlap section numbers (qNo) in one column
	if(type==3) {
		collapsedDF <- data.frame(IDs=as.vector(unlist(qlist)), qNo=rep(names(qlist), sapply(qlist, length)))
		collapsedDF <- collapsedDF[order(collapsedDF$IDs), ]
		rownames(collapsedDF) <- 1:length(collapsedDF[, 1])
		return(collapsedDF)
	}

	## Output data frame with overlap counts
	if(type==2) {
		qStat <- data.frame(count=sapply(qlist, length))
		return(qStat)
	}

	## Output presence-absence matrix
	if(type==4) {
		ids <- sort(unique(as.vector(unlist(qlist))))
		qDF <- matrix(ids, nrow=length(ids), ncol=length(qlist), dimnames=list(1:length(ids), names(qlist)))
		lqDF <- as.data.frame(lapply(names(qlist), function(x) qDF[,x] %in% qlist[[x]]))
		colnames(lqDF) <- colnames(qDF)
		lqDF <- as.matrix(lqDF)
		lqDF[lqDF] <- 1
		lqDF[!lqDF] <- 0
		rownames(lqDF) <- ids
		lqDF <- lqDF[names(rev(sort(rowSums(lqDF)))),]
		return(lqDF)
	}
}

