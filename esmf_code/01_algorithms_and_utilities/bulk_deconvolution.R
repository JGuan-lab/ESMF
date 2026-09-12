
deconv_lsc<-function(T, C, method, phenoDataC, P = NULL, elem = NULL, STRING = NULL, marker_distrib, refProfiles.var){ 
  
  bulk_methods = c("CIBERSORT","DeconRNASeq","OLS","nnls","FARDEEP","RLR","DCQ","elastic_net","lasso","ridge","EPIC","DSA","ssKL","ssFrobenius","dtangle")
  sc_methods = c("MuSiC","BisqueRNA","DWLS","deconvSeq","SCDC")
  
  ########## Using marker information for bulk_methods
  if(method %in% bulk_methods){
    
    C = C[rownames(C) %in% marker_distrib$gene,]
    T = T[rownames(T) %in% marker_distrib$gene,]
    refProfiles.var = refProfiles.var[rownames(refProfiles.var) %in% marker_distrib$gene,]
    
  } else { ### For scRNA-seq methods 
    
    #BisqueRNA requires "SubjectName" in phenoDataC
    if(length(grep("[N-n]ame",colnames(phenoDataC))) > 0){
      sample_column = grep("[N-n]ame",colnames(phenoDataC))
    } else {
      sample_column = grep("[S-s]ample|[S-s]ubject",colnames(phenoDataC))
    }
    colnames(phenoDataC)[sample_column] = "SubjectName"
    rownames(phenoDataC) = phenoDataC$cellID
    
    require(xbioc)
    C.eset <- Biobase::ExpressionSet(assayData = as.matrix(C),phenoData = Biobase::AnnotatedDataFrame(phenoDataC))
    T.eset <- Biobase::ExpressionSet(assayData = as.matrix(T))
    
  }
  
  ##########    MATRIX DIMENSION APPROPRIATENESS    ##########
  keep = intersect(rownames(C),rownames(T)) 
  C = C[keep,]
  T = T[keep,]
  
  ###################################
  if(method=="CIBERSORT"){ #without QN. By default, CIBERSORT performed QN (only) on the mixture.
    
    RESULTS = CIBERSORT(sig_matrix = C, mixture_file = T, QN = FALSE) 
    RESULTS = t(RESULTS[,1:(ncol(RESULTS)-3)])
    
  } else if(method=="DeconRNASeq"){ #nonnegative quadratic programming; lsei function (default: type=1, meaning lsei from quadprog)
    #datasets and reference matrix: signatures, need to be non-negative. 
    #"use.scale": whether the data should be centered or scaled, default = TRUE
    unloadNamespace("Seurat") #needed for PCA step
    library(pcaMethods) #needed for DeconRNASeq to work
    RESULTS = t(DeconRNASeq::DeconRNASeq(datasets = as.data.frame(T), signatures = as.data.frame(C), proportions = NULL, checksig = FALSE, known.prop = FALSE, use.scale = FALSE, fig = FALSE)$out.all)
    colnames(RESULTS) = colnames(T)
    require(Seurat)
    
  } else if (method=="OLS"){
    
    RESULTS = apply(T,2,function(x) lm(x ~ as.matrix(C))$coefficients[-1])
    RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
    RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
    rownames(RESULTS) <- unlist(lapply(strsplit(rownames(RESULTS),")"),function(x) x[2]))
    
  } else if (method=="nnls"){
    
    require(nnls)
    RESULTS = do.call(cbind.data.frame,lapply(apply(T,2,function(x) nnls::nnls(as.matrix(C),x)), function(y) y$x))
    RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
    rownames(RESULTS) <- colnames(C)
    
  } else if (method=="FARDEEP"){
    
    require(FARDEEP)
    RESULTS = t(FARDEEP::fardeep(C, T, nn = TRUE, intercept = TRUE, permn = 10, QN = FALSE)$abs.beta)
    RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
    
  } else if (method=="RLR"){ #RLR = robust linear regression
    
    require(MASS)
    RESULTS = do.call(cbind.data.frame,lapply(apply(T,2,function(x) MASS::rlm(x ~ as.matrix(C), maxit=100)), function(y) y$coefficients[-1]))
    RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
    RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
    rownames(RESULTS) <- unlist(lapply(strsplit(rownames(RESULTS),")"),function(x) x[2]))
    
  } else if (method=="DCQ"){#default: alpha = 0.05, lambda = 0.2. glmnet with standardize = TRUE by default
    
    require(ComICS)
    RESULTS = t(ComICS::dcq(reference_data = C, mix_data = T, marker_set = as.data.frame(row.names(C)) , alpha_used = 0.05, lambda_min = 0.2, number_of_repeats = 10)$average)
    RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
    RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
    
  } else if (method=="elastic_net"){#standardize = TRUE by default. lambda=NULL by default 
    
    require(glmnet)# gaussian is the default family option in the function glmnet. https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html
    RESULTS = apply(T, 2, function(z) coef(glmnet::glmnet(x = as.matrix(C), y = z, alpha = 0.2, standardize = TRUE, lambda = glmnet::cv.glmnet(as.matrix(C), z)$lambda.1se))[1:ncol(C)+1,])
    RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
    RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
    
  } else if (method=="ridge"){ #alpha=0
    
    require(glmnet)
    RESULTS = apply(T, 2, function(z) coef(glmnet::glmnet(x = as.matrix(C), y = z, alpha = 0, standardize = TRUE, lambda = glmnet::cv.glmnet(as.matrix(C), z)$lambda.1se))[1:ncol(C)+1,])
    RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
    RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
    
  } else if (method=="lasso"){ #alpha=1; shrinking some coefficients to 0. 
    
    require(glmnet)
    RESULTS = apply(T, 2, function(z) coef(glmnet::glmnet(x = as.matrix(C), y = z, alpha = 1, standardize = TRUE, lambda = glmnet::cv.glmnet(as.matrix(C), z)$lambda.1se))[1:ncol(C)+1,])
    RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
    RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
    
    RESULTS[is.na(RESULTS)] <- 0 #Needed for models where glmnet drops all terms of a model and fit an intercept-only model (very unlikely but possible).
    
  } else if (method=="EPIC"){
    
    require(EPIC)
    marker_distrib = marker_distrib[marker_distrib$gene %in% rownames(C),]
    markers = as.character(marker_distrib$gene)
    C_EPIC <- list()
    
    common_CTs <- intersect(colnames(C),colnames(refProfiles.var))
    
    C_EPIC[["sigGenes"]] <- rownames(C[markers,common_CTs])
    C_EPIC[["refProfiles"]] <- as.matrix(C[markers,common_CTs])
    C_EPIC[["refProfiles.var"]] <- refProfiles.var[markers,common_CTs]
    
    RESULTS <- t(EPIC::EPIC(bulk=as.matrix(T), reference=C_EPIC, withOtherCells=TRUE, scaleExprs=FALSE)$cellFractions) #scaleExprs=TRUE by default: only keep genes in common between matrices
    RESULTS = RESULTS[!rownames(RESULTS) %in% "otherCells",]
    
  } else if (method=="DSA"){ #DSA algorithm assumes that the input mixed data are in linear scale; If log = FALSE the data is left unchanged
    
    require(CellMix)
    md = marker_distrib
    ML = CellMix::MarkerList()
    ML@.Data <- tapply(as.character(md$gene),as.character(md$CT),list)
    RESULTS = CellMix::ged(as.matrix(T), ML, method = "DSA", log = FALSE)@fit@H
    RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
    RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
    
  } else if (method=="ssKL"){ 
    
    require(CellMix)
    md = marker_distrib #Full version, irrespective of C
    ML = CellMix::MarkerList()
    ML@.Data <- tapply(as.character(md$gene),as.character(md$CT),list)
    
    RESULTS <- CellMix::ged(as.matrix(T), ML, method = "ssKL", sscale = FALSE, maxIter=500, log = FALSE)@fit@H 
    
  } else if (method=="ssFrobenius"){
    
    require(CellMix)
    md = marker_distrib #Full version, irrespective of C
    ML = CellMix::MarkerList()
    ML@.Data <- tapply(as.character(md$gene),as.character(md$CT),list)
    
    RESULTS <- CellMix::ged(as.matrix(T), ML, method = "ssFrobenius", sscale = TRUE, maxIter = 500, log = FALSE)@fit@H #equivalent to coef(CellMix::ged(T,...)
    
  } else if(method=="dtangle"){#Only works if T & C are log-transformed
    
    require(dtangle)
    mixture_samples = t(T)
    reference_samples = t(C)
    
    marker_distrib <- tidyr::separate_rows(marker_distrib,"CT",sep="\\|")    
    marker_distrib = marker_distrib[marker_distrib$gene %in% rownames(C),] 
    MD <- tapply(marker_distrib$gene,marker_distrib$CT,list)
    MD <- lapply(MD,function(x) sapply(x, function(y) which(y==rownames(C))))
    
    RESULTS = t(dtangle::dtangle(Y=mixture_samples, reference=reference_samples, markers = MD)$estimates)
    
    ###################################
    ###################################
    
  } else if (method == "MuSiC"){
    
    require(MuSiC)
    RESULTS = t(MuSiC::music_prop(bulk.eset = T.eset, sc.eset = C.eset, clusters = 'cellType',
                                  markers = NULL, normalize = FALSE, samples = 'cellID', 
                                  verbose = F)$Est.prop.weighted)
    
  } else if (method == "DWLS"){
    
    require(DWLS)
    path=paste(getwd(),"/results_",STRING,sep="")
    
    if(! dir.exists(path)){ #to avoid repeating marker_selection step when removing cell types; Sig.RData automatically created
      
      dir.create(path)
      Signature <- DWLS::buildSignatureMatrixMAST(scdata = C, id = as.character(phenoDataC$cellType), path = path, diff.cutoff = 0.5, pval.cutoff = 0.01)
      
    } else {#re-load signature and remove CT column + its correspondent markers
      
      load(paste(path,"Sig.RData",sep="/"))
      Signature <- Sig
      
      if(!is.null(elem)){#to be able to deal with full C and with removed CT
        
        Signature = Signature[,!colnames(Signature) %in% elem]
        CT_to_read <- dir(path) %>% grep(paste(elem,".*RData",sep=""),.,value=TRUE)
        load(paste(path,CT_to_read,sep="/"))
        
        Signature <- Signature[!rownames(Signature) %in% cluster_lrTest.table$Gene,]
        
      }
      
    }
    
    RESULTS <- apply(T,2, function(x){
      b = setNames(x, rownames(T))
      tr <- DWLS::trimData(Signature, b)
      RES <- t(DWLS::solveDampenedWLS(tr$sig, tr$bulk))
    })
    
    rownames(RESULTS) <- as.character(unique(phenoDataC$cellType))
    RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
    RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
    
  } else if (method == "BisqueRNA"){#By default, BisqueRNA uses all genes for decomposition. However, you may supply a list of genes (such as marker genes) to be used with the markers parameter
    
    require(BisqueRNA)
    RESULTS <- BisqueRNA::ReferenceBasedDecomposition(T.eset, C.eset, markers=NULL, use.overlap=FALSE)$bulk.props #use.overlap is when there's both bulk and scRNA-seq for the same set of samples
    
  } else if (method == "deconvSeq"){
    
    singlecelldata = C.eset 
    celltypes.sc = as.character(phenoDataC$cellType) #To avoid "Design matrix not of full rank" when removing 1 CT 
    tissuedata = T.eset 
    
    design.singlecell = model.matrix(~ -1 + as.factor(celltypes.sc))
    colnames(design.singlecell) = levels(as.factor(celltypes.sc))
    rownames(design.singlecell) = colnames(singlecelldata)
    
    dge.singlecell = deconvSeq::getdge(singlecelldata,design.singlecell, ncpm.min = 1, nsamp.min = 4, method = "bin.loess")
    b0.singlecell = deconvSeq::getb0.rnaseq(dge.singlecell, design.singlecell, ncpm.min =1, nsamp.min = 4)
    dge.tissue = deconvSeq::getdge(tissuedata, NULL, ncpm.min = 1, nsamp.min = 4, method = "bin.loess")
    
    RESULTS = t(deconvSeq::getx1.rnaseq(NB0 = "top_fdr",b0.singlecell, dge.tissue)$x1) #genes with adjusted p-values <0.05 after FDR correction
    
  } else if (method == "SCDC"){ ##Proportion estimation with traditional deconvolution + >1 subject
    
    require(SCDC)
    RESULTS <- t(SCDC::SCDC_prop(bulk.eset = T.eset, sc.eset = C.eset, ct.varname = "cellType", sample = "SubjectName", ct.sub = unique(as.character(phenoDataC$cellType)), iter.max = 200)$prop.est.mvw)
    
  }
  
  return(RESULTS) 
  
}
