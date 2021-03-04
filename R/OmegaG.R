#' @name OmegaG
#' @title Coefficient Omega-Generic
#' @author Yujiao Mai, Deo Kumar Srivastava, and Kevin R. Krull
#' @description
#' This function is used to estimate the composite reliability coefficient Omega-generic (\cite{Mai, Srivastava, & Krull, 2021}), given factor loadings, factor covariance matrix, and covariance matrix of item residuals.
#' @param Lambda The input factor lading matrix. Each row contains the loadings of one item on factors. Each column includes the loadings of one factor. In the case of bi-factor structure, the first column of loadings is on the global factor.
#' @param Phi The input factor covariance matrix.
#' @param Psi The input covariance matrix of item residuals. Typically, \code{Psi} is a diagonal matrix.
#' @param items.index The vector indexing the items of which the composite reliability is being estimated. It is an optional argument. If it is specified, the argument \code{scale.structure} is not effective. If it is not specified, the argument \code{scale.structure} is expected to be specified and effective.
#' @param factor.index The vector indexing the factor(s)/construct(s) regarding which the composite reliability is being estimated. It is an optional argument. If it is not specified, the function will estimate the composite reliability regarding each factor/construct.
#' @param scale.structure The scale structure in a list or a Boolean matrix form. In a list form, each element is a vector of items (names) of a subscale. If in a Boolean form, the element on the i-th row and the j-th column indicates whether the i-th item is within the j-th subscale. If both the argument \code{scale.structure} and \code{Lambda} include colnames and rownames, the names must match with each other. This argument \code{scale.structure} is optional. This argument is effective only when \code{item-index} is not specified, 
#' @param modeltype The type of factor structure (\code{"correlated-factor"} or \code{"bi-factor"}). The current version supports correlated-factor structure and bi-factor structure. A bi-factor model typically assumes factors are uncorrelated. The default is correlated-factor structure.
#' @return The estimated composite reliability coefficient OmegaG. It is a list of CR estimates with headers. 
#' @references Mai, Y., Srivastava, D.K., & Krull, K.R. (2021). Estimating Composite reliability of Multidimensional Measurement with Overlapping Items. Present at the 2021 Eastern North American Region (ENAR) Spring Virtual Meeting.
#' @examples
#' \dontshow{
#'  #Toy example for check:
#'  Lambda <- cbind(c(0.7,0.65, 0.5, 0.0, 0.05, 0.15),
#'                  c(0.10, 0.05, 0.21, 0.72, 0.70, 0.65),
#'                  c(0.5, 0.6, 0.71, 0.65, 0.4, 0.5))
#'  Phi <- diag(rep(1,3))
#'  Psi <- diag(c(0.2500, 0.2150, 0.2018, 0.0591, 0.3475, 0.3050))
#'  varname <- paste('y',1:6,sep=''); factorname <- c(paste('fs',1:2,sep=''),'fg')
#'  colnames(Lambda) <- factorname; rownames(Lambda)<-varname;
#'  colnames(Phi)<-rownames(Phi)<-factorname; colnames(Psi)<-rownames(Psi)<-varname
#'
#'  OmegaG(Lambda, Phi, Psi,
#'       scale.structure = cbind(c(rep(TRUE,3),rep(FALSE,3)),
#'                               c(rep(FALSE,3), rep(TRUE,3))),
#'       modeltype='bi-factor')
#' # Omegaesem.s1  Omegaesem.s2 Omegaesem.hs1 Omegaesem.hs2  Omegaesem.hg
#' # 0.8950369     0.8303052     0.4107981     0.3055959     0.6024197
#'
#'
#' OmegaG(Lambda, Phi, Psi,items.index=1:3,factor.index=1,
#'       scale.structure = cbind(c(rep(TRUE,3),rep(FALSE,3)),
#'                               c(rep(FALSE,3), rep(TRUE,3))),
#'       modeltype='bi-factor')
#' OmegaG(Lambda, Phi, Psi,items.index=1:3,factor.index=1)
#' }
#' \donttest{
#' #### Example 1:
#' OmegaG(Lambda = PedsQLMFS$ESEM$Lambda,
#'                        Phi = PedsQLMFS$ESEM$Phi,
#'                         Psi = PedsQLMFS$ESEM$Psi,
#'                         modeltype = "correlated-factor",
#'                         scale.structure = PedsQLMFS$ScaleStructure
#'                         )
#'
#' #  Model type = correlated-factor
#' #
#' #  CR of each subscale:
#' #       GeneralFatigue :    0.770
#' #         SleepFatigue :    0.690
#' #     CognitiveFatigue :    0.777
#'
#'
#' #### Example 2:
#'  OmegaG(Lambda = PedsQLMFS$biESEM$Lambda,
#'                 Phi = PedsQLMFS$biESEM$Phi,
#'                 Psi = PedsQLMFS$biESEM$Psi,
#'                 modeltype = "bi-factor",
#'                 scale.structure = PedsQLMFS$ScaleStructure
#'  )
#'
#' # Model type = bi-factor
#' #
#' # Hierarchy and Hierarchical-subscale CR:
#' #                          GlobalFatigue :    0.806
#' #                         GeneralFatigue :    0.174
#' #                           SleepFatigue :    0.361
#' #                       CognitiveFatigue :    0.190
#' #
#' # Scale Total and Subscale CR:
#' #   GlobalFatigue + all sepcific factors :    0.926
#' #         GlobalFatigue + GeneralFatigue :    0.859
#' #           GlobalFatigue + SleepFatigue :    0.758
#' #       GlobalFatigue + CognitiveFatigue :    0.839
#'
#'
#' # Example 3:
#'  OmegaG::OmegaG(Lambda = PedsQLMFS$biESEM$Lambda,
#'        Phi = PedsQLMFS$biESEM$Phi,
#'        Psi = PedsQLMFS$biESEM$Psi,
#'        modeltype = "bi-factor",
#'        items.index = 1:6,factor.index = 2
#'  )
#'
#' # Model type = bi-factor
#' #
#' # CR of Items 1 2 3 4 5 6 regarding factor 2:
#' #                      GeneralFatigue :    0.174
#'
#'
#' # Example 4:
#'   OmegaG::OmegaG(Lambda = PedsQLMFS$ESEM$Lambda,
#'                  Phi = PedsQLMFS$ESEM$Phi,
#'                   Psi = PedsQLMFS$ESEM$Psi,
#'                  modeltype = "correlated-factor",
#'                   items.index = 7:12,factor.index = 2
#'    )
#'
#' # Model type = correlated-factor
#' #
#' # CR of Items 7 8 9 10 11 12 regarding factor 2:
#' #   SleepFatigue :    0.690
#'
#' }
#'
##' @export


OmegaG <- function(Lambda=NULL,Phi=NULL,Psi=NULL,
                   items.index=NULL, factor.index=NULL, scale.structure=NULL,
                   modeltype=c('correlated-factor','bi-factor') ){
  CRlist <- NULL
  CRset <- -1
  if(missing(Lambda)) stop('Missing Lambda');
  if(missing(Phi)) stop('Missing Phi');
  if(missing(Psi)) stop('Missing Psi');
  if( (missing(items.index) & missing(factor.index)) & missing(scale.structure)) {CRset <- 0 # scale total
    } else if (!missing(items.index) & missing(factor.index)){ CRset <- 1 ##  the items regarding each factor and all common factor
    } else if (!missing(items.index) & !missing(factor.index)) { CRset <- 2 ##  the items regarding the factor
    } else if (!missing(scale.structure) & missing(factor.index)) {CRset <- 3  ## use the scale.structure to output all CR
    } else if (!missing(scale.structure) & !missing(factor.index) ) {CRset <- 4  ## use factor.index and the scale.structure
    }
  if (!missing(modeltype)) {modeltype<-match.arg(modeltype)} else {modeltype<-'correlated-factor'}
  if(is.matrix(Lambda)){}else if ( class(Lambda) == 'data.frame') Lambda <- as.matrix(Lambda)
  if(is.matrix(Phi)){}else if ( class(Phi) == 'data.frame') Phi <- as.matrix(Phi)
  if(is.matrix(Psi)){}else if ( class(Psi) == 'data.frame') Psi <- as.matrix(Psi)
  if(CRset %in% c(3,4)){## transform the scale.structure into a matrix
    if(is.matrix(scale.structure)){
      if(nrow(Lambda)==nrow(scale.structure)) {} else stop("The number of rows dosnt match for Lambda and scale.structure as input")
      if(!is.null(rownames(scale.structure)) & !is.null(rownames(Lambda))) scale.structure <- scale.structure[rownames(Lambda),drop=F]
      if(!is.null(colnames(scale.structure)) & !is.null(colnames(Lambda))){
        if(modeltype=='correlated-factor')  scale.structure <- scale.structure[, colnames(Lambda), drop=F]
        if(modeltype=='bi-factor') scale.structure <- scale.structure[, colnames(Lambda)[-1], drop=F]
      }
    } else{
      if(is.list(scale.structure)){
        L <- length(scale.structure)
        p <- nrow(Lambda)
        tmatrix <- matrix(FALSE,nrow=p,ncol=L)
        rownames(tmatrix)<-rownames(Lambda)
        colnames(tmatrix)<-names(scale.structure)
        for(i in 1:p){  for(j in 1:L){   tmatrix[scale.structure[[j]],j]<- TRUE  }  }
        scale.structure<-tmatrix
      }else {stop("the argument scale.structure must be a matrix or a list of vectors") }
    }
  }

  ####### CR function
  ## Omega ESEM: Model-based estimation of composite reliability
  Omega.g <- function(Lambda=NULL,Phi=NULL,Psi=NULL,items.index=NULL, factor.index=NULL){

    if(missing(Lambda)) stop('Missing Lambda');
    if(missing(Phi)) stop('Missing Phi');
    if(missing(Psi)) stop('Missing Psi');
    if(missing(items.index)) stop('Missing item.index');
    if(missing(factor.index)) stop('Missing factor.index');

    ## get the dim of matrix Lambda
    dim.Lambda = dim(Lambda)
    ## get the dim of matrix Sigma.f
    dim.Phi = dim(Phi)
    ## get the dim of matrix Sigma.delta
    dim.Psi = dim(Psi)
    if (dim.Lambda[2] == dim.Phi[2]) {} else {stop("dim Lambda doesn match dim Phi");}
    if (dim.Lambda[1] == dim.Psi[2]) {} else {stop('dim Lambda doesn match Psi');}

    ## get the number of items
    pp = length(items.index)
    ## get the number of factors
    L = dim.Phi[2]
    ## the index of the intended/targeted factor
    h = factor.index

    if (dim.Lambda[1] >= pp ) {} else {stop('dim Lambda row is less than the number of items');}
    ## To check if the items.index is valid / in the Lambda matrix
    if (max(items.index)<=dim.Lambda[1]){}else{stop('items.index is out of bound.')}
    ## To check if the factor.index is valid
    if (max(factor.index) <= dim.Phi[2]) {} else {stop('factor.index is out of bound.');}
    if (length(factor.index)==1) {} else {stop('factor.index contains more than one number.');}
    ## To check if duplicate item index
    if ( sum( duplicated(items.index))>=1) {stop('Duplicated item index input.')}

    ## The total variance of the sum score:
    diag.Phi = as.vector(diag(Phi))
    vec.one.p = rep(1,pp)
    vec.sum.Lambda = as.vector(vec.one.p%*%Lambda[items.index,])
    #First term
    t1 = sum( diag(vec.sum.Lambda)%*%diag(diag.Phi)%*%t(diag(vec.sum.Lambda)) )
    #Second term
    t2 = sum(  (vec.sum.Lambda %*% t(vec.sum.Lambda)-diag(vec.sum.Lambda)%*%t(diag(vec.sum.Lambda))) * Phi ) #Nov30 2020
    Sigma.delta = diag(Psi)[items.index]
    t3= sum( Sigma.delta )

    var.c = t1+t2+t3


    ## The variance of the sum score that reflects the variance of the intended factor h:
    vec.Phi.h = as.vector( Phi[1:L,h] )
    #First term
    t1 = sum( diag(vec.sum.Lambda) %*% diag(vec.Phi.h)^2 %*% t(diag(vec.sum.Lambda)) * 1/Phi[h,h] )
    #Second term
    sum.Lambda.h = sum(Lambda[items.index,h])
    vec.Phi.h = Phi[1:L,h]
    f.index = 1:L;  f.index = f.index[-h];
    t2 =  2*sum.Lambda.h * sum (( vec.sum.Lambda * vec.Phi.h )[f.index])

    var.t = t1+t2

    var.t/var.c

  }

  ## Omega ESEM hierarchial: Model-based estimation of composite reliability
  Omega.h = function(Lambda=NULL,Phi=NULL,Psi=NULL,items.index=NULL, factor.index=NULL,
                          model=c('bi-factor-Orthogonal','bi-factor-Oblique','correlated-factor') ){

    if(missing(Lambda)) stop('Missing Lambda');
    if(missing(Phi)) stop('Missing Phi');
    if(missing(Psi)) stop('Missing Psi');
    if(missing(items.index)) stop('Missing item.index');
    if(missing(factor.index)) stop('Missing factor.index');
    if(missing(model)) model='bi-factor-Orthogonal';

    ## get the dim of matrix Lambda
    dim.Lambda = dim(Lambda)

    ## get the dim of matrix Sigma.f
    dim.Phi = dim(Phi)

    ## get the dim of matrix Sigma.delta
    dim.Psi = dim(Psi)

    if (dim.Lambda[2] == dim.Phi[2]) {} else {stop("dim Lambda doesn match dim Phi");}
    if (dim.Lambda[1] == dim.Psi[2]) {} else {stop('dim Lambda doesn match Psi');}

    ## get the number of items
    pp = length(items.index)
    ## get the number of factors
    L = dim.Phi[2]
    ## the index of the targeted factor
    h = factor.index

    if (dim.Lambda[1] >= pp ) {} else {stop('dim Lambda row is less than the number of items');}
    ## To check if the items.index is valid / in the Lambda matrix
    if (max(items.index)<=dim.Lambda[1]){}else{stop('items.index is out of bound.')}
    ## To check if the factor.index is valid
    if (max(factor.index) <= dim.Phi[2]) {} else {stop('factor.index is out of bound.');}

    vec.one.p = rep(1,pp)
    ## The total variance of the sum score:  using matrix operation
    t12 <- t(vec.one.p)%*%Lambda[items.index,]%*%Phi%*%t(Lambda[items.index,])%*%vec.one.p
    Sigma.delta = diag(Psi)[items.index]
    t3= sum( Sigma.delta )
    var.c = t12+t3
    ## The variance of the sum score that reflects the variance of the targeted factor:
    if(model=='bi-factor-Orthogonal') {
      Phi.d<-diag(diag(Phi))
      var.t <- t(vec.one.p)%*%Lambda[items.index,factor.index]%*%Phi.d[factor.index,factor.index]%*%t(Lambda[items.index,factor.index])%*%vec.one.p
    } else if (model=='bi-factor-Oblique'){ stop("The current version dosent support Oblique bi-factor")
    } else if (model=='correlated-factor'){
      var.t <- t(vec.one.p)%*%Lambda[items.index,factor.index]%*%Phi[factor.index,factor.index]%*%t(Lambda[items.index,factor.index])%*%vec.one.p
    }

    var.t/var.c

  }
  CRlist <- list()
  if(CRset ==0){ #scale total
    L <- ncol(Lambda)
    p <- nrow(Lambda)
    if(modeltype=='bi-factor') {
      CRlist[[1]]<- list(header='Scale total',  OmegaG.hat = Omega.h(Lambda, Phi, Psi, items.index =1:p, factor.index=1:L))
    } else if (modeltype=='correlated-factor'){ }
  } else if (CRset==1){ ##  the items regarding each factor and all common factor
    L <- ncol(Lambda)
    factornames <- colnames(Lambda); if(is.null(factornames)) factornames <- paste("F",1:ncol(scale.structure), sep='')
    if(modeltype=='bi-factor'){
      CRlist[[1]]<- list(header=paste('Items',paste(items.index,collapse = ' '), 'regarding all common factor', sep =' '),
                     OmegaG.hat=Omega.h(Lambda, Phi, Psi, items.index = items.index, factor.index=1:L) )
      for (h in 1:L){
        CRlist[[h+1]]<- list(header=paste('Items',paste(items.index,collapse = ' '), 'regarding factor',factornames[h], sep =' '),
                       OmegaG.hat = Omega.h(Lambda, Phi, Psi, items.index = items.index, factor.index=h))
      }
    } else if(modeltype=='correlated-factor'){
      # CRlist[[1]]<- list(header=paste('Items',items.index, 'regarding all common factor', collapse =' '),
      #                    OmegaG.hat=Omega.h(Lambda, Phi, Psi, items.index = items.index, factor.index=1:L) )
      for (h in 1:L){
        CRlist[[h+1]]<- list(header=paste('Items',paste(items.index,collapse = ' '), 'regarding factor',factornames[h], sep =' '),
                             OmegaG.hat = Omega.h(Lambda, Phi, Psi, items.index = items.index, factor.index=h))
      }
    }
  } else if (CRset==2){ ##  the items regarding the factor
    CR <- NULL
    if(modeltype=='correlated-factor'){
      subscale.name <- colnames(Lambda); if (is.null(subscale.name))subscale.name <- paste("F",1:ncol(scale.structure), sep='')
      if(length(factor.index) >=2){
        CR <- sapply(factor.index, function(f){Omega.g(Lambda,Phi,Psi,items.index, f)})
        names(CR) <- subscale.name[factor.index]
      } else { CR<-Omega.g(Lambda,Phi,Psi,items.index, factor.index); names(CR)<-subscale.name[factor.index]}
    } else if (modeltype=='bi-factor'){
        subscale.name <- colnames(Lambda); if (is.null(subscale.name))subscale.name <- paste("F",1:ncol(scale.structure), sep='')
        CR <- Omega.h(Lambda,Phi,Psi,items.index, factor.index)
        names(CR)<-subscale.name[factor.index]
    }
    CRlist[[1]]<- list(header=paste('CR of Items',paste(items.index,collapse = ' '), 'regarding factor', factor.index, sep =' '),
                   OmegaG.hat=CR )

  } else if (CRset==3){ ## use the scale.structure to output all CR
      subscale.name <- colnames(scale.structure); if (is.null(subscale.name))subscale.name <- paste("F",1:ncol(scale.structure), sep='')
      L <- ncol(scale.structure)
      p <- nrow(scale.structure)

      if(modeltype=='correlated-factor') {
        CR <- NULL
          for(i in 1:L){
            items.index<- which(scale.structure[,i])
            CR <- c(CR,Omega.g(Lambda,Phi,Psi,items.index, i))
          }
        names(CR) <- subscale.name
        CRlist[[1]]<- list(header="CR of each subscale", OmegaG.hat=CR )
      } else if(modeltype=='bi-factor') {
          subscale.name <- colnames(scale.structure); if (is.null(subscale.name))subscale.name <- paste("F",1:ncol(scale.structure), sep='')
          CR <- NULL
          items.index <- 1:p
          OmegaG.h <- Omega.h(Lambda,Phi,Psi,items.index, 1)
          CR <- c(CR, OmegaG.h)
          for(i in 1:L){
            items.index<- which(scale.structure[,i])
            OmegaG.hs <- Omega.h(Lambda,Phi,Psi,items.index, i+1)
            CR <- c(CR, OmegaG.hs)
          }
          tname <- colnames(Lambda); if(is.null(tname)){fgname<-'GlobalFactor'}else{fgname <- tname[1]}
          names(CR)<- c(fgname, subscale.name)
          CRlist[[1]]<- list(header="Hierarchy and Hierarchical-subscale CR", OmegaG.hat=CR )

          CR<-NULL
          items.index <- 1:p
          OmegaG.t <- Omega.h(Lambda,Phi,Psi,items.index, 1:(L+1))
          CR <- c(CR, OmegaG.t)
          for(i in 1:L){
            items.index <- which(scale.structure[,i])
            OmegaG.s <- Omega.h(Lambda,Phi,Psi,items.index, c(1,i+1))
            CR <- c(CR, OmegaG.s)
          }
          tname <- colnames(Lambda); if(is.null(tname)){fgname<-'GlobalFactor'}else{fgname <- tname[1]}
          names(CR)<- c(paste(fgname, "+ all sepcific factors"), paste(fgname, subscale.name, sep=' + ') )
          CRlist[[2]]<- list(header="Scale Total and Subscale CR", OmegaG.hat=CR )

      }

  } else if (CRset==4){ ## use factor.index and the scale.structure


  }
  str <- paste("Model type = ", modeltype, sep='')
  for(i in 1:length(CRlist) ){
    x<-CRlist[[i]];
    str <- paste(str, '\n\n    ', x$header, sep='')
    tname <- names(x$OmegaG.hat);
    astr <- paste(sapply(seq_len(length(x$OmegaG.hat)), function(j){ sprintf('\t %50s :    %0.3f', tname[j], x$OmegaG.hat[j])    }), collapse = '\n')
    str <- paste(str,':\n', astr, sep='')
  }
  cat(str)
  invisible(CRlist)
}



