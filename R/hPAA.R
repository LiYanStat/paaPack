##' Hierarchical Principal Amalgamation Analysis with/without Taxonomic Tree Guidance
##'
##' Conduct hierarchical amalgamation using diversity indexes under the
##' guidance of taxonomic tree (optional) and produce a set of objects which describe
##' the hierarchical tree for pairwise amalgamation of components of the compositional 
##' matrix until there is just a single component.
##'
##' @usage
##' hPAA(X, taxonomy = NULL, 
##'      method = c("Simpson", "Shannon", "Bray-Curtis", "Weighted UniFrac"), 
##'      W = NULL, strong = FALSE)
##'
##' @param X a data matrix of compositions (\eqn{n \times p}).
##' @param taxonomy a vector of taxonomic tree structure related to the compositions
##'        which is used to conduct tree-guided hierarchical analysis.
##'        Each element of the vector denotes the taxonomic ranks of the components, 
##'        separated by semicolon, in the format of mothur taxonomy output (\eqn{p \times 1}).
##' @param method the diversity index used to conduct the hierarchical amalgamation.
##' @param W a vector of weight for each sample in measure of diversity (\eqn{n \times 1}).
##' @param strong indicator of tree-guided hierarchical approach, default is weak hierarchy taxonomy
##'        method. The parameter is valid only if taxonomic tree matrix is provided.
##' @return
##' An object of class "hPAA" which describes the tree produced by the amalgamation process
##' \item{merge}{a \eqn{p-1} by 2 matrix. Row \eqn{i} of merge describes the merging of components at step \eqn{i} 
##'              of the amalgamation. If an element \eqn{j} in the row is negative, then observation \eqn{-j} was merged 
##'              at this stage. If \eqn{j} is positive then the merge was with the new composition formed at the (earlier) stage 
##'              \eqn{j} of the algorithm. Thus negative entries in merge indicate agglomerations of singletons, and positive 
##'              entries indicate agglomerations of non-singletons.}
##' \item{imerge}{a \eqn{p-1} by 2 matrix. Row \eqn{i} of merge describes the merging of compositions at step \eqn{i} 
##'               of the amalgamation, as the indexes of columns in the original compositional data matrix. }
##' \item{height}{a set of \eqn{p-1} real values. The value of the diversity measure at each amalgamation step.}
##' \item{crit}{a set of \eqn{p-1} real values. The decrease of diversity meaure after amalgamation at each step.}
##' \item{order}{a vector giving the permutation of the original compositions suitable for plotting, 
##'              in the sense that a hierarchical dendrogram using this ordering and matrix merge will not have 
##'              crossings of the branches.}
##' \item{labels}{labels for each of the components being amalgamated.}
##' \item{diversity.method}{the measure of diversity used for conducting the analysis.}
##' \item{hierarchical.method}{the approach used for conducting the analysis. The value can be 
##'                            "unconstrained", "weak" and "strong".}
##' \item{lines.Ind}{a set of numbers associated with the steps. The numbers are meaningful only for strong taxonomic hierarchy 
##'                  where the numbers indicate at which step the depth of taxonomic tree decreases by 1.}
##' \item{compDat}{The compositional data for performing the PAA.}
##' \item{taxonomy}{The corresponding taxonomic tree structure of the compositional data if available.}
##' 
##' @examples
##' data(HIV)
##' ## get the compositional data
##' compDat <- HIV$compDat[, 1:60]
##' ## taxonomic tree structure
##' taxonomy <- HIV$taxonomy
##' ## weak taxonomic hierarchy
##' weak <- hPAA(compDat, taxonomy = taxonomy, method = "Simpson", strong = FALSE)
##' ## strong taxonomic hierarchy
##' strong <- hPAA(compDat, taxonomy = taxonomy, method = "Simpson", strong = TRUE)
##' @export
hPAA <- function(X, taxonomy = NULL, method = c("Simpson", "Shannon", "Bray-Curtis", "Weighted UniFrac"),
                 W = NULL, strong = FALSE) {
  method <- match.arg(method)
  ## if taxonomic structure is not provided, the strong param is not used
  if(is.null(taxonomy)) {
    strong <- FALSE
  }
  X <- as.matrix(X)
  ## Weighted UniFrac is not considered in the current version
  if (method == "Weighted UniFrac") {
    stop("Weighted UniFrac to be added")
  }
  if(! (any(X < 0) | max(abs(rowSums(X) - 1)) < 1e-8)) {
    stop("Invalid compositional data matrix")
  }
  ## check if the taxonomy matrix is valid
  if(!is.null(taxonomy)) {
    if(is.matrix(taxonomy)) {
      if(!all(rownames(taxonomy) == colnames(X))) {
        stop("Missing taxonomy structure for species")
      }
    } else {
      stop("Invalid taxonomy matrix")
    }
    ## taxonomy structure, assign levels to each of the structure for
    ## constructing the lines for taxonomy level splitting in dendrogram
    ## Domain, Kingdom, Phylum, Class, Order, Family, Genus
    tax.level <- c("D", "K", "P", "C", "O", "F", "G")
    ## split the taxonomy matrix into list
    taxonomy.1 <- taxonomy
    taxonomy.1[, 1] <- paste("root", taxonomy, sep = ";")
    taxonomy.list <- strsplit(taxonomy.1, split = ";", fixed = TRUE)
    ## get the taxonomy level labels for each of the taxa
    tax.mat <- lapply(taxonomy.list,
                      function(x) {
                        tax.level[seq_len(length(x))]
                      })
    taxonomy.list <-
      lapply(1:length(taxonomy.list),
             function(x) {
               c(taxonomy.list[[x]], rownames(taxonomy.1)[x])
             })
    dtmp <- max(sapply(taxonomy.list, length)) - 1
    ## clean the taxonomy structure
    # taxonomy.list
    rmInd <- matrix(NA, nrow = length(taxonomy.list), ncol = dtmp)
    for(k in 1:dtmp) {
      ## indicate whether remove the root
      tmp <- t(sapply(taxonomy.list, function(x) {x[1:(k + 1)]}))
      tmp.u <- unique(tmp)[apply(unique(tmp), 1, function(x) all(!is.na(x))), ,
                           drop = FALSE]
      Ind <- which(apply(tmp[, 1:k, drop = FALSE], 1, paste, collapse = " ") %in%
                     names(which(table(
                       apply(tmp.u[, 1:k, drop = FALSE], 1, paste, collapse = " ")
                     ) == 1)))
      rmInd[Ind, k] <- Ind
    }
    taxonomy.list <-
      lapply(1:length(taxonomy.list),
             function(i) {
               output <- head(taxonomy.list[[i]], -1)
               output <- output[is.na(rmInd[i, ])]
               output[!is.na(output)]
             })
    tax.mat <-
      lapply(1:length(taxonomy.list),
             function(i) {
               output <- tax.mat[[i]]
               output <- output[is.na(rmInd[i, ])]
               output[!is.na(output)]
             })
    ## name the taxonomy.list
    names(taxonomy.list) <- names(tax.mat) <- rownames(taxonomy.1)
    depth.list <- tax.level[tax.level %in% unique(unlist(tax.mat))]
    depth <- tail(tax.level[tax.level %in% unique(unlist(tax.mat))], 1)
    ## depth <- max(sapply(taxonomy.list, length))
  }
  n <- dim(X)[1]  ## number of individuals
  p <- dim(X)[2]  ## number of Species
  label <- colnames(X)  ## labels of species
  ## weight matrix
  if(is.null(W)) {
    W <- diag(n) / n
  } else {
    W <- diag(W)
  }
  ## ia, ib for merge as function stats::hclust and crit stores the decrease of diversity in each step
  ia <- ib <- crit <- rep(0, length = p - 1)
  flag <- rep(1, length = p)  ## indicate whether we should include the species in combination
  #### the initial distance matrix
  ## for method Simpson, we can directly perform the reduction with the dissimilarity matrix
  # if(method == "Simpson") {
  diss <- divM(X, W, method = method)  ## the dissimilarity measure in terms of decrease in diversity
  
  disnn <- rep(Inf, length = p - 1)  ## store the dissimilarity for ith entity i = 1, 2, ..., 9
  nn <- rep(0, length = p - 1)  ## store the dissimilarity for ith entity i = 1, 2, ..., 9
  ## Carry out the first agglomeration
  for (i in 1:(p - 1)) {
    dmin <- Inf  ## initialize the minimum distance
    jm <- 0  ## initalize the jth entity such that for fixed i, the dissimiarity of (i, j) is minimized
    for (j in (i+1):p) {
      ## keep the taxonomy structure
      if(!is.null(taxonomy)) {
        if(! identical(taxonomy.list[[label[i]]], taxonomy.list[[label[j]]])) {
          # print(paste("(", i, j, "Not Count )"))
          next
        }
      }
      # print(paste("(", i, j, "Count )"))
      ind <- ioffM(p,i, j)
      if(diss[ind] < dmin) {
        dmin <- diss[ind]
        jm <- j
      }
    }
    disnn[i] <- dmin
    nn[i] <- jm
  }
  
  ## Repeat previous steps until p-1 agglomerations carried out.
  if(!is.null(taxonomy)) {
    lines.Ind <- depth
  }
  num.entity <- p
  Xnew <- X  ## initial Xnew
  while(num.entity > 1) {
    # print(taxonomy.list[as.logical(flag)])
    # print("+++++++++++++++++")
    ###### Next, determine least diss. using list of NNs
    dmin <- Inf
    for(i in 1:(p - 1)) {
      if(flag[i]) {
        if(strong) {
          ## handling weak hierarchical amalgamation
          # if(! length(taxonomy.list[[label[i]]]) == depth) {
          if(! tail(tax.mat[[label[i]]], 1) == depth) {
            # print(taxonomy.list[label[i]])
            next
          }
        }
        # print(i)
        if (disnn[i] < dmin) {
          dmin <-disnn[i]
          im <- i
          jm <- nn[i]
        }
      }
    }
    ## the number of entity decreases by 1
    num.entity <- num.entity - 1
    ###### This allows an agglomeration to be carried out.
    ###### At step n-ncl, we found dmin=d[i2,j2]
    i2 <- im; j2 <- jm
    ia[p - num.entity] <- i2; ib[p - num.entity] <- j2
    crit[p - num.entity] <- dmin
    ###### get the X matrix after combination
    Xnew[, i2] <- Xnew[, i2] + Xnew[, j2]
    ###### Update dissimilarities from new entity, 
    ###### we merge the j2 into i2, recalculate the dissimilarity and 
    ###### add to the diss of i2 and kth (k = 1, 2, ..., p)
    flag[j2] <- 0
    ## also update the taxonomy list according to j2 and i2
    if(!is.null(taxonomy)) {
      Upgrade <- sum(sapply(taxonomy.list[as.logical(flag)], function(x) {
        if(length(x) >= length(taxonomy.list[[label[i2]]])) {
          identical(x[1:length(taxonomy.list[[label[i2]]])], taxonomy.list[[label[i2]]])
        } else {
          return(FALSE)
        }
      })) == 1
      if(Upgrade) {
        taxonomy.list[[label[i2]]] <- head(taxonomy.list[[label[i2]]], -1)
        taxonomy.list[[label[j2]]] <- head(taxonomy.list[[label[j2]]], -1)
        tax.mat[[label[i2]]] <- head(tax.mat[[label[i2]]], -1)
        tax.mat[[label[j2]]] <- head(tax.mat[[label[j2]]], -1)
      }
      ## update the depth
      depth <- tail(tax.level[tax.level %in% unique(unlist(tax.mat[as.logical(flag)]))], 1)
      lines.Ind <- c(lines.Ind, depth)
      # depth <- max(sapply(taxonomy.list[as.logical(flag)], length))
    }
    dmin <- Inf
    jj <- 0
    ## update disnn and nn
    # print(diss)
    for(k in 1:p) {
      if(flag[k] & (k != i2)) {
        ind1 <- ioffM(p, min(i2, k), max(i2, k))
        ind2 <- ioffM(p, min(j2, k), max(j2, k))
        # print(paste(min(j2, k), max(j2, k), "/", min(i2, k), max(i2, k)))
        ## the dissimilarity of new entity to kth entity is just dis[i2, k] + dis[j2, k]
        ############ update the distance based on different method
        if(method == "Simpson") {
          diss[ind1] <- diss[ind1] + diss[ind2]  ## Simpson
        } else if(method == "Shannon") {
          diss[ind1] <- crossprod(Xnew[, i2] + Xnew[, k], W %*% log(Xnew[, i2] + Xnew[, k])) - 
            crossprod(Xnew[, i2], W %*% log(Xnew[, i2])) - crossprod(Xnew[, k], W %*% log(Xnew[, k]))
          # diss[ind1] <- SHdist(Xnew[, i2], Xnew[, k], W)
        } else if(method == "Bray-Curtis") {
          diss[ind1] <- BCdistcpp(Xnew[, i2], Xnew[, k])
        }
        # diss[ind1] <- diss[ind1] + diss[ind2]  ## Shannon
        ## continue here
        
        ## update disnn and nn
        if(!is.null(taxonomy)) { ## keep the taxonomy structure
          if(! identical(taxonomy.list[[label[i2]]], taxonomy.list[[label[k]]])) {
            # print(paste("(", i2, jk, "Not Count )"))
            next
          }
        }
        if((i2 < k) && (diss[ind1] < dmin)) {  ## update for the i2 row
          # print(k)
          dmin <- diss[ind1];
          jj <- k
        }
        if( (i2 > k) & (nn[k] != j2) & (diss[ind1] < disnn[k]) ) { ## update for other row
          disnn[k] <- diss[ind1]
          nn[k] <- i2
        }
      }
    } ## end for loop k
    disnn[i2] <- dmin
    nn[i2] <-jj
    
    ##### Update list of NNs insofar as this
    for (i in 1:(p-1)) {
      if( flag[i] & ((nn[i] == i2) | (nn[i] == j2))) {
        dmin <- Inf
        for ( j in (i+1):p) {
          if(!is.null(taxonomy)) { ## keep the taxonomy structure
            if(! identical(taxonomy.list[[label[i]]], taxonomy.list[[label[j]]])) {
              # print(paste("(", i, j, "Not Count )"))
              next
            }
          }
          ind <- ioffM(p, i, j)
          if(flag[j] & (i != j) & (diss[ind] < dmin)) {
            dmin <- diss[ind]
            jj <- j
          }
          nn[i] <- jj
          disnn[i] <- dmin
        }
      }
    }
    # for(d in depth:1) {
    #   if(all(table(sapply(taxonomy.list[as.logical(flag)], function(x) x[d])) == 1)) {
    #     taxonomy.list[as.logical(flag)] <- lapply(taxonomy.list[flag == 1], function(x) head(x, -1))
    #   } else {
    #     break
    #   }
    # }
    # depth <- max(sapply(taxonomy.list[as.logical(flag)], length))
    # print(num.entity)
  } ## end for while
  ##### prepare the merged information for plot dendrogram
  iia <- iib <- rep(0, length = p - 1)
  for (i in 1:(p-1)) {
    iia[i] <- - ia[i]
    iib[i] <- - ib[i]
  }
  for (i in 1:(p-2)) {
    k <- min(ia[i], ib[i])
    for (j in (i+1):(p-1)) {
      if(ia[j] == k) iia[j] <- i
      if(ib[j] == k) iib[j] <- i
    }
  }
  for (i in 1:(p-1)) {
    if((iia[i] > 0) & (iib[i] < 0)) {
      k <- iia[i];
      iia[i] <- iib[i]
      iib[i] <- k
    }
    if((iia[i] > 0) & (iib[i] > 0)){
      k1 <- min(iia[i], iib[i])
      k2 <- max(iia[i], iib[i])
      iia[i] <- k1
      iib[i] <- k2
    }
  }
  ## sort the order for dendrogram
  iorder <- rep(0, length = p)
  iorder[1] <- iia[p-1]
  iorder[2] <- iib[p-1]
  loc <- 2
  for (i in (p-2):1) {
    for (j in 1:loc) {
      if (iorder[j] == i) {
        iorder[j] <- iia[i]
        if (j == loc) {
          loc <- loc + 1
          iorder[loc] <- iib[i]
        } else {
          loc <- loc + 1
          for (k in loc:(j+2)){
            iorder[k] <- iorder[k-1]
          }
          iorder[j+1] <- iib[i]
        }
        break
      }
    }
  }
  iorder <- - iorder
  # }  ## end for Simpson method
  if(!is.null(taxonomy)) {
    if(strong) {
      h.method <- "Strong Hierarchical"
    } else {
      h.method <- "Weak Hierarchical"
    }
    ## get lines.Ind
    lines.Ind <- cumsum(rev(table(lines.Ind)[depth.list]))
  } else {
    h.method <- "Unsupervised Hierarchical"
    lines.Ind <- NULL
  }
  
  ## get the height for different method
  if(method == "Simpson") {
    height <- (1 - sum(diag(t(X) %*% W %*% X))) - cumsum(crit * 2)
  } else if(method == "Shannon") {
    height <- - sum(diag(t(X) %*% W %*% log(X))) - cumsum(crit)
  } else if(method == "Bray-Curtis") {
    height <- sum(crit) - cumsum(crit)
  }
  ## sort all the results
  tree <- list(merge = cbind(iia, iib),
               imerge = cbind(ia, ib),
               height = height,
               crit = crit,
               order = iorder,
               labels = label,
               diversity.method = method,
               hierarchical.method = h.method,
               lines.Ind = lines.Ind,
               compDat = X, 
               taxonomy = taxonomy
  )
  class(tree) <- "hPAA"
  tree
}