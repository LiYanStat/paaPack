####### internal functions for hPAA

## indentify the location of the measure of ith and jth species (i < j)
ioffM <- function(p, i, j) {
  (i - 1) * p + j - i * (i + 1) / 2
}

## the compute the similarity measures, return a vector with dimension p * (p - 1) / 2
## where the measure of ith and jth species (i < j) are stored in spot (i - 1) * p + j - i * (i + 1)
divM <- function(X, W, method = c("Simpson", "Shannon", "Bray-Curtis", "Weighted UniFrac")) {
  method <- match.arg(method)
  output <- matrix(0, nrow = ncol(X), ncol = ncol(X))
  if(method == "Simpson") {
    output <- t(X) %*% W %*% X
  } else if(method == "Shannon") {
    for(j in 1:(ncol(X) - 1)) {
      for(k in (j+1):ncol(X)) {
        output[j, k] <- output[k, j] <- 
          crossprod(X[, j] + X[, k], W %*% log(X[, j] + X[, k])) - crossprod(X[, j], W %*% log(X[, j])) - crossprod(X[, k], W %*% log(X[, k]))
      }
    }
  } else if(method == "Bray-Curtis") {
    for(j in 1:(ncol(X) - 1)) {
      for(k in (j+1):ncol(X)) {
        output[j, k] <- output[k, j] <- BCdistcpp(X[, j], X[, k])
      }
    }
  }
  output[lower.tri(output)]
}

SHdist <- function(X, Y, W) {
  crossprod(X + Y, W %*% log(X + Y)) - crossprod(X, W %*% log(X)) - crossprod(Y, W %*% log(Y))
}

####### internal functions for plots
max_labels_height <- utils::getFromNamespace("max_labels_height", "dendextend")
rescale <- utils::getFromNamespace("rescale", "dendextend")
rotated_str_dim <- utils::getFromNamespace("rotated_str_dim", "dendextend")


colored_bars_modified <- 
  function (colors, dend, rowLabels = NULL, cex.rowLabels = 0.9, 
            add = TRUE, y_scale, y_shift, logy = FALSE, text_shift = 1, sort_by_labels_order = TRUE, 
            horiz = FALSE, ...) {
  n_colors <- if (is.null(dim(colors))) 
    length(colors)
  else nrow(colors)
  n_groups <- if (is.null(dim(colors))) 
    1
  else ncol(colors)
  if (!missing(dend)) {
    if (is.hclust(dend)) 
      dend <- as.dendrogram(dend)
    if (!is.dendrogram(dend)) 
      stop("'dend' should be a dendrogram.")
    dend_labels <- labels(dend)
    dend_order <- order.dendrogram(dend)
  }
  else {
    dend_labels <- rep("W", n_colors)
    dend_order <- seq_len(n_colors)
  }
  if (!sort_by_labels_order) 
    dend_order <- seq_len(n_colors)
  if (!horiz) {
    if (missing(y_shift)) 
      y_shift <- -max_labels_height(dend_labels) + par("usr")[3L] - 
        strheight("X")
    if (missing(y_scale)) 
      y_scale <- strheight("X") * n_groups
  }
  else {
    if (missing(y_shift)) 
      y_shift <- -(min(strwidth(dend_labels)) + par("usr")[2L] + 
                     strwidth("X"))
    if (missing(y_scale)) 
      y_scale <- strwidth("X") * n_groups
  }
  y_shift <- y_shift - y_scale
  colors <- as.matrix(colors)
  dimC <- dim(colors)
  if (is.null(rowLabels) & (length(dimnames(colors)[[2]]) == 
                            dimC[2])) 
    rowLabels <- names(as.data.frame(colors))
  op <- options()
  pr <- par(no.readonly = TRUE)
  options(stringsAsFactors = FALSE)
  par(xpd = TRUE)
  if (length(dend_order) != dimC[1]) {
    stop("ERROR: length of colors vector not compatible with number of objects in the hierarchical tree.")
  }
  C <- colors[dend_order, ]
  C <- as.matrix(C)
  step <- 1/(n_colors - 1)
  ystep <- 1/n_groups
  if (!add) {
    barplot(height = 1, col = "white", border = FALSE, space = 0, 
            axes = FALSE, ...)
  }
  charWidth <- strwidth("W")/2
  charHeight <- strheight("W")/2
  if(logy) {
    if(y_shift <=0) {
      stop("invalid y_shift value")
    } else {
      ysteps <- (exp(seq(log(y_shift), log(y_shift + y_scale), length.out = n_groups + 1)) - y_shift) / y_scale
    }
  } else {
    ysteps <- 0:n_groups * ystep
  }
  for (j in 1:n_groups) {
    ind <- (1:n_colors)
    xl <- (ind - 1.5) * step
    xr <- (ind - 0.5) * step
    yb <- rep(ysteps[j], n_colors)
    yt <- rep(ysteps[j + 1], n_colors)
    # yb <- rep(ystep * (j - 1), n_colors)
    # yt <- rep(ystep * j, n_colors)
    if (add) {
      xl <- rescale(xl, to = c(1 - 0.5, n_colors - 0.5))
      xr <- rescale(xl, to = c(1 + 0.5, n_colors + 0.5))
      yb <- yb * y_scale + y_shift
      yt <- yt * y_scale + y_shift
    }
    if (horiz) {
      rect(-yb, xl, -yt, xr, col = as.character(C[, j]), 
           # border = as.character(C[, j]))
           border = NA)
      par(srt = 90)
      if (is.null(rowLabels)) {
        s <- as.character(j)
        text(s, pos = 1, offset = 0.5, y = charHeight * 
               text_shift - rotated_str_dim(s)[2]/2, x = -(ystep * 
                                                             (j) * y_scale + y_shift), cex = cex.rowLabels)
      }
      else {
        s <- rowLabels[j]
        text(s, pos = 1, offset = 0.5, y = charHeight * 
               text_shift - rotated_str_dim(s)[2]/2, x = -(ystep * 
                                                             (j) * y_scale + y_shift), cex = cex.rowLabels)
      }
    }
    else {
      rect(xl, yb, xr, yt, col = as.character(C[, j]), 
           # border = as.character(C[, j]))
           border = NA)
      if (is.null(rowLabels)) {
        text(as.character(j), pos = 2, x = charWidth * 
               text_shift, y = (ysteps[j] +  ysteps[j + 1]) / 2 * y_scale + 
               y_shift, cex = cex.rowLabels)
      }
      else {
        text(rowLabels[j], pos = 2, x = charWidth * text_shift, 
             y = (ysteps[j] +  ysteps[j + 1]) / 2 * y_scale + y_shift, 
             cex = cex.rowLabels)
      }
    }
  }
  # for (j in 0:n_groups) {
  #   the_x <- dendextend:::rescale(c(0, 1), to = c(1 - 0.5, n_colors +
  #                                      0.5))
  #   if (horiz) {
  #     lines(y = the_x, x = -(c(ystep * j, ystep * j) *
  #                              y_scale + y_shift), col = "white", lwd = 0.8)
  #   }
  #   else {
  #     lines(x = the_x, y = c(ystep * j, ystep * j) * y_scale +
  #             y_shift, col = "white", lwd = 0.8)
  #   }
  # }
  options(op)
  # par(pr)
  return(invisible(C))
  }

######### internal functions for the plot Lines
## taxon.matrix: the taxonomy matrix use in the tree guided amalgamation
pComp <- function(x) {
  compDat <- x$compDat
  method <- x$diversity.method
  set.seed(123)
  jit <- c(head(sort(jitter(rep(0, length(x$height)), factor = 1)), -1),rep(0, 1))
  x$height <- as.double(x$height) - jit
  if(method == "Simpson") {
    x$height <- (1 - x$height / (1 - sum(diag(t(compDat) %*% compDat)) / nrow(compDat)  - jit[1]*1.1)) * 100
  } else if(method == "Shannon") {
    x$height <- (1 - x$height / (- sum(diag(t(compDat) %*% log(compDat))) / nrow(compDat)  - jit[1]*1.1)) * 100
  } else if(method == "Bray-Curtis") {
    x$height <- (1 - x$height / (x$height[1]  - jit[1]*1.1)) * 100
  }
  return(x)
}


######### internal functions for the plot MDS

## logy: indicate whether set y axis to log scale
## taxon.matrix: the taxonomy matrix use in the tree guided amalgamation
cutoffComp <- function(x, compDat, cutoff, num.entities = FALSE, ...) {
  method <- x$diversity.method
  set.seed(123)
  jit <- c(head(sort(jitter(rep(0, length(x$height)), factor = 1)), -1),rep(0, 1))
  # jit <- rep(0, length(x$height))
  x$height <- as.double(x$height) - jit
  if(method == "Simpson") {
    x$height <- (1 - x$height / (1 - sum(diag(t(compDat) %*% compDat)) / nrow(compDat)  - jit[1]*1.1)) * 100
  } else if(method == "Shannon") {
    x$height <- (1 - x$height / (- sum(diag(t(compDat) %*% log(compDat))) / nrow(compDat)  - jit[1]*1.1)) * 100
  } else if(method == "Bray-Curtis") {
    x$height <- (1 - x$height / (x$height[1]  - jit[1]*1.1)) * 100
  }
  if(num.entities) {
    cut.ind <- 1:(nrow(x$merge) + 1 - cutoff)
  } else {
    cut.ind <- which(x$height < cutoff)
  }
  merge <- x$imerge[cut.ind, , drop = FALSE]
  output <- compDat
  for(i in 1:nrow(merge)) {
    output[, merge[i, 1]] <- output[, merge[i, 1]] + output[, merge[i, 2]]
  }
  output <- output[, - merge[, 2]]
  output
}