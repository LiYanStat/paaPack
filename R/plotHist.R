##' Plot Method for the Object with class "hPAA" from hPAA Analysis
##' 
##' Plot the results of hierarchical Principal Amalgamation Analysis (PAA) as 
##' dendrogram with the y-axis showing the y-axis shows the percentage decrease in total diversity as measured by 
##' the selected diversity index, optionally include color bars to indicate the taxonomic tree structure of the 
##' compositional data matrix
##'
##' @usage 
##' plotHPAA(object, taxonomy = TRUE, rank.line = FALSE, 
##' label = NULL, hang = -1, axes = TRUE, frame.plot = FALSE,
##' ann = TRUE, lab.cex = 0.8, main = "Dendrogram", sub = NULL, xlab = "",
##' ylab = "Percentage change (\%)", lim.y = -2, logy = TRUE,
##' jitters = FALSE, select.level = "all", ...)
##'
##' @param object the output with class "hPAA" from function hPAA().
##' @param taxonomy indicate if colored bars representing taxonomic tree structure is added to the dendrogram.
##' @param rank.line indicate if horizontal lines showing the taxonomic ranks is added to the dendrogram, 
##'                  only valid in the strong hierarchy.
##' @param label a character vector of labels for the leaves of the tree.
##'              By default the row names or row numbers of the original data are used.
##' @param hang the fraction of the plot height by which labels should hang below the rest of the plot. 
##'             A negative value will cause the labels to hang down from 0.
##' @param axes,frame.plot,ann logical flags as in plot.default.
##' @param lab.cex size of the lablels for the leaves of the dendrogram
##' @param main,sub,xlab,ylab character strings for title.
##' @param lim.y lower bound of y for the color plot
##' @param logy indicate whether plot the y-axis in log scale.
##' @param jitters indicate whether jitters need to be added to make the heights of dendrogram 
##'                more distinguising.
##' @param select.level indicate the levels of taxonomic ranks/depths to display on the dendrogram (as horizontal lines).
##' @param ... further graphical arguments for the colored_bars(). 
##'            E.g., y_shift controls the distance between the labels of dendrogram leaves and color bar.
##' 
##' @examples
##' data(HIV)
##' ## get the compositional data
##' compDat <- HIV$compDat[, 1:60]
##' ## taxonomic tree structure
##' taxonomy <- HIV$taxonomy
##' ## weak taxonomic hierarchy
##' weak <- hPAA(compDat, method = "Simpson", taxonomy = taxonomy, strong = FALSE)
##' plotHPAA(weak, lab.cex = 0.7, ylab = "Percentage change (%)", 
##'          y_scale = 1.1, y_shift = -1, cex.rowLabels = 1.1)
##' 
##' @import dendextend grDevices graphics utils
##' @importFrom stats as.dendrogram order.dendrogram
##' @importFrom RColorBrewer brewer.pal
##' @export
plotHPAA <- function(object, taxonomy = TRUE, rank.line = FALSE, 
                     label = NULL, hang = -1, axes = TRUE, frame.plot = FALSE, 
                     ann = TRUE, lab.cex = 0.8, main = "Dendrogram", sub = NULL, xlab = "", 
                     ylab = "Percentage change (%)", lim.y = -2, logy = TRUE, 
                     jitters = FALSE, select.level = "all", ...) {
  class(object) <- "hclust"
  compDat <- object$compDat
  taxon.matrix <- object$taxonomy
  ## check if the taxonomic tree structure is null
  if(is.null(taxon.matrix)) {
    taxonomy <- FALSE
  }
  ## rank line is only valid for strong hierarchy
  if(object$hierarchical.method != "Strong Hierarchical") {
    rank.line <- FALSE
  }
  ## check the label of the components in the compositional matrix
  if(is.null(label)) {
    label <- paste0("Taxon ", seq_len(length(object$labels)))
  }
  object$labels <- label
  ## get the diversity method
  method <- object$diversity.method
  if(jitters) {
    set.seed(123)
    jit <- c(head(sort(jitter(rep(0, length(object$height)), factor = 1)), -1),rep(0, 1))
    object$height <- as.double(object$height) - jit
  }
  if(logy) {
    if(jitters) {
      if(method == "Simpson") {
        object$height <- (1 - object$height / (1 - sum(diag(t(compDat) %*% compDat)) / nrow(compDat) - jit[1]*1.1)) * 100
      } else if(method == "Shannon") {
        object$height <- (1 - object$height / (- sum(diag(t(compDat) %*% log(compDat))) / nrow(compDat) - jit[1]*1.1)) * 100
      } else if(method == "Bray-Curtis") {
        object$height <- (1 - object$height / (sum(object$crit)  - jit[1]*1.1)) * 100
      }
    } else {
      if(method == "Simpson") {
        object$height <- (1 - object$height / (1 - sum(diag(t(compDat) %*% compDat)) / nrow(compDat))) * 100
      } else if(method == "Shannon") {
        object$height <- (1 - object$height / (- sum(diag(t(compDat) %*% log(compDat))) / nrow(compDat))) * 100
      } else if(method == "Bray-Curtis") {
        object$height <- (1 - object$height / sum(object$crit)) * 100
      }
    }
    ## adjust for nearly 0 entities
    if(any(object$height == 0)) {
      object$height[object$height == 0] <- min(object$height[object$height != 0]) * 1e-4
    }
    ## take the log
    object$height <- log10(object$height)
    thres <- abs(min(object$height)) + 0.1 ## adjust value under zeros
    object$height <- object$height + thres
    den <- as.dendrogram(object)
    if(rank.line) {
      plot(hang.dendrogram(den, hang = hang), nodePar = list(pch = NA, lab.cex = lab.cex),
           axes = axes, ylim = c(lim.y, max(object$height)), yaxt = "n",
           xlim = c(0.5, ncol(compDat) + 2),
           frame.plot = frame.plot, ann = ann, main = main, sub = sub)
    } else {
      plot(hang.dendrogram(den, hang = hang), nodePar = list(pch = NA, lab.cex = lab.cex),
           axes = axes, ylim = c(lim.y, max(object$height)), yaxt = "n",
           frame.plot = frame.plot, ann = ann, main = main, sub = sub)
    }
    axis(2, c(0, log10(c(0.5, 1, 5, 25, 50, 100)) + thres), labels = c("", c(0.5, 1, 5, 25, 50, 100)))
  } else {
    if(jitters) {
      if(method == "Simpson") {
        object$height <- (1 - object$height / (1 - sum(diag(t(compDat) %*% compDat)) / nrow(compDat) - jit[1]*1.1))
      } else if(method == "Shannon") {
        object$height <- (1 - object$height / (- sum(diag(t(compDat) %*% log(compDat))) / nrow(compDat) - jit[1]*1.1))
      } else if(method == "Bray-Curtis") {
        object$height <- (1 - object$height / (sum(object$crit) - jit[1]*1.1))
      }
    } else {
      if(method == "Simpson") {
        object$height <- (1 - object$height / (1 - sum(diag(t(compDat) %*% compDat)) / nrow(compDat)))
      } else if(method == "Shannon") {
        object$height <- (1 - object$height / (- sum(diag(t(compDat) %*% log(compDat))) / nrow(compDat)))
      } else if(method == "Bray-Curtis") {
        object$height <- (1 - object$height / sum(object$crit))
      }
    }
    den <- as.dendrogram(object)
    if(rank.line) {
      plot(hang.dendrogram(den, hang = hang), nodePar = list(pch = NA, lab.cex = lab.cex),
           axes = axes, ylim = c(lim.y, max(object$height)), 
           xlim = c(0.5, ncol(compDat) + 2),
           frame.plot = frame.plot, ann = ann, main = main, sub = sub, yaxt = "n")
    } else {
      plot(hang.dendrogram(den, hang = hang), nodePar = list(pch = NA, lab.cex = lab.cex),
           axes = axes, ylim = c(lim.y, max(object$height)),
           frame.plot = frame.plot, ann = ann, main = main, sub = sub, yaxt = "n")
    }
    axis(2, c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c(0, 20, 40, 60, 80, 100))
  }
  title(ylab = ylab, mgp=c(2.5,1,0), adj = 0.75)
  if(rank.line) {
    ## add the lines for different depth
    if(select.level == "all") {
      selected <- 1:(length(object$lines.Ind) - 1)
    } else {
      if(all(select.level %in% 1:(length(object$lines.Ind) - 1))) {
        selected <- (1:(length(object$lines.Ind) - 1))[select.level]
      } else {
        stop("invalid level for line")
      }
    }
    for(i in selected) {
      lines(x = c(-2, ncol(compDat) + 4), y = rep(object$height[object$lines.Ind[i]], 2), lty  = 2, col = "red")
      text(ncol(compDat) + 1, y = object$height[object$lines.Ind[i]], adj = c(0.2, 1.1),
           labels = paste0(c("Family", "Order", "Class", "Phylum")[i], 
                           # " (", ncol(compDat) - object$lines.Ind[i], ")"),
                           ""),
           # cex = 1)
           cex = 1.3)
      if(logy) {
        text(-1.2, y = object$height[object$lines.Ind[i]], adj = c(0.2, 1.1),
             labels = sprintf("%.1f", 10^(object$height[object$lines.Ind[i]] - thres)), cex = 1)
      } else {
        text(-1.2, y = object$height[object$lines.Ind[i]], adj = c(0.2, 1.1),
             labels = sprintf("%.1f", object$height[object$lines.Ind[i]] * 100), cex = 1)
      }
    }
  }
  ## plot the colored bar, only if taxonomic tree structure is provided
  if(taxonomy) {
    # taxon.matrix <- taxon.matrix[object$order, , drop = FALSE]
    taxonomy.list <- strsplit(taxon.matrix, split = ";", fixed = TRUE)
    names(taxonomy.list) <- rownames(taxon.matrix)
    depth <- max(sapply(taxonomy.list, length))
    # col.list <- c("YlOrRd", "Greens", "Reds", "Purples", "Oranges","Blues")
    col.list <- c("Set1", "Set2", "Pastel1", "Accent", "Set3", "Paired")  ## color of the bar
    the_bars <- NULL
    for(i in 1:depth) {
      classes <- sapply(taxonomy.list, function(x) x[i])
      cols <- colorRampPalette(brewer.pal(n=7, name=col.list[i]))(length(unique(classes)))
      cols <- cols[match(classes, unique(classes))]
      cols[is.na(classes)] <- NA
      the_bars <- cbind(the_bars, cols)
    }
    # if(is.null(y_shift)) {
    #   y_shift = min(object$height) - 0.2
    # }
    colored_bars_modified(colors = the_bars, dend = object, 
                          rowLabels = tail(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), depth),
                          ...)
  }
}