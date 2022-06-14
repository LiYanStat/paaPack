##' Plot Method for the Object with class "hPAA" from hPAA Analysis
##' 
##' Plot the results of hierarchical Principal Amalgamation Analysis (PAA) as 
##' NMDS ordination plot to visualize the changes in the between-sample distance
##' patterns before and after HPAA with any given number of principal compositions.
##'
##' @usage 
##' plotMDS(object, cutoff, ...)
##'
##' @param object the output with class "hPAA" from function hPAA().
##' @param cutoff number of principal compostions to cut off from the HPAA path.
##' @param ... further graphical arguments for the geom_points(). 
##'            E.g., size controls the size of points in the plot.
##' 
##' @examples
##' data(HIV)
##' ## get the compositional data
##' compDat <- HIV$compDat[, -61]
##' ## taxonomic tree structure
##' taxonomy <- HIV$taxonomy
##' ## weak taxonomic hierarchy
##' weak <- hPAA(compDat, method = "Simpson", taxonomy = taxonomy, strong = FALSE)
##' plotMDS(weak, cutoff = 20)
##' @import ggplot2 ggforce
##' @importFrom vegan metaMDS
##' @export
plotMDS <- function(object, cutoff, ...) {
  compDat <- object$compDat
  aggDat <- cutoffComp(object, compDat, cutoff = cutoff, num.entities = TRUE)
  origPlot <- metaMDS(compDat, trymax = 1000, trace = FALSE)
  aggPlot <- metaMDS(aggDat, trymax = 1000,trace = FALSE)
  
  c.points <- (aggPlot$points + origPlot$points) / 2
  distance <- sqrt(rowSums((aggPlot$points - origPlot$points)^2)) / 2
  plotCircle <- data.frame(cbind(c.points, distance))
  colnames(plotCircle) <- c("MDS1", "MDS2", "distance")
  
  comb.point <- 
    data.frame(rbind(cbind(aggPlot$points, "Principal compositions"),
                     cbind(origPlot$points, "Original compositions")))
  colnames(comb.point) <- c("MDS1", "MDS2", "Type")
  comb.point$MDS1 <- as.numeric(comb.point$MDS1)
  comb.point$MDS2 <- as.numeric(comb.point$MDS2)
  
  comb.point$Type <- factor(comb.point$Type, levels = c("Original compositions", "Principal compositions"))
  ## plot the results as ggplot2
  p <- ggplot() + 
    geom_circle(data = plotCircle, aes_string(x0="MDS1", y0="MDS2", r="distance"), color="darkgrey") + 
    coord_equal() +
    geom_point(data = comb.point, aes_string(x="MDS1", y="MDS2", color="Type"), ...) +
    scale_color_manual(values = c("red", "blue")) +
    ggtitle(paste0(object$diversity.method, " Index, ", "Number of PCs: ", ncol(aggDat))) + 
    theme_bw() + 
    # lims(x = range(comb.point$MDS1), y = range(comb.point$MDS2)) + 
    theme(legend.position = "bottom", 
          legend.direction = "horizontal", 
          legend.title = element_blank(), 
          legend.text = element_text(size = 12), 
          legend.background = element_rect(size=0.5, linetype="solid", colour ="black"),
          plot.title = element_text(hjust = 0.5))
  return(p)
}