##' Plot Method for the Object with class "hPAA" from hPAA Analysis
##' 
##' Plot the results of hierarchical Principal Amalgamation Analysis (PAA) as 
##' scree plot, which shows the percentage change in the diversity loss as a 
##' function of the number of principal compositions.
##'
##' @usage 
##' plotLine(object, 
##'          xlab = "Number of PCs", 
##'          ylab = "Percentage change (\%)", 
##'          main = paste0(object$diversity.method, " Index"), ...)
##'
##' @param object the output with class "hPAA" from function hPAA().
##' @param xlab a label for the x axis.
##' @param ylab	a label for the y axis, defaults to a description of y.
##' @param main a main title for the plot.
##' @param ... further graphical arguments as in the plot.default.
##' 
##' @examples
##' data(HIV)
##' ## get the compositional data
##' compDat <- HIV$compDat[, 1:60]
##' ## taxonomic tree structure
##' taxonomy <- HIV$taxonomy
##' ## weak taxonomic hierarchy
##' weak <- hPAA(compDat, method = "Simpson", taxonomy = taxonomy, strong = FALSE)
##' plotLine(weak)
##' 
##' @export
plotLine <- function(object, 
                     xlab = "Number of PCs", 
                     ylab = "Percentage change (%)", 
                     main = paste0(object$diversity.method, " Index"), ...) {
  plotData <- pComp(object)
  plot(1:length(plotData$height), plotData$height, xlab = xlab,
       ylab = ylab, main = main, type = "l", xaxt = 'n', ...)
  points(1:length(plotData$height), plotData$height, ...)
  axis(1, at=1:length(plotData$height), labels=factor(length(plotData$height):1))
}