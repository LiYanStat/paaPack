##' Microbiome, HIV infection and MSM factor
##'
##' A dataset containing the number of counts of 60 different genera in a group
##' of 128 HIV - infected samples.
##'
##' @format The \code{list} is composed by a data.frame of compositional data for 
##' 60 taxa at genus level and a classification variable, and a taxonomic matrix 
##' for corresponding taxnomic tree structure of the 60 genera.
##' \describe{
##'   \item{compDat}{The first 60 columns is the compositional data matrix for the 60 taxa, with 
##'   last column indicating if the individual is \code{MSM} (\emph{Men Sex with Men}) or not (\code{nonMSM}).}
##'   \item{taxonomy}{The taxonomic matrix for the 60 genera, from taxonomic rank Kingdom to Genus.}
##' }
##' @docType data
##' @name HIV
##' @references 
##' \itemize{
##'  \item Rivera-Pinto et al. (2018). Balances: a new perspective for microbiome analysis.
##'  \item Quinn and Erb (2020). Amalgams: data-driven amalgamation for the dimensionality reduction of compositional data.
##'  \item Chen et al. (2021). Principal Amalgamation Analysis.
##' }
##' @examples
##' data(HIV)
##' compDat <- HIV$compDat
##' taxonomy <- HIV$taxonomy
NULL