library(GenomicRanges)

#' Remove scaffolds. This funciton removes all the scaffold chromosomes from a genomic ranges object and then removes the sequence factors
#'
#' @param GR Genomic ranges object with scaffolds to remove
#'
#' @return Returns the same GR object with the scaffolds removed. Also sets the chromosome style to UCSC
#' @export
#'
#' @examples
rmScf <- function(GR){

  GR <- GR[seqnames(GR) == 'chr1' | seqnames(GR) == 'chr2' | seqnames(GR) == 'chr3' |
             seqnames(GR) == 'chr4' | seqnames(GR) == 'chr5' | seqnames(GR) == 'chr6' |
             seqnames(GR) == 'chr7' | seqnames(GR) == 'chr8' | seqnames(GR) == 'chr9' |
             seqnames(GR) == 'chr10' | seqnames(GR) == 'chr11' | seqnames(GR) == 'chr12' |
             seqnames(GR) == 'chr13' | seqnames(GR) == 'chr14' | seqnames(GR) == 'chr15' |
             seqnames(GR) == 'chr16' | seqnames(GR) == 'chr17' | seqnames(GR) == 'chr18' |
             seqnames(GR) == 'chr19' | seqnames(GR) == 'chr20' | seqnames(GR) == 'chr21' |
             seqnames(GR) == 'chr22' | seqnames(GR) == 'chrX' | seqnames(GR) == 'chrY'
           ]  # Remove peaks on unassigned scaffolds

  seqlevels(GR) <- seqlevelsInUse(GR) # remove sequence levels not in use
  seqlevelsStyle(GR) <- "UCSC" # change the style of the chromosome seqnames
  GR
}
