#' Title
#'
#' @param tpDrHealthy
#' @param tpDrExp
#' @param sampleID
#' @param fileIDExp
#' @param fileIDHealthy
#' @param bamFilesHealthInput
#' @param bamFilesExpInput
#' @param tissue
#' @param replicate
#' @param bedExt
#' @param bamExt
#' @param peakCaller
#' @param bamDr
#' @param bedDr
#' @param cond1
#' @param cond2
#' @param bamFiles
#'
#' @return
#' @export
#'
#' @examples
createDiffBindSampleTable <- function(tpDrHealthy, tpDrExp, sampleID, fileIDExp, fileIDHealthy, bamFiles = F, bamFilesHealthInput = F, bamFilesExpInput = F, tissue = NULL, replicate = NULL,
                                      bedExt = "_trim_bowtie2UPBlkLstRm_MACS2_peaks.narrowPeak", bamExt = "_trim_bowtie2UPBlkLstRm.sorted.bam", peakCaller = "narrow",
                                      bamDr = "alignments/", bedDr = "peakCalls/", cond1 = "Healthy", cond2 = "OC"){


  #if (length(sampleID) == length(tissue))

  df <- data.frame(SampleID = sampleID, Tissue = tissue, PeakCaller = peakCaller, Replicate = replicate)

  condition <- c(rep(cond1, length(fileIDHealthy)), rep(cond2, length(fileIDExp)))
  # Add condition --------------------------------------------------------------------------------------------------- #
  if (nrow(df) != length(condition)){
    print("Length of sampleID not equal to length of condition")
    return()

  }

  df <- cbind(df, data.frame(Condition = condition))

  # Add bedfiles --------------------------------------------------------------------------------------------------- #
  # Check if the bed files exist ctrl

  bedFls <- c()

  for (dx in fileIDHealthy){
    fn <- paste0(tpDrHealthy, bedDr,dx,  bedExt)
    if (file.exists(fn)){

      bedFls <- append(bedFls, fn)

    } else {

      print(paste0("File " , fn, " does not exist"))
      return()
    }
  }

  # Check if the bed files exist exp
  for (dx in fileIDExp){
    fn <- paste0(tpDrExp, bedDr,dx,  bedExt)
    if (file.exists(fn)){

      bedFls <- append(bedFls, fn)

    } else {

      print(paste0("File " , fn, " does not exist"))
      return()
    }
  }


  if (assertthat::are_equal(nrow(df), length(bedFls))){
    df <- cbind(df, data.frame(Peaks = bedFls))
  } else {

    print("Length of bedfiles not equal the sampleID")
    return()

  }


  # Bam files ============================================================================= #
  if (bamFiles){

    bamFls <- c()
    bamFls <- c()


    # loop through healthy
    for (dx in fileIDHealthy){
      fn <- paste0(tpDrHealthy, bamDr, dx, bamExt)
      if (file.exists(fn)){

        bamFls <- append(bamFls, fn)


      } else {
        print(paste0("File " , fn, " does not exist"))
        return()
      }

    }
    # loop through exp
    for (dx in fileIDExp){
      fn <- paste0(tpDrExp, bamDr, dx, bamExt)
      if(file.exists(fn)){
        bamFls <- append(bamFls, fn)
      }
      else {
        print(paste0("File " , fn, " does not exist"))
        return()

      }



    }
    df <- cbind(df, data.frame(bamReads = bamFls))
  }

  # Bam Files Input ======================================================================================
  bamFlsInput <- c()
  #browser()
  if (!isEmpty(bamFilesExpInput)){

    # healthy
    for (dx in bamFilesHealthInput){

      fn <- paste0(tpDrHealthy, bamDr, dx, bamExt)

      if (file.exists(fn)){

        bamFlsInput <- append(bamFlsInput, fn)

      } else {
        print(paste0("File " , fn, " does not exist"))
        return()
      }
    }


    # Exp

    for (dx in bamFilesExpInput){

      fn <- paste0(tpDrExp, bamDr, dx, bamExt)

      if (file.exists(fn)){

        bamFlsInput <- append(bamFlsInput, fn)

      } else {
        print(paste0("File " , fn, " does not exist"))
        return()

      }

    }

    df <- cbind(df, data.frame(bamControl = bamFlsInput))

  } # if (bamFilesExpInput)



  df
}

