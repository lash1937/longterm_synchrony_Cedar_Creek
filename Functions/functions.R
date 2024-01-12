

#' Smith and Wilson's Evenness Index (Evar)
#' 
#' This function computes the Evar index of evenness from a vector of species abundances
#' for a given plot or community unit
#'
#' @param abundances 
#'
#' @return Numeric index between 0 and 1
#' @export
#'
#' @examples
#' 
#' abundances <- rpois(10, lambda = 10)
#' Evar(abundances)
#' 
Evar <- function(abundances){
  
  if(sum(abundances == 0) > 0){
    warning(
      "Removing species with zero abundance"
    )
  }
  abund_sub <- abundances[abundances != 0]
  
  if(length(abund_sub) <= 1){
    warning(
      "Returning NA: At least two species must be present to compute evenness"
    )
    return(NA)
  }
  # log the abundances
  l_n <- log(abund_sub)
  n <- length(l_n)
  
  
  # compute variance
  vln <- var(l_n) * (n - 1) / n
  
  # transform to constrain to [0, 1]
  return(
    1 - (2 / pi) * atan(vln)
  )
  
}






#' Compute Box-Cox Transformation for a given lambda
#'
#' @param x Vector of strictly positive data.
#' @param lambda Chosen lambda based on profile likelihood plot
#'
#' @return Box-Cox transformed x.
#' @export
#'
#' @examples
#' 
boxcox_transform <- function(x, lambda){
  
  if(sum(x <= 0) > 0){
    stop("Box-Cox Transform only appropriate for strictly positive data.
         Some observations <= 0.")
  }
  if(lambda == 0){
    return(log(x))
  } else{
    return((x^lambda - 1) / lambda)
  }
  
}




#' Export semEff table
#' 
#' Specialized function to write output from summary(semEff) to a LaTeX table
#'
#' @param x semEff object
#' @param outfile Output table filename
#'
#' @return .tex file
#' 
export_sem_table <- function(x, outfile, longtable = FALSE){
  
  library(tidyr)
  library(stringr)
  # capture each line character vector
  out <- capture.output(
    summary(x)
  )
  
  # capture headings
  subnames <- stringr::str_subset(
    as.vector(capture.output(
      summary(x),
      type = "message"
    )),
    pattern = "\\("
  )
  
  # create cleaner subnames
  subnames <- stringr::str_extract_all(
    subnames,
    "T[:alpha:]+"
  ) %>% stringr::str_replace(
    .,
    "^T", ""
  )
  
  # replace any multi-dash lines with \hline
  out <- out[grep("---", out, invert = T)]
  
  # handle first line separately
  out[1] <- str_trim(out[1]) %>%
    str_replace_all(., "[\\s]+", " & ") %>%
    paste0(" & & ", .) %>%
    str_replace(., "Std. &", "Std. ") %>%
    str_replace(., "Lower &", "Lower ") %>%
    str_replace(., "Upper &", "Upper ")
  
  # remove blank line at the end, and add hline
  out <- c(
    out[1],
    "",
    "\\hline",
    out[2:(length(out) - 1)]
  )
  
  # repace the blank lines with subtable names
  ndivs <- str_count(out[1], "&") 
  subnames_full <- subnames
  for(i in 1:length(subnames_full)){
    subnames_full[i] <- paste0("\\textbf{", subnames[i], "}", paste0(rep(" & ", ndivs), collapse = ""))
  }
  
  out[str_which(out, "^$")] <- subnames_full
  
  out[str_which(out, "Effect")[-1]] <- "\\hline"
  
  out <- str_replace_all(out, "\\|", " ") %>%
    str_replace_all(., "[\\s]+\\*", "") %>%
    str_trim(., side = "right") %>%
    str_replace_all(., "[\\s]{3,}+", " & ") %>%
    str_replace(., "^$", " & & & & & & ") %>%
    str_replace(., "INDIRECT ", "INDIRECT &") %>%
    str_replace(., "MEDIATORS ", "MEDIATORS &")
  
  # add horizontal space for some lines
  for(i in 1:length(out)){
    cnt_i <- str_count(out[i], pattern = "&")
    if(cnt_i < ndivs & out[i] != "\\hline"){
      pre <- paste0(rep(" &", ndivs - cnt_i), collapse = "")
      out[i] <- paste0(pre, out[i])
    }
  }
  
  # add linebreak characters
  for(i in 1:length(out)){
    if(out[i] != "\\hline"){
      out[i] <- str_replace(out[i], "$", "\\\\\\\\")
    }
  }
  
  if(longtable){
    out <- c(
      paste0("\\begin{longtable}[c]{", paste0(rep("l ", ndivs + 1), collapse = ""), "}"),
      "\\caption{Cap goes here. \\label{tbl:table-name}}\\\\",
      out,
      "\\end{longtable}"
    )
  } else{
    out <- c(
      paste0("\\begin{tabular}{", paste0(rep("l ", ndivs + 1), collapse = ""), "}"),
      out,
      "\\end{tabular}"
    )
  }
  
  out_con <- file(outfile)
    writeLines(print(out), outfile)
  close(out_con)
}







