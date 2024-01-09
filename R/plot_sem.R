##################################################################
# This function takes a few arguments to create a nice path diagram
##################################################################

#' Plot a path diagram for a semeff object
#'
#' @param semeff SemEff object with bootstrapped effects
#' @param rowcounts Vector with the numer of nodes (variables) desired in each row of the
#' path diagram (see example)
#' @param order Order for the variables. These get filled from left to right, top to bottom based
#' on the rowcounts argument
#' @param height Total height of the plotting area
#' @param width Width of the plotting area
#'
#' @return ggplot object
#' 
plot_sem <- function(semeff, rowcounts, order, height = 20, width = 10){
  library(tidyverse)
  
  # store some handy variables
  p <- length(order)
  
  # create a list of matrices to pull from
  mats <- vector(mode = "list", length = 4)
  names(mats) <- c("estim", "se", "low", "high")
  
  # create epmty matrices
  for(i in 1:length(mats)){
    mats[[i]] <- matrix(nrow = p, ncol = p)
    rownames(mats[[i]]) <- order
    colnames(mats[[i]]) <- order
  }
  
  # fill in bootstrapped results
  for(i in 1:p){
    for(j in 1:p){
      if(
        order[j] %in% names(semeff[["Bootstrapped Effects"]]) &
        order[i] %in% colnames(semeff[["Bootstrapped Effects"]][[order[j]]][["Direct"]])
      ){
        # mean from bootstraps
        mats$estim[i, j] <- mean(
          semeff[["Bootstrapped Effects"]][[order[j]]][["Direct"]][, order[i]]
        )
        
        # sd from bootstraps
        mats$se[i, j] <- sd(
            semeff[["Bootstrapped Effects"]][[order[j]]][["Direct"]][, order[i]]
        )
        
        # lower and upper CI bounds
        mats$low[i, j] <- quantile(
          semeff[["Bootstrapped Effects"]][[order[j]]][["Direct"]][, order[i]],
          probs = 0.025
        )
        mats$high[i, j] <- quantile(
          semeff[["Bootstrapped Effects"]][[order[j]]][["Direct"]][, order[i]],
          probs = 0.975
        )
      }
    }
  }

  # creating the layout
  x <- as.vector(sapply(
    rowcounts,
    function(s, w){
      return(seq(w / (s + 1), w - w / (s + 1), by = w / (s + 1)))
    },
    w = width
  ) %>% unlist())
  
  # y coordinates
  y <- seq(
    height / (length(rowcounts) + 1), 
    height - height / (length(rowcounts) + 1), 
    by = height / (length(rowcounts) + 1)
  ) 
  y <- y[order(y, decreasing = T)]
  y <- rep(y, rowcounts)
  
  # create node dataframe
  nodes <- tibble(
    label = order,
    x = x, y = y,
    xmin = x - 1,
    xmax = x + 1,
    ymin = y - 1,
    ymax = y + 1
  )
  
  # create the edges dataframe
  edges <- tibble(
    from = rep(order, p),
    to = rep(order, each = p),
    estim = as.vector(mats$estim),
    se = as.vector(mats$se),
    low = as.vector(mats$low),
    high = as.vector(mats$high)
  )
  
  # remove NAs
  edges <- edges[complete.cases(edges), ]
  
  # add start and end coords
  edges <- left_join(
    edges, nodes, 
    by = join_by(from == label)
  ) %>% left_join(
    ., nodes,
    by = join_by(to == label)
  ) %>% dplyr::rename(
    x = x.x, y = y.x,
    xend = x.y, yend = y.y
  )
  
  # subtract 1 from the yends
  edges$yend <- edges$yend + 1
  
  # add aesthetics
  edges <- edges %>% mutate(
    weight = abs(estim) / (0.5 * max(abs(estim))),
    color = ifelse(estim < 0, "negative", "positive"),
    linetype = ifelse(
      low < 0 & high > 0, "CI overlaps 0", "No overlap"
    )
  )
  
  # create plot
  ggplot() + theme_void() +
    geom_segment(
      data = edges,
      aes(
        x = x, y = y, xend = xend, yend = yend, 
        color = color, linetype = linetype
      ),
      arrow = arrow(
        length = unit(0.02, units = "npc"),
        type = "closed"
      ),
      size = edges$weight
    ) +
    geom_rect(
      data = nodes, 
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "white",
      color = "black"
    ) +
    geom_text(data = nodes, aes(x = x, y = y, label = label)) +
    scale_color_manual(values = c("red", "black")) +
    scale_linetype_manual(values = c("dashed", "solid"))
  
}




