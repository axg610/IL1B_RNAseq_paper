cleanup_objects <- function(keep = readLines("../R/cleanup_ignore.txt")) {
  all_objs <- ls(envir = .GlobalEnv)
  to_remove <- setdiff(all_objs, keep)
  rm(list = to_remove, envir = .GlobalEnv)
  # gc()
}

plot_euler <- function(
    logicalMatrix, 
    plotType = "euler",
    shapeType = "circle", 
    plotTitle = "title here",
    quants = c("counts"),     # vector of metrics to show (e.g., counts, percent)
    cutoff = 0,               # intersections smaller than this number will be removed
    returnAsFunction = FALSE, # helps for patchworking multiple eulers together
    showLabels = T,
    fillColors = c("#FFFFFF", "#D9D9D9", "#A6CEE3"),
    aspectRatio = 1
){
  
  # a wrapper for the eulerr:euler function to create a Euler diagram given a
  # logical matrix specifying groups as colnames, presence/absence as 1/0, and
  # gene names as rownames.
  
  # prepare intersections
  obj <- logicalMatrix %>%
    mutate(
      across(
        everything(),
        ~ ifelse(. == 1, cur_column(), NA)
      )
    ) %>%
    rowwise() %>%
    mutate(group = paste(na.omit(c_across(cols = everything())), collapse = "&")) %>%
    count(group) %>%
    deframe()
  
  # impose cutoff
  obj <- obj[obj >= cutoff]
  
  # calculate and plot
  if(plotType == "euler"){
    obj <- euler(obj, shape = shapeType)
  }
  else if(plotType == "venn"){
    obj <- venn(obj)
  }
  
  if(returnAsFunction){
    function() {
      plot(
        obj,
        quantities = list(type = quants),
        main = plotTitle,
        asp = aspectRatio,
        labels = showLabels,
        fills = list(fill = fillColors)
      )
    }
  }
  
  else if(!returnAsFunction){
    plot(
      obj,
      quantities = list(type = quants),
      main = plotTitle,
      asp = aspectRatio,
      labels = showLabels,
      fills = list(fill = fillColors)
    )
    
  }
  
}

generic_l2f_heatmap_colors <- function(
    mat,
    colorscale = c(
      "steelblue4", "steelblue2" , "white", "firebrick2", "firebrick4")
){
  
  # return a pretty log2fold heatmap color palette given a matrix
  
  colors = circlize::colorRamp2(
    c(
      seq(
        quantile(mat, 0.02),
        -0.05,
        length = 75
      ),
      seq(
        -0.49,
        0.49,
        length = 50
      ),
      seq(
        0.5,
        quantile(mat, 0.98),
        length = 75
      )
    ),
    colorRampPalette(
      colorscale
    )(200)
  )
  
  colors
}
