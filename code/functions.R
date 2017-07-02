# Order substrings of a character vector
order_substr <- 
  function(x, split = "\\s*,\\s*|\\s*/\\s*|\\s*-\\s*|\\s+", 
           perl = TRUE, collapse = " ", reverse = FALSE) {
    # Split each element into its own character vector at user-defined split patterns
    tmp <- strsplit(x, split = split, perl = perl)
    # Sort the substrings alphabetically and collapse.
    tmp2 <- vapply(tmp, 
                   function(x) paste(sort(x, decreasing = reverse), 
                                     collapse = collapse
                   ),
                   FUN.VALUE = character(1)
    )
    # Trim leading and trailing whitespace
    tmp3 <- gsub(tmp2, patt = "^\\s|\\s$", rep = "")
    
    return(tmp3)
  }

# Save objects in a list as separate files
save_list <- 
  function(x, dir, ext = ".rds", names = NULL, 
           cores = 1, ...){
    if(!file.exists(dir)) dir.create(dir)
    
    if(is.null(names)) names <- names(x)
    require(parallel)
    
    dir2 <- ifelse(grepl(dir, patt = "/$"), dir, paste0(dir, "/"))
    
    if(grepl(ext, patt = "\\.?rds$", ignore.case = TRUE)){
      mclapply(names, function(i) saveRDS(x[[i]],
                                          paste0(dir2, i, ext),
                                          ...),
               mc.cores = cores
      )
    }
    
    if(grepl(ext, patt = "\\.?csv$", ignore.case = TRUE)){
      if(require(data.table)){
        mclapply(names, function(i) fwrite(x[[i]],
                                           file = paste0(dir2, i, ext),
                                           ...
        ),
        mc.cores = cores
        )
      }else{
        mclapply(names, function(i) write.csv(x[[i]],
                                              file = paste0(dir2, i, ext),
                                              ...),
                 mc.cores = cores
        )
      }
    }
    
    if(grepl(ext, patt = "\\.?shp$", ignore.case = TRUE)){
      if(require(maptools)){
        mclapply(names, function(i) writeSpatialShape(x[[i]],
                                                      paste0(dir2, i, ext),
                                                      ...),
                 mc.cores = cores
        )
      }else{ stop("Package 'maptools' must be installed to save spatial objects") }
    }
  }