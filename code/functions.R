# Cole Tanigawa-Lau
# Tue Jul  4 16:00:26 2017
# Description: Functions for approximate georeferencing guide.

# CONVENIENCE FUNCTIONS -----------------
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

zero_na <- function(x) if(length(x) == 0) NA else x

class_nrow <- function(x, class){
  if(class(x) == class){nrow(x)} else 0
}
# SPECIALIZED MATCHING FUNCTIONS ---------------------------------------------
# Using locations and sublocations
match_lsl <- function(tps, const, urb_const, polys){
  if( !require(purrr) ) stop("match_lsl requires the purrr package.")
  
  tconst <- ifelse(is.na(tps$constid07), tps$constid02, tps$constid07)
  
  if(!is.na(tconst)){
    # Load sublocations for the polling station's constituency
    tsl <- paste0("data/sl_constb/", tconst, ".rds") %>% 
      readRDS()
    # Load sublocations for the polling station's constituency
    tmp <- paste0("data/mp_constb/", tconst, ".rds") %>% 
      readRDS()
    
    # Hand check of threshold
    # sl_97 <- data.frame(ps_sl97 = tps$sl97, tsl$SUBLOCATIO, sd = stringdist(tps$sl97, tsl$SUBLOCATIO, method = "jw", p = 0.15)) %>% filter(sd < 0.3)
    # loc_97 <- data.frame(ps_loc97 = tps$loc97, tsl$LOCATION, sd = stringdist(tps$loc97, tsl$LOCATION, method = "jw", p = 0.15)) %>% filter(sd < 0.3)
    # sl_02 <- data.frame(ps_sl02 = tps$sl02, tsl$SUBLOCATIO, sd = stringdist(tps$sl02, tsl$SUBLOCATIO, method = "jw", p = 0.15)) %>% filter(sd < 0.3)
    # loc_07 <- data.frame(ps_loc07 = tps$loc07, tsl$LOCATION, sd = stringdist(tps$loc02, tsl$LOCATION, method = "jw", p = 0.15)) %>% filter(sd < 0.3)
    # sl_07 <- data.frame(ps_sl07 = tps$sl07, tsl$SUBLOCATIO, sd = stringdist(tps$loc02, tsl$SUBLOCATIO, method = "jw", p = 0.15)) %>% filter(sd < 0.3)
    
    ## Quality measures for location-sublocation information.
    
    # Calculate string and geographic distances between ps 1997 sublocation and the sublocation shapefile
    strdist_sl97 <- stringdist(a = tps$sl97, b = tsl$SUBLOCATIO, method = 'jw', p = 0.15)
    min_loc97 <- min(strdist_sl97)
    sl_ind97 <- which(strdist_sl97 == min_loc97 & strdist_sl97 < 0.117)
    
    strdist_loc97 <- stringdist(a = tps$loc97, b = tsl$LOCATION, method = 'jw', p = 0.15)
    min_loc97 <- min(strdist_loc97)
    loc_ind97 <- which(strdist_loc97 == min_loc97 & strdist_loc97 < 0.106)
    
    # Calculate string distances between ps 2002 sublocation and the sublocation shapefile
    strdist_sl02 <- stringdist(a = tps$sl02, b = tsl$SUBLOCATIO, method = 'jw', p = 0.15)
    min_loc02 <- min(strdist_sl02)
    sl_ind02 <- which(strdist_sl02 == min_loc02 & strdist_sl02 < 0.097)
    
    strdist_loc02 <- stringdist(a = tps$loc02, b = tsl$LOCATION, method = 'jw', p = 0.15)
    min_loc02 <- min(strdist_loc02)
    loc_ind02 <- which(strdist_loc02 == min_loc02 & strdist_loc02 < 0.109)
    
    # Calculate string distances between ps 2007 sublocation and the sublocation shapefile
    strdist_sl07 <- stringdist(a = tps$sl07, b = tsl$SUBLOCATIO, method = 'jw', p = 0.15)
    min_loc07 <- min(strdist_sl07)
    sl_ind07 <- which(strdist_sl07 == min_loc07 & strdist_sl07 < 0.0757)
    
    strdist_loc07 <- stringdist(a = tps$loc07, b = tsl$LOCATION, method = 'jw', p = 0.15)
    min_loc07 <- min(strdist_loc07)
    loc_ind07 <- which(strdist_loc07 == min_loc07 & strdist_loc07 < 0.095)
    
    # Calculate 1997 quality measures: distance to nearest sublocation, distance to nearest location, sublocation is nested in location.
    if(length(sl_ind97) > 0){
      sl97_dist <- gDistance(tsl[sl_ind97, ], tmp, byid = TRUE) %>% apply(1, min)
      match_sl97 <- tsl[sl_ind97, ]$SL_ID
      nest_97 <- tsl[sl_ind97[which(sl_ind97 %in% loc_ind97)], ]$SL_ID
      nest_97 <- ifelse(length(nest_97) == 0, NA, nest_97)
    }else{sl97_dist <- NA; match_sl97 <- NA; nest_97 <- NA}
    
    
    if(length(loc_ind97) > 0){
      loc97_dist <- gDistance(tsl[loc_ind97, ], tmp, byid = TRUE) %>% apply(1, min)
      match_loc97 <- tsl[loc_ind97, ]$LOC_ID
    }else{loc97_dist <- NA; match_loc97 <- NA}
    
    
    # Calculate 2002 quality measures: distance from sublocation, distance from location, sublocation is nested in location.
    if(length(sl_ind02) > 0){
      sl02_dist <- gDistance(tsl[sl_ind02, ], tmp, byid = TRUE) %>% apply(1, min)
      match_sl02 <- tsl[sl_ind02, ]$SL_ID
      nest_02 <- tsl[sl_ind02[which(sl_ind02 %in% loc_ind02)], ]$SL_ID
      nest_02 <- ifelse(length(nest_02) == 0, NA, nest_02)
    }else{sl02_dist <- NA; match_sl02 <- NA; nest_02 <- NA}
    
    
    if(length(loc_ind02) > 0){
      loc02_dist <- gDistance(tsl[loc_ind02, ], tmp, byid = TRUE) %>% apply(1, min)
      match_loc02 <- tsl[loc_ind02, ]$LOC_ID
    }else{loc02_dist <- NA; match_loc02 <- NA}
    
    if(length(sl_ind07) > 0){
      sl07_dist <- gDistance(tsl[sl_ind07, ], tmp, byid = TRUE) %>% apply(1, min)
      match_sl07 <- tsl[sl_ind07, ]$SL_ID
      nest_07 <- tsl[sl_ind07[which(sl_ind07 %in% loc_ind07)], ]$SL_ID
      nest_07 <- ifelse(length(nest_07) == 0, NA, nest_07)
    }else{sl07_dist <- NA; match_sl07 <- NA; nest_07 <- NA}
    
    
    if(length(loc_ind07) > 0){
      loc07_dist <- gDistance(tsl[loc_ind07, ], tmp, byid = TRUE) %>% apply(1, min)
      match_loc07 <- tsl[loc_ind07, ]$LOC_ID
    }else{loc07_dist <- NA; match_loc07 <- NA}
    
    # Poly match columns
    mpolys <- list(match_sl97, match_loc97, match_sl02, match_loc02, match_sl07, match_loc07)
    mpolys2 <- map2(polys, mpolys, paste, sep = "_") %>% unlist()
    
    strdists <- list(strdist_sl97[sl_ind97], strdist_loc97[loc_ind97], 
                     strdist_sl02[sl_ind02], strdist_loc02[loc_ind02],
                     strdist_sl02[sl_ind07], strdist_loc02[loc_ind07]) %>%
      lapply(zero_na) %>% unlist()
    
    if(sum(sapply(list(sl_ind97, loc_ind97, sl_ind02, loc_ind02, sl_ind07, loc_ind07), length)) > 0){
      
      # Match ps to potential mp using 2002 truncs, types, and subtypes.
      tdist <- stringdist(tps$trunc02, tmp$trunc, method = "jw", p = 0.15)
      t_ind <- which(tdist == min(tdist) & tdist < 0.07)
      m_trunc <- tmp[t_ind, ]
      
      type_dum <- ifelse(m_trunc$type == tps$type02, 1, 0)
      subtype_dum <- ifelse(m_trunc$subtype == tps$subtype02, 1, 0)
      trunc_dist <- tdist[t_ind]
      
      if(length(t_ind) > 0){
        m_trunc@data <- mutate(m_trunc@data, type_dum, subtype_dum, trunc_dist, 
                               nested97  = ifelse(SL_ID %in% nest_97, 1, 0), 
                               nested02  = ifelse(SL_ID %in% nest_02, 1, 0),
                               sl97_dist = sl97_dist[t_ind], loc97_dist = loc97_dist[t_ind],
                               sl02_dist = sl02_dist[t_ind], loc02_dist = loc02_dist[t_ind],
                               sl07_dist = sl07_dist[t_ind], loc07_dist = loc07_dist[t_ind],
                               ps_name   = tps$name02, ps_trunc = tps$trunc02)
        
        # Find shortest distance to a matched polygon and compare it to a threshold that is a function of the constituency size.
        min_dist <- min(m_trunc@data[, grep(names(m_trunc@data), patt = "[0-9]_dist")], na.rm = TRUE)
        bdist <- bbox(const[const$id == tconst, ]) %>% apply(1, diff) %>% min()
        # If a match is in a rural constituency, allow for a distance threshold of up to 5km.
        # If a match is in an urban constituency, set a threshold of 2km.
        dist_thresh <- ifelse(tconst %in% urb_const,
                              2000,
                              min(bdist * 0.05, 5000)
        )
        
        # If there are multiple matches, pick those with the high sum of quality scores.
        if(sum(m_trunc$trunc_dist == max(trunc_dist)) > 1){
          tqual <- m_trunc$type_dum + m_trunc$subtype_dum + m_trunc$nested97 + m_trunc$nested02
          m_trunc2 <- m_trunc[which((tqual) == max(tqual)), ]
        }else{
          m_trunc2 <- m_trunc
        }
        
        if(min_dist <= dist_thresh){
          out <- list(ps = as.data.frame(tps), match_mp = m_trunc2, 
                      match_polys = data.frame(poly = mpolys2, 
                                               strdist = strdists) %>% 
                        filter(., !duplicated(.))
          )
        }
        
        if(length(t_ind) == 0 || min_dist > dist_thresh){
          out <- list(ps = as.data.frame(tps), 
                      match_mp = "NO MP MATCH",
                      match_polys = data.frame(poly = mpolys2, 
                                               strdist = strdists) %>% 
                        filter(., !duplicated(.))
          )
        }
      }else{
        out <- list(ps = as.data.frame(tps), 
                    match_mp = "NO MP MATCH", 
                    match_polys = data.frame(poly = mpolys2, 
                                             strdist = strdists) %>% 
                      filter(., !duplicated(.))
        )
      }
      
    }else{
      out <- list(ps = as.data.frame(tps), 
                  match_mp = "NO MP MATCH", 
                  match_polys = "NO POLYGON MATCH")
    }
  }else{
    out <- list(ps = as.data.frame(tps),
                match_mp = "NO OLD CONSTITUENCY",
                match_polys = "NO OLD CONSTITUENCY")
  }
  
  return(out)
}

match_inf_loc13 <- function(tps, const, urb_const, polys){
  if( !require(purrr) ) stop("match_inf_loc13 requires the purrr package.")
  
  if(!is.na(tps$constid13)){
    # Load sublocation for the polling station's constituency
    tsl <- paste0("data/sl_const13b/", tps$constid13, ".rds") %>% 
      readRDS()
    # Load match points for the constituency
    tmp <- paste0("data/mp_const13b/", tps$constid13, ".rds") %>% 
      readRDS()
    
    # Calculate string and geographic distances between each potential ps 2013 sublocation and the sublocations from the sl shapefile
    sl_ind13 <- NULL
    sl_strdist <- NULL
    for(jj in 1:length(unlist(tps$sl13))){
      strdist_sl13 <- stringdist(a = tps$sl13[[1]][jj], b = tsl$SUBLOCATIO, method = "jw", p = 0.15)
      min_sl13 <- min(strdist_sl13)
      tind <- which(strdist_sl13 == min_sl13 & strdist_sl13 < .0814)
      sl_strdist <- c(sl_strdist, strdist_sl13[tind])
      sl_ind13 <- c(sl_ind13, tind)
    }
    
    # Calculate string and geographic distances between each potential ps 2013 location and the locations in the sl shapefile
    loc_ind13 <- NULL
    loc_strdist <- NULL
    for(jj in 1:length(unlist(tps$loc13))){
      strdist_loc13 <- stringdist(a = tps$loc13[[1]][jj], b = tsl$LOCATION, method = 'jw', p = 0.15)
      min_loc13 <- min(strdist_loc13)
      tind <- which(strdist_loc13 == min_loc13 & strdist_loc13 < 0.088)
      loc_strdist <- c(loc_strdist, strdist_loc13[tind])
      loc_ind13 <- c(loc_ind13, tind)
    }
    
    mpolys <- list(tsl$SL_ID[sl_ind13], tsl$LOC_ID[loc_ind13])
    mpolys2 <- map2(polys, mpolys, paste, sep = "_") %>% unlist()
    
    strdists <- list(strdist_sl13[sl_ind13], strdist_loc13[loc_ind13]) %>%
      lapply(zero_na) %>% unlist()
    
    # Quality measures
    if(length(c(sl_ind13, loc_ind13)) > 0){
      
      if(length(sl_ind13) > 0){
        match_sl13 <- tsl[unique(sl_ind13), ]$SL_ID
        dist_sl13 <- gDistance(tsl[sl_ind13, ], tmp, byid = TRUE) %>% apply(1, min)
      }else{match_sl13 <- NA; dist_sl13 <- NA}
      
      if(length(loc_ind13) > 0){
        match_loc13 <- tsl[unique(loc_ind13), ]$LOC_ID
        dist_loc13 <- gDistance(tsl[loc_ind13, ], tmp, byid = TRUE) %>% apply(1, min)
      }else{match_loc13 <- NA; dist_loc13 <- NA}
      
      tdist <- stringdist(tps$trunc13, tmp$trunc, method = "jw", p = 0.15)
      t_ind <- which(tdist == min(tdist) & tdist < 0.0808)
      m_trunc <- tmp[t_ind, ]
      
      type_dum <- ifelse(m_trunc$type == tps$type13, 1, 0)
      subtype_dum <- ifelse(m_trunc$subtype == tps$subtype13, 1, 0)
      trunc_dist <- tdist[t_ind]
      
      if(length(t_ind) > 0){
        m_trunc@data <- data.frame(m_trunc@data, type_dum, subtype_dum, trunc_dist, 
                                   loc13_dist = dist_loc13[t_ind], sl13_dist = dist_sl13[t_ind],
                                   ps_name = tps$name13, ps_trunc = tps$trunc13)
        
        # Find shortest distance to a matched polygon and compare it to a threshold that is a function of the constituency size.
        min_dist <- min(m_trunc@data[, grep(names(m_trunc@data), patt = "[0-9]_dist")], na.rm = TRUE)
        bdist <- bbox(const[const$constid == tps$constid13, ]) %>% apply(1, diff) %>% min()
        # If a match is in a rural constituency, allow for a distance threshold of up to 5km.
        # If a match is in an urban constituency, set a threshold of 2km.
        dist_thresh <- ifelse(tps$constid13 %in% urb_const,
                              2000,
                              min(bdist * 0.05, 5000)
        )
        
        # If more than one match, take those with the highest summed quality scores.
        if(sum(m_trunc$trunc_dist == max(trunc_dist)) > 1){
          m_trunc2 <- m_trunc[which((m_trunc$type_dum + m_trunc$subtype_dum) == max(m_trunc$type_dum + m_trunc$subtype_dum)), ]
        }else{
          m_trunc2 <- m_trunc
        }
        
        if(min_dist <= dist_thresh){
          out <- list(ps = as.data.frame(tps), 
                      match_mp = m_trunc2, 
                      match_polys = data.frame(poly = mpolys2, 
                                               strdist = strdists) %>% 
                        filter(., !duplicated(.))
          )
        }else{
          out <- list(ps = as.data.frame(tps), 
                      match_mp = "NO MP MATCH", 
                      match_polys = data.frame(poly = mpolys2, 
                                               strdist = strdists) %>% 
                        filter(., !duplicated(.))
          )
        }
      }else{
        out <- list(ps = as.data.frame(tps), 
                    match_mp = "NO MP MATCH", 
                    match_polys = data.frame(poly = mpolys2, 
                                             strdist = strdists) %>% 
                      filter(., !duplicated(.))
        )
      }
    }else{
      out <- list(ps = as.data.frame(tps),
                  match_mp = "NO MP MATCH",
                  match_polys = "NO LOCATION MATCH")
    }
  }else{
    out <- list(ps = as.data.frame(tps), 
                match_mp = "NO 2013 CONSTITUENCY", 
                match_polys = "NO 2013 CONSTITUENCY")
  }
  return(out)
}