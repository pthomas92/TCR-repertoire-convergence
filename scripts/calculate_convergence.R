
##### Load libraries #####

library(tidyverse)

setwd('~/OneDrive - University College London/_Leo Post Doc/_Human Challenge Repertoires/')

#####

##### define functions #####

mend_timepoints = function(x, tp, tol = 1, direction = c('earlier', 'later')){
  
  # function in case of misaligned timepoints
    # this function over-writes one column with the name to be tested (by default taking the closest timepoint).
    # tol sets the 'tolerance'. By default this is 1, meaning if there is more than 1 day to the next available timepoint
      # the function will error.
  
  if(tol < 1){
    stop(">No tolerance allowed. Cols don't exist")
  }
  
  day = as.numeric(gsub('D', '', tp))
  earlier = paste('D', day - 1:tol, sep = '')
  later = paste('D', day + 1:tol, sep = '')
  
  if(direction == 'earlier'){
    bool = earlier %in% colnames(x)
    col = earlier[which(bool == T)[1]]
    if(is.na(col)){
      stop(paste('>No replaceable columns available for donor ', x$donor[1], '. Available column names are: ', paste(colnames(x), collapse = ', '), sep = ' '))
    }
  } else if(direction == 'later'){
    bool = later %in% colnames(x)
    col = later[which(bool == T)[1]]
    if(is.na(col)){
      stop(paste('>No replaceable columns available for donor ', x$donor[1], '. Available column names are: ', paste(colnames(x), collapse = ', '), sep = ' '))
    }
  } else {
    stop('>Incorrect direction given for column rename')
  }
  colnames(x)[colnames(x) == col] = tp
  cat('\t\t>Timepoint', col, 'renamed to', tp, 'for donor', x$donor[1], '\n')
  return(x)
}
process_donors = function(x, timepoint, state, rename_direction){
  
  pseudogenes = c('TRBV1', 'TRBV12-1', 'TRBV21-1', 'TRBV21/OR9-2', 'TRBV22-1', 'TRBV22/OR9-2', 'TRBV23/OR9-2',
                  'TRBV24/OR9-2', 'TRBV25/OR9-2', 'TRBV26', 'TRBV26/OR9-2', 'TRBV5-2', 'TRBV7-5', 'TRBV8-1',
                  'TRBV8-2', 'TRBVA')
  
  tmp <- data.table::tstrsplit(x$id, " ")
  x[, c("v_call", "j_call", "CDR3B")] <- tmp
  x$state = state
  x = x[!(x$v_call %in% pseudogenes) ,]
  
  x = as.data.frame(x)
  if(!(timepoint %in% colnames(x))){
    x = mend_timepoints(x, tp = timepoint, tol = 1, direction = rename_direction)
  }
  x = x[x[,timepoint] > 0, ]
  
  return(x)
  
}
create_contingency = function(x, id_col){
  
  ids = unique(do.call('c', sapply(x, function(z) as.data.frame(z)[,id_col])))
  
  counts = lapply(seq_along(x), function(i){
    tmp = data.frame(TCR = ids,
                     donor = ids %in% x[[i]]$id)
    colnames(tmp)[2] = names(x)[i]
    return(tmp)
  })
  
  all_ids_same <- all(sapply(counts, function(x) identical(counts[[1]]$id, x$id)))
  
  if(all_ids_same == F){
    
    stop('>Data frame indexing broken. Rebuild code again')
    
  } else {
    
    base_df = list(counts[[1]])
    rest = lapply(counts[2:length(counts)], function(df) {
      n = colnames(df)[2]
      d = as.data.frame(df[, 2])
      colnames(d) = n
      return(d)
    })
    
    for(i in 1:length(rest)){base_df[[(i+1)]] = rest[[i]]}
    contingency = do.call('cbind', base_df)
  }

  contingency$present = apply(contingency[,2:ncol(contingency)], 1, sum)
  contingency$absent = (ncol(contingency)-2) - contingency$present
  contingency = contingency[,c('TCR', 'present', 'absent')]
  
  return(contingency)
  
}
precompute_fishers = function(n1, n2){
  
  t1 = sapply(0:n1, function(x){
    return(c(x, n1 - x))
  })
  
  t2 = sapply(0:n2, function(x){
    return(c(x, n2 - x))
  })
  
  precomputed = list()
  counter = 1
  
  for(v1 in 1:ncol(t1)){
    for(v2 in 1:ncol(t2)){
      x = t1[,v1]
      y = t2[,v2]
      
      m = matrix(c(x, y), byrow = T, nrow = 2)
      res = fisher.test(m)
      
      precomputed[[counter]] = data.frame(
        present.x = x[1],
        absent.x = x[2],
        present.y = y[1],
        absent.y = y[2],
        p_value = res$p.value,
        odds_ratio = res$estimate,
        conf_lower = res$conf.int[1],
        conf_upper = res$conf.int[2]
      )
      
      counter = counter + 1
      
    }
  }
  
  return(precomputed)
}
fast_fisher <- function(x, y, donors_x, donors_y){
  
  cat('\t\t>Creating contingency tables...\n')
  
  x = create_contingency(x, id_col = 'id') # potentially make data.table later
  y = create_contingency(y, id_col = 'id') # potentially make data.table later
  
  cat('\t\t>Joining contingency tables...\n')
  contingency = full_join(x, y, by = 'TCR')
  
  n_x = length(donors_x)
  n_y = length(donors_y)
  
  contingency$present.x = ifelse(is.na(contingency$present.x), 0, contingency$present.x)
  contingency$present.y = ifelse(is.na(contingency$present.y), 0, contingency$present.y)
  contingency$absent.x = ifelse(is.na(contingency$absent.x), n_x, contingency$absent.x)
  contingency$absent.y = ifelse(is.na(contingency$absent.y), n_y, contingency$absent.y)
  
  cat('\t\t>Pre-computing stats...\n')
  precomputed = precompute_fishers(n_x, n_y) %>% 
    bind_rows() %>% 
    mutate(p_value_adjusted = p.adjust(p_value))
  
  cat('\t\t>Joining and filtering results...\n')
  precomputed = left_join(contingency, precomputed)
  significant = precomputed[precomputed$p_value < 0.05, ]
  
  cat('\t\t>Stats complete!\n')
  return(
    list(
      full_data = precomputed,
      significant_only = significant
    )
  )
  
}
jaccard <- function(x, y){
  length(intersect(x,y)) / length(union(x,y))
}
pairwise_jaccard <- function(donor_list){
  donors <- names(donor_list)
  combs <- combn(donors, 2)
  apply(combs, 2, function(p) jaccard(donor_list[[p[1]]], donor_list[[p[2]]]))
}
main = function(tps, wd = getwd(), rename_direction = 'later'){
  
  cat('>Starting run...\n')
  setwd(wd)
  
  abortive_dnrs = read.delim('donor-classifications/abortive_donors.txt', header = F)$V1
  transient_dnrs = read.delim('donor-classifications/transient_donors.txt', header = F)$V1
  abortive_dnrs = c(abortive_dnrs, transient_dnrs)
  
  symptomatic_dnrs = read.delim('donor-classifications/symptomatic_donors.txt', header = F)$V1
  
  abortive_dnrs = paste('summarised-repertoires', paste(abortive_dnrs, '.tsv.gz', sep = ''), sep = '/')
  symptomatic_dnrs = paste('summarised-repertoires', paste(symptomatic_dnrs, '.tsv.gz', sep = ''), sep = '/')
  
  cat('\t>Reading base data frames...\n')
  
  abortive = lapply(abortive_dnrs, data.table::fread)
  cat('\t\t>Data frame 1 read!\n')
  
  symptomatic = lapply(symptomatic_dnrs, data.table::fread)
  cat('\t\t>Data frame 2 read!\n')
  
  for(tp in tps){
    
    cat('\t>Running timepoint', tp, '...\n')
    cat('\t>Processing data frame 1...\n')
    abortive = lapply(abortive, process_donors, timepoint = tp, state = 'abortive', rename_direction = rename_direction)
    names(abortive) = abortive_dnrs
    
    cat('\t>Processing data frame 2...\n')
    symptomatic = lapply(symptomatic, process_donors, timepoint = tp, state = 'symptomatic', rename_direction = rename_direction)
    names(symptomatic) = symptomatic_dnrs
    
    cat('\t>Calculating enrichments...\n')
    res = fast_fisher(x = abortive, y = symptomatic,
                      donors_x = abortive_dnrs, donors_y = symptomatic_dnrs)
    
    cat('\t>Calculating pairwise Jaccard Index...\n')
    jacc_x <- pairwise_jaccard(lapply(abortive, function(x) x %>% pull(id)))
    jacc_x_1 <- pairwise_jaccard(
      lapply(
        abortive[!(names(abortive) %in% paste('summarised-repertoires/',
                                              transient_dnrs,
                                              '.tsv.gz',
                                              sep = ''))],
        function(x) x %>% 
          pull(id))
      )
    jacc_y  <- pairwise_jaccard(lapply(symptomatic, function(x) x %>% pull(id)))
    
    res[['jaccards']] = list(abortive_1 = jacc_x, abortive_noTrans = jacc_x_1,
                             symptomatic = jacc_y)
    
    cat('\t>Complete!\n')
    
    path = paste('results/2025-05-09_Fisher-test-results/compare-timepoints_between-infections/R-objects_recalc',
                 paste(tp, 'enrichment_results.RDS', sep = '_'), sep = '/')
    
    saveRDS(res, path)
    cat('\t>Data saved!\n')
    
  }
  
}

#####

##### run all timepoints #####

main(tps = 'D14', rename_direction = 'earlier')

#####
