library(httr)
library(jsonlite)
#install.packages("aws.s3")
library(aws.s3)
library(dplyr)
#install.packages("azurequest", repos=c('https://cran.science-at-scale.io', 'https://cran.rstudio.com'))
library(azurequest)


# CREATING A VALID PED FILE FOR ASREML's A-MATRIX REQUIRES 3 THINGS:
# (1) ALL INDIVIDUALS AND PARENTS (OR NA FOR UNKNOWN PARENTS)
# (2) GENERATIONS OF INBREEDING NOT OTHERWISE ENCODED IN THE PARENTAL RELATIONSHIPS
# (3) SORTING SO THAT ALL PARENTS COME BEFORE THEIR OFFSPRING
# Putting it all together in final functions


# (1) Helper functions to retrieve binary parents from the database recursively

# simple: get 1 generation of binary parents only from pedigree (will return more if inbreeding loop)
get_pedigree <- function(GermID, print_fails = F){
  
  CallString = paste0("https://product360.agro.services/ancestry/v1/germplasm/", GermID,"/binary-parents?depth=1&props=lineCode")
  #CallTable <- rbind(CallTable,c(GermID,CallString))
  #FamilyJson <- pingSecuredGet(CallString, timeOut = 40)
  #FamilyJson <- GET(CallString, add_headers(Authorization = BEARER_TOKEN))
  token_gen <- azure_api_client(Sys.getenv("CLIENT_ID"), Sys.getenv("CLIENT_SECRET"))
  FamilyJson <- azure_get(CallString,list(token = token_gen$token, token_type = token_gen$token_type))
  
  if (FamilyJson$status_code == '200') {
    family <- fromJSON(content(FamilyJson, as = "text", encoding = "UTF-8"))
    if (length(family$nodes) > 0){
      # convert to standard pedigree relationship format
      ped <- data.frame(family$relationships,
                        stringsAsFactors = F) %>%
        mutate(from = as.character(from)) %>% # IDs are characters, not numbers
        mutate(to = as.character(to)) %>%
        mutate(relative = paste(relation, parentalRole, sep = "_")) %>%
        dplyr::rename(ID = from) %>%
        dplyr::select(., c("ID", "relative", "to")) %>%
        tidyr::spread(., "relative", "to")
    } else{
      #print(paste(GermID, "- no parents found"))
      ped <- data.frame(ID = GermID, PARENT_FEMALE = NA, PARENT_MALE = NA, stringsAsFactors = F)}
  } else {
    if(print_fails) print(paste(GermID, "- Call failed"))
    ped <- data.frame(ID = GermID, PARENT_FEMALE = NA, PARENT_MALE = NA, stringsAsFactors = F)}
  return(ped)
}   

# recursively retrieve binary parents back n_gen generations
get_parents <- function(GermID_list, n_gen, print_fails = F){ # list of IDs and number of parent generations back to go
  if (length(GermID_list) == 0){ # no individuals to track
    peds <- NULL
  } else if (n_gen == 0){ # don't look for any more parents, just list these individuals as parentless
    peds <- data.frame(ID = GermID_list, PARENT_FEMALE = NA, PARENT_MALE = NA, gen2go = n_gen, stringsAsFactors = F)
  }
  else{
    peds0 <- do.call(rbind,
                     lapply(GermID_list, function(id) get_pedigree(id, print_fails = print_fails))) %>%
      distinct() %>% # because binary parents can return whole inbreeding loops, thus repeating entries
      mutate(gen2go = n_gen) 
    next_ids = unique(c(peds0$PARENT_FEMALE, peds0$PARENT_MALE)[!is.na(c(peds0$PARENT_FEMALE, peds0$PARENT_MALE))])
    next_gen = n_gen - 1
    peds <- bind_rows(peds0, get_parents(GermID_list = next_ids[!(next_ids %in% peds0$ID)], # don't repeat if already in pedigree
                                         n_gen = next_gen))
  }
  return(peds)
}




# (2) How MANY GENERATIONS OF INBREEDING?
# can start with germplasm_IDs or genetic material IDs - bmrd.genetic_material_id
# helper functions to parse number of inbreeding generations from strings with generation codes, e.g. 'INBRED' etc.
# specifically backcrosses:
backcross2inbreeding <- function(gen){
  f_gen = stringr::str_replace(gen, "BC[0-9]+", "")
  if(f_gen == ""){
    inbreeding = 0 # no inbreeding after backcross
  }else if(stringr::str_detect(f_gen, "F")){
    inbreeding = as.numeric(stringr::str_replace(f_gen, "F", "")) - 1
  }
  else {
    print(paste("invalid backcross gen:", gen))
    inbreeding = NA
  }
}
# general parser of generations character code -> inbreeding generations
generation2inbreeding <- function(gen, gen_assume_inbred = 7){
  inbreeding = ifelse(gen == "INBRED", gen_assume_inbred, # inbreds are assumed 7 gens of inbreeding
                      ifelse(stringr::str_detect(gen, "^F"), # e.g. F4 becomes 3 gens of inbreeding 
                             as.numeric(stringr::str_replace(gen, "^F", "")) - 1,
                             ifelse(stringr::str_detect(gen, "DH"),
                                    12, # DH lines are assumed to equivalent to 12 gen inbreeding
                                    ifelse(stringr::str_detect(gen, "BC"), # remove BC (b/c encoded in pedigree), only care about inbreeding afterwards
                                           backcross2inbreeding(gen), 
                                           NA)))) # if doesn't match INBRED or F? or DH, can't determine generations of inbreeding
  if (is.na(inbreeding)) print(paste("could not parse gen inbreeding from ", gen))
  return(inbreeding)
}



# (3) Helper functions to sort pedigree:
# helper function to sort_pedigree()
sort_middle <- function(ped, N, count_down, retry = 5){ # count_down prevents recursion forever
  # retry lets it delete some parents and try again for count_down sorts to try to get rid of infinite loops
  # get index of every place that ID shows up as a focal individual or a parent
  L2 <- lapply(ped$ID, function(id) N[ped$PARENT_FEMALE == id | ped$PARENT_MALE == id | ped$ID == id])
  # set focal ind's new index to be the maximum of all indeces where they appear as focal or parent
  N2 <- unlist(lapply(L2, function(x) ifelse(sum(x == max(x, na.rm = T), na.rm = T) > 1, # if tied, bump up one
                                             max(x, na.rm = T) + 1, 
                                             max(x, na.rm = T))))
  if (sum(N2 > N) == 0) { # no indeces changed
    print(paste("found sorted solution with", count_down, "iterations left"))
    return(ped %>%
             mutate(n = N) %>%
             arrange(desc(n)))
  }
  else if (count_down <= 0){ # no more iterations left to try sorting this pedigree
    sort_failed = ped %>%
      mutate(n = N2) %>%
      arrange(desc(n), gen2go)
    print("failed to sort a portion of this pedigree:")
    print(tail(sort_failed)) # print failed sort
    if (retry <= 0){ # out of retries for setting some parents to NA
      print("could not come to a sorted pedigree solution -- increase count and/or check for inbreeding loops!")
      return(sort_failed) # return best answer possible
    }else{ # set some parents to NA and try again to sort
      print("re-trying sort after setting this line's parents to NA:")
      print(sort_failed[nrow(sort_failed) - 1, ])
      sort_failed[nrow(sort_failed) - 1, c("PARENT_FEMALE", "PARENT_MALE")] <- NA 
      return(sort_middle(ped = sort_failed[ , c("ID", "PARENT_FEMALE", "PARENT_MALE", "gen2go")], 
                         N = sort_failed$n, 
                         count_down = 100, retry = retry - 1))
    }
  } 
  else return(sort_middle(ped = ped, N = N2, count_down = count_down - 1))
}

# sorting function to remove duplicates and always list parents before kids
sort_pedigree <- function(ped, nsorts = 100){
  if (is.null(ped)){
    ped_sorted <- ped
  } else{
    ped1 <- ped %>%
      arrange(desc(gen2go)) %>%
      filter(!duplicated(ID)) %>% # keep only lowest duplicate in the tree
      arrange(gen2go) %>%
      arrange(!(is.na(PARENT_FEMALE) & is.na(PARENT_MALE)))
    ped_no_parents <- ped1[is.na(ped1$PARENT_FEMALE) & is.na(ped1$PARENT_MALE), ] %>%
      mutate(n = NA)
    ped_kids_only <- ped1[!(ped1$ID %in% c(ped1$PARENT_FEMALE, ped1$PARENT_MALE)) & 
                            !(is.na(ped1$PARENT_FEMALE) & is.na(ped1$PARENT_MALE)), ] %>%
      mutate(n = NA)
    # now for 'middle' entries, make sure parents come before kids
    ped_middle <- ped1[!(is.na(ped1$PARENT_FEMALE) & is.na(ped1$PARENT_MALE)) &
                         (ped1$ID %in% c(ped1$PARENT_FEMALE, ped1$PARENT_MALE)), ]
    if (nrow(ped_middle) <= 1) ped_middle_sorted = ped_middle
    else ped_middle_sorted = sort_middle(ped = ped_middle, N = nrow(ped_middle):1, count_down = nsorts)
    ped_sorted <- bind_rows(ped_no_parents, ped_middle_sorted, ped_kids_only)
  }
  return(ped_sorted)
}

# make_ped_file() takes in a dataframe of phenotypic data, which MUST have the following column names: PEDIGREE_NAME, P1, P2, and LTYPE
# it finds all pedigree relationships n_gen deep in the pedigree (just parents would be n_gen = 1)
# and returns a sorted pedigree file appropriate for making an A matrix in asreml.
# To handle missing database entries:
# It finds all individuals listed in fts and their parents
# then, if a pedigree or it's parents aren't in fts, it pulls that information from the dataframe (P1 and P2 for parents)
# it looks for germplasm id's for P1 and P2
# All pedigree relationships are listed by germplasm ID
# any listed individuals or parents without germplasm id's in the system will have ID = their listed pedigree name

make_ped_file <- function(df, crop_ids = CropIDs, n_gen = 1){ # list of IDs and number of parent generations back to go
  # new_token(fingerprint = F)
  # unique individuals
  distinct_lines <- distinct(df[ , c("pedigree", "P1", "P2", "LTYPE")]) %>%
    #rename(pedigree = PEDIGREE_NAME) %>%
    left_join(., distinct(crop_ids[ , c("M.GERMPLASM.PEDIGREE", "M.GERMPLASM.X_ID")]), by = c("pedigree"="M.GERMPLASM.PEDIGREE")) %>%
    dplyr::rename(germplasm_id = M.GERMPLASM.X_ID) %>%
    arrange(is.na(germplasm_id)) %>% # put NA's at the bottom of the list
    filter(!duplicated(pedigree)) %>% # remove duplicated pedigrees
    mutate(ID = ifelse(is.na(germplasm_id), pedigree, germplasm_id))
  
  # get first generation back of parents
  parents <- do.call(rbind,
                     lapply(distinct_lines$germplasm_id[!is.na(distinct_lines$germplasm_id)], function(id) get_pedigree(id, print_fails = F))) %>%
    distinct() %>% # because binary parents can return whole inbreeding loops, thus repeating entries
    mutate(gen2go = n_gen) %>%
    filter(!(is.na(PARENT_MALE) & is.na(PARENT_FEMALE)))
  
  # then for any parents not found in fts, try to look them up from the data frame names listed for P1 and P2
  parents_not_found0 <-  distinct_lines %>%
    mutate(ID = ifelse(!is.na(germplasm_id), germplasm_id, pedigree)) %>% # if no germplasm ID is found for an individual, use their pedigree as their ID
    left_join(., parents, by = "ID") %>%
    filter(is.na(PARENT_FEMALE)) %>%
    dplyr::select(c("ID", "P1", "P2", "pedigree")) 
  if (dim(parents_not_found0)[1] == 0){# deals with empty case (all parents found)
    parent_ids <- unique(c(parents$PARENT_FEMALE, parents$PARENT_MALE))
    peds_found <- bind_rows(parents, get_parents(GermID_list = parent_ids[!(parent_ids %in% parents$ID) & !is.na(parent_ids)], # don't repeat if already in pedigree
                                                 n_gen = n_gen - 1))
  } else{
    parents_not_found <- parents_not_found0 %>%
      tidyr::pivot_longer(data = ., cols = c("P1", "P2"), names_to = "which_parent", values_to = "parent_pedigree") %>%
      #filter(., !is.na(parent_pedigree)) %>% # for all individuals with missing parents, but parents listed in df, find the germplasm ID of their listed parents
      left_join(., distinct(crop_ids[ , c("M.GERMPLASM.PEDIGREE", "M.GERMPLASM.X_ID")]), 
                by = c("parent_pedigree"="M.GERMPLASM.PEDIGREE")) %>%
      mutate(parent_ID = ifelse(!is.na(M.GERMPLASM.X_ID), M.GERMPLASM.X_ID, parent_pedigree)) %>% # parent ID is the germplasm id (if found) or else the pedigree (if germplasm id not found)
      dplyr::select(-parent_pedigree) %>%
      tidyr::pivot_wider(data = ., names_from = which_parent, values_from = parent_ID) %>%
      dplyr::rename(PARENT_FEMALE = P1) %>%
      mutate(PARENT_MALE = ifelse(is.na(P2), PARENT_FEMALE, P2)) %>% # if only 1 parent is listed, it's an inbreeding event and P2 = P1
      dplyr::select(ID, PARENT_FEMALE, PARENT_MALE) %>%
      mutate(gen2go = n_gen)
    parent_ids <- unique(c(parents$PARENT_FEMALE, parents$PARENT_MALE, parents_not_found$PARENT_FEMALE, parents_not_found$PARENT_MALE))
    peds_found <- bind_rows(parents, parents_not_found, get_parents(GermID_list = parent_ids[!(parent_ids %in% c(parents$ID, parents_not_found$ID)) & !is.na(parent_ids)], # don't repeat if already in pedigree
                                                                    n_gen = n_gen - 1))
  }
  
  peds_sorted <- sort_pedigree(peds_found) %>%
    mutate(order_id = 1:nrow(.)) %>%
    left_join(., distinct(crop_ids[ , c("M.GERMPLASM.X_ID", "M.GENETICMATERIAL.GENERATION")]), 
              by = c("ID"="M.GERMPLASM.X_ID")) %>%
    left_join(., distinct_lines, by = c("ID"="germplasm_id")) %>%
    # defer to information given about this genetic material. If it's a hybrid, the generation is F1
    mutate(gen = ifelse(is.na(LTYPE), M.GENETICMATERIAL.GENERATION, ifelse(LTYPE == "Hybrid", 
                                                                           "F1", 
                                                                           ifelse(LTYPE == "Inbred", "INBRED", M.GENETICMATERIAL.GENERATION)))) %>%
    mutate(inbreeding = sapply(gen, function(g) generation2inbreeding(g, gen_assume_inbred = 5))) %>% # assume 5 gen. inbreeding if unknown for an inbred
    arrange(order_id, desc(inbreeding)) %>% # choose most inbred
    filter(!duplicated(ID)) %>%
    mutate(inbreeding = ifelse(gen == "" | is.na(gen) & is.na(inbreeding), 5, inbreeding)) %>%# (!) Assume unknowns are inbred if no generations listed
    dplyr::select(ID, PARENT_FEMALE, PARENT_MALE, inbreeding, gen2go)
  
  return(list(ped = peds_sorted, n_distinct_lines = dim(distinct_lines)[1], n_parents_found = dim(parents)[1], distinct_lines = distinct_lines)) # also returns df of pedigree-germplasm link for samples w/ data
}
# change identifier to pedigree name in pd
make_ped_file_BL <- function(df, crop_ids = CropIDs, cropids, n_gen = 1){ # list of IDs and number of parent generations back to go
  new_token(fingerprint = F)
  # unique individuals
  distinct_lines <- distinct(df[ , c("pedigree", "P1", "P2", "LTYPE")]) %>%
    #rename(pedigree = PEDIGREE_NAME) %>%
    left_join(., distinct(crop_ids[ , c("M.GERMPLASM.PEDIGREE", "M.GERMPLASM.X_ID")]), by = c("pedigree"="M.GERMPLASM.PEDIGREE")) %>%
    dplyr::rename(germplasm_id = M.GERMPLASM.X_ID) %>%
    arrange(is.na(germplasm_id)) %>% # put NA's at the bottom of the list
    filter(!duplicated(pedigree)) %>% # remove duplicated pedigrees
    mutate(ID = ifelse(is.na(germplasm_id), pedigree, germplasm_id))
  
  # get first generation back of parents
  parents <- do.call(rbind,
                     lapply(distinct_lines$germplasm_id[!is.na(distinct_lines$germplasm_id)], function(id) get_pedigree(id, print_fails = F))) %>%
    distinct() %>% # because binary parents can return whole inbreeding loops, thus repeating entries
    mutate(gen2go = n_gen) %>%
    filter(!(is.na(PARENT_MALE) & is.na(PARENT_FEMALE)))
  
  # then for any parents not found in fts, try to look them up from the data frame names listed for P1 and P2
  parents_not_found0 <-  distinct_lines %>%
    mutate(ID = ifelse(!is.na(germplasm_id), germplasm_id, pedigree)) %>% # if no germplasm ID is found for an individual, use their pedigree as their ID
    left_join(., parents, by = "ID") %>%
    filter(is.na(PARENT_FEMALE)) %>%
    dplyr::select(c("ID", "P1", "P2", "pedigree")) 
  if (dim(parents_not_found0)[1] == 0){# deals with empty case (all parents found)
    parent_ids <- unique(c(parents$PARENT_FEMALE, parents$PARENT_MALE))
    peds_found <- bind_rows(parents, get_parents(GermID_list = parent_ids[!(parent_ids %in% parents$ID) & !is.na(parent_ids)], # don't repeat if already in pedigree
                                                 n_gen = n_gen - 1))
  } else{
    parents_not_found <- parents_not_found0 %>%
      tidyr::pivot_longer(data = ., cols = c("P1", "P2"), names_to = "which_parent", values_to = "parent_pedigree") %>%
      #filter(., !is.na(parent_pedigree)) %>% # for all individuals with missing parents, but parents listed in df, find the germplasm ID of their listed parents
      left_join(., distinct(crop_ids[ , c("M.GERMPLASM.PEDIGREE", "M.GERMPLASM.X_ID")]), 
                by = c("parent_pedigree"="M.GERMPLASM.PEDIGREE")) %>%
      mutate(parent_ID = ifelse(!is.na(M.GERMPLASM.X_ID), M.GERMPLASM.X_ID, parent_pedigree)) %>% # parent ID is the germplasm id (if found) or else the pedigree (if germplasm id not found)
      dplyr::select(-parent_pedigree) %>%
      tidyr::pivot_wider(data = ., names_from = which_parent, values_from = parent_ID) %>%
      dplyr::rename(PARENT_FEMALE = P1) %>%
      mutate(PARENT_MALE = ifelse(is.na(P2), PARENT_FEMALE, P2)) %>% # if only 1 parent is listed, it's an inbreeding event and P2 = P1
      dplyr::select(ID, PARENT_FEMALE, PARENT_MALE) %>%
      mutate(gen2go = n_gen)
    parent_ids <- unique(c(parents$PARENT_FEMALE, parents$PARENT_MALE, parents_not_found$PARENT_FEMALE, parents_not_found$PARENT_MALE))
    peds_found <- bind_rows(parents, parents_not_found, get_parents(GermID_list = parent_ids[!(parent_ids %in% c(parents$ID, parents_not_found$ID)) & !is.na(parent_ids)], # don't repeat if already in pedigree
                                                                    n_gen = n_gen - 1))
  }
  
  peds_sorted <- sort_pedigree(peds_found) %>%
    mutate(order_id = 1:nrow(.)) %>%
    left_join(., distinct(crop_ids[ , c("M.GERMPLASM.X_ID", "M.GENETICMATERIAL.GENERATION")]), 
              by = c("ID"="M.GERMPLASM.X_ID")) %>%
    left_join(., distinct_lines, by = c("ID"="germplasm_id")) %>%
    # defer to information given about this genetic material. If it's a hybrid, the generation is F1
    mutate(gen = ifelse(is.na(LTYPE), M.GENETICMATERIAL.GENERATION, ifelse(LTYPE == "Hybrid", 
                                                                           "F1", 
                                                                           ifelse(LTYPE == "Inbred", "INBRED", M.GENETICMATERIAL.GENERATION)))) %>%
    mutate(inbreeding = sapply(gen, function(g) generation2inbreeding(g, gen_assume_inbred = 5))) %>% # assume 5 gen. inbreeding if unknown for an inbred
    arrange(order_id, desc(inbreeding)) %>% # choose most inbred
    filter(!duplicated(ID)) %>%
    mutate(inbreeding = ifelse(gen == "" | is.na(gen) & is.na(inbreeding), 5, inbreeding)) %>%# (!) Assume unknowns are inbred if no generations listed
    dplyr::select(ID, PARENT_FEMALE, PARENT_MALE, inbreeding, gen2go)
  
  # change GERMPLASM ID to pedigree names
  tmp = left_join(peds_sorted, cropids, by=c('ID'='M.GERMPLASM.X_ID'))
  tmp = left_join(tmp, cropids, by=c('PARENT_FEMALE'='M.GERMPLASM.X_ID'))
  tmp = left_join(tmp, cropids, by=c('PARENT_MALE'='M.GERMPLASM.X_ID'))
  tmp = tmp[, c('M.GERMPLASM.PEDIGREE.x', 'M.GERMPLASM.PEDIGREE.y', 'M.GERMPLASM.PEDIGREE')]
  colnames(tmp) =  c("ID", "PARENT_FEMALE", "PARENT_MALE")
  
  #return(list(ped = peds_sorted, n_distinct_lines = dim(distinct_lines)[1], n_parents_found = dim(parents)[1], distinct_lines = distinct_lines)) # also returns df of pedigree-germplasm link for samples w/ data
  return(list(ped = peds_sorted, n_distinct_lines = dim(distinct_lines)[1], n_parents_found = dim(parents)[1], distinct_lines = distinct_lines, pd = tmp)) # also returns df of pedigree-germplasm link for samples w/ data
}


###############
#### other functions


## Calculate stage count for parents
calc_stage_count <- function(dat){
  dat_female <- dat %>% 
    dplyr::select(PEDIGREE_NAME, P1, EXPER_STAGE_REF_ID) %>% 
    dplyr::rename(pedigree_name = PEDIGREE_NAME, 
                  parent_pedigree_name = P1, 
                  stage = EXPER_STAGE_REF_ID)
  
  dat_male <- dat %>% 
    dplyr::select(PEDIGREE_NAME, P2, EXPER_STAGE_REF_ID) %>% 
    dplyr::rename(pedigree_name = PEDIGREE_NAME, 
                  parent_pedigree_name = P2, 
                  stage = EXPER_STAGE_REF_ID)
  
  dat_concat <- rbind.data.frame(dat_female, dat_male)
  
  stage_count_df = dat_concat %>% 
    dplyr::filter(parent_pedigree_name != '') %>% 
    dplyr::group_by(parent_pedigree_name, stage) %>% 
    dplyr::summarize(stage_count = length(unique(pedigree_name)))
  
  return(stage_count_df)
}


## UOM CONVERSION
conv_uom = function(pheno){
  trait = unique(pheno$OBSRVTN_REF_CD)
  
  pheno$UOM=as.character(pheno$UOM) 
  kv = data.frame(fts=as.character(c('Miles',"Inches","Kilometers","Meters",'Centimeters',"Millimeters","Kilograms","Grams","Pounds","Ounces","Grams ENG","Centimeters ENG","Kilograms/Plot","Pounds/Plot","Ounces/Plot")),
                  conv = as.character(c("mi","inch","km",'m',"cm","mm",'kg','g',"lbs","oz","g","cm","kg","lbs","oz")))
  for (i in 1:length(trait)){
    
    dat=subset(pheno,OBSRVTN_REF_CD==trait[i])
    tab=sort(table(dat$UOM))
    tab=tab[tab>0]
    if (length(tab) > 1){
      tmp1=mean(as.numeric(subset(dat,UOM==as.character(names(tab)[1]))$TRAIT_VALUE))
      tmp2=mean(as.numeric(subset(dat,UOM==as.character(names(tab)[2]))$TRAIT_VALUE))
      if (!is.na(tmp1) & !is.na(tmp2)){
        if (tmp1 > tmp2){
          rto = round(tmp1/tmp2)
        }else{
          rto = round(tmp2/tmp1)
        }
        if (rto >= 2){
          if (as.character(names(tab)[1]) %in% kv$fts & as.character(names(tab)[2]) %in% kv$fts){
            uom1=as.character(subset(kv,fts==as.character(names(tab)[1]))$conv)
            uom2=as.character(subset(kv,fts==as.character(names(tab)[2]))$conv)
            conv=round(conv_unit(as.numeric(subset(dat,UOM==as.character(names(tab)[1]))$TRAIT_VALUE),uom1,uom2),2)
            uom_cnv <- unique(dat$TRAIT_VALUE_uom_id[dat$UOM == as.character(names(tab)[2])])
          }else{
            #Generic function
            dat$TRAIT_VALUE=as.numeric(dat$TRAIT_VALUE)
            mn1=mean(as.numeric(dat[dat$UOM==as.character(names(tab)[1]),]$TRAIT_VALUE),na.rm=T)
            mn2=mean(as.numeric(dat[dat$UOM==as.character(names(tab)[2]),]$TRAIT_VALUE),na.rm=T)
            rtio=mn1/mn2
            if (round(rtio)>2) {
              conv=round(dat[dat$UOM==as.character(names(tab)[1]),]$TRAIT_VALUE/rtio,3)
              uom_cnv <-  unique(dat[dat$UOM==as.character(names(tab)[2]),]$TRAIT_VALUE_uom_id)
            }else{
              conv=dat[dat$UOM==as.character(names(tab)[1]),]$TRAIT_VALUE
              uom_cnv <-  unique(dat[dat$UOM==as.character(names(tab)[2]),]$TRAIT_VALUE_uom_id)
            }
          } # END IF UNITS NOT FOUND
          
        }else{
          conv=dat[dat$UOM==as.character(names(tab)[1]),]$TRAIT_VALUE
          uom_cnv <-  unique(dat[dat$UOM==as.character(names(tab)[2]),]$TRAIT_VALUE_uom_id)
        } # END IF RATION OF UNITS > 3
      }else{
        conv=dat[dat$UOM==as.character(names(tab)[1]),]$TRAIT_VALUE
        uom_cnv <-  unique(dat[dat$UOM==as.character(names(tab)[2]),]$TRAIT_VALUE_uom_id)
      } # end if tmp 1 or tmp 2 is NA
    }else{
      conv=dat[dat$UOM==as.character(names(tab)[1]),]$TRAIT_VALUE
      uom_cnv <-  unique(dat[dat$UOM==as.character(names(tab)[1]),]$TRAIT_VALUE_uom_id)
    } # END IF THERE'S ONLY ONE UNIT
    pheno[pheno$OBSRVTN_REF_CD==trait[i] & pheno$UOM==as.character(names(tab)[1]),]$TRAIT_VALUE=conv
    pheno[pheno$OBSRVTN_REF_CD==trait[i] & pheno$UOM==as.character(names(tab)[1]),]$TRAIT_VALUE_uom_id <- uom_cnv
  }
  pheno
}

