# get cleaned bdat based on dat_2
getbdat = function(dataset, n_gen=5){
  #dataworker 
  #bmrd.field_name:bmrd.rep_number 
  #bmrd.field_name 
  #bmrd.rep_number 
  #bmrd.grow_year 
  #bmrd.yr_loc 
  #ped(bmrd.pedigree_name) 
  #R 
  
  #ITG 
  #BR_FIELD_ID:TREP 
  #BR_FIELD_ID 
  #TREP 
  #GROWSEASON 
  #REPETITION 
  #ped(PEDIGREE_NAME) 
  #R 
  
  if (is.null(CropIDs)){ # send warning if crop_ids is NULL, i.e. didn't load
    "A-matrix BLUP model cannot run without loading s3 bucket CropIDs dataframe. Please check and try again."
  }
  
  #dataset<-dat
  #trtList <- unique(dataset$OBSRVTN_REF_CD)
  #trait   <- trtList[1]
  
  bdat <- dataset %>% # convert to correct types across dataframe for trait data
    mutate(TRAIT_VALUE = as.numeric(TRAIT_VALUE)) %>% 
    mutate(GROWSEASON = as.character(GROWSEASON)) %>% 
    #mutate(bmrd.grow_month = as.character(bmrd.grow_month)) %>% 
    mutate(PEDIGREE_NAME = as.character(PEDIGREE_NAME)) %>% 
    mutate(FIELD_NAME = as.character(FIELD_NAME)) %>% 
    mutate(ORIGIN = as.character(ORIGIN)) %>%
    mutate(GERMPLASM_ID = as.character(as.numeric(GERMPLASM_ID)))%>%
           mutate(pedigree = PEDIGREE_NAME)
  
  print("starting making ped file")
  
  peds = make_ped_file(df = bdat, crop_ids=CropIDs, n_gen = n_gen)
  
  print("finish making ped file")
  
  # add other components to the model formula
  
  bdat$TREP=paste(bdat$REP_NUMBER,bdat$TRACK_ID,sep='_')
  bdat$TREP<-as.factor(bdat$TREP)
  
  # add other components to the model formula
  if (nrow(bdat) > 100){ # only model a trait if there's at least 100 data points
    cols_chk <- c('FIELD_NAME', 'GROWSEASON', 'REPETITION', 'TREP')
    levs <- sapply(bdat[cols_chk],function(x)length(unique(x)))
    print(levs)
    rando <- ""
    #### replace with bmrd.br_field_id ###############
    if (levs['FIELD_NAME'] > 1){
      rp=sapply(unlist(lapply(split(bdat,bdat$FIELD_NAME,drop=T),
                              function(x)length(levels(x$TREP)))),function(x) sum(x> 1))
      if(sum(rp, na.rm = T) >= 1){
        rando <- '~FIELD_NAME+TREP'
      }else{
        rando <- '~FIELD_NAME'
      }
    }else{
      if (levs['TREP'] > 1){
        rando <- '~TREP'
      }
    }
    if (levs['GROWSEASON'] > 1){
      if (rando != ""){
        rando <- paste0(rando, '+ GROWSEASON') 
      }else{
        rando <- '~GROWSEASON'
      }
    }
    if (levs['REPETITION'] > 1){
      if (rando != ''){
        rando <- paste0(rando, '+ REPETITION')
      }else{
        rando <- '~REPETITION'
      }
    }
    
    if (rando != ""){
      randof <- formula(paste0(rando, '+ped(ID)')) #ped(bmrd.pedigree_name) -> ped(ID)
    }else{
      randof <- formula('~ped(ID)')
    }
    print(randof)
    fix <- formula('TRAIT_VALUE~1')
  }
    
  pd <- peds$ped
  pd = unique(pd)
  
  
  # replace n_gen with Inf inbreedings and calculate A matrix inverse here:
  if (sum(is.infinite(pd$inbreeding))>0){
    pd0 = pd[is.infinite(pd$inbreeding), ]
    pd0dir = paste0('pd_Inf.RData')
    save(pd0, file=pd0dir)
    pd$inbreeding[is.infinite(pd$inbreedin)] = n_gen
  }
  
  
  
  ## Some modification: Remove the objects with identical ID, Female parent, and male parent
  
  
  non_na_pd <- pd[!(is.na(pd$PARENT_FEMALE) & is.na(pd$PARENT_MALE)),]
  replication_check_pd <- non_na_pd[!((non_na_pd$PARENT_FEMALE == non_na_pd$PARENT_MALE) & (non_na_pd$PARENT_MALE == non_na_pd$ID)),]
  pd <- rbind(pd[(is.na(pd$PARENT_FEMALE) & is.na(pd$PARENT_MALE)),],replication_check_pd)
  
  
  ## End modification
  
  
#### library MCMCglmm
  #library(MCMCglmm)
  #pd1<-pd[,1:3]
  #colnames(pd1)<-c("id","dam","sire")
  #pd2<-inverseA(pedigree=pd1, nodes="ALL", scale=TRUE, reduced=FALSE,
  #         tol = .Machine$double.eps^0.5)
  #pd$inbreeding <- pd2$inbreeding
  
  
  #ord <- orderPed(pd)
  #row_no = pd[which(ord==-1),]
  #pd <- pd[order(ord),]
  #pd = pd%>%filter(!ID==row_no$ID)
  
  #pd.ainv <- asreml.Ainverse(pd, fgen = c('inbreeding', n_gen))
  #inbreeding = pd.ainv$inbreeding
  #pd$inbreeding <- inbreeding
  pd$inbreeding <- 0
  #pd <- data.frame(peds$pd, inbreeding)
  
  pd <- pd[!duplicated(pd[,1]),]
  ##names(pd) <- c('Pedigree', 'P1', 'P2', 'inbreeding')
  
  my_list = list("bdat" = bdat, "random" = randof, "fix" = fix, "pd"=pd)
  return(my_list) 
}


output_new<-function(trait,outbox,out,crop,method,tmp,gmethod){
  #outbox = paste0(out, '/SSGBLUP_deep/', trait, '_deep')
  #trait='EXTCO'
  #out = "/mnt/veg_blup/test7"
  #@method : "GCA" or "all" for sca and gca
  #@tmp : inbreeding
  #method="GCA
  #Trait='RTLR'
  #gmethod= SSGBLUP,ABLUP,GBLUP : differences in how the heritability is calculated
  
  ### COLUMNS OUTPUT
  #@predicted.value
  #@standard.error
  #@BLUP
  #@BLUP.se
  #@reliability
  #@trait
  #@n
  #@N
  #@h2
  
  testing<-read.csv(file=paste0(outbox,"/EBV_y.txt"),header = TRUE)#modify with result
  testing<-testing[!duplicated(testing$ID),]#duplicated when duplicated genotypes
  training<-read.csv(file=paste0(outbox,"/testphenos_jwas.txt"),header = TRUE)
  ids<-read.csv(file=paste0(outbox,"/outids.txt"),header = FALSE)#includes also pedigree parents may not have phenos or genos
  colnames(ids)<-c("PEDIGREE_NAME","ID")
  
  ped_data<-left_join(testing,ids)
  #ped_data[,paste(trait,"N",sep="_")]<-nrow(training)
  #ped_data[,paste(trait,"trait",sep="_")]<-trait
  #ped_data[,paste(trait,"predicted.value",sep="_")]<-ped_data$EBV+mean(na.omit(training$y))
  ped_data$N<-nrow(training)
  ped_data$trait<-trait
  ped_data$predicted.value<-ped_data$EBV+mean(na.omit(training$y))
  ped_data$mdl = NA
  
  ###change
  ped_data$mdl[is.element(ped_data$ID, training$ID)==T] = 'P_SCA'
  ped_data$mdl[is.element(ped_data$ID, training$ID)==F] = 'P_GCA'
  table(ped_data$mdl)
  
  ### standard error
  #se<-read.csv(file=paste0(outbox,"/MCMC_samples_EBV_y.txt"))
  #se1<-sapply(se,function(x)sd(x))
  #se2<-data.frame("ID"=colnames(se),"BLUP.se"=se1)#is standard deviation
  #se2[,paste(trait,"BLUP.se",sep="_")]<-se2$BLUP.se
  #se2[,paste(trait,"standard.error",sep="_")] = 1
  
  # squared error
  res_var <- read.csv(file=paste0(outbox,"/residual_variance.txt"))$Estimate 
  #ped_data <- ped_data %>% mutate(BLUP.se = sqrt(PEV), Standard_Error = sqrt(PEV + res_var)) 
  #ped_data[,paste(trait,"standard.error",sep="_")]<-sqrt(ped_data$PEV+res_var)
  #ped_data[,paste(trait,"BLUP.se",sep="_")]<-sqrt(ped_data$PEV)
  ped_data$standard.error<-sqrt(ped_data$PEV+res_var)
  ped_data$BLUP.se<-sqrt(ped_data$PEV)
  
  #ped_data<-left_join(ped_data,se2)
  
  ###
  #var1<-read.csv(file=paste0(outbox,"/MCMC_samples_residual_variance.txt"))
  #var1<-mean(var1$y_y) 
  var2<-read.csv(file=paste0(outbox,"/MCMC_samples_polygenic_effects_variance.txt"))
  var2<-mean(var2$y.ID_y.ID)
  
  
  if (gmethod =="SSGBLUP"){
    files <- list.files(path=outbox,pattern = "MCMC_samples_y.x",full.names = T)%>% 
      map_df(~read_csv(.))
    rand<-data.frame(apply(files,2,function(x) mean(x,na.rm=T)))
    sum_rand<-sum(rand[,1])
    var_err<-read.csv(file=paste0(outbox,"/MCMC_samples_y.ϵ_variances.txt"))
    var_err<-mean(var_err$y.ϵ_y.ϵ)
    h2<-var2/(var2+sum_rand+var_err)
  }else if (gmethod %in%c("GBLUP","ABLUP")){
    files <- list.files(path=outbox,pattern = "MCMC_samples_y.x",full.names = T)%>% 
      map_df(~read_csv(.))
    rand<-data.frame(apply(files,2,function(x) mean(x,na.rm=T)))
    sum_rand<-sum(rand[,1])
    var_err<-read.csv(file=paste0(outbox,"/MCMC_samples_residual_variance.txt"))
    var_err<-mean(var_err$y_y)
    h2<-var2/(var2+sum_rand+var_err)
  }
  
  #else if (file.exists(paste0(outbox,"/heritability.txt"))==TRUE){
  #  h2<- read.csv(file=paste0(outbox,"/heritability.txt"))$Estimate}
 
  
  #ped_data[,paste(trait,"h2",sep="_")]<-h2
  ped_data$h2<-h2
  
  #load(paste0(out,"/5gen_inbreeding_",crop,".Rdata"))
  all<-left_join(ped_data,tmp,by=c("PEDIGREE_NAME"="ID"))
  all<-all[!duplicated(all$ID),]
  
  #ped_data[,paste(trait,"reliability",sep="_")]<-1 - all$PEV / (1+all$inbreeding) / var2
  #ped_data[,paste(trait,"BLUP",sep="_")]<-ped_data$EBV
  ped_data$reliability<-1 - all$PEV / (1+all$inbreeding) / var2
  ped_data$BLUP<-ped_data$EBV
  
  n = table(training$ID)
  n = data.frame(n)
  colnames(n) = c('ID', "n")
  ped_data<-left_join(ped_data,n)
  
  #if (method == "GCA"){
  #  ped_data<-ped_data%>%filter(mdl=="P_GCA")
  #  ped_data1<-ped_data%>%select(PEDIGREE_NAME,matches("_predicted.value"),matches("_standard.error"),matches("_BLUP"),matches("_BLUP.se"),matches("_reliability"),matches("_trait"),matches("_n"),matches("_N"),matches("_h2"))
    
  #}else{
  #  ped_data<-ped_data
  #  ped_data1<-ped_data%>%select(PEDIGREE_NAME,matches("_predicted.value"),matches("_standard.error"),matches("_BLUP"),matches("_BLUP.se"),matches("_reliability"),matches("_trait"),matches("_n"),matches("_N"),matches("_h2"),matches("mdl"))
    
  #}
  
  if (method == "GCA"){
    ped_data<-ped_data%>%filter(mdl=="P_GCA")
    ped_data1<-ped_data%>%dplyr:select(PEDIGREE_NAME,matches("predicted.value"),matches("standard.error"),matches("BLUP"),matches("BLUP.se"),matches("reliability"),matches("trait"),matches("n"),matches("N"),matches("h2"))
    #ped_data1<-ped_data%>%dplyr:select(PEDIGREE_NAME,matches("predicted.value"),matches("standard.error"),matches("BLUP"),matches("BLUP.se"),matches("trait"),matches("n"),matches("N"),matches("h2"))
    
  }else{
    ped_data<-ped_data
    ped_data1<-ped_data%>%dplyr::select(PEDIGREE_NAME,matches("predicted.value"),matches("standard.error"),matches("BLUP"),matches("BLUP.se"),matches("reliability"),matches("trait"),matches("n"),matches("N"),matches("h2"),matches("mdl"))
    #ped_data1<-ped_data%>%dplyr::select(PEDIGREE_NAME,matches("predicted.value"),matches("standard.error"),matches("BLUP"),matches("BLUP.se"),matches("trait"),matches("n"),matches("N"),matches("h2"),matches("mdl"))
    
  }
  
  #[1] "PEDIGREE_NAME"          "ASI50D_predicted.value" "ASI50D_standard.error"  "ASI50D_BLUP"            "ASI50D_BLUP.se"         "ASI50D_reliability"     "ASI50D_trait"          
  #[8] "ASI50D_n"               "ASI50D_N"               "ASI50D_h2" 
  return(ped_data1)
}


output2<-function(trait,outbox,out,crop,method,tmp){
  
  ### column names according to AOB ## coded by OS 
  
  #outbox = paste0(out, '/SSGBLUP_deep/', trait, '_deep')
  #trait='RTLR'
  #out = "/mnt/veg_blup/test7"
  #@method : "GCA" or "all" for sca and gca
  #@tmp : inbreeding
  #method="GCA
  #Trait='RTLR'
  #AOB column names
  
  testing<-read.csv(file=paste0(outbox,"/EBV_y.txt"),header = TRUE)#modify with result
  testing<-testing[!duplicated(testing$ID),]
  training<-read.csv(file=paste0(outbox,"/testphenos_jwas.txt"),header = TRUE)# find hybrids of the season keep this file here that has NA
  ids<-read.csv(file=paste0(outbox,"/outids.txt"),header = FALSE)
  colnames(ids)<-c("PEDIGREE_NAME","ID")
  
  ped_data<-left_join(testing,ids)
  ped_data$N<-nrow(training)
  ped_data$trait<-trait
  ped_data<-ped_data %>% 
    mutate(BLUP = ped_data$EBV,
           predicted.value = BLUP + mean(na.omit(training$y)),
           mdl = NA)
  
  ###change
  ped_data$mdl[is.element(ped_data$ID, training$ID)==T] = 'P_SCA'
  ped_data$mdl[is.element(ped_data$ID, training$ID)==F] = 'P_GCA'
  table(ped_data$mdl)
  
  
  # squared error
  var2<-read.csv(file=paste0(outbox,"/MCMC_samples_polygenic_effects_variance.txt"))
  var2<-mean(var2$y.ID_y.ID)
  res_var <- read.csv(file=paste0(outbox,"/residual_variance.txt"))$Estimate 
  all<-left_join(ped_data,tmp,by=c("PEDIGREE_NAME"="ID"))
  all<-all[!duplicated(all$ID),]
  ped_data$reliability <-1 - all$PEV / (1+all$inbreeding) / var2
  ped_data <- ped_data %>%
    mutate(BLUP.se = sqrt(PEV),
           standard.error = sqrt(PEV + res_var))
  
  n = table(training$ID)
  n = data.frame(n)
  colnames(n) = c('ID', "n")
  ped_data<-left_join(ped_data,n)
  
  ped_data<-ped_data  %>%
    dplyr::rename(lineName = PEDIGREE_NAME)
  
  
  if (file.exists(paste0(outbox,"/heritability.txt"))==TRUE){
    h2<- read.csv(file=paste0(outbox,"/heritability.txt"))$Estimate
  }else if(method=="SSGBLUP"){
    files <- list.files(path=outbox,pattern = "MCMC_samples_y.x",full.names = T)%>% 
      map_df(~read_csv(.))
    rand<-data.frame(apply(files,2,function(x) mean(x,na.rm=T)))
    sum_rand<-sum(rand[,1])
    var_err<-read.csv(file=paste0(outbox,"/MCMC_samples_y.ϵ_variances.txt"))
    var_err<-mean(var_err$y.ϵ_y.ϵ)
    h2<-var2/(var2+sum_rand+var_err)
  }else if(method == "GBLUP"){
    files <- list.files(path=outbox,pattern = "MCMC_samples_y.x",full.names = T)%>% 
      map_df(~read_csv(.))
    rand<-data.frame(apply(files,2,function(x) mean(x,na.rm=T)))
    sum_rand<-sum(rand[,1])
    var_err<-read.csv(file=paste0(outbox,"/MCMC_samples_residual_variance.txt"))
    var_err<-mean(var_err$var_err$y_y)
    h2<-var2/(var2+sum_rand+var_err)
  }
  
  ped_data$h2<-h2
  
  
  ###################
  ### new columns format
  
  GCA_result <- subset(ped_data, mdl == 'P_GCA')
  SCA_result <- subset(ped_data, mdl == 'P_SCA')
  
  GCA_result <- GCA_result %>% 
    dplyr::select(lineName, predicted.value, standard.error, BLUP, BLUP.se, reliability, trait, n, h2)
  colnames(GCA_result) <- c('lineName', 'PGCA_PRED', 'PGCA_PRED_STD_ERR', 'PGCA_BLUP', 'PGCA_STD_ERR', 'PGCA_RELIABILITY', 'TRAIT', 'PGCA_X', 'PGCA_HERITABILITY')
  
  
  SCA_result <- SCA_result %>% 
    dplyr::select(lineName, predicted.value, standard.error, BLUP, BLUP.se, reliability, trait, n, h2)
  colnames(SCA_result) <- c('lineName', 'PSCA_PRED', 'PSCA_PRED_STD_ERR', 'PSCA_BLUP', 'PSCA_STD_ERR', 'PSCA_RELIABILITY', 'TRAIT', 'PSCA_N', 'PSCA_HERITABILITY')
  
  gca_long <- melt(GCA_result, id.vars = c('lineName', 'TRAIT'))
  order_exp <- c('PGCA_PRED', 'PGCA_PRED_STD_ERR', 'PGCA_RELIABILITY', 'PGCA_X', 'PGCA_HERITABILITY','PGCA_BLUP','PGCA_STD_ERR')
  gca_long <- gca_long[order(gca_long$lineName, gca_long$TRAIT, match(gca_long$variable, order_exp)),]
  gca_long <- gca_long %>% 
    mutate(traitRefId  = paste(TRAIT, variable, sep = '_')) %>% 
    dplyr::select(lineName, traitRefId, value)
  
  sca_long <- melt(SCA_result, id.vars = c('lineName', 'TRAIT'))
  order_exp_sca <- c('PSCA_PRED', 'PSCA_PRED_STD_ERR', 'PSCA_RELIABILITY', 'PSCA_N', 'PSCA_HERITABILITY','PSCA_BLUP','PSCA_STD_ERR')
  sca_long <- sca_long[order(sca_long$lineName, sca_long$TRAIT, match(sca_long$variable, order_exp_sca)),]
  sca_long <- sca_long %>% 
    mutate(traitRefId  = paste(TRAIT, variable, sep = '_')) %>% 
    dplyr::select(lineName, traitRefId, value)
  
  
  if (method == "GCA"){
    blup_df <-gca_long
    blup_df <- blup_df %>% 
      mutate(value = as.character(round(value,4)))
  }else{
    blup_df <- rbind(gca_long, sca_long)
    blup_df <- blup_df %>% 
      mutate(value = as.character(round(value,4)))
  }
  
  
  number_uomid <- 125239296
  ratio_uomid <- 179175424
  blup_df$numUOMId[grepl('_STD_ERR', blup_df$traitRefId)] <- number_uomid
  blup_df$numUOMId[grepl('_BLUP', blup_df$traitRefId)] <- number_uomid
  blup_df$numUOMId[grepl('_N', blup_df$traitRefId)] <- number_uomid
  blup_df$numUOMId[grepl('_X', blup_df$traitRefId)] <- number_uomid
  
  blup_df$numUOMId[grepl('_HERITABILITY', blup_df$traitRefId)] <- ratio_uomid
  blup_df$numUOMId[grepl('_RELIABILITY', blup_df$traitRefId)] <- ratio_uomid
  
  #if (method == "GCA"){
  #  ped_data<-ped_data%>%filter(mdl=="P_GCA")
  #  ped_data1<-ped_data%>%select(PEDIGREE_NAME,matches("predicted.value"),matches("standard.error"),matches("BLUP"),matches("BLUP.se"),matches("reliability"),matches("trait"),matches("n"),matches("N"),matches("h2"))
  
  #}else{
  #  ped_data<-ped_data
  #  ped_data1<-ped_data%>%select(PEDIGREE_NAME,matches("predicted.value"),matches("standard.error"),matches("BLUP"),matches("BLUP.se"),matches("reliability"),matches("trait"),matches("n"),matches("N"),matches("h2"),matches("mdl"))
  
  #}
  
  #[1] "PEDIGREE_NAME"          "ASI50D_predicted.value" "ASI50D_standard.error"  "ASI50D_BLUP"            "ASI50D_BLUP.se"         "ASI50D_reliability"     "ASI50D_trait"          
  #[8] "ASI50D_n"               "ASI50D_N"               "ASI50D_h2" 
  return(blup_df)
}



inbred_vals<-function(pheno_1,crop){
#library(asreml)


#pheno_1 <- s3read_using(readRDS, bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace', 
#                        object = 'ycao1/ITG/Corn/pheno.rds')

VaultToken <- content(POST(url = paste0(VaultURL, "auth/approle/login"), body= AppRoleCredsJson, encode = 'json'))
print(VaultToken)
AWSCreds <- content(GET(url = paste0(VaultURL, VaultSecretPath), httr::add_headers('X-Vault-Token' = VaultToken$auth$client_token)))
print(AWSCreds)
Sys.setenv("AWS_ACCESS_KEY_ID" = AWSCreds$data$AWS_ACCESS_KEY_ID,
           "AWS_SECRET_ACCESS_KEY" = AWSCreds$data$AWS_SECRET_ACCESS_KEY,
           "AWS_DEFAULT_REGION" = "us-east-1")
Sys.setenv("CLIENT_ID" = CLIENT_ID,
           "CLIENT_SECRET" = CLIENT_SECRET)

#crop       = 'Corn'
IDTableName <- paste0(crop, '/', crop, '_IDs.RData')
print(IDTableName)
s3load(object = IDTableName, bucket = "veg-apd-sdi-predictiveanalytcs-prod-reference-data")
cropids    = CropIDs[, c('M.GERMPLASM.X_ID','M.GERMPLASM.PEDIGREE')]
cropids    = cropids[duplicated(cropids)==F, ]

raw_dat_sub <- pheno_1 %>% 
  mutate(P1 = as.character(P1), 
         P2 = as.character(P2),
         trait = OBSRVTN_REF_CD, 
         TRAIT_VALUE = as.numeric(TRAIT_VALUE), 
         GROWSEASON = as.character(GROWSEASON),
         pedigree = as.character(PEDIGREE_NAME))
raw_dat_sub[is.na(raw_dat_sub$P1), 'P2'] <- NA
raw_dat_sub[is.na(raw_dat_sub$P2), 'P1'] <- NA
raw_dat_sub <- unique(raw_dat_sub[, c('pedigree', 'P1', 'P2', 'LTYPE')])


gen5 <- make_ped_file(raw_dat_sub, crop_ids = CropIDs, n_gen = 5)
#A_matrix5 <- asreml.Ainverse(gen5$ped, fgen = list("inbreeding", 5))

pd <- gen5$ped
pd = unique(pd)
### replace pd$ID with pedigree name from cropID
# change GERMPLASM ID to pedigree names
tmp = left_join(pd, cropids, by=c('ID'='M.GERMPLASM.X_ID'))
tmp = left_join(tmp, cropids, by=c('PARENT_FEMALE'='M.GERMPLASM.X_ID'))
tmp = left_join(tmp, cropids, by=c('PARENT_MALE'='M.GERMPLASM.X_ID'))
tmp = tmp[, c('M.GERMPLASM.PEDIGREE.x', 'M.GERMPLASM.PEDIGREE.y', 'M.GERMPLASM.PEDIGREE')]
colnames(tmp) =  c("ID", "PARENT_FEMALE", "PARENT_MALE")
tmp$inbreeding = 0
file1<-paste0("/mnt/5gen_inbreeding_",crop,".Rdata")
save(tmp, file=file1)
#save(tmp,file="5gen_inbreeding_SweetCorn.Rdata")
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

