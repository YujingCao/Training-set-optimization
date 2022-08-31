##########
# @ ENVIRONMENT : Any environment without ASREML
########################################### Input parameters
args = commandArgs(trailingOnly = TRUE)


#For the split of pedigree
repeat_times <- as.numeric(args[1])
percent <- sapply(unlist(lapply(strsplit(args[2], ','), trimws, 'both')),as.numeric)
BLUPs <- unlist(lapply(strsplit(args[3], ','), trimws, 'both'))
out0_input <- args[4] 
TRAITS <- unlist(lapply(strsplit(args[5], ','), trimws, 'both'))
GROWSEASONS_TRAINING <- unlist(lapply(strsplit(args[6], ','), trimws, 'both'))
GROWSEASONS_TEST <- unlist(lapply(strsplit(args[7], ','), trimws, 'both'))
STAGE_TRAINING <- unlist(lapply(strsplit(args[8], ','), trimws, 'both'))
STAGE_TEST <- unlist(lapply(strsplit(args[9], ','), trimws, 'both'))
seed <- as.numeric(args[10])

set.seed(seed)


########################################### Parameters for testing


#repeat_times <- 1
#percent <- 1
#BLUPs <- c('GBLUP')#, 'ABLUP')
#out0_input <- '/mnt/veg_blup/train_spring_test_code/'
#TRAITS <- c("FRLGT")#,"FRLGT","AFW_C","FRNMK","SHAPE")
#GROWSEASONS_TEST <- c('2019:01', '2020:01', '2021:01')
#GROWSEASONS_TRAINING <- c('2018:01','2019:01', '2020:01', '2021:01')
#STAGE_TRAINING <- 'P1'
#STAGE_TEST <- c('P2', 'P3')
#set.seed(19)




install.packages('azurequest',repos=c('https://cran.science-at-scale.io', 'https://cran.rstudio.com'))
install.packages('pedigree',repos=c('https://cran.science-at-scale.io', 'https://cran.rstudio.com'))


# This is an example to run AoB with data pulled from H2H

# Dynamic way of checking required libraries in R
list.of.packages <- c("aws.s3", "httr", "tidyverse", "asreml", "asremlPlus", "measurements", "reshape", "tidyr", 
                      "jsonlite", "dplyr", "readr", "foreach", "parallel", "doParallel","data.table","reshape2","azurequest","pedigree")
dynamic_require <- function(package){
  if(eval(parse(text=paste("require(",package,")")))){return(TRUE)}
  install.packages(package)
  return(eval(parse(text=paste("require(",package,")"))))
}
sapply(list.of.packages, dynamic_require)
rename = dplyr::rename

# Source functions
source("/mnt/deep_ped_functions_ODB.R")
source('/mnt/vaultCredentials.R')
source('/mnt/ValutCreds_ancestry.R')
source('/mnt/getbdat_AOB_V3_git.R')
source('/mnt/JWAS_Model_git.R')
options(scipen = 999)

#### Vault credentials
VaultToken <- content(POST(url = paste0(VaultURL, "auth/approle/login"), body= AppRoleCredsJson, encode = 'json'))
print(VaultToken)
AWSCreds <- content(GET(url = paste0(VaultURL, VaultSecretPath), httr::add_headers('X-Vault-Token' = VaultToken$auth$client_token)))
print(AWSCreds)
Sys.setenv("AWS_ACCESS_KEY_ID" = AWSCreds$data$AWS_ACCESS_KEY_ID,
           "AWS_SECRET_ACCESS_KEY" = AWSCreds$data$AWS_SECRET_ACCESS_KEY,
           "AWS_DEFAULT_REGION" = "us-east-1")

CLIENT_ID <- "599849ba-355c-48aa-8154-994ab3dd79b1"
CLIENT_SECRET <- "bXpU8e0mN~-.cR9.3-14XV.9X3pjpc69pl"

Sys.setenv("CLIENT_ID" = CLIENT_ID,
           "CLIENT_SECRET" = CLIENT_SECRET)

################################################
#### read cropIDs
crop       = 'Cucumber'

IDTableName <- paste0(crop, '/', crop, '_IDs.RData')
print(IDTableName)
s3load(object = IDTableName, bucket = "veg-apd-sdi-predictiveanalytcs-prod-reference-data")
cropids    = CropIDs[, c('M.GERMPLASM.X_ID','M.GERMPLASM.PEDIGREE')]
cropids    = cropids[duplicated(cropids)==F, ]

################################################
#### read phenos
## example cucumber example AOB cucumber.R
## clean up with calculated traits and outlier detection of only test sets for genotyping see phenos_cucumber.R
#load("/mnt/Cln_cucumber_spring04272022.Rdata",verbose = T)

s3load(object = 'Cucumber/Cucumber_Parentals.RData', bucket = "veg-apd-sdi-predictiveanalytcs-prod-reference-data")

#if (season == 'spring'){
#if (season == 'spring'){
s3load(object = '/ycao1/ITG/Cucumber/Summer/Cln_cucumber_summer18022022_new.Rdata', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

dataset<-pheno

## UOM Conversion
UOMTable <- read.csv("/mnt/MIDAS_UOM_Table.csv", colClasses = rep(times=5,x="character"))
dataset$UOM <- UOMTable$NAME[match(dataset$UOM,UOMTable$UOM_ID)]
dataset_uom <- conv_uom(dataset) 
dataset_uom_copy <- dataset_uom
################################################
## load genotypes

## training set 2015-2019 might need reprocessing from the genomics team if names had been changed due to coded
## training set 2015-2019 might need reprocessing from the genomics team if names had been changed
load("/mnt/hybrid_geno_cucumber.Rdata")
g<-hybrid_wide_ped
gg = g[,2:ncol(g)]
colnames(gg)[1] = 'PEDIGREE'

## training set 2020
load("/mnt/summer_2022_DHxTsr/summer_20_traincuc.Rdata")
g1<-summer_20_train
gg1 = g1[,2:ncol(g1)]
colnames(gg1)[1] = 'PEDIGREE'

## training set 2021
load("/mnt/summer_2022_DHxTsr/summer_21_traincuc.Rdata")
g2<-summer_21_train
gg2 = g2[,2:ncol(g2)]
colnames(gg2)[1] = 'PEDIGREE'

## training set 2022
load("/mnt/summer_2022_DHxTsr/summer_22_traincuc.Rdata")
g3<-summer_22_train
gg3 = g3[,2:ncol(g3)]
colnames(gg3)[1] = 'PEDIGREE'


## all genos
gg4=rbind(gg,gg1,gg2,gg3)

tst <- gg4[, 3:ncol(gg4)]
tst <- as.matrix(tst)
for (i in 1:ncol(tst)){
  tst[,i][which(is.na(tst[,i]))] <- ceiling(mean(tst[,i], na.rm = T))
}

tst <- as.data.frame(tst)
gg4[,2:ncol(gg4)] <- tst
table(is.na(gg4))


#gg4$PEDIGREE = as.character(gg4$PEDIGREE)
#stopifnot(nrow(gg4) == length(unique(gg4$PEDIGREE)))
#table(gg4$PEDIGREE%in%dataset_uom$PEDIGREE_NAME)
#FALSE  TRUE 
#2694  1447 

################################################

#}


################################################ Filter the PEDIGREE existed in the genotype data.

dataset_uom <- dataset_uom %>%filter(PEDIGREE_NAME %in% gg4$PEDIGREE)

dataset_uom_copy_spring <- dataset_uom







################################################ Seperate out the test set.

load("/mnt/5gen_inbreeding_Cucumber.Rdata")


dataset_uom_test <- dataset_uom[dataset_uom$GROWSEASON %in% GROWSEASONS_TEST & 
                                  dataset_uom$EXPER_STAGE_REF_ID %in% STAGE_TEST & 
                                  dataset_uom$LTYPE == 'Hybrid' &
                                  dataset_uom$OBSRVTN_REF_CD %in% TRAITS,]

## Save the crucial information of the test set.

True_values <- dataset_uom_test[,c('PEDIGREE_NAME', 
                                   'OBSRVTN_REF_CD', 
                                   "TRAIT_VALUE",
                                   'EXPER_STAGE_REF_ID',
                                   'P1', 'P2', 'LTYPE', 'ORIGIN', 'FIELD_NAME', 'IS_CHK')]#if(file.exists(out0)==F) {dir.create(out0, recursive = T)}
#repeat_times <- 1
dataset_uom_test$TRAIT_VALUE <- NA
#percent <- c(0.2, 0.5, 0.8)
#BLUPs <- c('SSGBLUP', 'ABLUP')
#out0_input <- '/mnt/veg_blup/train_spring_rolling_exp/'


for (replication in 1:repeat_times){
  
  
  for(perc in 1:length(percent)){
    for (j in 1:length(BLUPs)){
      
      
      #out0 = "/mnt/veg_blup/train_spring_Dunia_2/"  
      out0 <- paste0(out0_input, percent[perc]*100, '/',replication, '/')
      if(file.exists(out0)==F) {dir.create(out0, recursive = T)}
      
      
      ################ Construct training set
      
      dataset_uom <- dataset_uom_copy_spring 
      dataset_uom<-dataset_uom%>%dplyr::filter(OBSRVTN_REF_CD%in%TRAITS)
      dataset_uom<-dataset_uom%>%filter(GROWSEASON%in%GROWSEASONS_TRAINING&
                                          EXPER_STAGE_REF_ID %in% STAGE_TRAINING &
                                          LTYPE == 'Hybrid')
      dataset_uom[dataset_uom$PEDIGREE_NAME %in% dataset_uom_test$PEDIGREE_NAME, 'TRAIT_VALUE'] <- NA
      candidate_set <- dataset_uom%>%filter(!is.na(TRAIT_VALUE))
      
      
      ########## Filter the unique pedigrees
      
      
      
      uni_ped <- unique(candidate_set$PEDIGREE_NAME)
      
      sample_size <- ceiling(length(uni_ped)*percent[perc])
      sample_training <- sample(uni_ped,size = sample_size)
      dataset_uom_training <- dataset_uom%>%filter(PEDIGREE_NAME%in%sample_training)
      
      
      ######## Save the crucial information of Training set.
      
      
      Training_main_info <- dataset_uom_training[, c('PEDIGREE_NAME', 
                                                     'OBSRVTN_REF_CD', 
                                                     "TRAIT_VALUE",
                                                     'EXPER_STAGE_REF_ID',
                                                     'P1', 'P2', 'LTYPE', 'ORIGIN', 'FIELD_NAME', 'IS_CHK')]
      
      dataset_uom <- rbind(dataset_uom_training, dataset_uom_test)
      
      ### QA/QC case 3
      #remove filler
      
      #dataset_uom<-dataset_uom%>%filter(!PEDIGREE_NAME=="FILLER")
      
      #### Run once the inbreeding code
      
      #inbred_vals(pheno_1=dataset_uom,crop='Cucumber')
      
      
      
      
      
      
      
      
      
      
      
      results <- data.frame()
      
      for (i in 1:length(TRAITS)){
        #for (i in 1:3){
        trait_choosen = TRAITS[i]
        print(paste(" ################--- ", trait_choosen,"__is completed ############"))
        method1 = BLUPs[j]
        if(method1 == 'SSGBLUP'){
          res = getoutput_new(out=out0,dataset=dataset_uom,trait = trait_choosen,gg=gg4,n_gen=5,seas=NULL,method=method1)
        }
        if(method1 == 'ABLUP'){
          res = getoutput_new(out=out0,dataset=dataset_uom,trait = trait_choosen,gg=NULL,n_gen=5,seas=NULL,method=method1)
        }
        if(method1 == 'GBLUP'){
          res = getoutput_new(out=out0,dataset=dataset_uom,trait = trait_choosen,gg=gg4,n_gen=5,seas=NULL,method=method1)
        }
        outbox<-paste0(out0, method1, '_deep', '/',trait=trait_choosen, '_deep')
        
        trait_GCA<-output_new(trait=trait_choosen,outbox=outbox,out=out0,crop=crop,method="SCA",tmp=tmp,gmethod=method1)
        results = rbind(results,trait_GCA)
        save(results,file="results_asreml_zixuan.Rdata")
        
        
      }
      save(Training_main_info, file = 'Training_data.Rdata')
      save(True_values, file="original_asreml_zixuan.Rdata")
    }
  }}