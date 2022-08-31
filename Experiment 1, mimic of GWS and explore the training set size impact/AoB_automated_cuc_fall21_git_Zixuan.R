##########
# @ ENVIRONMENT : Duplicate of R 3.5 and Python 3.7 -- EVA and Asreml

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
source('/mnt/getbdat_AOB_V3.R')
source('/mnt/JWAS_Model.R')
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

## example corn
#dataset <- s3read_using(readRDS, bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace', 
#                        object = 'ycao1/ITG/Corn/pheno.rds')

## example cucumber example AOB cucumber.R
## clean up with calculated traits and outlier detection of only test sets for genotyping see phenos_cucumber.R

s3load(object = '/ycao1/ITG/Cucumber/Fall/Cln_cucumber_fall12102021_new.Rdata', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

#load("/mnt/Cln_cucumber_fall12102021_new.Rdata")
dataset<-pheno

## UOM Conversion
UOMTable <- read.csv("/mnt/MIDAS_UOM_Table.csv", colClasses = rep(times=5,x="character"))
dataset$UOM <- UOMTable$NAME[match(dataset$UOM,UOMTable$UOM_ID)]
dataset_uom <- conv_uom(dataset) ####Look at the details about this function
dataset_uom_copy <- dataset_uom
################################################
## load genotypes


s3load(object = '/ycao1/ITG/Cucumber/hybrids/hybrid_geno_cucumber.Rdata', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
#load("/mnt/hybrid_geno_cucumber.Rdata")
g<-hybrid_wide_ped
gg = g[,2:ncol(g)]
colnames(gg)[1] = 'PEDIGREE'


## training set 2020 GROWSEASON=="2020:08"
# DH lines processing_2020_fall
s3load(object = '/ycao1/ITG/Cucumber/hybrids/Fall_2020_geno_traincuc.Rdata', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

#load("/mnt/DH_processing/Fall_2020_geno_traincuc.Rdata")

g1<-Fall_20_train
gg1 = g1[,2:ncol(g1)]
colnames(gg1)[1] = 'PEDIGREE'

## training set 2021 GROWSEASON=="2021:08"
# DH lines processing_2021_fall
s3load(object = '/ycao1/ITG/Cucumber/hybrids/Fall_2021_geno_traincuc.Rdata', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

load("/mnt/DH_processing/Fall_2021_geno_traincuc.Rdata")

g2<-Fall_21_train
gg2 = g2[,2:ncol(g2)]
colnames(gg2)[1] = 'PEDIGREE'

####### Try to find out how the genotype data is linked with the training set.


## all genos
gg4 <- rbind(gg, gg1, gg2)
tst <- gg4[, 3:ncol(gg4)]
tst <- as.matrix(tst)
for (i in 1:ncol(tst)){
  tst[,i][which(is.na(tst[,i]))] <- ceiling(mean(tst[,i], na.rm = T))
}

tst <- as.data.frame(tst)
gg4[,2:ncol(gg4)] <- tst
table(is.na(gg4))

gg4=gg4%>%distinct()

#table(gg$PEDIGREE%in%dataset_uom$PEDIGREE_NAME)
#FALSE  TRUE 
#2694  1447 


################################################
# create test folder to store intermediate files for JWAS to run
#out0 = "/mnt/veg_blup/swtest"
#if(file.exists(out0)==F) {dir.create(out0, recursive = T)}

out0 = "/mnt/veg_blup/train_fall_2021_1yr/Prediction_2017to2020_to_2021SSGBLUP_75%"
if(file.exists(out0)==F) {dir.create(out0, recursive = T)}

### QA/QC case 3
#remove filler

dataset_uom<-dataset_uom%>%filter(!PEDIGREE_NAME=="FILLER")

#### Run once the inbreeding code

inbred_vals(pheno_1=dataset_uom,crop='Cucumber')###Need to change the environment

#load("/mnt/5gen_inbreeding_Cucumber.Rdata")

################################################
## Calculate stage count for parents
## needed to run the master table needs to be run with the prediction code file6

#stage_count_dat <- calc_stage_count(dataset_uom)

#writing_csv <- function(dat, filename){
#  write.csv(dat, file = filename, row.names = F)
#}

#s3write_using(stage_count_dat, FUN = writing_csv,
#               object = '/Cucumber/Cucumber_STAGE_COUNTxxx.csv',
#               bucket = "veg-apd-sdi-predictiveanalytcs-prod-reference-data")

###################
## function SSGBLUP
###################

TRAITS<-c('FRLGT','EXTCO',"AFW_C","FRNMK","SHAPE",
          "WFSA_C","NFSA_C","FRNMKAvg_1-57","WFSA_CAvg_1-57",
          "NFSA_CAvg_1-57","WFSA_CSum_1-57","NFSA_CSum_1-57","NETWTSum_1-57" )

### scenario 1 reduce dataset 3 years data
dataset_uom<-dataset_uom%>%filter(!GROWSEASON%in%c("2015:08","2016:08","2017:08"))

### scenario 2 training set 2021
dataset_uom<-dataset_uom%>%filter(GROWSEASON%in%c("2021:08"))

dataset_uom<-dataset_uom%>%filter(OBSRVTN_REF_CD%in%TRAITS)

#@method1 = SSGBLUP or ABLUP or GBLUP   ## Check GBLUP definition.
#@Also try different methods. 
#@Start with the continuous variables. E.g., EXTCO, SHAPE.

#@Can also try different seasons (Or years).
#@seas=NULL is prediction with phenotype in model, seas="2022:08" masks the phenotype of the season
#@ gg = NULL for ABLUP
#@Cross-Validation to get more accuracy.

results <- data.frame()
for (i in 1:length(TRAITS)){
  trait_choosen = TRAITS[i]
  print(paste(" ################--- ", trait_choosen,"__is completed ############"))
  method1='SSGBLUP'
  #res = getoutput_new(out=out0,dataset=dataset_uom_new,trait = trait_choosen,gg=gg2,n_gen=5,seas=NULL,method=method1)
  res = getoutput_new(out=out0,dataset=dataset_uom_new,trait = trait_choosen,gg=gg,n_gen=5,seas=NULL,method=method1)
  outbox<- paste0(out0, method1,'_deep/', trait = trait_choosen, '_deep')
  trait_GCA<-output_new(trait=trait_choosen,outbox=outbox,out=out0,crop=crop,method="SCA",tmp=tmp,gmethod=method1)
  results = rbind(results,trait_GCA)
  save(results,file="results_asreml_ABLUP_with_gg_Zixuan.Rdata")
  #save(results,file="results_aBLUP.Rdata")
}

save(True_values_2021, file = 'Original_phenotypic_data.Rdata')


set.seed(54321)

### scenario 3: Use 2020:08 P1 to predict 2021:08 P2
year_train_test <- '2020-2021'
dataset_2020_fall_P1 <- dataset_uom%>%filter(GROWSEASON%in%c("2020:08") & EXPER_STAGE_REF_ID == 'P1')
dataset_2021_fall_P2 <- dataset_uom%>%filter(GROWSEASON%in%c("2021:08") & EXPER_STAGE_REF_ID == 'P2')
True_values_2021 <- dataset_2021_fall_P2[,c('PEDIGREE_NAME', 'OBSRVTN_REF_CD', "TRAIT_VALUE", 'EXPER_STAGE_REF_ID')]
dataset_2021_fall_P2$TRAIT_VALUE <- NA
dataset_uom_new <- rbind(dataset_2020_fall_P1, dataset_2021_fall_P2)
#sum(is.na(dataset_uom_new$TRAIT_VALUE))


### scenario 4: Use 2017:08 P1 to predict 2018:08 P2
set.seed(321)
year_train_test <- '2017-2018'
dataset_2017_fall_P1 <- dataset_uom%>%filter(GROWSEASON%in%c("2017:08") & EXPER_STAGE_REF_ID == 'P1')
dataset_2018_fall_P2 <- dataset_uom%>%filter(GROWSEASON%in%c("2018:08") & EXPER_STAGE_REF_ID == 'P2')
True_values_2018 <- dataset_2018_fall_P2[,c('PEDIGREE_NAME', 'OBSRVTN_REF_CD', "TRAIT_VALUE", 'EXPER_STAGE_REF_ID')]
dataset_2018_fall_P2$TRAIT_VALUE <- NA
dataset_uom_new <- rbind(dataset_2017_fall_P1, dataset_2018_fall_P2)
#sum(is.na(dataset_uom_new$TRAIT_VALUE))




### scenario 5: Use 75% of 2020:08 P1 to predict 2021:08 P2
set.seed(32)
year_train_test <- '2020-2021'
dataset_2020_fall_P1 <- dataset_uom%>%filter(GROWSEASON%in%c("2020:08") & EXPER_STAGE_REF_ID == 'P1' & OBSRVTN_REF_CD %in% TRAITS)
new_trainingset <- data.frame()
for (i in 1:length(TRAITS)){
  trait_chosen_set <- dataset_2020_fall_P1[dataset_2020_fall_P1$OBSRVTN_REF_CD == TRAITS[i],]
  number_to_chosen <- floor(nrow(trait_chosen_set)*0.75)
  new_training_chosen <- sample(nrow(trait_chosen_set), size = number_to_chosen)
  trait_chosen_set <- trait_chosen_set[new_training_chosen,]
  new_trainingset <- rbind(new_trainingset, trait_chosen_set)
}
dataset_2020_fall_P1 <- new_trainingset
dataset_2021_fall_P2 <- dataset_uom%>%filter(GROWSEASON%in%c("2021:08") & EXPER_STAGE_REF_ID == 'P2' & OBSRVTN_REF_CD %in% TRAITS)
True_values_2021 <- dataset_2021_fall_P2[,c('PEDIGREE_NAME', 'OBSRVTN_REF_CD', "TRAIT_VALUE", 'EXPER_STAGE_REF_ID')]
dataset_2021_fall_P2$TRAIT_VALUE <- NA
dataset_uom_new <- rbind(dataset_2020_fall_P1, dataset_2021_fall_P2)
#sum(is.na(dataset_uom_new$TRAIT_VALUE))


set.seed(12345)

### scenario 6: Use 75% of 2020:08 P1 to predict 2021:08 P2 with SSGBLUP
year_train_test <- '2020-2021'
dataset_2020_fall_P1 <- dataset_uom%>%filter(GROWSEASON%in%c("2020:08") & EXPER_STAGE_REF_ID == 'P1' & OBSRVTN_REF_CD %in% TRAITS)
new_trainingset <- data.frame()
for (i in 1:length(TRAITS)){
  trait_chosen_set <- dataset_2020_fall_P1[dataset_2020_fall_P1$OBSRVTN_REF_CD == TRAITS[i],]
  number_to_chosen <- floor(nrow(trait_chosen_set)*0.75)
  new_training_chosen <- sample(nrow(trait_chosen_set), size = number_to_chosen)
  trait_chosen_set <- trait_chosen_set[new_training_chosen,]
  new_trainingset <- rbind(new_trainingset, trait_chosen_set)
}
dataset_2020_fall_P1 <- new_trainingset
dataset_2021_fall_P2 <- dataset_uom%>%filter(GROWSEASON%in%c("2021:08") & EXPER_STAGE_REF_ID == 'P2' & OBSRVTN_REF_CD %in% TRAITS)
True_values_2021 <- dataset_2021_fall_P2[,c('PEDIGREE_NAME', 'OBSRVTN_REF_CD', "TRAIT_VALUE", 'EXPER_STAGE_REF_ID')]
dataset_2021_fall_P2$TRAIT_VALUE <- NA
dataset_uom_new <- rbind(dataset_2020_fall_P1, dataset_2021_fall_P2)
#sum(is.na(dataset_uom_new$TRAIT_VALUE))



### scenario 6: Use all of 2020:08, 2019:08, 2018:08, 2017:08 to predict 2021:08 with SSGBLUP
#year_train_test <- '2020-2021'
dataset_2017_to_2020_fall <- dataset_uom%>%filter(GROWSEASON%in%c("2020:08","2019:08","2018:08","2017:08") & OBSRVTN_REF_CD %in% TRAITS)
#new_trainingset <- data.frame()
#for (i in 1:length(TRAITS)){
#  trait_chosen_set <- dataset_2020_fall_P1[dataset_2020_fall_P1$OBSRVTN_REF_CD == TRAITS[i],]
#  number_to_chosen <- floor(nrow(trait_chosen_set)*0.75)
#  new_training_chosen <- sample(nrow(trait_chosen_set), size = number_to_chosen)
#  trait_chosen_set <- trait_chosen_set[new_training_chosen,]
#  new_trainingset <- rbind(new_trainingset, trait_chosen_set)
#}
#dataset_2020_fall_P1 <- new_trainingset
dataset_2021_fall <- dataset_uom%>%filter(GROWSEASON%in%c("2021:08") & OBSRVTN_REF_CD %in% TRAITS)
True_values_2021 <- dataset_2021_fall[,c('PEDIGREE_NAME', 'OBSRVTN_REF_CD', "TRAIT_VALUE", 'EXPER_STAGE_REF_ID', 'P1', 'P2', 'LTYPE')]
dataset_2017_to_2020_fall <- dataset_2017_to_2020_fall[!dataset_2017_to_2020_fall$PEDIGREE_NAME %in% unique(dataset_2021_fall$PEDIGREE_NAME),]
dataset_2021_fall$TRAIT_VALUE <- NA
dataset_uom_new <- rbind(dataset_2017_to_2020_fall, dataset_2021_fall)
#sum(is.na(dataset_uom_new$TRAIT_VALUE))





### scenario 7: DH0 optimization for 2021:08

dataset_2019_fall_P1 <- dataset_uom%>%filter(GROWSEASON%in%c("2019:08") & EXPER_STAGE_REF_ID == 'P1' & OBSRVTN_REF_CD %in% TRAITS)

