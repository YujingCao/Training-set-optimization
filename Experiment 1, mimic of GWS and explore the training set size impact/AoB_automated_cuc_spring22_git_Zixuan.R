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
s3load(object = '/ycao1/ITG/Cucumber/Spring/Cln_cucumber_spring04272022.Rdata', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

dataset<-pheno

## UOM Conversion
UOMTable <- read.csv("/mnt/MIDAS_UOM_Table.csv", colClasses = rep(times=5,x="character"))
dataset$UOM <- UOMTable$NAME[match(dataset$UOM,UOMTable$UOM_ID)]
dataset_uom <- conv_uom(dataset) 
dataset_uom_copy_spring <- dataset_uom
################################################
## load genotypes

## training set 2015-2019 might need reprocessing from the genomics team if names had been changed due to coded
load("/mnt/hybrid_geno_cucumber.Rdata")
g<-hybrid_wide_ped
gg = g[,2:ncol(g)]
colnames(gg)[1] = 'PEDIGREE'

## training set 2020
load("/mnt/spring_2022_PCM/spring_20_traincuc.Rdata")
g1<-spring_20_train
gg1 = g1[,2:ncol(g1)]
colnames(gg1)[1] = 'PEDIGREE'

## training set 2021
load("/mnt/spring_2022_PCM/spring_21_traincuc.Rdata")
g2<-spring_21_train
gg2 = g2[,2:ncol(g2)]
colnames(gg2)[1] = 'PEDIGREE'

## training set 2022
load("/mnt/spring_2022_PCM/spring_22_traincuc.Rdata")
g3<-spring_22_train
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

gg4=gg4%>%distinct()
#table(gg$PEDIGREE%in%dataset_uom$PEDIGREE_NAME)
#FALSE  TRUE 
#2694  1447 

################################################


#args = commandArgs(trailingOnly = TRUE)


set.seed(9883)#For the split within pedigree
set.seed(9887)#For the split of pedigree
repeat_times <- as.numeric(args[1])
percent <- sapply(unlist(lapply(strsplit(args[2], ','), trimws, 'both')),as.numeric)
BLUPs <- unlist(lapply(strsplit(args[3], ','), trimws, 'both'))
out0_input <- args[4] 
#if(file.exists(out0)==F) {dir.create(out0, recursive = T)}
#repeat_times <- 1
#percent <- c(0.2, 0.5, 0.8)
#BLUPs <- c('SSGBLUP', 'ABLUP')
#out0_input <- '/mnt/veg_blup/train_spring_rolling_exp/'
seeds <- sample(10000, size = length(percent))
#True_values_spring <- dataset_uom_copy_spring[dataset_uom_copy_spring$GROWSEASON == '2022:01', 
#                                              c('PEDIGREE_NAME', 
#                                                'OBSRVTN_REF_CD', 
#                                                "TRAIT_VALUE",
#                                                'EXPER_STAGE_REF_ID',
#                                                'P1', 'P2', 'LTYPE')]
#dataset_uom_test <- dataset_uom_copy_spring[dataset_uom_copy_spring$GROWSEASON == '2022:01',]
#dataset_uom_test$TRAIT_VALUE <- NA
for (replication in 1:repeat_times){
  

for(perc in 1:length(percent)){
for (j in 1:length(BLUPs)){

  
#out0 = "/mnt/veg_blup/train_spring_Dunia_2/"  
  out0 <- paste0(out0_input, percent[perc]*100, '/',replication, '/')
  if(file.exists(out0)==F) {dir.create(out0, recursive = T)}
##Dunia:


#Dunia_training_set <- dataset_uom_copy_spring%>%
#  filter(GROWSEASON%in%c("2017:01", "2018:01","2016:01") & LTYPE == 'Hybrid')
#Dunia_test_set_first_step <- dataset_uom_copy_spring%>%
#  filter(GROWSEASON%in%c("2019:01", "2020:01","2021:01") & 
#           LTYPE == 'Hybrid' & EXPER_STAGE_REF_ID %in% c('P2', 'P3'))
#Dunia_test_set <- Dunia_test_set_first_step%>%
#  filter(!PEDIGREE_NAME %in% Dunia_training_set$PEDIGREE_NAME)
#
#dataset_uom_test <- Dunia_test_set
#dataset_uom_training <- Dunia_training_set
#dataset_uom_test$TRAIT_VALUE <- NA
#dataset_uom <- rbind(dataset_uom_training, dataset_uom_test)
#True_values_spring <- Dunia_test_set[,c('PEDIGREE_NAME', 
#                                        'OBSRVTN_REF_CD', 
#                                        "TRAIT_VALUE",
#                                        'EXPER_STAGE_REF_ID',
#                                        'P1', 'P2', 'LTYPE')]
  TRAITS<-c("EXTCO","FRLGT","AFW_C","FRNMK","SHAPE")
#### Dunia done
  
  load("/mnt/spring_2022_PCM/5gen_inbreeding_Cucumber.Rdata")
  dataset_uom <- dataset_uom_copy_spring 
  ##Dunia:
  dataset_uom %>% filter(OBSRVTN_REF_CD %in% TRAITS)
  
  Dunia_training_set <- dataset_uom_copy_spring%>%
    filter(GROWSEASON%in%c("2018:01","2019:01", '2020:01', '2021:01') & LTYPE == 'Hybrid')
  Dunia_test_set <- dataset_uom_copy_spring%>%
    filter(GROWSEASON%in%c("2019:01", "2020:01","2021:01") & 
             LTYPE == 'Hybrid' & EXPER_STAGE_REF_ID %in% c('P2', 'P3'))
  Dunia_training_set[Dunia_training_set$PEDIGREE_NAME %in% Dunia_test_set$PEDIGREE_NAME,
                     "TRAIT_VALUE"] <- NA
  
  dataset_uom_test <- Dunia_test_set
  dataset_uom_training <- Dunia_training_set
  dataset_uom_test$TRAIT_VALUE <- NA
  dataset_uom <- rbind(dataset_uom_training, dataset_uom_test)

  #TRAITS<-c("EXTCO","FRLGT","AFW_C","FRNMK","SHAPE")
  #### Dunia done
  
  



#unique(dataset_uom$OBSRVTN_REF_CD)
#TRAITS<-c("EXTCO","FRLGT","AFW_C","FRNMK","SHAPE")#,"WFSA_C","NFSA_C")

### scenario 1 training set 2015-2018 + 2020-2021
#  dataset_uom<-dataset_uom%>%dplyr::filter(OBSRVTN_REF_CD%in%TRAITS)
#  dataset_uom<-dataset_uom%>%filter(GROWSEASON%in%c("2022:01"))
#  uni_ped <- unique(dataset_uom$PEDIGREE_NAME)
#  set.seed(seeds[perc])
#  sample_size <- ceiling(length(uni_ped)*percent[perc])
#  sample_training <- sample(uni_ped,size = sample_size)
#  dataset_uom_training <- dataset_uom%>%filter(PEDIGREE_NAME%in%sample_training)
#  dataset_uom_test <- dataset_uom[!dataset_uom$PEDIGREE_NAME %in% dataset_uom_training$PEDIGREE_NAME,]

          #,"FRNMKAvg_1-84","WFSA_CAvg_1-84","NFSA_CAvg_1-84","NETWTAvg_1-84","NFSA_CSum_1-84","WFSA_CSum_1-84", "NETWTSum_1-84")
  
#  True_values_spring <- dataset_uom_test[,c('PEDIGREE_NAME', 
#                                          'OBSRVTN_REF_CD', 
#                                          "TRAIT_VALUE",
#                                          'EXPER_STAGE_REF_ID',
#                                          'P1', 'P2', 'LTYPE')]
#True_values_spring <- dataset_uom_test[, c('PEDIGREE_NAME', 
#                                           'OBSRVTN_REF_CD', 
#                                            "TRAIT_VALUE",
#                                            'EXPER_STAGE_REF_ID',
#                                            'P1', 'P2', 'LTYPE')]

#  dataset_uom_test$TRAIT_VALUE <- NA
#  dataset_uom <- rbind(dataset_uom_training, dataset_uom_test)

### QA/QC case 3
#remove filler

#dataset_uom<-dataset_uom%>%filter(!PEDIGREE_NAME=="FILLER")

#### Run once the inbreeding code

#inbred_vals(pheno_1=dataset_uom,crop='Cucumber')



###################
## function SSGBLUP
###################

# loop function
#check traits

### scenario 2 training set 1 yr BLUPs
#dataset_uom<-dataset_uom%>%filter(GROWSEASON%in%c("2021:01") & TEST_SET_NAME%in%c('21EURP1B2', '21EURP1B4'))

#write_csv2(dataset_uom, 
#                path = "/mnt/veg_blup/traintest_spring_2022yr/dataset_uom.csv")
#write_csv2(dataset_uom_test, 
#                path = "/mnt/veg_blup/traintest_spring_2022yr/results_spring.csv")





#for (i in 1:length(uni_ped)){
#  sub_dataset_ped <- dataset_uom[dataset_uom$PEDIGREE_NAME == uni_ped[i],]
#  sample_size <- ceiling(nrow(sub_dataset_ped)*percent[perc])
#  sample_sub_train <- sample(nrow(sub_dataset_ped), size = sample_size)
#  sample_sub_test <- seq(1:nrow(sub_dataset_ped))[!seq(1:nrow(sub_dataset_ped)) %in% sample_sub_train]
#  sub_dataset_train <- sub_dataset_ped[sample_sub_train,]
#  sub_dataset_test <- sub_dataset_ped[sample_sub_test,]
#  True_values_spring_temp <- sub_dataset_test[,c('PEDIGREE_NAME', 
#                                                 'OBSRVTN_REF_CD', 
#                                                 "TRAIT_VALUE",
#                                                 'EXPER_STAGE_REF_ID',
#                                                 'P1', 'P2', 'LTYPE')]
#  True_values_spring <- rbind(True_values_spring, True_values_spring_temp)
#  sub_dataset_test$TRAIT_VALUE <- NA
#  new_dataset_uom <- rbind(new_dataset_uom, sub_dataset_train, sub_dataset_test)
#}
#dataset_uom <- new_dataset_uom








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
  outbox<-paste0(out0, method1, '_deep', '/',trait=trait_choosen, '_deep')
  #outbox<-paste0(out0, 'ABLUP_deep/', trait=trait_choosen, '_deep')
  trait_GCA<-output_new(trait=trait_choosen,outbox=outbox,out=out0,crop=crop,method="SCA",tmp=tmp,gmethod=method1)
  results = rbind(results,trait_GCA)
  save(results,file="results_asreml_zixuan.Rdata")
  #save(results,file="results_aBLUP.Rdata")
  
}
save(True_values_spring, file="original_asreml_zixuan.Rdata")
}
}}
