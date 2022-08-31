##########
# ANY Environment without ASREML

install.packages('azurequest',repos=c('https://cran.science-at-scale.io', 'https://cran.rstudio.com'))
install.packages('pedigree',repos=c('https://cran.science-at-scale.io', 'https://cran.rstudio.com'))


# This is an example to run AoB with data pulled from H2H

# Dynamic way of checking required libraries in R
list.of.packages <- c("aws.s3", "httr", "tidyverse", "measurements", "reshape", "tidyr", 
                      "jsonlite", "dplyr", "readr","data.table","reshape2","azurequest","pedigree")
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


############################################## Input parameters

args = commandArgs(trailingOnly = TRUE)


repeat_times <- as.numeric(args[1])
percent <- sapply(unlist(lapply(strsplit(args[2], ','), trimws, 'both')),as.numeric)
BLUPs <- unlist(lapply(strsplit(args[3], ','), trimws, 'both'))
out0_input <- args[4] 
TRAITS <- unlist(lapply(strsplit(args[5], ','), trimws, 'both'))
seed <- as.numeric(args[6])
clusters <- as.numeric(args[7])


#repeat_times <- 1
#percent <- 0.5
#BLUPs <- 'SSGBLUP'
#out0_input <- '/mnt/PCA/Spring_test/'
#TRAITS <- c("AFW_C")
#seed <- 129
#clusters <- 4


###########################################
## read genotypes 


load("/mnt/5gen_inbreeding_Cucumber.Rdata")

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

genotype_data <- rbind(gg1, gg2, gg3)


#######################################
## read phenotypes and filter out the pedigree without genotypes

s3load(object = '/ycao1/ITG/Cucumber/Spring/Cln_cucumber_spring04272022.Rdata', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')


genotype_data <- genotype_data%>%
  drop_na() %>%
  distinct() %>%
  dplyr::filter(!duplicated(PEDIGREE))



dataset<-pheno

## UOM Conversion
UOMTable <- read.csv("/mnt/MIDAS_UOM_Table.csv", colClasses = rep(times=5,x="character"))
dataset$UOM <- UOMTable$NAME[match(dataset$UOM,UOMTable$UOM_ID)]
dataset_uom <- conv_uom(dataset) 
dataset_2020_to_2022_with_geno <- dataset_uom%>%dplyr::filter(GROWSEASON %in% 
                                                                c('2020:01', 
                                                                  '2021:01', '2022:01') &
                                                                PEDIGREE_NAME %in% 
                                                                genotype_data$PEDIGREE)
phenotype_data <- dataset_2020_to_2022_with_geno
genotype_data <- genotype_data %>% dplyr::filter(PEDIGREE %in% phenotype_data$PEDIGREE_NAME)
genotype_data_copy <- genotype_data








out0 <- out0_input
if(file.exists(out0)==F) {dir.create(out0, recursive = T)}

#### Perform PCA to visualize the genotypes

pca_fit <-  genotype_data %>% 
  select(-c("PEDIGREE")) %>% 
  prcomp(rank. = 3)  # compute just 3 components, no need for more
save(pca_fit, file = paste0(out0_input, 'PCA_FIT.Rdata'))




row.names(genotype_data) <- genotype_data$PEDIGREE
genotype_data <- genotype_data %>% select(-c("PEDIGREE"))

#### Apply K-means to cluster the pedigrees

set.seed(135)

kmeans_geno <- genotype_data %>%   
  kmeans(centers = clusters, nstart = 5)



genotype_data$cluster <- as.factor(kmeans_geno$cluster)
genotype_data <- cbind(genotype_data, pca_fit$x)


genotype_data_copy <- genotype_data

#### Visualize the genotypes

genotype_data %>%
  ggplot(mapping = aes(x = PC1, y = PC2, color = cluster)) + geom_point(alpha = 0.7)



genotype_data %>%
  ggplot(mapping = aes(x = PC2, y = PC3, color = cluster)) + geom_point(alpha = 0.7)




genotype_data %>%
  ggplot(mapping = aes(x = PC1, y = PC3, color = cluster)) + geom_point(alpha = 0.7)

#### All the possible combinations.

combination_matrix <- combn(1:clusters, m = 2)


#pca_fit %>% broom::tidy(matrix = 'd')

set.seed(seed)
for(replication in 1:repeat_times){
  for (perc in 1:length(percent)){
    for (j in 1:length(BLUPs)){
      for (class in 1:ncol(combination_matrix)){
        
        
        
        out0 <- paste0(out0_input, combination_matrix[1,class], '_to_predict_', 
                       combination_matrix[2,class], '/', percent[perc]*100, '/',replication, '/')
        if(file.exists(out0)==F) {dir.create(out0, recursive = T)}
        
        
        genotype_data <- genotype_data_copy # Refresh with the original genotypic data.
        genotype_data$PEDIGREE <- row.names(genotype_data)
        dataset_uom <- phenotype_data # Refresh with the original phenotypic data.
        
        #### Use combination_matrix[1,class] to predict combination_matrix[2,class]
        
        between_cluster <- genotype_data[genotype_data$cluster == combination_matrix[1,class],]
        Test_cluster <- genotype_data[genotype_data$cluster == combination_matrix[2,class],]
        
        #### Construct training set
        
        uni_pedi <- unique(between_cluster$PEDIGREE)
        candidate_set <- dataset_uom[dataset_uom$PEDIGREE_NAME %in% uni_pedi,]
        sample_size <- ceiling(length(uni_pedi) * percent[perc])
        sample_training <- sample(uni_pedi, sample_size)
        
        #### Save useful information of Test set.
        
        dataset_uom_test <- dataset_uom[dataset_uom$OBSRVTN_REF_CD %in% TRAITS &
                                          dataset_uom$PEDIGREE_NAME %in% Test_cluster$PEDIGREE,]
        True_values <- dataset_uom_test[,c('PEDIGREE_NAME', 
                                           'OBSRVTN_REF_CD', 
                                           "TRAIT_VALUE",
                                           'EXPER_STAGE_REF_ID',
                                           'P1', 'P2', 'LTYPE', 'ORIGIN', 'FIELD_NAME', 'IS_CHK')]
        
        dataset_uom_test$TRAIT_VALUE <- NA
        
        #### Save useful information of Training set.
        
        dataset_uom_training <- candidate_set[candidate_set$OBSRVTN_REF_CD %in% TRAITS &
                                                candidate_set$PEDIGREE_NAME %in% sample_training,]
        Training_main_info <- dataset_uom_training[, c('PEDIGREE_NAME', 
                                                       'OBSRVTN_REF_CD', 
                                                       "TRAIT_VALUE",
                                                       'EXPER_STAGE_REF_ID',
                                                       'P1', 'P2', 'LTYPE', 'ORIGIN', 'FIELD_NAME', 'IS_CHK')]
        dataset_uom <- rbind(dataset_uom_training, dataset_uom_test)
        
        
        
        
        
        results <- data.frame()
        
        for (i in 1:length(TRAITS)){
          #for (i in 1:3){
          trait_choosen = TRAITS[i]
          print(paste(" ################--- ", trait_choosen,"__is completed ############"))
          method1 = 'SSGBLUP'
          if(method1 == 'SSGBLUP'){
            res = getoutput_new(out=out0,dataset=dataset_uom,
                                trait = trait_choosen,gg=genotype_data,n_gen=5,seas=NULL,method=method1)
          }
          if(method1 == 'ABLUP'){
            res = getoutput_new(out=out0,dataset=dataset_uom,
                                trait = trait_choosen,gg=NULL,n_gen=5,seas=NULL,method=method1)
          }
          if(method1 == 'GBLUP'){
            res = getoutput_new(out=out0,dataset=dataset_uom,
                                trait = trait_choosen,gg=genotype_data,n_gen=5,seas=NULL,method=method1)
          }
          outbox<-paste0(out0, method1, '_deep', '/',trait=trait_choosen, '_deep')
          #outbox<-paste0(out0, 'ABLUP_deep/', trait=trait_choosen, '_deep')
          trait_GCA<-output_new(trait=trait_choosen,outbox=outbox,out=out0,
                                crop=crop,method="SCA",tmp=tmp,gmethod=method1)
          results = rbind(results,trait_GCA)
          #
          #save(results,file="results_aBLUP.Rdata")
          
        }
        
        #### Delete processing results
        unlink(paste0(out0, method1,'_deep/'), recursive = TRUE)
        
        ## Save results
        
        out_results <- paste0(out0, 'results/')
        if(file.exists(out_results)==F) {dir.create(out_results, recursive = T)}
        setwd(out_results)
        
        save(genotype_data, file = 'genotypic_data.Rdata')
        save(results,file="results_asreml_zixuan.Rdata")
        save(Training_main_info, file = 'Training_data.Rdata')
        save(True_values, file = "original_asreml_zixuan.Rdata")

        
        
        
        ############################## Perform the same thing just inverse the Training cluster and test cluster.
        
        
        
        out0 <- paste0(out0_input, combination_matrix[2,class], '_to_predict_', 
                       combination_matrix[1,class], '/', percent[perc]*100, '/',replication, '/')
        if(file.exists(out0)==F) {dir.create(out0, recursive = T)}
        
        #TRAITS <- c("EXTCO","FRLGT")#,"AFW_C","FRNMK","SHAPE")
        genotype_data <- genotype_data_copy
        genotype_data$PEDIGREE <- row.names(genotype_data)
        dataset_uom <- phenotype_data
        between_cluster <- genotype_data[genotype_data$cluster == combination_matrix[2,class],]
        Test_cluster <- genotype_data[genotype_data$cluster == combination_matrix[1,class],]
        uni_pedi <- unique(between_cluster$PEDIGREE)
        candidate_set <- dataset_uom[dataset_uom$PEDIGREE_NAME %in% uni_pedi,]
        sample_size <- ceiling(length(uni_pedi) * percent[perc])
        sample_training <- sample(uni_pedi, sample_size)
        dataset_uom_test <- dataset_uom[dataset_uom$OBSRVTN_REF_CD %in% TRAITS &
                                          dataset_uom$PEDIGREE_NAME %in% Test_cluster$PEDIGREE,]
        True_values <- dataset_uom_test[,c('PEDIGREE_NAME', 
                                           'OBSRVTN_REF_CD', 
                                           "TRAIT_VALUE",
                                           'EXPER_STAGE_REF_ID',
                                           'P1', 'P2', 'LTYPE', 'ORIGIN', 'FIELD_NAME', 'IS_CHK')]
        
        dataset_uom_test$TRAIT_VALUE <- NA
        
        dataset_uom_training <- candidate_set[candidate_set$OBSRVTN_REF_CD %in% TRAITS &
                                                candidate_set$PEDIGREE_NAME %in% sample_training,]
        Training_main_info <- dataset_uom_training[, c('PEDIGREE_NAME', 
                                                       'OBSRVTN_REF_CD', 
                                                       "TRAIT_VALUE",
                                                       'EXPER_STAGE_REF_ID',
                                                       'P1', 'P2', 'LTYPE', 'ORIGIN', 'FIELD_NAME', 'IS_CHK')]
        dataset_uom <- rbind(dataset_uom_training, dataset_uom_test)
        
        
        
        
        
        results <- data.frame()
        
        for (i in 1:length(TRAITS)){
          #for (i in 1:3){
          trait_choosen = TRAITS[i]
          print(paste(" ################--- ", trait_choosen,"__is completed ############"))
          method1 = 'SSGBLUP'
          if(method1 == 'SSGBLUP'){
            res = getoutput_new(out=out0,dataset=dataset_uom,
                                trait = trait_choosen,gg=genotype_data,n_gen=5,seas=NULL,method=method1)
          }
          if(method1 == 'ABLUP'){
            res = getoutput_new(out=out0,dataset=dataset_uom,
                                trait = trait_choosen,gg=NULL,n_gen=5,seas=NULL,method=method1)
          }
          if(method1 == 'GBLUP'){
            res = getoutput_new(out=out0,dataset=dataset_uom,
                                trait = trait_choosen,gg=genotype_data,n_gen=5,seas=NULL,method=method1)
          }
          outbox<-paste0(out0, method1, '_deep', '/',trait=trait_choosen, '_deep')
          #outbox<-paste0(out0, 'ABLUP_deep/', trait=trait_choosen, '_deep')
          trait_GCA<-output_new(trait=trait_choosen,outbox=outbox,out=out0,
                                crop=crop,method="SCA",tmp=tmp,gmethod=method1)
          results = rbind(results,trait_GCA)
          #
          #save(results,file="results_aBLUP.Rdata")
          
        }
        unlink(paste0(out0, method1,'_deep/'), recursive = TRUE)
        
        out_results <- paste0(out0, 'results/')
        if(file.exists(out_results)==F) {dir.create(out_results, recursive = T)}
        setwd(out_results)
        
        save(genotype_data, file = 'genotypic_data.Rdata')
        save(results,file="results_asreml_zixuan.Rdata")
        save(Training_main_info, file = 'Training_data.Rdata')
        save(True_values, file = "original_asreml_zixuan.Rdata")        
        
      }
    }
  }
}
