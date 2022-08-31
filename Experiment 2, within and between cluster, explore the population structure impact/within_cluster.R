##########
# @ ENVIRONMENT : Any environment without ASREML

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

 
########################### Input parameters

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


###########################
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


###########################
## load phenotypes and filter the pedigree without genotypes

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


#########################
## Perform PCA to visualize the genotype data.


pca_fit <-  genotype_data %>% 
  select(-c("PEDIGREE")) %>% 
  prcomp(rank. = 3)  # compute just 3 components, no need for more
save(pca_fit, file = paste0(out0_input, 'PCA_FIT.Rdata'))




row.names(genotype_data) <- genotype_data$PEDIGREE
genotype_data <- genotype_data %>% select(-c("PEDIGREE"))

set.seed(135)

#####################
## Clustering the genotype with K-means

kmeans_geno <- genotype_data %>%   
  kmeans(centers = clusters, nstart = 5)



genotype_data$cluster <- as.factor(kmeans_geno$cluster)
genotype_data <- cbind(genotype_data, pca_fit$x)

####################### 
## PCA details
explanation_pca_fit <- pca_fit %>% broom::tidy(matrix = 'd')



######################
## Visualization of genotypes based on clusters
genotype_data %>%
  ggplot(mapping = aes(x = PC1, y = PC2, color = cluster)) + geom_point(alpha = 0.7) + 
  labs(x = paste0('PC1(', 100*explanation_pca_fit$cumulative[1], '%)'),
       y = paste0('PC2(', 100*(explanation_pca_fit$cumulative[2] - explanation_pca_fit$cumulative[1]), '%)'))



genotype_data %>%
  ggplot(mapping = aes(x = PC2, y = PC3, color = cluster)) + geom_point(alpha = 0.7) +
  labs(x = paste0('PC2(', 100*(explanation_pca_fit$cumulative[2] - explanation_pca_fit$cumulative[1]), '%)'),
       y = paste0('PC3(', 100*(explanation_pca_fit$cumulative[3] - explanation_pca_fit$cumulative[2]), '%)'))



genotype_data %>%
  ggplot(mapping = aes(x = PC1, y = PC3, color = cluster)) + geom_point(alpha = 0.7) + 
  labs(x = paste0('PC1(', 100*(explanation_pca_fit$cumulative[1]), '%)'),
       y = paste0('PC3(', 100*(explanation_pca_fit$cumulative[3] - explanation_pca_fit$cumulative[2]), '%)'))


genotype_data$PEDIGREE <- row.names(genotype_data)
genotype_data_tester <- left_join(genotype_data, unique(phenotype_data[,c('PEDIGREE_NAME', 'P2')]), by = c('PEDIGREE' = 'PEDIGREE_NAME'))
genotype_data_tester_visual <- genotype_data_tester[genotype_data_tester$P2 %in% names(table(genotype_data_tester$P2)[table(genotype_data_tester$P2) > 10]),]



genotype_data_tester_visual %>%
  ggplot(mapping = aes(x = PC1, y = PC2, color = P2)) + geom_point(alpha = 0.7) + 
  labs(x = paste0('PC1(', 100*explanation_pca_fit$cumulative[1], '%)'),
       y = paste0('PC2(', 100*(explanation_pca_fit$cumulative[2] - explanation_pca_fit$cumulative[1]), '%)'))






set.seed(seed)


for(replication in 1:repeat_times){
  for (perc in 1:length(percent)){
    for (j in 1:length(BLUPs)){
      for (class in 1:clusters){
        
        
        
        out0 <- paste0(out0_input, class, '/', percent[perc]*100, '/',replication, '/')
        if(file.exists(out0)==F) {dir.create(out0, recursive = T)}
        
        
        genotype_data$PEDIGREE <- row.names(genotype_data)
        dataset_uom <- phenotype_data
        within_cluster <- genotype_data[genotype_data$cluster == class,]
        
        #### Randomly split within cluster based on Pedigrees
        
        uni_pedi <- unique(within_cluster$PEDIGREE)
        dataset_uom <- dataset_uom[dataset_uom$PEDIGREE_NAME %in% uni_pedi,]
        sample_size <- ceiling(length(uni_pedi) * percent[perc])
        sample_training <- sample(uni_pedi, sample_size)
        dataset_uom_test <- dataset_uom[dataset_uom$OBSRVTN_REF_CD %in% TRAITS &
                                          !dataset_uom$PEDIGREE_NAME %in% sample_training,]
        
        #### Save the useful information of test set
        
        True_values <- dataset_uom_test[,c('PEDIGREE_NAME', 
                                           'OBSRVTN_REF_CD', 
                                           "TRAIT_VALUE",
                                           'EXPER_STAGE_REF_ID',
                                           'P1', 'P2', 'LTYPE', 'ORIGIN', 'FIELD_NAME', 'IS_CHK')]
        
        dataset_uom_test$TRAIT_VALUE <- NA
        
        dataset_uom_training <- dataset_uom[dataset_uom$OBSRVTN_REF_CD %in% TRAITS &
                                              dataset_uom$PEDIGREE_NAME %in% sample_training,]
        
        #### Save the useful information of training set 
        
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
        
        #### Delete processing files
        
        unlink(paste0(out0, method1,'_deep/'), recursive = TRUE)
        
        #### Save results
        
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

