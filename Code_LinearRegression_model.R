
##########################
library(stringr); library(reshape2); library(data.table); library(dplyr);
library(viridis);


memory.limit(size = 24000)

setwd("~/Documents/LLR/TEST_output/")

`%+%` = function(a, b) paste0(a, b)
clipOut = function( x ) { write.table(x,"clipboard-16384", sep="\t", row.names=F, quote=F) }
clipIn = function( ) { read.table("clipboard", header=TRUE, na.strings=c("NA","NaN","?"), sep="\t"); }
numNas = function( dt ) { sapply( pcaPointsForIndiv, function(x){sum(is.na(x))} )   }
cosineDist <- function(x){ as.dist(1 - x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))) }
mbind<-function(...){
    Reduce( function(x,y){cbind(x,y[match(row.names(x),row.names(y)),])}, list(...) )
}

# ============================================================

loadInputData = function( fn  ) {
    
    inData = fread(fn, data.table = T)  # glimpse(inData)
    
    inData = inData[, colnames(inData) %>% str_subset( "(Mut|NoMut)_(Del|Loss|WT|Gain|Amp)$") %>% c("Gene","Tissue", .), with=F ]
    
    inData =
    inData %>% melt( id.vars=c("Gene", "Tissue"), variable.name="variantType", value.name="N" ) %>%
    mutate( mut = c("Mut"=1, "NoMut"=0)[ str_match( variantType, "([^_]+)_[^_]+" )[,2] ] ) %>%
    mutate( cna = c("Del"=-2, "Loss"=-1, "WT"=0, "Gain"=1, "Amp"=2)[ str_match( variantType, "[^_]+_([^_]+)")[,2] ] ) %>%
    select(-variantType)
    return(inData)
}



Input.dataLong = loadInputData( "~/Documents/LLR/TEST_input" %+%
"/Gain_TCGA2018_TSG_OG_Value_2_vs_5_table_Targetingtissues_KICH_non.txt")

Input.data = Input.dataLong %>% reshape2::dcast( Gene+Tissue ~ mut+cna, value.var="N") ## updated because of error message

####SK: =============================================================
Input.data <- fread("~/Documents/LLR/TEST_input/Gain_TCGA2018_TSG_OG_Value_2_vs_5_table_Targetingtissues_KICH_non.txt", data.table = F)


colnames(Input.data)
Input.data <- Input.data[,c(1:2,12:15)]

names(Input.data) <- gsub(x = names(Input.data), pattern = "NoMut", replacement = 0) %>%
gsub(pattern = "Mut", replacement = 1) %>%
gsub(pattern = "Loss", replacement = -1) %>%
gsub(pattern = "WT", replacement = 0) %>%
gsub(pattern = "Gain", replacement = 1) %>%
gsub(pattern = "Tquery", replacement = 1) %>%
gsub(pattern = "nonTquery", replacement = 0)

Input.data <- as.data.table(Input.data)

Input.dataLong[]

Input.dataLlm = copy(Input.dataLong)
Input.dataLlm <- as.data.table(Input.dataLlm)

Input.dataLlm[, mut:=as.factor(mut)]
Input.dataLlm[, cna:=as.ordered(cna)]
str(Input.dataLlm)
glimpse(Input.dataLlm)

#########################################
outputDirectory='~/Documents/LLR/TEST_output/'

loadInputData = function( fn  ) {
    
    inData = fread(fn, data.table = T)  # glimpse(inData)
    
    inData = inData[, colnames(inData) %>% str_subset( "(Mut|NoMut)_(Del|Loss|WT|Gain|Amp)$") %>% c("Gene","Tissue", .), with=F ]
    
    inData =
    inData %>% melt( id.vars=c("Gene", "Tissue"), variable.name="variantType", value.name="N" ) %>%
    mutate( mut = c("Mut"=1, "NoMut"=0)[ str_match( variantType, "([^_]+)_[^_]+" )[,2] ] ) %>%
    mutate( cna = c("Del"=-2, "Loss"=-1, "WT"=0, "Gain"=1, "Amp"=2)[ str_match( variantType, "[^_]+_([^_]+)")[,2] ] ) %>%
    select(-variantType)
    
    return(inData)
}



#Test1 :2-way interactions between mutation and CNAs across gene-tissue pairs

Input_Input_direc <- paste("~/Documents/LLR/TEST_input/CNV_Loss_LLR_input_type1.txt", sep="")

Input.dataLong = loadInputData( Input_Input_direc)
Input.data = Input.dataLong %>% reshape2::dcast( Gene+Tissue ~ mut+cna, value.var="N") ## updated because of error message

cat(sprintf( "%s\n", Input_Input_direc ));

Input.data <- fread(file = Input_Input_direc, data.table = F)

colnames(Input.data)
Input.data <- Input.data[,c(1:2,12:15)]

names(Input.data) <- gsub(x = names(Input.data), pattern = "NoMut", replacement = 0) %>%
gsub(pattern = "Mut", replacement = 1) %>%
gsub(pattern = "Loss", replacement = -1) %>%
gsub(pattern = "WT", replacement = 0) %>%
gsub(pattern = "Gain", replacement = 1) %>%
gsub(pattern = "Tquery", replacement = 1) %>%
gsub(pattern = "nonTquery", replacement = 0)

Input.data <- as.data.table(Input.data)

Input.dataLong[]

Input.dataLlm = copy(Input.dataLong)
Input.dataLlm <- as.data.table(Input.dataLlm)

Input.dataLlm[, mut:=as.factor(mut)]
Input.dataLlm[, cna:=as.ordered(cna)]
str(Input.dataLlm)
glimpse(Input.dataLlm)

Input_results <- c()

#####################
cat(sprintf( "Tissue %s\n", tumor ));
for ( aGene in Input.dataLlm[, unique(Gene)] %>% setdiff("_allCancerGenes") ) {
    
    
    aLlm <- glm(N ~ mut + cna + mut:cna, family=poisson(link="log"),
    control=glm.control(epsilon = 1e-6,maxit=100), data=Input.dataLlm[Gene == aGene])
    
    aLlCoef = summary(aLlm) %>% coef
    row.names(aLlCoef)[1] <- paste('intercept_',aGene,sep='')
    Input_results <-Input_results %>%rbind(aLlCoef)
}

Input_results <- as.data.frame(Input_results) #make data frame out of it
write.table(Input_results %>% mutate(Tissue=tumor,variable=row.names(Input_results)),
file= paste(outputDirectory,"CNV_Loss_LLR_output_type1.txt",sep=''),
quote=FALSE, sep='\t',row.names=F,col.names = T)


#Test2 :tissue-specific interactions between mutation and CNAs of a cancer gene between two cancer types

Input_Input_direc <- paste("~/Documents/LLR/TEST_input/CNV_Loss_LLR_input_type2.txt", sep="")

Input.dataLong = loadInputData( Input_Input_direc)
Input.data = Input.dataLong %>% reshape2::dcast( Gene+Tissue ~ mut+cna, value.var="N") ## updated because of error message

cat(sprintf( "%s\n", Input_Input_direc ));

Input.data <- fread(file = Input_Input_direc, data.table = F)

colnames(Input.data)
Input.data <- Input.data[,c(1:2,3:6)]

names(Input.data) <- gsub(x = names(Input.data), pattern = "NoMut", replacement = 0) %>%
gsub(pattern = "Mut", replacement = 1) %>%
gsub(pattern = "Loss", replacement = -1) %>%
gsub(pattern = "WT", replacement = 0) %>%
gsub(pattern = "Gain", replacement = 1) %>%
gsub(pattern = "Tquery", replacement = 1) %>%
gsub(pattern = "nonTquery", replacement = 0)

Input.data <- as.data.table(Input.data)

Input.dataLong[]
Input.dataLlm = copy(Input.dataLong)
Input.dataLlm <- as.data.table(Input.dataLlm)

Input.dataLlm[, mut:=as.factor(mut)]
Input.dataLlm[, cna:=as.ordered(cna)]
str(Input.dataLlm)
glimpse(Input.dataLlm)

Input_results <- c()

#####################
cat(sprintf( "Tissue %s\n", tumor ));
for ( aGene in Input.dataLlm[, unique(Gene)] %>% setdiff("_allCancerGenes") ) {
    
    
    aLlm <- glm(N ~ mut + cna + Tissue+ mut:cna + mut:Tissue + cna:Tissue + mut:cna:Tissue, family=poisson(link="log"),
    control=glm.control(epsilon = 1e-6,maxit=100), data=solip.dataLlm[Gene == aGene])
    
    aLlCoef = summary(aLlm) %>% coef
    row.names(aLlCoef)[1] <- paste('intercept_',aGene,sep='')
    Input_results <-Input_results %>%rbind(aLlCoef)
}

Input_results <- as.data.frame(Input_results) #make data frame out of it
write.table(Input_results %>% mutate(Tissue=tumor,variable=row.names(Input_results)),
file= paste(outputDirectory,"CNV_Loss_LLR_output_type2.txt",sep=''),
quote=FALSE, sep='\t',row.names=F,col.names = T)


#Test3 :3-way interactions between two genes with three genomics alterations

Input_Input_direc <- paste("~/Documents/LLR/TEST_input/CNV_Loss_LLR_input_type3.txt", sep="")

Input.dataLong = loadInputData( Input_Input_direc)
Input.data = Input.dataLong %>% reshape2::dcast( Gene+Tissue ~ mut+cna, value.var="N") ## updated because of error message

cat(sprintf( "%s\n", Input_Input_direc ));

Input.data <- fread(file = Input_Input_direc, data.table = F)

colnames(Input.data)
Input.data <- Input.data[,c(1:2,3:6)]


names(Input.data) <- gsub(x = names(Input.data), pattern = "NoMut", replacement = 0) %>%
gsub(pattern = "Mut", replacement = 1) %>%
gsub(pattern = "Loss", replacement = -1) %>%
gsub(pattern = "WT", replacement = 0) %>%
gsub(pattern = "Gain", replacement = 1) %>%
gsub(pattern = "TidsA", replacement = 1) %>%
gsub(pattern = "TidsB", replacement = 0) %>%

Input.data <- as.data.table(Input.data)

Input.dataLong[]

Input.dataLlm = copy(Input.dataLong)
Input.dataLlm <- as.data.table(Input.dataLlm)

Input.dataLlm[, mut:=as.factor(mut)]
Input.dataLlm[, cna:=as.ordered(cna)]
str(Input.dataLlm)
glimpse(Input.dataLlm)

Input_results <- c()

#####################
cat(sprintf( "Tissue %s\n", tumor ));
for ( aGene in Input.dataLlm[, unique(Gene)] %>% setdiff("_allCancerGenes") ) {
    
    
    aLlm <- glm(N ~ mut +cna + Target + cna:mut + cna:Target + mut:Target + cna:mut:Target, family=poisson(link="log"),
    control=glm.control(epsilon = 1e-6,maxit=100), data=solip.dataLlm[Gene == aGene])
    
    aLlCoef = summary(aLlm) %>% coef
    row.names(aLlCoef)[1] <- paste('intercept_',aGene,sep='')
    Input_results <-Input_results %>%rbind(aLlCoef)
}

Input_results <- as.data.frame(Input_results) #make data frame out of it
write.table(Input_results %>% mutate(Tissue=tumor,variable=row.names(Input_results)),
file= paste(outputDirectory,"CNV_Loss_LLR_output_type3.txt",sep=''),
quote=FALSE, sep='\t',row.names=F,col.names = T)


