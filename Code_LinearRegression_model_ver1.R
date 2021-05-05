
##########################
library(stringr); library(reshape2); library(data.table); library(dplyr);
library(viridis);


memory.limit(size = 24000)

setwd("~/Documents/input/")
outputDirectory='~/Documents/output/'

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
    mutate( cna = c("Loss"=-1, "WT"=0, "Gain"=1)[ str_match( variantType, "[^_]+_([^_]+)")[,2] ] ) %>%
    select(-variantType)
    return(inData)
}




#Model1: 2-way interactions between mutation and CNAs within a cancer gene across cancer types

for(tumor in c('BRCA')){
    
    Input_Input_direc <- paste("~/Documents/input/Model1_MutCNV_int/Loss_Model1_",tumor,"_SingleCancer.txt", sep="")
    
    Input.dataLong = loadInputData( Input_Input_direc)
    Input.data = Input.dataLong %>% reshape2::dcast( Gene+Tissue ~ mut+cna, value.var="N") ## updated because of error message
    
    cat(sprintf( "%s\n", Input_Input_direc ));
    
    Input.data <- fread(file = Input_Input_direc, data.table = F)
    Input.data <- Input.data[,c(1:2,3:6)]
    
    names(Input.data) <- gsub(x = names(Input.data), pattern = "NoMut", replacement = 0) %>%
    gsub(pattern = "Mut", replacement = 1) %>%
    gsub(pattern = "Loss", replacement = -1) %>%
    gsub(pattern = "WT", replacement = 0) %>%
    gsub(pattern = "Gain", replacement = 1) %>%
    gsub(pattern = "Tquery", replacement = 1) %>%
    gsub(pattern = "nonTquery", replacement = 0)
    
    Input.data <- as.data.table(Input.data)
    
    Input.dataLlm = copy(Input.dataLong)
    Input.dataLlm <- as.data.table(Input.dataLlm)
    
    Input.dataLlm[, mut:=as.factor(mut)]
    Input.dataLlm[, cna:=as.ordered(cna)]
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
    
    Input_results <- as.data.frame(Input_results)
    write.table(Input_results %>% mutate(Tissue=tumor,variable=row.names(Input_results)),
    file= paste(outputDirectory,"Test2_Model1_Loss_","_Gene","_",tumor,".txt",sep=''),
    quote=FALSE, sep='\t',row.names=F,col.names = T)
    
}#end of loop for tumor/tissues

##################################################################################################
#TEST2:tissue-specific interactions between mutation and CNAs of a cancer gene between two cancer types
#CancerA = detected cancer type with FDR < 0.2, CancerB = other signficantly mutated cancer type(s)
#Number of CancerB can be different depending on cancer genes

for(tumor in c('BRCA')){
    
    Input_Input_direc <- paste("~/Documents/input/Model2_TissueSpecific_int/Loss_Model2_",tumor,"_vs_Others.txt", sep="")
    
    Input.dataLong = loadInputData(Input_Input_direc)
    Input.data = Input.dataLong %>% reshape2::dcast( Gene+Tissue ~ mut+cna, value.var="N") ## updated because of error message
    
    cat(sprintf( "%s\n", Input_Input_direc ));
    
    Input.data <- fread(file = Input_Input_direc, data.table = F)
    Input.data <- Input.data[,c(1:2,3:6)]
    
    names(Input.data) <- gsub(x = names(Input.data), pattern = "NoMut", replacement = 0) %>%
    gsub(pattern = "Mut", replacement = 1) %>%
    gsub(pattern = "Loss", replacement = -1) %>%
    gsub(pattern = "WT", replacement = 0) %>%
    gsub(pattern = "Gain", replacement = 1) %>%
    gsub(pattern = "Tquery", replacement = 1) %>%
    gsub(pattern = "nonTquery", replacement = 0)
    
    Input.data <- as.data.table(Input.data)
    
    Input.dataLlm = copy(Input.dataLong)
    Input.dataLlm <- as.data.table(Input.dataLlm)
    
    Input.dataLlm[, mut:=as.factor(mut)]
    Input.dataLlm[, cna:=as.ordered(cna)]
    
    Input_results <- c()
    
    #####################
    cat(sprintf( "Tissue %s\n", tumor ));
    for ( aGene in Input.dataLlm[, unique(Gene)] %>% setdiff("_allCancerGenes") ) {
        
        aLlm <- glm(N ~ mut + cna + Tissue+ mut:cna + mut:Tissue + cna:Tissue + mut:cna:Tissue, family=poisson(link="log"),
        control=glm.control(epsilon = 1e-6,maxit=100), data=Input.dataLlm[Gene == aGene])
        
        aLlCoef = summary(aLlm) %>% coef
        row.names(aLlCoef)[1] <- paste('intercept_',aGene,sep='')
        Input_results <-Input_results %>%rbind(aLlCoef)
    }
    
    Input_results <- as.data.frame(Input_results)
    write.table(Input_results %>% mutate(Tissue=tumor,variable=row.names(Input_results)),
    file= paste(outputDirectory,"Model2_Loss_","_Gene","_",tumor,".txt",sep=''),
    quote=FALSE, sep='\t',row.names=F,col.names = T)
    
}


##################################################################################################
#Test3 :3-way interactions between two genes with three genomics alterations across cancer type ####


loadInputData = function( fn  ) {
    
    inData = fread(fn, data.table = T)  # glimpse(inData)
    
    inData = inData[, colnames(inData) %>% str_subset( "(Mut|NoMut)_(Del|Loss|WT|Gain|Amp)$") %>% c("Gene","Target", .), with=F ]
    
    inData =
    inData %>% melt( id.vars=c("Gene", "Target"), variable.name="variantType", value.name="N" ) %>%
    mutate( mut = c("Mut"=1, "NoMut"=0)[ str_match( variantType, "([^_]+)_[^_]+" )[,2] ] ) %>%
    mutate( cna = c("Loss"=-1, "WT"=0, "Gain"=1)[ str_match( variantType, "[^_]+_([^_]+)")[,2] ] ) %>%
    select(-variantType)
    
    return(inData)
}


for(tumor in c('BRCA')){
    
    Input_Input_direc <- paste("~/Documents/input/Model3_ThreeWay_int/Loss_Model3_",tumor,"_SingleCancer.txt", sep="")
    Input.dataLong = loadInputData(Input_Input_direc)
    Input.data = Input.dataLong %>% reshape2::dcast( Gene+Target ~ mut+cna, value.var="N")
    
    cat(sprintf( "%s\n", Input_Input_direc ));
    Input.data <- fread(file = Input_Input_direc, data.table = F)
    

    Input.data <- Input.data[,c(1:2,3:6)]
    
    names(Input.data) <- gsub(x = names(Input.data), pattern = "NoMut", replacement = 0) %>%
    gsub(pattern = "Mut", replacement = 1) %>%
    gsub(pattern = "Loss", replacement = -1) %>%
    gsub(pattern = "WT", replacement = 0) %>%
    gsub(pattern = "Gain", replacement = 1) %>%
    gsub(pattern = "YesSecondG", replacement = 1) %>%
    gsub(pattern = "NoSecondG", replacement = 0)
    
    Input.data <- as.data.table(Input.data)
    
    Input.dataLlm = copy(Input.dataLong)
    Input.dataLlm <- as.data.table(Input.dataLlm)
    
    Input.dataLlm[, mut:=as.factor(mut)]
    Input.dataLlm[, cna:=as.ordered(cna)]
    Input_results <- c()
    
    #####################
    cat(sprintf( "Tissue %s\n", tumor ));
    for ( aGene in Input.dataLlm[, unique(Gene)] %>% setdiff("_allCancerGenes") ) {
        
        #cat(sprintf( "Gene %s\n", aGene ));
        
        aLlm <- glm(N ~ mut +cna + Target + cna:mut + cna:Target + mut:Target + cna:mut:Target, family=poisson(link="log"),
        control=glm.control(epsilon = 1e-6,maxit=100), data=Input.dataLlm[Gene == aGene])
        
        aLlCoef = summary(aLlm) %>% coef
        row.names(aLlCoef)[1] <- paste('intercept_',aGene,sep='')
        Input_results <-Input_results %>%rbind(aLlCoef)
    }
    
    Input_results <- as.data.frame(Input_results)
    write.table(Input_results %>% mutate(Tissue=tumor,variable=row.names(Input_results)),
    file= paste(outputDirectory,"Model3_Loss_","_Gene","_",tumor,".txt",sep=''),
    quote=FALSE, sep='\t',row.names=F,col.names = T)
}


