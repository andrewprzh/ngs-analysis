#!/bin/R
# By Anoushka Joglekar 11.2022
# Based on old code but tweaks implemented for new data format

## Takes in a list of exon inclusion - exclusion counts and performs
## chi-sq tests by specified groups 

# Setup -----
library(data.table)
library(dplyr)
library(tidyr)
library(parallel)

# Reading in arguments and input ------
## Argument list:
# 1. path to dataframe of exon counts
# 2. level names in "cell type" string. Should use "::" as delimiter
# 3. Either the name of level for which you want 1 vs. all comparison, or path to config file
# 4. path to output folder wrt current directory
# 5. number of threads to use when parallelizing

args <- commandArgs(trailing=TRUE)
## one vs. all example
args <- c('CROPseq.ONT.exon_counts.formatted.tsv','levelNames','GeneName','oneVSee','6')
## two sample example
args <- c('CROPseq.ONT.exon_counts.formatted.tsv','levelNames','vsControl.config','allVSnegative','6') 

levelNames <- as.vector(read.table(args[2],header = FALSE)$V1)

countData <- fread(args[1])
colnames(countData) <- c("GE","Inclusion","Exclusion","Total",levelNames)


if(file.exists(args[3])){
	configFile <- read.table(args[3],sep = "\t")
	colnames(configFile) <- c("groupingFactor","group1","group2")
	comparisonType <- "twoConditions"
} else {
	groupingFactor <- args[3]
	comparisonType <- "oneVsAll"
}

outputDir <- args[4]
if(!dir.exists(outputDir)){dir.create(outputDir)}

numThreads <- as.integer(args[5])



# Define functions -------
getProcessedDF <- function(countDF_perCond,gName){ ## one vs. all
  threshold <- 10
	countDF_perCond <- countDF_perCond %>% mutate_at(vars(groupingFactor), 
		~ case_when(. == gName ~ "Group1", TRUE ~ "Group2" ) )
	processedDF_perCond <- countDF_perCond %>% select(toInclude[1:4],all_of(toInclude[-c(1:4)])) %>% 
		group_by(GE,!!sym(groupingFactor)) %>%
		summarise(Inclusion = sum(Inclusion),Exclusion = sum(Exclusion), 
			Total = sum(Total), .groups = "keep") %>% 
		filter(sum(Total) >= threshold) %>% ## mostly not required. Can uncomment to minimize number of tests
		group_by(GE) %>% filter(n() > 1 ) %>% ungroup()
	return(processedDF_perCond)
}

getProcessedDF_twoSample <- function(countDF_perCond,cond1,cond2){
  threshold <- 10
	cond2 = unlist(strsplit(cond2,","))
	countDF_perCond <- countData %>% mutate_at(vars(groupingFactor), ~ case_when(. == cond1 ~ "Group1", 
					. %in% cond2 ~ "Group2", TRUE ~ "Discard" ) ) %>% 
		filter(!!sym(groupingFactor) != "Discard")
        processedDF_perCond <- countDF_perCond %>% select(toInclude[1:4],all_of(toInclude[-c(1:4)])) %>%
               	group_by(GE,!!sym(groupingFactor)) %>%
                summarise(Inclusion = sum(Inclusion),Exclusion = sum(Exclusion),
                        Total = sum(Total), .groups = "keep") %>%
                filter(sum(Total) >= threshold) %>% ## mostly not required. Can uncomment to minimize number of tests
                group_by(GE) %>% filter(n() > 1 ) %>% ungroup()
        return(processedDF_perCond)
}

checkAndCompute <- function(inputMat,ix){	## common for both
	mat <- inputMat %>% filter(GE == ix) %>%
		select(Inclusion,Exclusion) %>% as.matrix()
	if (sum(mat) > 0 & max(colSums(mat))/sum(mat) <= 0.9 &
	min(rowSums(mat))*min(colSums(mat))/sum(mat) > 5){
		exonL <- as.matrix(inputMat %>% filter(GE == ix) %>% select(GE))[1]
		psis <- mat[,1]/rowSums(mat)
		psi1 <- c(psi1,psis[1])
		psi2 <- c(psi2,psis[2])
		pval <- c(pval,chisq.test(mat)$p.value)
		dpsi <- psis[1]-psis[2]
	}
	return(list(exonL,pval,dpsi,psi1,psi2))}

computeResPerCond <- function(countDF,gName){	## one vs. all 
	processedDF_perCond <- getProcessedDF(countDF,gName)
	exonList <- unique(processedDF_perCond$GE)
	res <- lapply(exonList, function(ix) checkAndCompute(processedDF_perCond,ix))
	res <- unlist(plyr::compact(res))
	res <- as.data.frame(matrix(res, ncol = 5,  byrow = TRUE), stringsAsFactors = FALSE)
	colnames(res) <- c("GE","Pval","dPSI","psi1","psi2")
	res$FDR <- p.adjust(res$Pval,method="BY")
	res[[groupingFactor]] <- gName
	return(res)
}

computeRes_twoSample <- function(countDF,cond1,cond2){	## similar as above
        processedDF_perCond <- getProcessedDF_twoSample(countDF,cond1,cond2)
        exonList <- unique(processedDF_perCond$GE)
        res <- mclapply(exonList, function(ix) checkAndCompute(processedDF_perCond,ix),	mc.cores = numThreads)
        res <- unlist(plyr::compact(res))
        res <- as.data.frame(matrix(res, ncol = 5,  byrow = TRUE), stringsAsFactors = FALSE)
        colnames(res) <- c("GE","Pval","dPSI","psi1","psi2")
        res$FDR <- p.adjust(res$Pval,method="BY")
	res[[groupingFactor]] <- paste(cond1,cond2,sep = "-")
        return(res)
}


## one vs. all configuration
if(comparisonType == "oneVsAll"){
	toInclude <- c(colnames(countData)[1:4],groupingFactor)
	print(paste0("Calculating one vs all differential exon usage for",groupingFactor))
	conditionList <- unique(countData[[groupingFactor]])
	pval <- exonL <- dpsi <- psi1 <- psi2 <- c()
	print(paste("Starting",length(conditionList),"comparisons"))
	fullResList <- mclapply(conditionList, function(condName) computeResPerCond(countData,condName),
				mc.cores = numThreads)
	fullResDF <- do.call('rbind',fullResList)
	write.table(fullResDF,file.path(outputDir,paste0(comparisonType,"_per_",groupingFactor,"_diffExonUsage")),
		row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}

# Re-formatting input for testing conditions ------

#countData <- countData %>% separate(CT, into = levelNames, sep = "::",remove = FALSE)

## two sample comparison

if(comparisonType == "twoConditions"){
	for(i in 1:nrow(configFile)){
		groupingFactor <- as.vector(configFile$groupingFactor[i])
		cond1 <- as.vector(configFile$group1[i])
		cond2 <- as.vector(configFile$group2[i])
		pval <- exonL <- dpsi <- psi1 <- psi2 <- c()
		print(paste0("Calculating two sample differential exon usage for ",cond1," vs ",cond2))
		toInclude <- c(colnames(countData)[1:4],groupingFactor)
		resDF <- computeRes_twoSample(countData,cond1,cond2)
		write.table(resDF, file.path(outputDir, paste(comparisonType,groupingFactor,cond1,substr(cond2, 1, 4),"diffExonUsage",sep="_")),
			row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
	}
}
