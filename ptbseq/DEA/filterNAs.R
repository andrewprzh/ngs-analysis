subs_m <- as.matrix(subs) ## Main matrix

DR <- dist(subs_m) ## Get the euclidean distance
DR <- as.data.frame(as.matrix(DR))


DR2 <- DR %>% mutate(numNA = ncol(DR) - rowMeans(is.na(.))*ncol(DR)) ## find number of NAs per row
tooManyNAs_R <- unname(which(DR2$numNA <= ncol(DR)-200 )) ## I set it to 200 for a 6000 row matrix but you will have to adjust the parameter
length(tooManyNAs_R) ## was like 1200 


subs_m <- as.matrix(subs)
subs_m <- subs_m[-c(tooManyNAs_R),] ## Remove rows with too many NAs
dim(subs_m) ## was like 4500