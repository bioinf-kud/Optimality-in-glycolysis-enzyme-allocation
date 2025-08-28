library(ggplot2)
library(ggrepel)
dir<-"-path to folder-/figures_and_source_data"
head<-"/figureS5/data/calculated/"
tail1<-".glycolysis.tsv"
tail2<-".glycolysis.a_permuted.tsv"
tail3<-".glycolysis.e_permuted.tsv"
tail4<-".glycolysis.both_permuted.tsv"
data5 <- read.table(paste0(dir,head,"TCGA-COAD",tail1), sep = "\t")
data6 <- read.table(paste0(dir,head,"TCGA-COAD",tail2), sep = "\t")
data7 <- read.table(paste0(dir,head,"TCGA-COAD",tail3), sep = "\t")
data8 <- read.table(paste0(dir,head,"TCGA-COAD",tail4), sep = "\t")
data9 <- read.table(paste0(dir,head,"TCGA-BRCA",tail1),sep='\t')
data10 <- read.table(paste0(dir,head,"TCGA-BRCA",tail2),sep='\t')
data11 <- read.table(paste0(dir,head,"TCGA-BRCA",tail3),sep='\t')
data12 <- read.table(paste0(dir,head,"TCGA-BRCA",tail4),sep='\t')
data13 <- read.table(paste0(dir,head,"TCGA-HNSC",tail1),sep='\t')
data14 <- read.table(paste0(dir,head,"TCGA-HNSC",tail2),sep='\t')
data15 <- read.table(paste0(dir,head,"TCGA-HNSC",tail3),sep='\t')
data16 <- read.table(paste0(dir,head,"TCGA-HNSC",tail4),sep='\t')
data17 <- read.table(paste0(dir,head,"TCGA-KIRC",tail1),sep='\t')
data18 <- read.table(paste0(dir,head,"TCGA-KIRC",tail2),sep='\t')
data19 <- read.table(paste0(dir,head,"TCGA-KIRC",tail3),sep='\t')
data20 <- read.table(paste0(dir,head,"TCGA-KIRC",tail4),sep='\t')
data21 <- read.table(paste0(dir,head,"TCGA-LUAD",tail1),sep='\t')
data22 <- read.table(paste0(dir,head,"TCGA-LUAD",tail2),sep='\t')
data23 <- read.table(paste0(dir,head,"TCGA-LUAD",tail3),sep='\t')
data24 <- read.table(paste0(dir,head,"TCGA-LUAD",tail4),sep='\t')
data25 <- read.table(paste0(dir,head,"TCGA-LUSC",tail1),sep='\t')
data26 <- read.table(paste0(dir,head,"TCGA-LUSC",tail2),sep='\t')
data27 <- read.table(paste0(dir,head,"TCGA-LUSC",tail3),sep='\t')
data28 <- read.table(paste0(dir,head,"TCGA-LUSC",tail4),sep='\t')
data29 <- read.table(paste0(dir,head,"TCGA-PRAD",tail1),sep='\t')
data30 <- read.table(paste0(dir,head,"TCGA-PRAD",tail2),sep='\t')
data31 <- read.table(paste0(dir,head,"TCGA-PRAD",tail3),sep='\t')
data32 <- read.table(paste0(dir,head,"TCGA-PRAD",tail4),sep='\t')
data33 <- read.table(paste0(dir,head,"TCGA-THCA",tail1),sep='\t')
data34 <- read.table(paste0(dir,head,"TCGA-THCA",tail2),sep='\t')
data35 <- read.table(paste0(dir,head,"TCGA-THCA",tail3),sep='\t')
data36 <- read.table(paste0(dir,head,"TCGA-THCA",tail4),sep='\t')

data5 <- t(data5)
colnames(data5)<-make.unique(data5[1, ])
data5 <- data5[-1, ]
data5 <- as.data.frame(data5)
data6 <- t(data6)
colnames(data6)<-make.unique(data6[1, ])
data6 <- data6[-1, ]
data6 <- as.data.frame(data6)
data7 <- t(data7)
colnames(data7)<-make.unique(data7[1, ])
data7 <- data7[-1, ]
data7 <- as.data.frame(data7)
data8 <- t(data8)
colnames(data8)<-make.unique(data8[1, ])
data8 <- data8[-1, ]
data8 <- as.data.frame(data8)




data9 <- data9[, -10]
data9 <- t(data9)
colnames(data9) <- data9[1, ]
data9 <- data9[-1, ]
data9 <- as.data.frame(data9)
data10 <- data10[, -10]
data10 <- t(data10)
colnames(data10) <- data10[1, ]
data10 <- data10[-1, ]
data10 <- as.data.frame(data10)
data11 <- data11[, -10]
data11 <- t(data11)
colnames(data11) <- data11[1, ]
data11 <- data11[-1, ]
data11 <- as.data.frame(data11)
data12 <- data12[, -10]
data12 <- t(data12)
colnames(data12) <- data12[1, ]
data12 <- data12[-1, ]
data12 <- as.data.frame(data12)
data13 <- data13[, -10]
data13 <- t(data13)
colnames(data13) <- data13[1, ]
data13 <- data13[-1, ]
data13 <- as.data.frame(data13)
data14 <- data14[, -10]
data14 <- t(data14)
colnames(data14) <- data14[1, ]
data14 <- data14[-1, ]
data14 <- as.data.frame(data14)
data15 <- data15[, -10]
data15 <- t(data15)
colnames(data15) <- data15[1, ]
data15 <- data15[-1, ]
data15 <- as.data.frame(data15)
data16 <- data16[, -10]
data16 <- t(data16)
colnames(data16) <- data16[1, ]
data16 <- data16[-1, ]
data16 <- as.data.frame(data16)
data17 <- data17[, -10]
data17 <- t(data17)
colnames(data17) <- data17[1, ]
data17 <- data17[-1, ]
data17 <- as.data.frame(data17)
data18 <- data18[, -10]
data18 <- t(data18)
colnames(data18) <- data18[1, ]
data18 <- data18[-1, ]
data18 <- as.data.frame(data18)
data19 <- data19[, -10]
data19 <- t(data19)
colnames(data19) <- data19[1, ]
data19 <- data19[-1, ]
data19 <- as.data.frame(data19)
data20 <- data20[, -10]
data20 <- t(data20)
colnames(data20) <- data20[1, ]
data20 <- data20[-1, ]
data20 <- as.data.frame(data20)
data21 <- data21[, -10]
data21 <- t(data21)
colnames(data21) <- data21[1, ]
data21 <- data21[-1, ]
data21 <- as.data.frame(data21)
data22 <- data22[, -10]
data22 <- t(data22)
colnames(data22) <- data22[1, ]
data22 <- data22[-1, ]
data22 <- as.data.frame(data22)
data23 <- data23[, -10]
data23 <- t(data23)
colnames(data23) <- data23[1, ]
data23 <- data23[-1, ]
data23 <- as.data.frame(data23)
data24 <- data24[, -10]
data24 <- t(data24)
colnames(data24) <- data24[1, ]
data24 <- data24[-1, ]
data24 <- as.data.frame(data24)
data25 <- data25[, -10]
data25 <- t(data25)
colnames(data25) <- data25[1, ]
data25 <- data25[-1, ]
data25 <- as.data.frame(data25)
data26 <- data26[, -10]
data26 <- t(data26)
colnames(data26) <- data26[1, ]
data26 <- data26[-1, ]
data26 <- as.data.frame(data26)
data27 <- data27[, -10]
data27 <- t(data27)
colnames(data27) <- data27[1, ]
data27 <- data27[-1, ]
data27 <- as.data.frame(data27)
data28 <- data28[, -10]
data28 <- t(data28)
colnames(data28) <- data28[1, ]
data28 <- data28[-1, ]
data28 <- as.data.frame(data28)
data29 <- data29[, -10]
data29 <- t(data29)
colnames(data29) <- data29[1, ]
data29 <- data29[-1, ]
data29 <- as.data.frame(data29)
data30 <- data30[, -10]
data30 <- t(data30)
colnames(data30) <- data30[1, ]
data30 <- data30[-1, ]
data30 <- as.data.frame(data30)
data31 <- data31[, -10]
data31 <- t(data31)
colnames(data31) <- data31[1, ]
data31 <- data31[-1, ]
data31 <- as.data.frame(data31)
data32 <- data32[, -10]
data32 <- t(data32)
colnames(data32) <- data32[1, ]
data32 <- data32[-1, ]
data32 <- as.data.frame(data32)
data33 <- data33[, -10]
data33 <- t(data33)
colnames(data33) <- data33[1, ]
data33 <- data33[-1, ]
data33 <- as.data.frame(data33)
data34 <- data34[, -10]
data34 <- t(data34)
colnames(data34) <- data34[1, ]
data34 <- data34[-1, ]
data34 <- as.data.frame(data34)
data35 <- data35[, -10]
data35 <- t(data35)
colnames(data35) <- data35[1, ]
data35 <- data35[-1, ]
data35 <- as.data.frame(data35)
data36 <- data36[, -10]
data36 <- t(data36)
colnames(data36) <- data36[1, ]
data36 <- data36[-1, ]
data36 <- as.data.frame(data36)


vector <- c(8.841552744, -0.6857715, 7.833065244, -9.641140506, -6.172457509, 5.478969765, -1.774938001, 1.290864001)


data<-data5
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column<-as.numeric(current_column)
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}

correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"original"
correlations$dat <- "TCGA-COAD"
dataf<-correlations
data<-data6
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column<-as.numeric(current_column)
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}

correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"a_permutated"
correlations$dat <- "TCGA-COAD"
dataf<-rbind(dataf,correlations)

data<-data7
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column<-as.numeric(current_column)
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}

correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"e_permutated"
correlations$dat <- "TCGA-COAD"
dataf<-rbind(dataf,correlations)

data<-data8
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column<-as.numeric(current_column)
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}

correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"both_permutated"
correlations$dat <- "TCGA-COAD"
dataf<-rbind(dataf,correlations)

data<-data9
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column <- as.numeric(current_column)
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}

correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"original"
correlations$dat <- "TCGA-BRCA"
dataf<-rbind(dataf,correlations)
data<-data10
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column <- as.numeric(current_column)
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}

correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"a_permutated"
correlations$dat <- "TCGA-BRCA"
dataf<-rbind(dataf,correlations)

data<-data11
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column <- as.numeric(current_column)
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}

correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"e_permutated"
correlations$dat <- "TCGA-BRCA"
dataf<-rbind(dataf,correlations)

data<-data12
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column <- as.numeric(current_column)
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}

correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"both_permutated"
correlations$dat <- "TCGA-BRCA"
dataf<-rbind(dataf,correlations)

data<-data13
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column <- as.numeric(current_column)
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}

correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"original"
correlations$dat <- "TCGA-HNSC"
dataf<-rbind(dataf,correlations)
data<-data14
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column <- as.numeric(current_column)
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}

correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"a_permutated"
correlations$dat <- "TCGA-HNSC"
dataf<-rbind(dataf,correlations)

data<-data15
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column <- as.numeric(current_column)
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}

correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"e_permutated"
correlations$dat <- "TCGA-HNSC"
dataf<-rbind(dataf,correlations)

data<-data16
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column <- as.numeric(current_column)
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}

correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"both_permutated"
correlations$dat <- "TCGA-HNSC"
dataf<-rbind(dataf,correlations)

data<-data17
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column <- as.numeric(current_column)
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}

correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"original"
correlations$dat <- "TCGA-KIRC"
dataf<-rbind(dataf,correlations)
data<-data18
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column <- as.numeric(current_column)
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}

correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"a_permutated"
correlations$dat <- "TCGA-KIRC"
dataf<-rbind(dataf,correlations)

data<-data19
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column <- as.numeric(current_column)
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}

correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"e_permutated"
correlations$dat <- "TCGA-KIRC"
dataf<-rbind(dataf,correlations)

data<-data20
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column <- as.numeric(current_column)
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}

correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"both_permutated"
correlations$dat <- "TCGA-KIRC"
dataf<-rbind(dataf,correlations)

data<-data21
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column <- as.numeric(current_column)
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}

correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"original"
correlations$dat <- "TCGA-LUAD"
dataf<-rbind(dataf,correlations)
data<-data22
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column <- as.numeric(current_column)
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}

correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"a_permutated"
correlations$dat <- "TCGA-LUAD"
dataf<-rbind(dataf,correlations)

data<-data23
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column <- as.numeric(current_column)
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}

correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"e_permutated"
correlations$dat <- "TCGA-LUAD"
dataf<-rbind(dataf,correlations)

data<-data24
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column <- as.numeric(current_column)
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}

correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"both_permutated"
correlations$dat <- "TCGA-LUAD"
dataf<-rbind(dataf,correlations)

data<-data25
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column <- as.numeric(current_column)
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}

correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"original"
correlations$dat <- "TCGA-LUSC"
dataf<-rbind(dataf,correlations)
data<-data26
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column <- as.numeric(current_column)
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}

correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"a_permutated"
correlations$dat <- "TCGA-LUSC"
dataf<-rbind(dataf,correlations)

data<-data27
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column <- as.numeric(current_column)
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}

correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"e_permutated"
correlations$dat <- "TCGA-LUSC"
dataf<-rbind(dataf,correlations)

data<-data28
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column <- as.numeric(current_column)
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}

correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"both_permutated"
correlations$dat <- "TCGA-LUSC"
dataf<-rbind(dataf,correlations)

data<-data29
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column <- as.numeric(current_column)
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}

correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"original"
correlations$dat <- "TCGA-PRAD"
dataf<-rbind(dataf,correlations)
data<-data30
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column <- as.numeric(current_column)
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}

correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"a_permutated"
correlations$dat <- "TCGA-PRAD"
dataf<-rbind(dataf,correlations)

data<-data31
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column <- as.numeric(current_column)
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}

correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"e_permutated"
correlations$dat <- "TCGA-PRAD"
dataf<-rbind(dataf,correlations)

data<-data32
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column <- as.numeric(current_column)
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}

correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"both_permutated"
correlations$dat <- "TCGA-PRAD"
dataf<-rbind(dataf,correlations)

data<-data33
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column <- as.numeric(current_column)
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}

correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"original"
correlations$dat <- "TCGA-THCA"
dataf<-rbind(dataf,correlations)
data<-data34
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column <- as.numeric(current_column)
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}

correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"a_permutated"
correlations$dat <- "TCGA-THCA"
dataf<-rbind(dataf,correlations)

data<-data35
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column <- as.numeric(current_column)
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}

correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"e_permutated"
correlations$dat <- "TCGA-THCA"
dataf<-rbind(dataf,correlations)

data<-data36
correlations <- numeric(ncol(data))
datanum<- numeric(ncol(data))
for (i in 1:ncol(data)) {
  current_column <- data[[i]]
  current_column <- as.numeric(current_column)
  current_column<- log(current_column)
  valid_indices <- which(current_column != -Inf)
  filtered_column <- current_column[valid_indices]
  filtered_vector <- vector[valid_indices]
  valid_indices <- which(current_column != Inf)
  filtered_column <- filtered_column[valid_indices]
  filtered_vector <- filtered_vector[valid_indices]
  if (length(valid_indices) > 1) {
    datanum[i]<- length(valid_indices)
    correlations[i] <- cor(filtered_column, filtered_vector, method = "spearman")
  } else {
    correlations[i] <- NA
  }
}

correlations <- as.data.frame(correlations)
names(correlations)<-c("value")
row.names(correlations) <- names(data)
correlations$datanum <-datanum
correlations$per<-"both_permutated"
correlations$dat <- "TCGA-THCA"
dataf<-rbind(dataf,correlations)
dataf$type=row.names(dataf)
dataf$type <- gsub("\\..*", "", dataf$type)

dataf$per <- factor(dataf$per, levels = c("original", "a_permutated", "e_permutated", "both_permutated"))
dataf_sub<-dataf[dataf$type=="Tumor",]
dataf_sub<-dataf_sub[dataf_sub$per=="original",]
group_counts <- dataf_sub %>%
  group_by(dat) %>%
  summarise(count = n()) %>%
  mutate(new_label = paste0(dat, " (n=", count, ")"))
label_map <- setNames(group_counts$new_label, group_counts$dat)
dataf_sub<-dataf[dataf$type=="Tumor",]
ggplot(dataf_sub) + 
  geom_boxplot(aes(y = value, x = per, fill = per, colour = per), alpha = 0.5) +  
  scale_fill_manual(
    values = c("#824098","#70AD47","#ED7D31","#4472C4"),  
    labels = c("original data", "a permuted", "[E] permuted", "both permuted") 
  ) +
  scale_colour_manual(
    values = c("#824098","#70AD47","#ED7D31","#4472C4")  
  ) +
  labs(
    title = expression(
      atop(
        "Correlation between " * a[i] * "[" * E[i] * "]"^2 / (a[i+1] * "[" * E[i+1] * "]"^2) * " and " * K[i],
        "for glycolysis enzymes in TCGA tumor samples"
      )
    ),
    y = expression("Spearman's " ~ rho),  
    col = "method"
  ) +
  facet_wrap(~ dat, labeller = labeller(dat = label_map)) +  
  theme(
    axis.text.x = element_blank(),  
    axis.title.x = element_blank(),  
    axis.text.y = element_text(size = 20),  
    axis.title = element_text(size = 20),  
    plot.title = element_text(size = 22, hjust = 0.5, vjust = 1), 
    plot.title.position = "plot",  
    panel.grid = element_blank(),  
    panel.background = element_rect(fill = "white", color = NA),  
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    legend.position = c(0.85, 0.17), 
    legend.background = element_rect(fill = "white", color = NA), 
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 16) 
  ) +
  guides(fill = guide_legend(title = NULL), colour = "none")  


dataf_sub<-dataf[dataf$type=="Normal",]
dataf_sub<-dataf_sub[dataf_sub$per=="original",]
group_counts <- dataf_sub %>%
  group_by(dat) %>%
  summarise(count = n()) %>%
  mutate(new_label = paste0(dat, " (n=", count, ")"))
label_map <- setNames(group_counts$new_label, group_counts$dat)
dataf_sub<-dataf[dataf$type=="Normal",]
ggplot(dataf_sub) + 
  geom_boxplot(aes(y = value, x = per, fill = per, colour = per), alpha = 0.5) +  
  scale_fill_manual(
    values = c("#824098","#70AD47","#ED7D31","#4472C4"),  
    labels = c("original data", "a permuted", "[E] permuted", "both permuted") 
  ) +
  scale_colour_manual(
    values = c("#824098","#70AD47","#ED7D31","#4472C4")  
  ) +
  labs(
    title = expression(
      atop(
        "Correlation between " * a[i] * "[" * E[i] * "]"^2 / (a[i+1] * "[" * E[i+1] * "]"^2) * " and " * K[i],
        "for glycolysis enzymes in TCGA normal samples"
      )
    ),
    y = expression("Spearman's " ~ rho),  
    col = "method"
  ) +
  facet_wrap(~ dat, labeller = labeller(dat = label_map)) +  
  theme(
    axis.text.x = element_blank(),  
    axis.title.x = element_blank(),  
    axis.text.y = element_text(size = 20),  
    axis.title = element_text(size = 20),  
    plot.title = element_text(size = 22, hjust = 0.5, vjust = 1), 
    plot.title.position = "plot",  
    panel.grid = element_blank(),  
    panel.background = element_rect(fill = "white", color = NA),  
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    legend.position = c(0.85, 0.17), 
    legend.background = element_rect(fill = "white", color = NA), 
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 16)   
  ) +
  guides(fill = guide_legend(title = NULL), colour = "none")  
