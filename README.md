# serum-iron-GWAS-meta-analysis
R script for serum iron concentration GWAS meta-analysis

#######iron meta-analysis#######
##########Before editing######## 
Iron meta-analysis data cleaning R script


###set working directory###
setwd("C:/Users/DELL/Documents/iron meta analysis")
###uploading datasets from my documents in my desktop###
Denis_Iron = fread("C:/Users/DELL/Documents/iron meta analysis/Denis.txt")
head(Denis_Iron)
dim(Denis_Iron)
Thareja_Iron = fread("C:/Users/DELL/Documents/iron meta analysis/Thareja.txt")
head(Thareja_Iron)
dim(Thareja_Iron)
Benyamin_Iron = fread("C:/Users/DELL/Documents/iron meta analysis/ Benyamin.txt")
head(Benyamin_Iron)
dim(Benyamin_Iron)

####################

###calculating SE of Thareja dataset###

Thareja_Iron$p_value <- as.numeric(Thareja_Iron$p_value)
Thareja_Iron$BETA <- as.numeric(Thareja_Iron$BETA)

calculate_se <- function(beta, p_value) {
  if (is.na(beta) | is.na(p_value) | beta == 0 | p_value == 0) {
    return(NA)
  }
  Z <- qnorm(1 - p_value / 2)
  SE <- abs(beta / Z)
  return(SE)
}
#Apply the function to the data
Thareja_Iron <- Thareja_Iron %>%
  mutate(SE = mapply(calculate_se, BETA, p_value))
###adding N and SNP columns to Thareja dataset 
Thareja_Iron <- Thareja_Iron %>%
  mutate(SNP = paste0(CHR, ':', BP, '_', A1, '/', A2))

Thareja_Iron$N <- 6010

###checking columns names of Thareja dataset###

colnames(Thareja_Iron)
colnames(Thareja_Iron) <- c("CHR", "BP", "PVALUE", "Allele1", "Allele2", "FREQ", "EFFECT", "STDERR", "SNP",  "N")

###downloading the edited dataset into my desktop### 

write.table (x=Thareja_Iron, file='Thareja.txt', sep=' ', row.names=FALSE, quote=FALSE)


###adding SNP column into Dennis dataset###

Denis_Iron <- Denis_Iron %>%
  mutate(SNP = paste0(CHR, ':', BP, '_', Allele1, '/', Allele2))
#editing columns names
colnames(Denis_Iron)
colnames(Denis_Iron) <- c("CHR", "SNPrsid", "BP", "Allele1", "Allele2", "N", "FREQ", "EFFECT", "STDERR", "PVALUE", "SNP")
#download the edtied dataset
write.table (x=Denis_Iron, file='Denis.txt', sep=' ', row.names=FALSE, quote=FALSE)

###adding N and SNP columns to benyamin’s dataset###

benyaminiron$N <- 23986

benyaminiron <- benyaminiron %>% mutate(MarkerName = paste0(CHR, ':', BP, '_', A1, '/', A2))
head(benyaminiron)
#edit columns names
colnames(benyaminiron) <- c("SNPrsid", "Chromosome", "Position", "Allele1", "Allele2", "FREQ", "Beta", "StdErr", "P-value", "N", "MarkerName")
#download the edited dataset
write.table (x=benyaminiron, file='Benyamin.txt', sep=' ', row.names=FALSE, quote=FALSE)

###############liftover of benyamin###########

setwd("C:/Users/DELL/Documents/iron meta analysis")

benyaminiron <- fread("C:/Users/DELL/Documents/iron meta analysis/Benyamin.txt")
bed_iron <- benyaminiron[, c("Chromosome", "Position")]
head(bed_iron)
bed_iron$end <- bed_iron$Position + 1
bed_iron$end <- as.integer(bed_iron$end)
bed_iron$Chromosome <- paste0("chr", bed_iron$Chromosome)

formatted_rows <- paste0(bed_iron$Chromosome, ":", bed_iron$Position, "-", bed_iron$end)

# Define the path to save the file, e.g., to the Desktop
file_path <- "C:/Users/DELL/Documents/iron meta analysis/bed_iron.txt"

# Write the formatted rows to the txt file
writeLines(formatted_rows, file_path)

######bed format
 

#############edditing Benyamin file after the liftover#####

setwd("C:/Users/DELL/Documents/iron meta analysis")

Benyamin_LO <- fread("C:/Users/DELL/Documents/iron meta analysis/outbenyamin_liftover.txt")

head(Benyamin_LO)

###the data was Downloaded from ucsc as a one column, i splited it into Three columns, as integers
Benyamin_LO <- Benyamin_LO %>%
  separate(v1, into = c("Chromosome", "Position_end"), sep = ":")%>%
  separate(Position_end, into = c("Position", "End"), sep = "-") %>%
  mutate(
    Chromosome = as.integer(Chromosome),
    Position = as.integer(Position),
    End = as.integer(End)
  )

Benyamin_LO$Chromosome <- gsub("^chr", "", Benyamin_LO$Chromosome)

Benyamin_LO$Chromosome <- as.integer(Benyamin_LO$Chromosome)

################
Benyamin_deleted_SNPs <- fread("C:/Users/DELL/Documents/iron meta analysis/Benyamin_deleted.txt") 
head(Benyamin_deleted_SNPs)

###the data was Downloaded from ucsc as a one column, i splited it into Three columns, as integers
Benyamin_deleted_SNPs <- Benyamin_deleted_SNPs %>%
  separate(v1, into = c("Chromosome", "Position_end"), sep = ":")%>%
  separate(Position_end, into = c("Position", "End"), sep = "-") %>%
  mutate(
    Chromosome = as.integer(Chromosome),
    Position = as.integer(Position),
    End = as.integer(End)
  )

Benyamin_deleted_SNPs$Chromosome <- gsub("^chr", "", Benyamin_deleted_SNPs$Chromosome)

Benyamin_deleted_SNPs$Chromosome <- as.integer(Benyamin_deleted_SNPs$Chromosome)

Benyamin_deleted_SNPs <- Benyamin_deleted_SNPs %>%
  mutate(chrBP = paste(Chromosome, Position, sep = ":"))

######filtering benyamin Iron###

head(benyaminiron)
benyaminiron <- benyaminiron %>%
  mutate(chrBP = paste(Chromosome, Position, sep = ":"))

# Extract the Position values from Benyamin_deleted_SNPs
positions_to_remove <- Benyamin_deleted_SNPs$chrBP

# Filter out rows from benyaminiron where Position is in positions_to_remove
benyaminiron_filtered <- benyaminiron %>%
  filter(!(chrBP %in% positions_to_remove))

# Print the filtered dataset to check the result
print(benyaminiron_filtered)

###edit columns names 
colnames(benyaminiron_filtered) <- c("SNPrsid", "Chromosome", "Position36", "Allele1", "Allele2", "FREQ", "Beta", "StdErr", "P-value", "N", "MarkerName", "chrBP")
head(benyaminiron_filtered)
benyaminiron_filtered$Position <- Benyamin_LO$Position
head(benyaminiron_filtered)
benyaminiron_filtered <- benyaminiron_filtered %>%
  mutate(MarkerName = paste0(Chromosome, ':', Position, '_', Allele1, '/', Allele2))

write.table (x=benyaminiron_filtered, file='Benyamin37.txt', sep=' ', row.names=FALSE, quote=FALSE)
###Meta analysis has conducted using METAL package###
###configeration file###
# analysis by uncommenting the following line:
SCHEME  SAMPLESIZE
GENOMICCONTROL OFF
SEPARATOR  WHITESPACE
AVERAGEFREQ ON 

#LOAD THE TWO INPUT FILES
PROCESS Denis.txt
PROCESS Thareja.txt
PROCESS Benyamin37.txt


# === DESCRIBE AND PROCESS THE INPUT FILES ===
MARKERLABEL   MARKER
ALLELELABELS  Allele1 Allele2
PVALUELABEL   PVALUE
EFFECTLABEL EFFECT 
STDERR STDERR
A

#OUTPUT FILE 
OUTFILE iron_meta .tbl

# Execute meta-analysis

ANALYZE 

QUIT


###meta data then re uploaded on R for further cleaning###
setwd("C:/Users/DELL/Documents/iron meta analysis")
ironmeta <- fread("C:/Users/DELL/Downloads/Meta_Ikhlas_Iron1.tbl")
####creat CHR column the dataset### 

ironmeta <- ironmeta %>%
  mutate(Chromosome = sapply(str_split(MarkerName, ":"), `[`, 1))

####Add BP column to the dataset
ironmeta <- ironmeta %>%
  mutate(Position = sapply(str_split(MarkerName, "[:_]"), `[`, 2))
#view columns names
head(ironmeta)
##change columns names
colnames(ironmeta)
colnames(ironmeta) <- c("MarkerName", "Allele1", "allele2", "Freq1", "FreqSE", "Beta", "StdErr", "P-value", "Direction", "N", "Chromosome", "Position")

write.table (x=ironmeta, file='ironmeta.txt', sep=' ', row.names=FALSE, quote=FALSE)

###calculating N for the metadata from direction of effect of two datasets###
ironmeta$N <- NA
n_denis <- 15274
n_thareja <- 6010

summary(as.factor(ironmeta$Direction))
variants_in_all <- which(ironmeta$Direction %in% c("--", "-+", "+-", "++", "0-", "0+", "0"))
variants_denis_only <- which(ironmeta$Direction %in% c("-?", "+?", "0?"))
variants_thareja_only <- which(ironmeta$Direction %in% c("?+","?-","?0"))

ironmeta$N <- NA
ironmeta$N[variants_in_all] <- n_denis+n_thareja
ironmeta$N[variants_denis_only] <- n_denis
ironmeta$N[variants_thareja_only] <- n_thareja

summary(as.factor(ironmeta$N))
summary(as.factor(ironmeta$N))
6010   15274   21284 
3718261 1803163 4437446 

head(ironmeta)
colnames(ironmeta) <- c("MarkerName", "Allele1", "Allele2", "Freq1", "FreqSE", "Beta", "StdErr", "P-value", "Direction", "N", "Chromosome", "Position")

write.table (x=ironmeta, file='ironmeta.txt', sep=' ', row.names=FALSE, quote=FALSE)

###calculating N for three datasets in meta 
Iron_3_meta <- fread("C:/Users/DELL/Downloads/Meta_Ikhlas_Iron21.tbl")

n_denis <- 15274
n_thareja <- 6010
n_benyamin <- 23986

variants_in_all <- which(Iron_3_meta$Direction %in% c("---", "--+", "--0", "-+-", "-++", "-+0", "+--", 
                                               "+-+", "+-0", "++-", "+++", "++0"))
variants_Denis_only <- which(Iron_3_meta$Direction %in% c("-??", "+??", "0??"))
variants_Thareja_only <- which(Iron_3_meta$Direction %in% c("?+?","?-?","?0?"))
variants_benyamin_only <- which(Iron_3_meta$Direction %in% c("??+","??-","??0"))
variants_Denis_Thareja <- which(Iron_3_meta$Direction %in% c("-+?", "+-?", "++?", "--?"))
variants_Denis_benyamin <- which(Iron_3_meta$Direction %in% c("-?-", "-?+", "+?-", "+?+"))
variants_Thareja_benyamin <- which(Iron_3_meta$Direction %in% c("?-+", "?+-", "?++", "?--", "?+0"))

Iron_3_meta$N <- NA
Iron_3_meta$N[variants_in_all] <- n_denis+n_thareja+n_benyamin 
Iron_3_meta$N[variants_Denis_only] <- n_denis 
Iron_3_meta$N[variants_Thareja_only] <- n_thareja
Iron_3_meta$N[variants_benyamin_only] <- n_benyamin
Iron_3_meta$N[variants_Denis_Thareja] <- n_denis+n_thareja
Iron_3_meta$N[variants_Denis_benyamin] <- n_denis+n_benyamin 
Iron_3_meta$N[variants_Thareja_benyamin] <- n_thareja+n_benyamin

summary(as.factor(Iron_3_meta$N))
head (Iron_3_meta)

####creat CHR column the dataset 

Iron_3_meta <- Iron_3_meta %>%
  mutate(Chromosome = sapply(str_split(MarkerName, ":"), `[`, 1))

####Add BP column to the dataset

Iron_3_meta <- Iron_3_meta %>%
  mutate(Position = sapply(str_split(MarkerName, "[:_]"), `[`, 2))


###editing position and chromosome as integars for FUMA uploading file###

setwd("C:/Users/DELL/Documents/iron meta analysis")

ironmeta <- fread("C:/Users/DELL/Documents/iron meta analysis/ironmeta.txt")
head(ironmeta)
ironmeta$Position <- as.integer(ironmeta$Position)
ironmeta$Chromosome <- as.integer(ironmeta$Chromosome)

write.table (x=ironmeta, file='ironmeta.txt', sep=' ', row.names=FALSE, quote=FALSE)
head(ironmeta)

#Perform joins
common_data <- Tharejairon %>%
  +   inner_join(Denisiron, by = "Position")%>%
  +   inner_join(Benyamin_LO, by = "Position")


####rewrite the columns names for metal meta analysis file 
setwd("C:/Users/DELL/Documents/METAL")

colnames(benyaminiron) <- c("SNPrsid", "CHR", "BP", "Allele1", "Allele2", "FREQ", "EFFECT", "STDERR", "PVALUE", "N", "MARKER")
       
write.table (x=benyaminiron, file='Benyamin.txt', sep=' ', row.names=FALSE, quote=FALSE)

##########benferoni correction####

setwd("C:/Users/DELL/Documents/iron meta analysis")

ironmetta_chrom <- fread("C:/Users/DELL/Documents/iron meta analysis/IRONCHROM.txt")
head(ironmetta_chrom)
colnames(ironmetta_chrom)[10] <- "pvalue"

ironmetta_chrom$bonferroni_pvalue <- p.adjust(ironmetta_chrom$pvalue, method = "bonferroni")

print(ironmetta_chrom)

write.table (x=ironmetta_chrom, file='iRONCHROMBONFERONI.txt', sep=' ', row.names=FALSE, quote=FALSE)


######editing mapped genes regional plot 

# Load the image
image <- image_read("C:/Users/DELL/Downloads/annotPlot_FUMA_jobs511242.png")

# Change the background color to white
image <- image_background(image, "white")

# Save the modified image
image_write(image, "C:/Users/DELL/Downloads/annotPlot_FUMA_jobs511242.png")

####edit chromosome and position as integers 
Iron_3_meta$Position <- as.integer(Iron_3_meta$Position)
Iron_3_meta$Chromosome <- as.integer(Iron_3_meta$Chromosome)


colnames(Iron_3_meta)
colnames(Iron_3_meta) <- c("MarkerName", "Allele1", "Allele2", "FREQ", "FreqSE", "Weight", "Zscore", "P-value", "Direction", "Chromosome", "Position")

summary(as.factor(Iron_3_meta$Weight))
###Manhattan-Plot on R
install.packages("qqman")
library(qqman)
manhattan(gwas_data, chr="CHR", bp="BP", snp="SNP", p="P", 
          col = c("blue4", "orange3"), 
          genomewideline = -log10(5e-8), 
          suggestiveline = -log10(1e-5), 
          main = "GWAS Manhattan Plot")
png("manhattan_plot.png", width=800, height=600)
manhattan(gwas_data, chr="CHR", bp="BP", snp="SNP", p="P")
dev.off()

###Q-Q plot on R
install.packages("qqman")
library(qqman)
qq(gwas_data$P, main = "Q-Q Plot of GWAS P-values")
observed <- -log10(sort(gwas_data$P))
expected <- -log10(ppoints(length(gwas_data$P)))
plot(expected, observed, 
     xlab = "Expected -log10(p)", 
     ylab = "Observed -log10(p)", 
     main = "Q-Q Plot of GWAS P-values",
     pch = 19, col = "blue")
abline(0, 1, col = "red")
observed <- -log10(sort(gwas_data$P))
expected <- -log10(ppoints(length(gwas_data$P)))

plot(expected, observed, 
     xlab = "Expected -log10(p)", 
     ylab = "Observed -log10(p)", 
     main = "Q-Q Plot of GWAS P-values",
     pch = 19, col = "blue")
abline(0, 1, col = "red")
png("qq_plot.png", width=800, height=600)
qq(gwas_data$P, main = "Q-Q Plot of GWAS P-values")
dev.off()


#######iron meta-analysis#######
##########Before editing######## 
setwd("C:/Users/DELL/Documents/iron meta analysis")
Denis_Iron = fread("C:/Users/DELL/Documents/iron meta analysis/Denis.txt")
head(Denis_Iron)
summary(Denis_Iron$N)
Denis_Iron = fread("C:/Users/DELL/Documents/iron meta analysis/serum_iron_Dennis_37.fastGWA")
dim(Denis_Iron)
head(Denis_Iron)
colnames(Denis_Iron)
colnames(Denis_Iron) <- c("CHR", "SNPrsid", "BP", "Allele1", "Allele2", "N", "FREQ", "EFFECT", "STDERR", "PVALUE", "MARKER")


Denis_Iron <- Denis_Iron %>%
  mutate(SNP = paste0(CHR, ':', BP, '_', Allele1, '/', Allele2))

write.table (x=Denis_Iron, file='Denis.txt', sep=' ', row.names=FALSE, quote=FALSE)

Thareja_Iron = fread("C:/Users/DELL/Documents/iron meta analysis/serum_iron_Thareja_37.tsv")
head(Thareja_Iron)
colnames(Thareja_Iron)
colnames(Thareja_Iron) <- c("CHR", "BP", "p_value", "A1", "A2", "AF1", "BETA")

Thareja_Iron <- Thareja_Iron %>%
  mutate(SNP = paste0(CHR, ':', BP, '_', A1, '/', A2))


#######calculate SE#####

Thareja_Iron$p_value <- as.numeric(Thareja_Iron$p_value)
Thareja_Iron$BETA <- as.numeric(Thareja_Iron$BETA)

calculate_se <- function(beta, p_value) {
  if (is.na(beta) | is.na(p_value) | beta == 0 | p_value == 0) {
    return(NA)
  }
  Z <- qnorm(1 - p_value / 2)
  SE <- abs(beta / Z)
  return(SE)
}

# Apply the function to the data

Thareja_Iron <- Thareja_Iron %>%
  mutate(SE = mapply(calculate_se, BETA, p_value))

Thareja_Iron$N <- 6010

head(Thareja_Iron)
colnames(Thareja_Iron)
colnames(Thareja_Iron) <- c("CHR", "BP", "PVALUE", "Allele1", "Allele2", "FREQ", "EFFECT", "MARKER", "STDERR", "N")

write.table (x=Thareja_Iron, file='Thareja.txt', sep=' ', row.names=FALSE, quote=FALSE)


################


setwd("C:/Users/DELL/Documents/iron meta analysis")


ironmeta <- fread("C:/Users/DELL/Downloads/Meta_Ikhlas_Iron1.tbl")

ironmeta$N <- NA


####creat CHR column the dataset 

ironmeta <- ironmeta %>%
  mutate(Chromosome = sapply(str_split(MarkerName, ":"), `[`, 1))

####Add BP column to the dataset

ironmeta <- ironmeta %>%
  mutate(Position = sapply(str_split(MarkerName, "[:_]"), `[`, 2))

head(ironmeta)

colnames(ironmeta)
colnames(ironmeta) <- c("MarkerName", "Allele1", "allele2", "Freq1", "FreqSE", "Beta", "StdErr", "P-value", "Direction", "N", "Chromosome", "Position")

write.table (x=ironmeta, file='ironmeta.txt', sep=' ', row.names=FALSE, quote=FALSE)



#######solve the direction of effect to find N (sample size)



n_denis <- 15274
n_thareja <- 6010

summary(as.factor(ironmeta$Direction))
variants_in_all <- which(ironmeta$Direction %in% c("--", "-+", "+-", "++", "0-", "0+", "0"))
variants_denis_only <- which(ironmeta$Direction %in% c("-?", "+?", "0?"))
variants_thareja_only <- which(ironmeta$Direction %in% c("?+","?-","?0"))

ironmeta$N <- NA
ironmeta$N[variants_in_all] <- n_denis+n_thareja
ironmeta$N[variants_denis_only] <- n_denis
ironmeta$N[variants_thareja_only] <- n_thareja


summary(as.factor(ironmeta$N))

summary(as.factor(ironmeta$N))
6010   15274   21284 
3718261 1803163 4437446 

head(ironmeta)
colnames(ironmeta) <- c("MarkerName", "Allele1", "Allele2", "Freq1", "FreqSE", "Beta", "StdErr", "P-value", "Direction", "N", "Chromosome", "Position")

write.table (x=ironmeta, file='ironmeta.txt', sep=' ', row.names=FALSE, quote=FALSE)

###########

setwd("C:/Users/DELL/Documents/iron meta analysis")

ironmeta <- fread("C:/Users/DELL/Documents/iron meta analysis/ironmeta.txt")
head(ironmeta)
ironmeta$Position <- as.integer(ironmeta$Position)
ironmeta$Chromosome <- as.integer(ironmeta$Chromosome)

write.table (x=ironmeta, file='ironmeta.txt', sep=' ', row.names=FALSE, quote=FALSE)
head(ironmeta)

###########
setwd("C:/Users/DELL/Documents/iron meta analysis")

Denisiron <- fread("C:/Users/DELL/Documents/iron meta analysis/Denis.txt")
Tharejairon <- fread("C:/Users/DELL/Documents/iron meta analysis/Thareja.txt")

head(Denisiron)
colnames(Denisiron)
colnames(Denisiron) <- c("Chromosome", "SNPrsid", "Position", "Allele1",  "Allele2", "N", "FREQ", "Beta", "StdErr", "P-value", "MarkerName")

write.table (x=Denisiron, file='Denis.txt', sep=' ', row.names=FALSE, quote=FALSE)

head(Tharejairon)

Tharejairon$Position <- as.integer(Tharejairon$Position)
Tharejairon$Chromosome <- as.integer(Tharejairon$Chromosome)
colnames(Tharejairon) <- c("Chromosome", "Position", "P-value", "Allele1",  "Allele2", "FREQ", "Beta", "MarkerName", "StdErr", "N")

write.table (x=Tharejairon, file='Thareja.txt', sep=' ', row.names=FALSE, quote=FALSE)


#########3 datasets meta results#####

iron_3_meta <- fread("C:/Users/DELL/Downloads/Meta_Ikhlas_Iron21.tbl")
head(iron_3_meta)


iron_3_meta  <- iron_3_meta  %>%
       mutate(Chromosome = sapply(str_split(MarkerName, ":"), `[`, 1))
iron_3_meta  <- iron_3_meta %>%
       mutate(Position = sapply(str_split(MarkerName, "[:_]"), `[`, 2))


head(iron_3_meta)




colnames(iron_3_meta)

colnames(iron_3_meta) <- c("MarkerName", "Allele1", "Allele2", "Freq1", "FreqSE", "N", "Zscore", "P-value", "Direction", "Chromosome", "Position")

iron_3_meta$Beta <- iron_3_meta$Zscore * sqrt(1 / (2 * pi))  # Assuming SE is 1 for simplicity, adjust based on context


iron_3_meta$StdErr <- sqrt(1 / (2 * pi))

setwd("C:/Users/DELL/Documents/iron meta analysis")
write.table (x=iron_3_meta, file='iron_3_meta.txt', sep=' ', row.names=FALSE, quote=FALSE)


martairon$N <- 236612



#Perform joins
common_data <- Tharejairon %>%
  +   inner_join(Denisiron, by = "Position")%>%
  +   inner_join(Benyamin_LO, by = "Position")


common_data <- Tharejairon %>%
  inner_join(Denisiron, by = "MarkerName")
  
colnames(yangiron)
dim(common_data)
[1]  0 28

###########

benyaminiron <- fread("C:/Users/DELL/Documents/iron meta analysis/Benyamin.txt")
Denisiron <- fread("C:/Users/DELL/Documents/iron meta analysis/Denis.txt")
Tharejairon <- fread("C:/Users/DELL/Documents/iron meta analysis/Thareja.txt")

head(benyaminiron)
dim(benyaminiron)

benyaminiron$N <- 23986
colnames(benyaminiron)


benyaminiron <- benyaminiron %>% mutate(MarkerName = paste0(CHR, ':', BP, '_', A1, '/', A2))

head(benyaminiron)

colnames(benyaminiron) <- c("SNPrsid", "Chromosome", "Position", "Allele1", "Allele2", "FREQ", "Beta", "StdErr", "P-value", "N", "MarkerName")

write.table (x=benyaminiron, file='Benyamin.txt', sep=' ', row.names=FALSE, quote=FALSE)


####rewrite the columns names for metal meta analysis file 
setwd("C:/Users/DELL/Documents/METAL")

colnames(benyaminiron) <- c("SNPrsid", "CHR", "BP", "Allele1", "Allele2", "FREQ", "EFFECT", "STDERR", "PVALUE", "N", "MARKER")
       
write.table (x=benyaminiron, file='Benyamin.txt', sep=' ', row.names=FALSE, quote=FALSE)


#############
yangiron <- fread("C:/Users/DELL/Documents/iron meta analysis/serum_iron_Yang_37.tsv")
head(yangiron)

yangiron$N <- 1794
colnames(yangiron)

yangiron <- yangiron %>% mutate(MarkerName = paste0(chromosome, ':', base_pair_location, '_', effect_allele, '/', other_allele))


colnames(yangiron) <- c("Chromosome", "Position", "SNPrsid", "Allele1", "Allele2", "FREQ", "sample_size", "Beta", "StdErr", "P-value", "N", "MarkerName")

setwd("C:/Users/DELL/Documents/iron meta analysis")
write.table (x=yangiron, file='yang.txt', sep=' ', row.names=FALSE, quote=FALSE)







############
Iron_3_meta <- fread("C:/Users/DELL/Downloads/Meta_Ikhlas_Iron21.tbl")

n_denis <- 15274
n_thareja <- 6010
n_benyamin <- 23986



variants_in_all <- which(Iron_3_meta$Direction %in% c("---", "--+", "--0", "-+-", "-++", "-+0", "+--", 
                                               "+-+", "+-0", "++-", "+++", "++0"))
variants_Denis_only <- which(Iron_3_meta$Direction %in% c("-??", "+??", "0??"))
variants_Thareja_only <- which(Iron_3_meta$Direction %in% c("?+?","?-?","?0?"))
variants_benyamin_only <- which(Iron_3_meta$Direction %in% c("??+","??-","??0"))
variants_Denis_Thareja <- which(Iron_3_meta$Direction %in% c("-+?", "+-?", "++?", "--?"))
variants_Denis_benyamin <- which(Iron_3_meta$Direction %in% c("-?-", "-?+", "+?-", "+?+"))
variants_Thareja_benyamin <- which(Iron_3_meta$Direction %in% c("?-+", "?+-", "?++", "?--", "?+0"))


Iron_3_meta$N <- NA
Iron_3_meta$N[variants_in_all] <- n_denis+n_thareja+n_benyamin 
Iron_3_meta$N[variants_Denis_only] <- n_denis 
Iron_3_meta$N[variants_Thareja_only] <- n_thareja
Iron_3_meta$N[variants_benyamin_only] <- n_benyamin
Iron_3_meta$N[variants_Denis_Thareja] <- n_denis+n_thareja
Iron_3_meta$N[variants_Denis_benyamin] <- n_denis+n_benyamin 
Iron_3_meta$N[variants_Thareja_benyamin] <- n_thareja+n_benyamin

summary(as.factor(Iron_3_meta$N))
head (Iron_3_meta)


####creat CHR column the dataset 

Iron_3_meta <- Iron_3_meta %>%
  mutate(Chromosome = sapply(str_split(MarkerName, ":"), `[`, 1))

####Add BP column to the dataset

Iron_3_meta <- Iron_3_meta %>%
  mutate(Position = sapply(str_split(MarkerName, "[:_]"), `[`, 2))

Iron_3_meta$Position <- as.integer(Iron_3_meta$Position)
Iron_3_meta$Chromosome <- as.integer(Iron_3_meta$Chromosome)


colnames(Iron_3_meta)
colnames(Iron_3_meta) <- c("MarkerName", "Allele1", "Allele2", "FREQ", "FreqSE", "Weight", "Zscore", "P-value", "Direction", "Chromosome", "Position")

summary(as.factor(Iron_3_meta$Weight))

###############liftover of benyamin###########

setwd("C:/Users/DELL/Documents/iron meta analysis")


benyaminiron <- fread("C:/Users/DELL/Documents/iron meta analysis/Benyamin.txt")
bed_iron <- benyaminiron[, c("Chromosome", "Position")]
head(bed_iron)
bed_iron$end <- bed_iron$Position + 1
bed_iron$end <- as.integer(bed_iron$end)
bed_iron$Chromosome <- paste0("chr", bed_iron$Chromosome)


formatted_rows <- paste0(bed_iron$Chromosome, ":", bed_iron$Position, "-", bed_iron$end)

# Define the path to save the file, e.g., to the Desktop
file_path <- "C:/Users/DELL/Documents/iron meta analysis/bed_iron.txt"

# Write the formatted rows to the txt file
writeLines(formatted_rows, file_path)


#############edditing Benyamin file after the liftover#####

setwd("C:/Users/DELL/Documents/iron meta analysis")

Benyamin_LO <- fread("C:/Users/DELL/Documents/iron meta analysis/outbenyamin_liftover.txt")

head(Benyamin_LO)

###the data was Downloaded from ucsc as a one column, i splited it into Three columns, as integers
Benyamin_LO <- Benyamin_LO %>%
  separate(v1, into = c("Chromosome", "Position_end"), sep = ":")%>%
  separate(Position_end, into = c("Position", "End"), sep = "-") %>%
  mutate(
    Chromosome = as.integer(Chromosome),
    Position = as.integer(Position),
    End = as.integer(End)
  )

Benyamin_LO$Chromosome <- gsub("^chr", "", Benyamin_LO$Chromosome)

Benyamin_LO$Chromosome <- as.integer(Benyamin_LO$Chromosome)

################
Benyamin_deleted_SNPs <- fread("C:/Users/DELL/Documents/iron meta analysis/Benyamin_deleted.txt") 
head(Benyamin_deleted_SNPs)

###the data was Downloaded from ucsc as a one column, i splited it into Three columns, as integers
Benyamin_deleted_SNPs <- Benyamin_deleted_SNPs %>%
  separate(v1, into = c("Chromosome", "Position_end"), sep = ":")%>%
  separate(Position_end, into = c("Position", "End"), sep = "-") %>%
  mutate(
    Chromosome = as.integer(Chromosome),
    Position = as.integer(Position),
    End = as.integer(End)
  )

Benyamin_deleted_SNPs$Chromosome <- gsub("^chr", "", Benyamin_deleted_SNPs$Chromosome)

Benyamin_deleted_SNPs$Chromosome <- as.integer(Benyamin_deleted_SNPs$Chromosome)

Benyamin_deleted_SNPs <- Benyamin_deleted_SNPs %>%
  mutate(chrBP = paste(Chromosome, Position, sep = ":"))



######filtering benyamin Iron###

head(benyaminiron)

benyaminiron <- benyaminiron %>%
  mutate(chrBP = paste(Chromosome, Position, sep = ":"))

# Extract the Position values from Benyamin_deleted_SNPs
positions_to_remove <- Benyamin_deleted_SNPs$chrBP

# Filter out rows from benyaminiron where Position is in positions_to_remove
benyaminiron_filtered <- benyaminiron %>%
  filter(!(chrBP %in% positions_to_remove))

# Print the filtered dataset to check the result
print(benyaminiron_filtered)



###########


colnames(benyaminiron_filtered) <- c("SNPrsid", "Chromosome", "Position36", "Allele1", "Allele2", "FREQ", "Beta", "StdErr", "P-value", "N", "MarkerName", "chrBP")
head(benyaminiron_filtered)

benyaminiron_filtered$Position <- Benyamin_LO$Position

head(benyaminiron_filtered)

benyaminiron_filtered <- benyaminiron_filtered %>%
  mutate(MarkerName = paste0(Chromosome, ':', Position, '_', Allele1, '/', Allele2))

write.table (x=benyaminiron_filtered, file='Benyamin37.txt', sep=' ', row.names=FALSE, quote=FALSE)




merged_df <- benyaminiron_filtered %>%
  inner_join(Denisiron, by = "MarkerName")
  
##########benferoni correction####

setwd("C:/Users/DELL/Documents/iron meta analysis")

ironmetta_chrom <- fread("C:/Users/DELL/Documents/iron meta analysis/IRONCHROM.txt")
head(ironmetta_chrom)
colnames(ironmetta_chrom)[10] <- "pvalue"

ironmetta_chrom$bonferroni_pvalue <- p.adjust(ironmetta_chrom$pvalue, method = "bonferroni")

print(ironmetta_chrom)

write.table (x=ironmetta_chrom, file='iRONCHROMBONFERONI.txt', sep=' ', row.names=FALSE, quote=FALSE)
