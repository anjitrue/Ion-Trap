library(tidyr)
library(reshape2)
library(ggplot2)
library(ggExtra)
library(data.table)
library(dplyr)
library(plyr)


######### Read in Data #########
# Optimize max injection time
Turbo <- read.csv("E:/Projects/Proteomics/IontrapScanRange/Turbo_scans_perMaxInjectionTime.csv", header = TRUE, sep = ",")
Rapid <- read.csv("E:/Projects/Proteomics/IontrapScanRange/Rapid_scans_perMaxInjectionTime.csv", header = TRUE, sep = ",")
Normal <- read.csv("E:/Projects/Proteomics/IontrapScanRange/Normal_scans_perMaxInjectionTime.csv", header = TRUE, sep = ",")

Allspeeds_1000_optimalInjection <- read.csv("E:/Projects/Proteomics/IontrapScanRange/OptimalMaxInjection/OptimalInjection_allSpeeds.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

##### Slice #### 
#Turbo
NumPeptides_2000 <- read.csv("E:/Projects/Proteomics/IontrapScanRange/NumPeptides_slicefrom2000.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
NumPeptides_100 <- read.csv("E:/Projects/Proteomics/IontrapScanRange/NumPeptides_slicefrom100.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

#Rapid
rapid_numPeptides_100 <- read.csv("E:/Projects/Proteomics/IontrapScanRange/Slice_Rapid/Slice_from100end/FDR summary_100-1700.csv", header = T, sep = ",", stringsAsFactors = FALSE)

rapid_numPeptides_2000 <- read.csv("E:/Projects/Proteomics/IontrapScanRange/Slice_Rapid/Slice_from2000end/FDR summary_250-2000.csv", header = T, sep = ",", stringsAsFactors = FALSE)

#Normal
normal_numPeptides_100 <- read.csv("E:/Projects/Proteomics/IontrapScanRange/Slice_normal/From_100end/FDR summary_100_1450.csv", header = T, sep = ",", stringsAsFactors = F)

normal_numPeptides_2000 <- read.csv("E:/Projects/Proteomics/IontrapScanRange/Slice_normal/From_2000end/FDR summary_300_2000.csv", header = T, sep = ",", stringsAsFactors = F)

#Removed - TURBO
removed_turbo_100 <- read.csv("E:/Projects/Proteomics/IontrapScanRange/Slice_Turbo/Slice25_from 100end/DTA_files/outfile.csv", header = T, sep = ",", stringsAsFactors = F)

removed_turbo_2000 <- read.csv("E:/Projects/Proteomics/IontrapScanRange/Slice_Turbo/Slice25_from2000end/DTA_files/outfile.csv", header = T, sep = ",", stringsAsFactors = F)

#Removed - RAPID
removed_rapid_100 <- read.csv("E:/Projects/Proteomics/IontrapScanRange/Slice_Rapid/Slice_from100end/DTA_files/outfile.csv", header = T, sep = ",", stringsAsFactors = F)

removed_rapid_2000 <- read.csv("E:/Projects/Proteomics/IontrapScanRange/Slice_Rapid/Slice_from2000end/DTA_files/outfile.csv", header = T, sep = ",", stringsAsFactors = F)

#Removed - NORMAL
removed_normal_100 <- read.csv("E:/Projects/Proteomics/IontrapScanRange/Slice_normal/From_100end/DTA_files/outfile.csv", header = T, sep = ",", stringsAsFactors = F)

removed_normal_2000 <- read.csv("E:/Projects/Proteomics/IontrapScanRange/Slice_normal/From_2000end/DTA_files/outfile.csv", header = T, sep = ",", stringsAsFactors = F)

#Charges - Turbo
charges_turbo_100 <- read.csv("E:/Projects/Proteomics/IontrapScanRange/Slice_Turbo/Slice25_from 100end/peptides/ChargesOutfile.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

charges_turbo_2000 <- read.csv("E:/Projects/Proteomics/IontrapScanRange/Slice_Turbo/Slice25_from2000end/peptides/ChargesOutfile.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

##### Ranges #####
turbo_numPeptides <- read.csv("E:/Projects/Proteomics/IontrapScanRange/Ranges_turbo/FDR summary_turboRanges.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

rapid_numPeptides <- read.csv("E:/Projects/Proteomics/IontrapScanRange/Ranges_rapid/FDR summary_rapidRanges.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

old_new <- read.csv("E:/Projects/Proteomics/IontrapScanRange/old_vs_new_ranges_FDRoutput.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

###### Product Tolerances #####
PSMs_productTolerance <- read.csv("E:/Projects/Proteomics/IontrapScanRange/PMS_ProductIonTolerance_20181221.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

######### Format Data #########

###### Optimize max injection time  ######
#- TURBO
injectionTimes_turbo <- c("m/z",4:22) #injection range used
colnames(Turbo) <- injectionTimes_turbo #set column names to injectio, time
# [1] "m/z" "4"   "5"   "6"   "7"   "8"   "9"   "10"  "11"  "12"  "13"  "14"  "15"  "16"  "17"  "18"  "19"  "20"  "21" 
# [20] "22" 

turbo_long <- reshape(Turbo, direction = "long", #Reshape data to long format
                      varying = list(names(Turbo)[2:ncol(Turbo)]), 
                      v.names = "NumberOfScans", idvar = c("m/z"), 
                      timevar = "Max_InjectionTime", times = 4:22) 

turbo_long <- turbo_long[complete.cases(turbo_long),] #Remove rows with NA's

# - RAPID
injectionTimes_rapid <- c("m/z", 7:36)
colnames(Rapid) <- injectionTimes_rapid

rapid_long <- reshape(Rapid, direction = "long", #Reshape data to long format
                      varying = list(names(Rapid)[2:ncol(Rapid)]), 
                      v.names = "NumberOfScans", idvar = c("m/z"), 
                      timevar = "Max_InjectionTime", times = 7:36)

rapid_long <- rapid_long[complete.cases(rapid_long),]

# - NORMAL
injectionTimes_normal <- c("m/z", 9:66)
colnames(Normal) <- injectionTimes_normal

normal_long <- reshape(Normal, direction = "long", #Reshape data to long format
                       varying = list(names(Normal)[2:ncol(Normal)]), 
                       v.names = "NumberOfScans", idvar = c("m/z"), 
                       timevar = "Max_InjectionTime", times = 9:66)

normal_long <- normal_long[complete.cases(normal_long),]

#1000 m/z
injectionTimes_1000 <- c("scanSpeed", 6:42)
colnames(Allspeeds_1000_optimalInjection) <- injectionTimes_1000

long_1000 <- reshape(Allspeeds_1000_optimalInjection, direction = "long",
                     varying = list(names(Allspeeds_1000_optimalInjection)[2:ncol(Allspeeds_1000_optimalInjection)]),
                     v.names = "NumberOfScans", idvar = c("scanSpeed"),
                     timevar = "Max_InjectionTime", times = 6:42)

long_1000 <- long_1000[complete.cases(long_1000),]

# Create a smaller data frame with the maxInjection time based on max number of scans for the m/z 

long_format <- function(x){
  group <- data.table(x)
  
  group = group[,list(maxNumberOfScans=max(NumberOfScans)), by = group$`m/z`]
  
  colnames(group) <- c("Range", "Max_NumScans")
  
  max_injection <- data.table(x[which(x$NumberOfScans %in% group$Max_NumScans),])
  
  max_injection <- max_injection[,list(Max_InjectionTime = max(Max_InjectionTime)), by = max_injection$`m/z`]
  
  colnames(max_injection) <- c("Range", "Max_InjTime")
  
  x <- right_join(group, max_injection, by = "Range")
  
  return(x)
}

turbo_long <- long_format(turbo_long)
turbo_long$Speed <- rep("turbo", nrow(turbo_long))

rapid_long <- long_format(rapid_long)
rapid_long$Speed <-rep("rapid", nrow(rapid_long))

normal_long <- long_format(normal_long)
normal_long$Speed <-rep("normal", nrow(normal_long))

maxInjection_all <- rbind(turbo_long, rapid_long, normal_long)

maxInjection_all$Hz <- maxInjection_all$Max_NumScans/15

###### Slice data  ######
# Slice data - TURBO
#slice data from 2000 end
NumPeptides_2000$Peptides <- gsub(",", "",NumPeptides_2000$Peptides)
NumPeptides_2000$Peptides <- as.numeric(NumPeptides_2000$Peptides)

#slice data from 100 end
NumPeptides_100$Peptides <- gsub(",", "",NumPeptides_100$Peptides)
NumPeptides_100$Peptides <- as.numeric(NumPeptides_100$Peptides)

# Slice reformat from FDR output
# take FDR summary output for Rapid search (slice 25 m/z) and reformat to contain the m/z and PSMs/peptides (or any other relevant columns)

z <- vector()
q <- vector()
v <- vector()
y <- vector()
x <- data.frame(MaxRange = y, PSMs = q, Peptides = v)

reformat_fdr_names_100 <- function(r){
  
  y <- gsub(".*ITMS_", "", unlist(r[[1]]))
  y <- gsub(".*20sec_", "", y)
  y <- gsub(".*35sec_", "", y)
  y <- gsub("_2000.*", "", y)
  z <- as.numeric(y)
  q <- as.numeric(r[[14]])
  v <- as.numeric(r[[17]])
  
  
  x <- data.frame(z, q, v)
  colnames(x) <- c("MinRange", "PSMs" , "Peptides")
  # Note that you must first convert to character, then to numeric. 
  #Otherwise you may end up with the order according to the factors.
  x <- x[order(as.numeric(as.character(x$MinRange))),]
  
  return(x)
}

reformat_fdr_names_2000 <- function(r){
  
  y <- gsub(".*100_", "", r[,1])
  y <- gsub("_rapid.*", "", y)
  y <- gsub("_normal.*", "", y)
  y <- gsub(".csv.*", "", y)
  z <- as.numeric(y)
  q <- as.numeric(r[[14]])
  v <- as.numeric(r[[17]])
  
  
  x <- data.frame(z, q, v)
  colnames(x) <- c("MaxRange", "PSMs" , "Peptides")
  return(output)
}

# Slice data - RAPID
reformated_rapid_numPeptides_100 <- reformat_fdr_names_100(rapid_numPeptides_100)
reformated_rapid_numPeptides_2000 <- reformat_fdr_names_2000(rapid_numPeptides_2000)

# Slice data - NORMAL
reformated_normal_numPeptides_100 <- reformat_fdr_names_100(normal_numPeptides_100)
reformated_normal_numPeptides_2000 <- reformat_fdr_names_2000(normal_numPeptides_2000)

format_charge_100 <- function(r){
  y <- gsub("_2000.*", "", r[[1]])
  y <- gsub(".*ITMS_", "", y)
  z <- as.numeric(y)
  
  x<- data.frame(z,r[,2:5])
  
  colnames(x) <- c("Min","2","3","4","5")
  
  return(x)
}


# Charges - turbo
RFcharges_100_t <- format_charge_100(charges_turbo_100)
RFcharges_100_t <- RFcharges_100_t[order(as.numeric(as.character(RFcharges_100_t$Min))),]
charges_t_100_wide <- melt(RFcharges_100_t, 'Min')
colnames(charges_t_100_wide) <- c("Min", "Charge", "Count")

typeof(RFcharges_100_t$Min)

reformat_removedlines <- function(r){
  y <- gsub(".*ITMS_", "", r[,1])
  y <- gsub("_2000.*", "", y)
  
  z <- as.numeric(y)
  x <- data.frame(z,r[,2])
  colnames(x) <- c("Min", "NumberOfLinesRemoved")
  
  return(x)
}

reformat_removedlines_2000 <- function(r){
  y <- gsub(".*100_", "", r[,1])
  y <- gsub(".txt.*", "", y)
  
  z <- as.numeric(y)
  x <- data.frame(z,r[,2])
  colnames(x) <- c("Max", "NumberOfLinesRemoved")
  return(x)
}

#reformat remove lines for turbo
reformated_removed_turbo_100 <- reformat_removedlines(removed_turbo_100)
reformated_removed_turbo_100 <- reformated_removed_turbo_100[order(as.numeric(as.character(reformated_removed_turbo_100$NumberOfLinesRemoved))),]

reformated_removed_turbo_2000 <- reformat_removedlines_2000(removed_turbo_2000)
reformated_removed_turbo_2000 <- reformated_removed_turbo_2000[order(as.numeric(as.character(reformated_removed_turbo_2000$NumberOfLinesRemoved))),]

# reformat remove lines for rapid
reformated_removed_rapid_100 <- reformat_removedlines(removed_rapid_100)
reformated_removed_rapid_100 <- reformated_removed_rapid_100[order(as.numeric(as.character(reformated_removed_rapid_100$NumberOfLinesRemoved))),]

reformated_removed_rapid_2000 <- reformat_removedlines_2000(removed_rapid_2000)
reformated_removed_rapid_2000 <- reformated_removed_rapid_2000[order(as.numeric(as.character(reformated_removed_rapid_2000$NumberOfLinesRemoved))),]

# reforamt remove lines from normal
reformated_removed_normal_100 <- reformat_removedlines(removed_normal_100)
reformated_removed_normal_100 <- reformated_removed_normal_100[order(as.numeric(as.character(reformated_removed_normal_100$NumberOfLinesRemoved))),]

reformated_removed_normal_2000 <- reformat_removedlines_2000(removed_normal_2000)
reformated_removed_normal_2000 <- reformated_removed_normal_2000[order(as.numeric(as.character(reformated_removed_normal_2000$NumberOfLinesRemoved))),]

######### Product Ion Tolerance #########
PSMs_productTolerance <- data.frame(PSMs_productTolerance)
PSMs_productTolerance$ProductIonTolerance <- as.factor(PSMs_productTolerance$ProductIonTolerance)
PSMs_productTolerance$Mean <- rowMeans(PSMs_productTolerance[,3:5])
PSMs_productTolerance$sd <- apply(PSMs_productTolerance[,3:5],1, sd, na.rm = TRUE)

tolerance <- as.factor(PSMs_productTolerance$ProductIonTolerance)

#mean_PSMs <- data.frame(ProductIonTolerance = PSMs_productTolerance[,1], Means = rowMeans(PSMs_productTolerance[,3:5]))

turbo_tolerance <- PSMs_productTolerance[which(PSMs_productTolerance$Speed == "turbo"),]
rapid_tolerance <- PSMs_productTolerance[which(PSMs_productTolerance$Speed == "rapid"),]
normal_tolerance <- PSMs_productTolerance[which(PSMs_productTolerance$Speed == "normal"),]

###### Ranges  ######
ranges_reformated <- function (r){
  y <- gsub(".*sec_", "", r[,1])
  y <- gsub("_range.*", "", y)
  
  v <- as.numeric(gsub("_.*", "", y))
  q <- as.numeric(gsub(".*_", "", y))
  
  z <- q-v
  z[is.na(z)] = 600
  v[is.na(v)] = "auto"
  output <- data.frame(v,q,z, r[[4]], r[[14]], r[[17]])
  
  colnames(output) <- c("Min", "Max", "Range", "MSMS","PSMs", "Peptides")
  return(output)
}

#not working yet
# cast_ranges_PSms <- function(r){
#   r$PSMs_run <- PSMs_run
#   
#   r$Range <- as.factor(r$Range)
#   output <- dcast(r, Range ~ PSMs_run, value.var = "PSMs")
#   
#   r$Mean <- rowMeans(output[,2:4], na.rm = T)
#   r$sd <- apply(output[,2:4], 1 , sd, na.rm = TRUE)
#   
#   return(output)
# }

# Ranges - TURBO
ranges_turbo <- ranges_reformated(turbo_numPeptides)
ranges_turbo <- ranges_turbo[order(as.numeric(as.character(ranges_turbo$Range))),]

Peptide_run = c(rep(c("Peptide1","Peptide2", "Peptide3"),3), "Peptide1", "Peptide2", rep(c("Peptide1","Peptide2", "Peptide3"),2), "Peptide1")
PSMs_run = c(rep(c("PSMs1","PSMs2", "PSMs3"),3), "PSMs1", "PSMs2", rep(c("PSMs1","PSMs2", "PSMs3"),2), "PSMs1")
MSMS_run = c(rep(c("MS/MS1","MS/MS2", "MS/MS3"),3), "MS/MS1", "MS/MS2", rep(c("MS/MS1","MS/MS2", "MS/MS3"),2), "MS/MS1")

ranges_turbo$PSMs_run <- PSMs_run
ranges_turbo$Peptide_run <- Peptide_run
ranges_turbo$MSMS_run <- MSMS_run

ranges_turbo$Range <- as.factor(ranges_turbo$Range)
ranges_turbo_wide <- dcast(ranges_turbo, Range ~ Peptide_run, value.var= "Peptides")


ranges_turbo_wide$Mean <- rowMeans(ranges_turbo_wide[,2:4], na.rm = TRUE)
ranges_turbo_wide$sd <- apply(ranges_turbo_wide[,2:4],1, sd, na.rm = TRUE)


# Ranges - RAPID
ranges_rapid <- ranges_reformated(rapid_numPeptides)
ranges_rapid <- ranges_rapid[order(as.numeric(as.character(ranges_rapid$Range))),]

Peptide_run= c(rep(c("Peptide1","Peptide2", "Peptide3"),6), "Peptide1")
PSMs_run = c(rep(c("PSMs1","PSMs2", "PSMs3"),6), "PSMs1")
MSMS_run = c(rep(c("MS/MS1","MS/MS2", "MS/MS3"),6), "MS/MS1")

ranges_rapid$PSMs_run <- PSMs_run
ranges_rapid$Peptide_run <- Peptide_run
ranges_rapid$MSMS_run <- MSMS_run

ranges_rapid$Range <- as.factor(ranges_rapid$Range)
ranges_rapid_wide <- dcast(ranges_rapid, Range ~ MSMS_run, value.var= "MSMS")


ranges_rapid_wide$Mean <- rowMeans(ranges_rapid_wide[,2:4], na.rm = TRUE)
ranges_rapid_wide$sd <- apply(ranges_rapid_wide[,2:4],1, sd, na.rm = TRUE)

#ranges_rapid_wide_minus <- ranges_rapid_wide[-7,]


# old and new
old_new_df <- data.frame(old_new$Raw.File, old_new$Total.MS.MS.Spectra, old_new$PSMs, old_new$Peptides)
colnames(old_new_df) <- c("File", "MSMS", "PSMs", "Peptides")

Peptide_run = c(rep(c("Peptide1","Peptide2"),3))
PSMs_run <- c(rep(c("PSMs1","PSMs2"),3))
MSMS_run <- c(rep(c("MS/MS1","MS/MS2"),3))

old_new_wide <- dcast(old_new_df, File ~ MSMS_run, value.var= "MSMS")
old_new_wide$Mean <- rowMeans(old_new_wide[,2:3], na.rm = TRUE)
old_new_wide$sd <- apply(old_new_wide[,2:3], 1, sd, na.rm = TRUE)



######### Plotting #########
palette(c("steelblue2", "olivedrab3", "mediumorchid2","darkgoldenrod1", 
          "hotpink3", "red2", "sienna2","slategray4","mediumturquoise", 
          "deepskyblue", "orangered", "midnightblue"))

#turbo_long$Max_InjectionTime = as.factor(turbo_long$Max_injectionTime)
colnames(maxInjection_all) #[1] "m/z"          "Max_NumScans" "Max_InjTime" "Range" "Max_NumScans" "Max_InjTime" "Speed" 

######### Max Injection (ms) For Various Scan Speeds #########
speed.colors = as.numeric(factor(maxInjection_all$Speed))
p <- ggplot(maxInjection_all, aes(x = Range, y = Max_InjTime, 
                                  group = factor(Speed), color = Speed)) +
  scale_x_continuous(breaks = seq(200, 2000, 200)) +
  scale_y_continuous(breaks = seq(4,68, 8)) +
  #xlim(200, 2000) +
  geom_point(size = 5) +
  geom_line(linetype = "dashed")+
  theme_light() +
  scale_color_manual(values=c("steelblue2", "olivedrab3", "mediumorchid2"))


p + labs(title = "Max Injection at m/z Range in MS/MS", x = "m/z", y = "Maximum Injection Time (ms)") 

speed.colors = as.numeric(factor(long_1000$scanSpeed))
p <- ggplot(long_1000, aes(x = Max_InjectionTime, y = NumberOfScans, group = factor(scanSpeed))) + #, color = speed.colors)) +
  scale_x_continuous(breaks = seq(2, 44, 2)) +
  scale_y_continuous(breaks = seq(200,600, 100)) +
  geom_point(size = 5) +
  geom_line(linetype = "dashed")+
  theme_light()
#scale_color_manual(values=c("steelblue2", "olivedrab3", "mediumorchid2"))

p

######### Number of Scans For Various Scan Speeds at Maximum Injection #########

q <- ggplot(maxInjection_all, aes(x = Range, y = Hz, group = factor(Speed), color = Speed)) +
  scale_x_continuous(breaks = seq(200, 2000, 200)) +
  scale_y_continuous(breaks = seq(10,55, 5)) +
  geom_point(size = 4, color = speed.colors) +
  geom_line(linetype = "dashed") +
  theme_light() +
  scale_color_manual(values=c("steelblue2", "olivedrab3", "mediumorchid2"))

q + labs(title = "Scan Speed", x = "m/z", y = "Hz") 

######### SLICED - For different scan speeds Peptides #########
#Turbo - sliced from 2000
colnames(NumPeptides_2000) #[1] "MaxRange" "Peptides"
m <- ggplot(NumPeptides_2000, aes(x = MaxRange, y = Peptides)) +
  scale_x_continuous(breaks = seq(100, 2000, 200)) +
  scale_y_continuous(breaks = seq(0,40000, 5000)) +
  geom_point(size = 5, color = "mediumorchid2") +
  geom_point(data = reformated_rapid_numPeptides_2000, size = 5, color = "olivedrab3")+
  geom_point(data = reformated_normal_numPeptides_2000, size = 5, color = "steelblue2")+
  geom_line(linetype = "dashed")+
  geom_line(data = reformated_rapid_numPeptides_2000, linetype = "dashed")+
  geom_line(data = reformated_normal_numPeptides_2000, linetype = "dashed")+
  theme_light()

m + labs(title = "Number of Peptides from Spectra sliced from 2000", x = "m/z", y = "Number of Peptides") 

ggMarginal(m, type = "histogram")

#Turbo - sliced from 100 
colnames(NumPeptides_100) #[1] "MinRange" "Peptides"
# scatter plot with marginal histogram
m <- ggplot(NumPeptides_100, aes(x = MinRange, y = Peptides)) +
  scale_x_continuous(breaks = seq(100, 2000, 200)) +
  scale_y_continuous(breaks = seq(0,40000, 5000)) +
  geom_point(size = 5, color = "mediumorchid2") +
  geom_point(data = reformated_rapid_numPeptides_100, size = 5, color = "olivedrab3")+
  geom_point(data = reformated_normal_numPeptides_100, size = 5, color = "steelblue2" )+
  geom_line(linetype = "dashed")+
  geom_line(data = reformated_rapid_numPeptides_100, linetype = "dashed")+
  geom_line(data = reformated_normal_numPeptides_100, linetype = "dashed")+
  theme_light()

m + labs(title = "Number of Peptides from Spectra sliced from 100", x = "m/z", y = "Number of Peptides")

ggMarginal(m, type = "histogram", xlab = "m/z", ylab = "Number of Peptides")

# turbo histogram of peptides 
m <- ggplot(NumPeptides_100, aes(x=Peptides)) +
  geom_histogram(bins = 100)

m <- ggplot(reformated_removed_turbo_100, aes(x=NumberOfLinesRemoved)) +
  geom_histogram(bins = 100)

df_3 <- data.frame(charges_turbo_100$Min,charges_turbo_100$charge3)
colnames(df_3) <- c("Min", "charge3")

colnames(charges_turbo_100)

m <- ggplot(charges_turbo_100, aes(x = Min, y = charge2)) +
  scale_x_continuous(breaks = seq(100, 2000, 200)) +
  scale_y_continuous(breaks = seq(0,25000, 5000)) +
  geom_point(size = 5, color = "mediumorchid2") + 
  geom_point(data = df_3, size = 5, color = "olivedrab3")


geom_point(data = reformated_normal_numPeptides_100, size = 5, color = "steelblue2" )+
  geom_line(linetype = "dashed")+
  geom_line(data = reformated_rapid_numPeptides_100, linetype = "dashed")+
  geom_line(data = reformated_normal_numPeptides_100, linetype = "dashed")+
  theme_light()

m
#Rapid - sliced from 100
colnames (reformated_rapid_numPeptides_2000) #[1] "MaxRange" "PSMs"     "Peptides"

m <- ggplot(reformated_rapid_numPeptides_2000, aes(x = MaxRange, y = Peptides)) +
  scale_x_continuous(breaks = seq(100, 2000, 200)) +
  scale_y_continuous(breaks = seq(0,40000, 5000)) +
  geom_point(size = 4, color = "olivedrab3") +
  theme_light()

m + labs(title = "Rapid - Number of Peptides from Spectra Sliced starting at 2000 m/z", x = "m/z", y = "Number of Peptides")

m <- ggplot(reformated_removed_rapid_100, aes(Min,NumberOfLinesRemoved)) +
  scale_x_continuous(breaks = seq(100, 2000, 200)) +
  geom_point(size = 5, color = "olivedrab3") +
  theme_light()
m + labs(title = "Number of Product Ions Removed - Rapid ", x = "m/z", y = "Number of Product Ions")


#Rapid - sliced from 2000
colnames (reformated_rapid_numPeptides_100) #[1] "MinRange" "PSMs"     "Peptides"

m <- ggplot(reformated_rapid_numPeptides_100, aes(x = MinRange, y = Peptides)) +
  scale_x_continuous(breaks = seq(100, 2000, 200)) +
  scale_y_continuous(breaks = seq(0,40000, 5000)) +
  geom_point(size = 4, color = "olivedrab3") +
  theme_light()

m + labs(title = "Rapid - Number of Peptides from Spectra Sliced starting at 100 m/z", x = "m/z", y = "Number of Peptides")


#Normal - sliced from 100
colnames (reformated_normal_numPeptides_100) #[1] "MinRange" "PSMs"     "Peptides"

m <- ggplot(reformated_normal_numPeptides_100, aes(x = MinRange, y = Peptides)) +
  scale_x_continuous(breaks = seq(100, 2000, 200)) +
  scale_y_continuous(breaks = seq(0,40000, 5000)) +
  geom_point(size = 4, color = "steelblue2") +
  theme_light()

m + labs(title = "Normal - Number of Peptides from Spectra Sliced starting at 100 m/z", x = "m/z", y = "Number of Peptides")

#Normal - sliced from 2000

m <- ggplot(reformated_normal_numPeptides_2000, aes(x = MaxRange, y = Peptides)) +
  scale_x_continuous(breaks = seq(100, 2000, 200)) +
  scale_y_continuous(breaks = seq(0,22000, 5000)) +
  geom_point(size = 5, color = "steelblue2") +
  theme_light()

m + labs(title = "Normal - Number of Peptides from Spectra Sliced starting at 2000 m/z", x = "m/z", y = "Number of Peptides")

######### Charges ############
charge.colors = as.numeric(factor(charges_t_100_wide$Charge))
charges_t_100_wide$Charge <- as.factor(charges_t_100_wide$Charge)

m <- ggplot(charges_t_100_wide, aes(x = Min, y = Count, group = factor(Charge), color = Charge)) +
  scale_x_continuous(breaks = seq(100, 2000, 200)) +
  #xlim(200, 2000) +
  geom_point(size = 4) +
  geom_line(linetype = "dashed")+
  theme_light() 
#scale_color_manual(values=c("steelblue2", "olivedrab3", "mediumorchid2"))

m + labs(title = "Turbo - Charge distrubution", x = "m/z", y = "Count")


######### SLICED - Lines removed #########
colnames(reformated_removed_turbo_100)
m <- ggplot(reformated_removed_turbo_100, aes(x = Min, y = NumberOfLinesRemoved)) +
  scale_x_continuous(breaks = seq(100, 2000, 200)) +
  scale_y_continuous(breaks = seq(0,78000000, 5000000)) +
  geom_point(size = 4, color = "mediumorchid2") +
  theme_light()

m + labs(title = "Rapid - Number of Product Ions Removed", x = "m/z", y = "Number Product Ions Removed")


# turbo removed product ions - 100
m <- ggplot(reformated_removed_turbo_100, aes(Min,NumberOfLinesRemoved)) +
  scale_x_continuous(breaks = seq(100, 2000, 200)) +
  geom_line(linetype = "dashed", color = "mediumorchid2") +
  geom_line(data = reformated_removed_rapid_100, linetype = "dashed", color = "olivedrab3")+
  geom_line(data = reformated_removed_normal_100, linetype = "dashed" , color = "steelblue2")+
  theme_light()
m + labs(title = "Number of Product Ions Removed ", x = "m/z", y = "Number of Product Ions")

# turbo removed product ions - 2000
m <- ggplot(reformated_removed_turbo_2000, aes(Max,NumberOfLinesRemoved)) +
  scale_x_continuous(breaks = seq(100, 2000, 200)) +
  geom_line(linetype = "dashed", color = "mediumorchid2") +
  geom_line(data = reformated_removed_rapid_2000, linetype = "dashed" , color = "olivedrab3")+
  geom_line(data = reformated_removed_normal_2000, linetype= "dashed", color = "steelblue2") +
  theme_light()
m + labs(title = "Number of Product Ions Removed ", x = "m/z", y = "Number of Product Ions")


######### Product Ion Tolerance #########
#combined bar plot with all scan speeds
p <- ggplot(turbo_tolerance, aes(x = ProductIonTolerance, y = Mean, group = 1)) +
  geom_bar(stat = "identity", color = "#6C686E" , fill = "#89469B", position = position_dodge())+
  geom_bar(data = rapid_tolerance,stat = "identity", color = "#6D696F" ,  fill = "#7C9E3D", position = position_dodge())+
  geom_bar(data = normal_tolerance, stat = "identity", color = "#6D696F" , fill = "#6A9ABA", position = position_dodge())+
  scale_y_continuous(breaks = seq(0,57000, 10000)) +
  geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd) , width = .2, position = position_dodge(0.052))+
  geom_errorbar(data = rapid_tolerance, aes(ymin=Mean-sd, ymax=Mean+sd) , width = .2, position = position_dodge(0.052))+
  geom_errorbar(data = normal_tolerance, aes(ymin=Mean-sd, ymax=Mean+sd) , width = .2, position = position_dodge(0.052))+
  theme_light()
p +labs(title = "Product Ion Tolerance - By Ion Trap Speed", x = "Product Ion Tolerance in Dalton", y = "Number of PSMs") 

# Rapid tolerance bar graph
m <- ggplot(rapid_tolerance, aes(x = ProductIonTolerance, y = Mean, group = 1)) +
  geom_bar(stat = "identity", color = "#6D696F" ,  fill = "#7C9E3D", position = position_dodge())+
  scale_y_continuous(breaks = seq(0,50000, 10000)) +
  geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd) , width = .2, position = position_dodge(0.052))+
  theme_light()
m +labs(title = "Product Ion Tolerance - Rapid Speed", x = "Product Ion Tolerance in Dalton", y = "Number of PSMs")

# Normal tolerance bar graph
p <- ggplot(normal_tolerance, aes(x = ProductIonTolerance, y = Mean, group = 1)) +
  geom_bar(stat = "identity", color = "#6D696F" , fill = "#6A9ABA", position = position_dodge())+
  scale_y_continuous(breaks = seq(0,50000, 10000)) +
  geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd) , width = .2, position = position_dodge(0.052))+
  p +labs(title = "Product Ion Tolerance - Normal Speed", x = "Product Ion Tolerance in Dalton", y = "Number of PSMs")

# Normal tolerance scatter plot with dashed line
p <- ggplot(normal_tolerance, aes(x = ProductIonTolerance, y = Mean, group = 1)) +
  geom_line(linetype = "dashed")+
  geom_point(size = 3, color = "steelblue2")+
  geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd) , width = .2, position = position_dodge(0.052))+
  theme_light()
p

######### Ranges #########
#Turbo MSMS bar graph
p <- ggplot(ranges_turbo_wide[-7,], aes(x = Range, y = Mean, group = 1)) +
  geom_bar(stat = "identity", color = "#6D696F" , fill = "#89469B", position = position_dodge())+
  scale_y_continuous(breaks = seq(0,130000, 20000)) +
  scale_x_continuous(breaks = seq(500, 1300, 200))+
  geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd) , width = .2, position = position_dodge(0.052))+
  theme_light()

p +labs(title = "Total Number of MS/MS - TURBO", x = "Mass Range (m/z)", y = "Number of MS/MS")

p <- ggplot(ranges_turbo_wide[-7,], aes(x = Range, y = Mean, group = 1)) +
  scale_y_continuous(breaks = seq(40000,46000, 1000)) +
  scale_x_continuous(breaks = seq(500, 1300, 200))+
  geom_point(size = 5 , color = "mediumorchid2")+
  geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd) , width = .2, position = position_dodge(0.052))+
  theme_light()

p +labs(title = "Total Number of Peptides - TURBO", x = "Mass Range (m/z)", y = "Number of Uniqe Peptides")  


#Rapid MSMS bar graph
p <- ggplot(ranges_rapid_wide[-7,], aes(x = Range, y = Mean, group = 1)) +
  geom_bar(stat = "identity", color = "#6D696F" , fill = "#7C9E3D", position = position_dodge())+
  scale_y_continuous(breaks = seq(0,120000, 20000)) +
  geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd) , width = .2, position = position_dodge(0.052))+
  theme_light()

p +labs(title = "Total Number of MS/MS - rapid", x = "Mass Range (m/z)", y = "Number of MS/MS")

p <- ggplot(ranges_rapid_wide[-7,], aes(x = Range, y = Mean, group = 1)) +
  scale_y_continuous(breaks = seq(30000,42000, 1000)) +
  #scale_x_continuous(breaks = seq(500, 1300, 200))+
  geom_point(size = 5 , color = "olivedrab3")+
  geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd) , width = .2, position = position_dodge(0.052))+
  theme_light()

p +labs(title = "Total Number of Peptides - Rapid", x = "Mass Range (m/z)", y = "Number of Uniqe Peptides")

old_new_wide <- old_new_wide[c(3,2,1),]
p <- ggplot(old_new_wide, aes(x = File, y = Mean, group = 1)) +
  scale_y_continuous(breaks = seq(43000,46000, 500)) +
  #scale_x_continuous(breaks = seq(500, 1300, 200))+
  geom_point(size = 5 , color = "mediumorchid2")+
  geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd) , width = .2, position = position_dodge(0.052))+
  theme_light()

p <- ggplot(old_new_wide, aes(x = File, y = Mean, group = 1)) +
  scale_y_continuous(breaks = seq(110000,120000, 1000)) +
  #scale_x_continuous(breaks = seq(500, 1300, 200))+
  geom_point(size = 5 , color = "mediumorchid2")+
  #geom_bar(stat = "identity", color = "#6D696F" , fill = "#89469B", position = position_dodge())+
  geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd) , width = .2, position = position_dodge(0.052))+
  theme_light()

p +labs(title = "Total Number of MSMS - Turbo", x = "Mass Range (m/z)", y = "Total Number of MSMS")
