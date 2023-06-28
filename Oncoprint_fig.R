## Oncoprint map for Yvonne

### NOTE: make sure to load "tcga analysis 9-24-18.RData" before running

setwd("/Volumes/Big Backup/Pre-2022/TCGA LUAD matched 2021/Yvonne NCI resub 2022/")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
library("ComplexHeatmap")

###############################

library(dplyr)

# generate GLRX z-score data

GLRX <- read.csv("GLRX_507_z-scores.csv")
rownames(GLRX) = GLRX[, 1]
GLRX = GLRX[, c(-1:-2)]
MA.1 <- as.matrix(GLRX)
MA.2 <- t(na.omit(t(MA.1)))
test <- colnames(MA.2)

# grab variant data for 507 cases

PD <- read.csv("PATIENT_DATA_oncoprint 2.csv")
PD[PD == "homdel_rec"] <- "DEL"
PD[PD == "Deep Deletion"] <- "DEL"
PD[PD == "amp_rec"] <- "AMP"
PD[PD == "Amplification"] <- "AMP"
PD[PD == "Missense Mutation (putative driver)"] <- "MUT"
PD[PD == "Inframe Mutation (putative driver)"] <- "MUT"
PD[PD == "Truncating mutation (putative driver)"] <- "TRUNC"
PD[PD == "Truncating mutation (putative passenger)"] <- "TRUNC"
PD[PD == "Missense Mutation (putative passenger)"] <- "VUS"
PD[PD == "splice_rec"] <- "SPLICE"
PD[PD == "splice"]<- "SPLICE"

PD.1 <- PD[c(-25:-30),]

PD.1[is.na(PD.1)] = ""

compress <- PD.1 %>%
  group_by(track_name) %>%
  summarise_all(~ toString(na.omit(.)))

compress.1 <- compress[,-2]
compress.1 <- data.frame(compress.1)
rownames(compress.1) = compress.1[, 1]
compress.2 <- compress.1[,-1]
compress.3 <- as.matrix(compress.2)
compress.3[is.na(compress.3)] = ""

test <- colnames(MA.2)

compress.4 <- compress.3[,test]
compress.4 <- as.matrix(compress.4)

# generate and combine oncoprint and heatmap data

col = c("DEL" = "black", "MUT" = "red", "TRUNC" = "#008000", "AMP" = "blue", "VUS" = "yellow", "SPLICE" = "orange")

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(4, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  
  
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(4, "pt"), h*0.8, 
              gp = gpar(fill = col["AMP"], col = NA))
  },
  
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w-unit(4, "pt"), h*0.6, 
              gp = gpar(fill = col["MUT"], col = NA))
  },
  
  TRUNC = function(x, y, w, h) {
    grid.rect(x, y, w-unit(4, "pt"), h*0.4, 
              gp = gpar(fill = col["TRUNC"], col = NA))
  },
  
  VUS = function(x, y, w, h) {
    grid.rect(x, y, w-unit(4, "pt"), h*0.3, 
              gp = gpar(fill = col["VUS"], col = NA))
  },
  
  DEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(4, "pt"), h*0.2, 
              gp = gpar(fill = col["DEL"], col = NA))
  },
  
  SPLICE = function(x, y, w, h) {
    grid.rect(x, y, w-unit(4, "pt"), h*0.2, 
              gp = gpar(fill = col["SPLICE"], col = NA))
  }
  
)

column_title = "Lung Adenocarcinomas (TCGA)"
heatmap_legend_param = list(title = "Alternations", at = c("VUS", "AMP", "MUT", "TRUNC", "DEL", "SPLICE"), 
                            labels = c("VUS", "Amplification", "Mutation", "Truncation", "Deep Del", "Splice Variant"))


p1 <- oncoPrint(compress.4,
                alter_fun = alter_fun,
                col = col,
                show_column_names = FALSE,
                pct_gp = gpar(fontsize = 12, 
                              fontface = 2),
                column_title = column_title, 
                heatmap_legend_param = heatmap_legend_param) %v%
  Heatmap(MA.2, 
          split=split, 
          show_row_dend = FALSE,
          
          row_title = NULL,
          row_title_gp = gpar(fontsize = 12,
                              fontface = 2),
          row_title_rot = 0,
          row_names_gp = gpar(fontsize = 12,
                              fontface = 2),
          show_column_names = FALSE,
          heatmap_legend_param = list(title = "mRNA z-score"))

p1

####################







matt = read.table("alterations_across_samples.tsv", 
                 header = TRUE,stringsAsFactors=FALSE, sep = "\t")

matt[is.na(matt)] = ""
matt = matt[, -1]
rownames(matt) = matt[, 1]
matt = matt[, c(-1:-3)]
matt = matt[,c(-7:-30)]
matt = t(as.matrix(matt))
matt[1:3, 1:3]


test.1 <- data.frame(row.names(LUAD))
SW <- data.frame(gsub("\\.", "\\-", (test.1$row.names.LUAD.)))

matt.2 <- matt[, colnames(matt) %in% SW$gsub................test.1.row.names.LUAD...]

matt.3 <- data.frame(matt.2)

write.csv(matt.3, "test mut file.csv")

matt.4 <- read.csv("test mut file edit.csv", header = TRUE, row.names = 1)
matt.5 <- as.matrix(matt.4)
matt.5[is.na(matt.5)] = ""

### all cases

all_cases = read.table("cbio_adeno.csv", 
                  header = TRUE,stringsAsFactors=FALSE, sep = ",")

all_cases[is.na(all_cases)] = ""
all_cases = all_cases[, c(-1,-3:-10)]
rownames(all_cases) = all_cases[, 1]



all_cases = t(as.matrix(all_cases))
all_cases[1:3, 1:3]
####
# expression data for GLRX and other genes SEE cbioportal_2022.R

#####
col = c("DEL" = "black", "MUT" = "red", "TRUNC" = "#008000", "AMP" = "blue", "INDEL" = "yellow")

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
 
 
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h*0.8, 
              gp = gpar(fill = col["AMP"], col = NA))
  },
  
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h*0.6, 
              gp = gpar(fill = col["MUT"], col = NA))
  },
  
  TRUNC = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h*0.4, 
              gp = gpar(fill = col["TRUNC"], col = NA))
  },
  
  INDEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h*0.3, 
              gp = gpar(fill = col["INDEL"], col = NA))
  },
  
  DEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h*0.2, 
              gp = gpar(fill = col["DEL"], col = NA))
  }
  
)

column_title = "OncoPrint for TCGA Lung Adenocarcinoma, sub genes"
heatmap_legend_param = list(title = "Alternations", at = c("INDEL", "AMP", "MUT", "TRUNC", "DEL"), 
                            labels = c("INDEL", "Amplification", "Mutation", "Truncation", "Deep Del"))
oncoPrint(matt.5,
          alter_fun = alter_fun, col = col,
          show_column_names = TRUE,
          column_title = column_title, heatmap_legend_param = heatmap_legend_param)
###

### Now need to incorporate GLRX expression (or other genes); perhaps as a heatmap underneath???



### table manipulations from new cBioprotal data structures: 507 cases for Yvonne

library(dplyr)

PD <- read.csv("PATIENT_DATA_oncoprint 2.csv")
PD[PD == "homdel_rec"] <- "DEL"
PD[PD == "Deep Deletion"] <- "DEL"
PD[PD == "amp_rec"] <- "AMP"
PD[PD == "Amplification"] <- "AMP"
PD[PD == "Missense Mutation (putative driver)"] <- "MUT"
PD[PD == "Inframe Mutation (putative driver)"] <- "MUT"
PD[PD == "Truncating mutation (putative driver)"] <- "TRUNC"
PD[PD == "Truncating mutation (putative passenger)"] <- "TRUNC"
PD[PD == "Missense Mutation (putative passenger)"] <- "VUS"
PD[PD == "splice_rec"] <- "SPLICE"
PD[PD == "splice"]<- "SPLICE"

PD.1 <- PD[c(-25:-30),]

PD.1[is.na(PD.1)] = ""

compress <- PD.1 %>%
    group_by(track_name) %>%
    summarise_all(~ toString(na.omit(.)))

compress.1 <- compress[,-2]
compress.1 <- data.frame(compress.1)
rownames(compress.1) = compress.1[, 1]
compress.2 <- compress.1[,-1]
compress.3 <- as.matrix(compress.2)
compress.3[is.na(compress.3)] = ""

test <- colnames(MA.2)

compress.4 <- compress.3[,test]
compress.4 <- as.matrix(compress.4)

col = c("DEL" = "black", "MUT" = "red", "TRUNC" = "#008000", "AMP" = "blue", "VUS" = "yellow", "SPLICE" = "orange")

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(4, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  
  
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(4, "pt"), h*0.8, 
              gp = gpar(fill = col["AMP"], col = NA))
  },
  
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w-unit(4, "pt"), h*0.6, 
              gp = gpar(fill = col["MUT"], col = NA))
  },
  
  TRUNC = function(x, y, w, h) {
    grid.rect(x, y, w-unit(4, "pt"), h*0.4, 
              gp = gpar(fill = col["TRUNC"], col = NA))
  },
  
  VUS = function(x, y, w, h) {
    grid.rect(x, y, w-unit(4, "pt"), h*0.3, 
              gp = gpar(fill = col["VUS"], col = NA))
  },
  
  DEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(4, "pt"), h*0.2, 
              gp = gpar(fill = col["DEL"], col = NA))
  },
  
  SPLICE = function(x, y, w, h) {
    grid.rect(x, y, w-unit(4, "pt"), h*0.2, 
              gp = gpar(fill = col["SPLICE"], col = NA))
  }
  
)

column_title = "Lung Adenocarcinomas (TCGA)"
heatmap_legend_param = list(title = "Alternations", at = c("VUS", "AMP", "MUT", "TRUNC", "DEL", "SPLICE"), 
                            labels = c("VUS", "Amplification", "Mutation", "Truncation", "Deep Del", "Splice Variant"))
p1 <- oncoPrint(compress.4,
          alter_fun = alter_fun, 
          col = col,
          column_order = p2,
          show_column_names = TRUE,
          column_title = column_title, 
          heatmap_legend_param = heatmap_legend_param)

### add in GLRX expression with fixed order at heatmap of z-scores

GLRX <- read.csv("GLRX_507_z-scores.csv")
rownames(GLRX) = GLRX[, 1]
GLRX = GLRX[, c(-1:-2)]


MA.1 <- as.matrix(GLRX)
MA.2 <- t(na.omit(t(MA.1)))
kclus <- kmeans(MA.2,1)
kclus$cluster
pdf(file='MA.pdf', width=8.5, height=22)  

split <- paste0("Cluster\n", kclus$cluster)
default.hmap <- Heatmap(MA.2, 
                        split=split, 
                        cluster_columns = FALSE, 
                        cluster_rows = FALSE,
                        row_title = "mRNA",
                        #column_order = p1,
                        heatmap_legend_param = list(title = "")
                        )
draw(default.hmap, newpage=FALSE)
dev.off()

p1 <- oncoPrint(compress.4,
          alter_fun = alter_fun, 
          col = col,
          show_column_names = TRUE,
          column_title = column_title, 
          heatmap_legend_param = heatmap_legend_param) 

p2<-Heatmap(N.2, 
        split=split, 
        cluster_rows = FALSE,
        show_column_dend = FALSE,
        #column_order = p1,
        row_title = "mRNA",
        heatmap_legend_param = list(title = "")
)

p2
#combine figures using pushViewport

pushViewport(viewport(layout=grid.layout(nr=2, nc=1)))

pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))

draw(p1, newpage=FALSE)
upViewport()

pushViewport(viewport(layout.pos.row=2, layout.pos.col=1))
draw(p2, newpage=FALSE)
upViewport()

#### FINALLY

col = c("DEL" = "black", "MUT" = "red", "TRUNC" = "#008000", "AMP" = "blue", "VUS" = "yellow", "SPLICE" = "orange")

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(4, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  
  
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(4, "pt"), h*0.8, 
              gp = gpar(fill = col["AMP"], col = NA))
  },
  
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w-unit(4, "pt"), h*0.6, 
              gp = gpar(fill = col["MUT"], col = NA))
  },
  
  TRUNC = function(x, y, w, h) {
    grid.rect(x, y, w-unit(4, "pt"), h*0.4, 
              gp = gpar(fill = col["TRUNC"], col = NA))
  },
  
  VUS = function(x, y, w, h) {
    grid.rect(x, y, w-unit(4, "pt"), h*0.3, 
              gp = gpar(fill = col["VUS"], col = NA))
  },
  
  DEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(4, "pt"), h*0.2, 
              gp = gpar(fill = col["DEL"], col = NA))
  },
  
  SPLICE = function(x, y, w, h) {
    grid.rect(x, y, w-unit(4, "pt"), h*0.2, 
              gp = gpar(fill = col["SPLICE"], col = NA))
  }
  
)

column_title = "Lung Adenocarcinomas (TCGA)"
heatmap_legend_param = list(title = "Alternations", at = c("VUS", "AMP", "MUT", "TRUNC", "DEL", "SPLICE"), 
                            labels = c("VUS", "Amplification", "Mutation", "Truncation", "Deep Del", "Splice Variant"))


p1 <- oncoPrint(compress.4,
                alter_fun = alter_fun, 
                col = col,
                show_column_names = FALSE,
                pct_gp = gpar(fontsize = 12, 
                              fontface = 2),
                column_title = column_title, 
                heatmap_legend_param = heatmap_legend_param) %v%
Heatmap(N.2, 
            split=split, 
        show_row_dend = FALSE,
       
            row_title = NULL,
        row_title_gp = gpar(fontsize = 12,
                            fontface = 2),
        row_title_rot = 0,
        row_names_gp = gpar(fontsize = 12,
                            fontface = 2),
        show_column_names = FALSE,
            heatmap_legend_param = list(title = "mRNA z-score"))

p1

p2 <- oncoPrint(compress.4,
                   alter_fun = alter_fun, 
                   col = col,
                   show_column_names = FALSE,
                   column_title = column_title, 
                   heatmap_legend_param = heatmap_legend_param) %v%
  Heatmap(N.2, 
          split=split, 
          row_title = "mRNA",
          show_column_names = FALSE,
          heatmap_legend_param = list(title = "mRNA z-score")
  )
p2

