#!/software/mambaforge/bin/Rscript

rm(list = ls())

### Check package availability ###

list.of.packages <- c("ape", "ggtree","tidyverse","gggenes","tidytree","dplyr","ggnewscale","ggstance")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

if("ggtree" %in% new.packages){
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("ggtreeExtra")
}

### Load packages ###

suppressMessages(suppressWarnings(library("ape")))
suppressMessages(suppressWarnings(library("ggtree")))
suppressMessages(suppressWarnings(library("tidyverse")))
suppressMessages(suppressWarnings(library("gggenes")))
suppressMessages(suppressWarnings(library("tidytree")))
suppressMessages(suppressWarnings(library("dplyr")))
suppressMessages(suppressWarnings(library("ggnewscale")))
suppressMessages(suppressWarnings(library("ggstance")))

### Parse Arguments ###

args <- commandArgs(trailingOnly = TRUE)

# REP

args <- c("--tree","/Users/muratbuyukyoruk/Documents/DRT2/tree_output.nexus","--genemap","/Users/muratbuyukyoruk/Documents/DRT2/DRT2_refined_motif_positions_p2_redundant_added_promoter_cleaned_Clade7_added.txt","--clade","/Users/muratbuyukyoruk/Documents/DRT2/DRT2_clade_info_v3.txt","--flip","n","--anchor","3'Loop","--output","outfile.txt","--interactive","y")

filename <- args[which(args=="--tree")+1]
genemap_in <-  args[which(args=="--genemap")+1]
anc_point <-  args[which(args=="--anchor")+1]
clade_in <-  args[which(args=="--clade")+1]
file.out <- args[which(args=="--output")+1]
FLIP <- args[which(args=="--flip")+1]
interactive_plot <- args[which(args=="--interactive")+1]

color_val <- c("#ffffff",
               "#4c87c6",
               "#db2726",
               "#008b46",
               "#b91c8d",
               "#ff00ff",
               "#008b46",
               "#b91c8d",
               "#f15f22",
               "#2ac0d5",
               "#4c87c6",
               "#4c87c6",
               "#db2726",
               "#4c87c6",
               "#4c87c6",
               "#f0e904",
               "#008b46",
               "#008b46",
               "#4c87c6",
               "#2ac0d5",
               "#4c87c6",
               "#ffffff")

res <- c("Start",
         "A",
         "R",
         "N",
         "D",
         "C",
         "Q",
         "E",
         "G",
         "H",
         "I",
         "L",
         "K",
         "M",
         "F",
         "P",
         "S",
         "T",
         "W",
         "Y",
         "V",
         "-")

### Import Tree ###

treein <- read.tree(filename)
x <- treein

p1 <- ggtree(x) 

if(length(clade_in)!=0){
  quants <- list()
  
  clade_df <- read.delim(clade_in, header = TRUE)
  clade_append <- NULL
  
  for(i in 1:ncol(clade_df)){
    clade_info <- paste(as.list(clade_df[,i][clade_df[,i]!=""]),collapse = "|")
    assign(paste("clade",i,sep = "_"),grep(clade_info,treein$tip.label))
    if(i==1){
      clade_append <- clade_info
    }else{
      clade_append <- paste(clade_append,clade_info,sep = "|")
    }
    if(i==ncol(clade_df)){
      uncladed <- grep(clade_append,treein$tip.label,invert = TRUE)
    }
  }
  for(i in 1:ncol(clade_df)){
    if (length(uncladed)!=0){
      quants[["Rest"]] <- uncladed
    }
    quants[[paste0(colnames(clade_df)[i])]] <- get(paste("clade",i,sep = "_"))
  }
  
  p1 <- groupOTU(p1, quants, 'Clade') + aes(color=Clade)
}

p1 <- p1 + guides(color=guide_legend(ncol=1)) + new_scale_color()

### Import genemap dataframe PART I ###

genemap <- read.delim(genemap_in,header = TRUE)

### Decide to flip or not ###

if(FLIP=="y"){
  flip_check <- data.frame(unique(genemap$molecule))
  flip_check$flip <- NA
  
  colnames(flip_check) <- c("molecule","flip")
  
  for (i in 1:nrow(genemap)){
    if (genemap$gene[i]==anc_point){
      if (genemap$orientation[i]==0){
        flip_check[which(grepl(genemap$molecule[i], flip_check$molecule)),]$flip <- "yes"
      } else{
        flip_check[which(grepl(genemap$molecule[i], flip_check$molecule)),]$flip <- "no"
      }
    }
  }
  
  genemap$start_new <- NA
  genemap$end_new <- NA
  
  for (i in 1:nrow(genemap)){
    if (flip_check[which(grepl(genemap$molecule[i], flip_check$molecule)),]$flip=="yes"){
      genemap$start_new[i] <- ( 10000 - genemap$end[i]) + 1
      genemap$end_new[i] <- ( 10000 - genemap$start[i]) + 1
      if(genemap$orientation[i] == 1){
        genemap$orientation[i] <-  0
      }else{
        genemap$orientation[i] <-  1
      }
    } else{
      genemap$start_new[i] <- genemap$start[i]
      genemap$end_new[i] <- genemap$end[i]
    }
  }
  
  genemap$start <- genemap$start_new
  genemap$end <- genemap$end_new
}

### Map domains to ggtree plot ###

genemap <- subset(genemap,gene!="GAP")

tip_count <- length(x$tip.label)

psub <-  p1 + 
  geom_facet(mapping = aes(xmin = start, xmax = end, fill = gene),
             data = subset(genemap,gene!="Length" | gene == anc_point), geom = geom_motif, panel = 'ncRNA',
             on = "5'Loop"  , arrowhead_height = unit(0, "mm"), arrowhead_width = unit(0, "mm"), arrow_body_height = unit(200/tip_count,"mm"),size = 0) + theme_tree2(legend.key.size = unit(0.2, "cm"))

### Import genemap dataframe PART II ###

genemap <- read.delim("/Users/muratbuyukyoruk/Documents/DRT2/DRT2_refined_motif_positions_p1_redundant_added_Clade7_added.txt",header = TRUE)

### Decide to flip or not ###

if(FLIP=="y"){
  flip_check <- data.frame(unique(genemap$molecule))
  flip_check$flip <- NA
  
  colnames(flip_check) <- c("molecule","flip")
  
  for (i in 1:nrow(genemap)){
    if (genemap$gene[i]==anc_point){
      if (genemap$orientation[i]==0){
        flip_check[which(grepl(genemap$molecule[i], flip_check$molecule)),]$flip <- "yes"
      } else{
        flip_check[which(grepl(genemap$molecule[i], flip_check$molecule)),]$flip <- "no"
      }
    }
  }
  
  genemap$start_new <- NA
  genemap$end_new <- NA
  
  for (i in 1:nrow(genemap)){
    if (flip_check[which(grepl(genemap$molecule[i], flip_check$molecule)),]$flip=="yes"){
      genemap$start_new[i] <- ( 10000 - genemap$end[i]) + 1
      genemap$end_new[i] <- ( 10000 - genemap$start[i]) + 1
      if(genemap$orientation[i] == 1){
        genemap$orientation[i] <-  0
      }else{
        genemap$orientation[i] <-  1
      }
    } else{
      genemap$start_new[i] <- genemap$start[i]
      genemap$end_new[i] <- genemap$end[i]
    }
  }
  
  genemap$start <- genemap$start_new
  genemap$end <- genemap$end_new
}

### Map domains to ggtree plot ###

genemap <- subset(genemap,gene!="GAP")

tip_count <- length(x$tip.label)


psub2 <-  psub + 
  geom_facet(mapping = aes(xmin = start, xmax = end, fill = gene),
             data = subset(genemap,gene!="Length" | gene == anc_point), geom = geom_motif, panel = 'ncRNA2',
             on = "3'Loop"  , arrowhead_height = unit(0, "mm"), arrowhead_width = unit(0, "mm"), arrow_body_height = unit(200/tip_count,"mm"),size = 0) + theme_tree2(legend.key.size = unit(0.2, "cm"))

### Import genemap dataframe PART III ###

genemap <- read.delim("/Users/muratbuyukyoruk/Documents/DRT2/DRT2_redundant_pro_genemap_YXDD.txt",header = TRUE)

p <- psub2 + geom_facet(mapping = aes(xmin = start, xmax = end, fill = gene),
                        data = subset(subset(genemap,start>0), gene == 'Start' | gene == 'N-domain' | gene == 'IS' | gene == 'C-domain'), geom = geom_motif, panel = 'Alignment',
                        on = "Start" , arrowhead_height = unit(0, "mm"), arrowhead_width = unit(0, "mm"), arrow_body_height = unit(200/tip_count,"mm"),size = 0) + new_scale_fill() +
  geom_facet(mapping = aes(xmin = start, xmax = end, fill = gene),
             data = subset(subset(genemap,start>0), gene != 'N-domain' & gene != 'IS' & gene != 'C-domain'), geom = geom_motif, panel = 'Alignment',
             on = "Start" , arrowhead_height = unit(0, "mm"), arrowhead_width = unit(0, "mm"), arrow_body_height = unit(200/tip_count,"mm"),size = 0) + scale_fill_manual(values=color_val,breaks = res) 
#facet_widths(p, c(1,3,5))

### Import genemap dataframe PART IV ###

promoter <-  read.delim("/Users/muratbuyukyoruk/Documents/DRT2/DRT2_refined_motif_positions_cleaned_for_cargo_promoter_dist_redundant_added.txt",header = TRUE)

colnames(promoter) <- c("Accession","Position","Distance")

p_plot <- p + geom_facet(panel = "Cargo Length (nt)", data = promoter, geom = geom_col, 
                         aes(x = Distance, fill = Position) , orientation = 'y', width = .6)

facet_widths(p_plot, c(2,4,4,3,2))
