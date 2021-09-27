# Christopher Hall, Babraham Institute, hallc@babraham.ac.uk
# EBI-EMBL Introduction to Flow Cytomtry DAta Anlaysi Using R 2021
# https://www.babraham.ac.uk/science-services/flow-cytometry
# https://github.com/hally166/Cytometry-R-scripts 

# set working directory
setwd('C:/Users/hallc/Documents/FCSData/FlowRepositoryData')

# install the packages
install.packages('BioCManager')
BiocManager::install('flowCore')
BiocManager::install('ggcyto')
BiocManager::install('flowAI')

# load packages
library('flowCore')
library('ggcyto')
library('flowAI')

# List files to be analysed
files_list <- list.files('C:/Users/hallc/Data Transfers/FlowRepository_FR-FCM-Z3WR_files',full.names = TRUE, patter = 'fcs')
files_list

# flowFrame - a single file
ff <- read.FCS(files_list[1]) # takes the fisrt file
ff <- read.FCS('C:/Users/hallc/Data Transfers/FlowRepository_FR-FCM-Z3WR_files/180E1A_0.fcs', truncate_max_range = FALSE)
keyword(ff)$'$CYT'
ff
markernames(ff)
keyword(ff)
t<-keyword(ff)
ff@exprs[,1]
ff@parameters@data[1:2]
ff@exprs

# flowSet - a group of files from the same experiment 
fs<-read.flowSet(files_list, truncate_max_range = FALSE)
fs[[1]]
markernames(fs)
fs[[2]]@exprs
ff
??read.flowSet
fsApply(fs, each_col, max)
fsApply(fs[1:3],function(x)autoplot(x))
fsApply(fs[1:3],function(x)median(x@exprs[,1]))

# Load data, compensate, clean, transform 
# compensation
spillover(ff)
ff <- read.FCS('C:/Users/hallc/Documents/FCSData/FlowRepositoryData/FlowRepository_FR-FCM-Z3GH_files/20P09_CFSE_PREGNANT_Experiment_Group_PREG 17_CFSE NS.fcs',emptyValue = FALSE)
ff@description$`$CYT`
ff_comp <- compensate(ff, spillover(ff)$`$SPILLOVER`)

# flowAI file cleaning
??flow_auto_qc
fs_clean<-flow_auto_qc(fs, remove_from = "FR_FS", fcs_QC = FALSE)
fs_clean
fs_clean[[1]]
fs[[1]]

#transform 
??estimateLogicle
ff@parameters@data[1:2]
ff@description$`$SPILLOVER`
colnames(ff@description$`$SPILLOVER`)
trans<-estimateLogicle(ff, channels = colnames(ff[,7:42]), t = 4194304)
trans<-estimateLogicle(ff, channels = colnames(ff@description$`$SPILLOVER`))
ff_trans<-transform(ff_clean,trans)
ff_trans<-transform(ff,trans)
autoplot(ff_trans)
autoplot(ff)

fs[[1]]@parameters@data[1:2]
trans<-estimateLogicle(fs[[1]], channels = colnames(fs[[1]][,7:42]), t = 4194304)
fs_trans<-transform(fs_clean,trans)
autoplot(fs_trans[[1]])

fs_trans
ggcyto(fs_trans, aes(x = CD4, y=CD8)) + geom_hex(bins=256)

# FSC SSC Live Singlets - We can gate using flowCore, but do we want to?
# Some people prepare their data in other softwares for inmput into R
# We can do three easy gates - none debris, singlet, live - using the gatingset

library(flowWorkspace)
gs <- GatingSet(fs_trans)

# None debris
autoplot(fs_trans[[1]],x = 'FSC-H', y = 'SSC-H', bins = 256)

rg1 <- rectangleGate("FSC-H"=c(500000, Inf), filterId="NonDebris")
gs_pop_add(gs, rg1, parent = "root")
gs_get_pop_paths(gs)
recompute(gs)
autoplot(gs,x = 'FSC-H', y = 'SSC-H', "NonDebris", bins = 256)

# singlet
autoplot(fs_trans[[1]],x = 'FSC-H', y = 'FSC-A', bins = 256)
rg2 <- rectangleGate("FSC-H"=c(500000, 2000000),"FSC-A"=c(5000, 2500000))
gs_pop_add(gs, rg2, parent = "NonDebris", name = "singlets")
gs_get_pop_paths(gs)
recompute(gs)
autoplot(gs,x = 'FSC-H', y = 'FSC-A', "singlets", bins = 256)

gs_pop_remove(gs,'singlets')

# singlet with openCyto
library(openCyto)

import_data <- gs_pop_get_data(gs, "NonDebris") #get parent data
singlet_gate <- fsApply(import_data, function(x) openCyto:::.singletGate(x, channels =c("FSC-H", "FSC-A")))
gs_pop_add(gs, singlet_gate, parent = "NonDebris", name = "singlets")
recompute(gs)
autoplot(gs,x = 'FSC-H', y = 'FSC-A', "singlets", bins = 256) 

#live cells
fs[[1]]@parameters@data[1:2]
autoplot(fs_trans[[1]],x = 'LIVE DEAD Blue-A', bins = 256)

import_data <- gs_pop_get_data(gs, "singlets") #get parent data
live_gate <- fsApply(import_data, function(x) openCyto:::.mindensity(x, channels ='LIVE DEAD Blue-A',positive = FALSE))
gs_pop_add(gs, singlet_gate, parent = "singlets", name = "live_cells")
recompute(gs)
gs_pop_remove(gs,'Live')
autoplot(fs_trans[[1]], 'LIVE DEAD Blue-A') + geom_gate(live_gate)
autoplot(gs[1], x = 'LIVE DEAD Blue-A', "live_cells", bins = 256) 

ggcyto(gs[1], aes(x ='LIVE DEAD Blue-A'), subset ='live_cells') + geom_histogram(bins=256)
gs_get_leaf_nodes(gs)
gs_pop_get_stats(gs)

import_data <- gs_pop_get_data(gs, "singlets") #get parent data
live_gate2 <- rectangleGate("LIVE DEAD Blue-A"=c(0, 1.75),filterId="Live")
gs_pop_add(gs, live_gate2, parent = "singlets", name = "Live")
recompute(gs)
autoplot(gs, x = 'LIVE DEAD Blue-A', "Live", bins = 256) 

# exporting from gs
autoplot(gs[[1]], bins=512)
plot(gs)
gs_pop_get_stats(gs)
gs_pop_get_stats(gs, "singlets", "percent")
k<-gs_pop_get_stats(gs, "singlets", "percent")

#move back tot he flowset for other packages
autoplot(gs[[1]], bins=512)
autoplot(gs, "singlets",bins=512)
autoplot(gs, "Live", bins=512)
autoplot(gs, "NonDebris",bins=512)
fs_gated <- gs_pop_get_data(gs, "Live")
fs_gated[[1]]
fs_trans[[1]]

# ggyto
ggplot(fs_gated)
p <- ggcyto(fs_gated, aes(x = `CD4`))
p1 <- p + geom_histogram(bins = 512, fill = 'red') 
p1

p <- ggcyto(fs_gated, aes(x = `CD4`, y=`CD8`))
p1 <- p + geom_hex(bins = 512, fill = 'red') 
p1

p <- ggcyto(fs_gated, aes(x = `CD4`, y=`CD8`))
p1 <- p + geom_hex(bins = 128)
p1

p <- ggcyto(gs, aes(x = `CD4`, y=`CD8`),subset = "Live")
p1 <- p + geom_hex(bins = 128)
p1

p <- ggcyto(gs, aes(x = `LIVE DEAD Blue-A`),subset = "singlets")+ geom_histogram(bins = 128)+geom_gate("Live")
p1 <- p + geom_histogram(bins = 128)
p2 <- p1 +geom_gate("Live") 
p3 <- p2 + xlim(c(1,3))
p3


p <- gs_pop_get_data(gs, "Live")
p1 + geom_overlay(data = fs_trans, size = 0.3, alpha = 0.7)

# add overlay layer by gate name
p + geom_overlay(data = "DNT", size = 0.3, alpha = 0.7)


p <- ggcyto(gs, aes(x = `LIVE DEAD Blue-A`), subset = "singlets") + geom_density(aes(y = ..count..))
p + geom_overlay("Live", aes(y = ..count..), fill = "red")

#useful functions
autovect_verbose<- function(ff){
  c<- data.frame(ff@parameters@data)
  d<- grep("FSC|SSC|TIME|-W|-H", c$name, invert = TRUE, value = TRUE)
  return(unname(d))
}
loadNtransform <- function(ff){
  ff<-tryCatch(read.FCS(ff), error = function(e){read.FCS(ff)})
  x<-autovect_verbose(ff)
  biexp  <- biexponentialTransform("myTransform")
  ff<-transform(ff,transformList(x,biexp))
  return(ff)
}
autovect_verbose(fs[[1]])
loadNtransform(files_list[1])


exprs(fs_gated[[1]])
