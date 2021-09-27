# FR-FCM-Z3WR 
# compensate, clean and transform
# Show me the change in CD8 alphabeta T cells overt time in the three patients
# Cells -> Single -> Live -> CD3+ -> TCR- -> CD8+ -> -> CD45RA- -> CD161+
# You will find it easier to look up how to make a gating template in openCyto to do this

#load required packages
library('flowCore')
library('ggcyto')
library('flowAI')
library('flowWorkspace')
library('openCyto')

#load data
files_list <- list.files('C:/Users/hallc/Data Transfers/FlowRepository_FR-FCM-Z3WR_files',full.names = TRUE, patter = 'fcs')
fs<-read.flowSet(files_list, truncate_max_range = FALSE)

#compensate (is it needed)
spillover(fs[[1]])
fs[[1]]@description$`$CYT` # Aurora so not needed, these are the instrument unmixed files

#transform the data
fs[[1]]@parameters@data[1:2]
trans<-estimateLogicle(fs[[1]], channels = colnames(fs[[1]][,7:42]), t = 4194304)
fs_trans<-transform(fs,trans)
autoplot(fs_trans[[1]])

#clean the data
fs_clean<-flow_auto_qc(fs_trans, remove_from = "FR_FS", fcs_QC = FALSE)

#gating the cells ## Usle plotGates()
autoplot(fs_clean[[1]]) # lets look at the dimensions
autoplot(fs_clean[[1]], x="TCR gd")

gs <- GatingSet(fs_clean) # make gatingset

library(data.table)
gt <- gatingTemplate('C:/Users/hallc/Desktop/gatingtemplate.csv') #laod template
plot(gt)
gt_gating(gt, gs)
autoplot(gs[[1]])

#export the stats
gs_pop_get_stats(gs, "CD161+", type = pop.MFI)
t<-gs_pop_get_stats(gs, type="percent")
write.csv(t,"stats.csv")

#plot the output
autoplot(gs[[1]], bins=256) #plot 1
lapply(gs, function(x)autoplot(x)) # plot them all

#output for high dimension tools
fs_gated<-gs_pop_get_data(gs, "CD161+")
fs_gated<-cytoset_to_flowSet(fs_gated)

fs[[1]]
fs_gated[[1]]

