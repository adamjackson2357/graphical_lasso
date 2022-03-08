rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

LoadPackages=function(packages){
  for (i in 1:length(packages)){
    suppressPackageStartupMessages(library(packages[i], character.only=TRUE))
  }
}

LoadPackages(c("pheatmap","corpcor","abind","parallel",
               "RColorBrewer","igraph","ppcor","mvtnorm",
               "pROC","glasso","stabs","huge","pulsar",
               "glassoFast","colorspace","glmnet"))
``
source("penalisation_functions.R")
source("QUIC.R")

covars=readRDS("Data/Covariates_selected_proteins.rds")
proteins=readRDS("Data/Proteins_selected_denoised_re.rds")
transcripts=readRDS("Data/Transcripts_log_transformed.rds")

ids=intersect(rownames(proteins), rownames(transcripts))
covars=covars[ids,]
transcripts=transcripts[ids,]
proteins=proteins[ids,]
print(all(rownames(proteins)==rownames(transcripts)))
print(all(rownames(proteins)==rownames(covars)))


#### LASSO REGRESSION ####

# Complete cases
ids_to_exclude=unique(c(which(is.na(covars$packyears)), 
                        which(is.na(covars$age.sample)), 
                        which(is.na(covars$bmi))))
if (length(ids_to_exclude)>0){
  covars=covars[-ids_to_exclude,]
}
proteins=proteins[rownames(covars),]
y=covars$packyears
x=cbind(proteins, age=covars$age.sample, bmi=covars$bmi)

# create the penalty factor
penalty.factor <- c(rep(1, ncol(x) - 2), rep(0, 2))

# run the cv glmnet
model.lasso = cv.glmnet(x=x, y=y, penalty.factor = penalty.factor, type.measure="mse", seed=42)

plot(model.lasso)

# list the variables
lassomodels.1se = coef(model.lasso, s = "lambda.1se")
lasso_vars <- rownames(which(lassomodels.1se!=0, arr.ind = TRUE))[-1]
lasso_vars

#### STABILITY SELECTION ####

# Running stability selection
t0=Sys.time()
out=CalibrateRegression(xdata=x, ydata=y, K=100, tau=0.5, verbose=FALSE, 
                        penalty.factor=c(rep(1,ncol(proteins)),0,0), family="gaussian")
t1=Sys.time()
print(t1-t0)

CalibrationPlot(out)

# Calibrated selection proportions 
selprop=CalibratedSelectionProportionsRegression(out, plot=FALSE)
selprop=selprop[-which(names(selprop)%in%c("age","bmi"))]
print(selprop)

# Calibrated parameters
hat_params=GetArgmax(out)
print(hat_params)

# Visualisation of selection proportions
par(mar=c(10,5,1,1))
plot(selprop, type="h", lwd=3, las=1, xlab="", ylab="Selection Proportion", xaxt="n",
     col=ifelse(selprop>=hat_params[2], yes="red", no="grey"), cex.lab=1.5)
abline(h=hat_params[2], lty=2, col="darkred")
for (i in 1:length(selprop)){
  axis(side=1, at=i, labels=names(selprop)[i], las=2, 
       col=ifelse(selprop[i]>=hat_params[2], yes="red", no="grey"),
       col.axis=ifelse(selprop[i]>=hat_params[2], yes="red", no="grey"))
}

stab_vars <- names(which(selprop>=hat_params[2], arr.ind = TRUE))
stab_vars
lasso_vars

par(mar=c(10,5,1,1))
plot(selprop, type="h", lwd=3, las=1, xlab="", ylab="Selection Proportion", xaxt="n",
     col=ifelse(names(selprop)%in%lasso_vars, yes="blue", no="grey"), cex.lab=1.5)
abline(h=hat_params[2], lty=2, col="darkred")
for (i in 1:length(selprop)){
  axis(side=1, at=i, labels=names(selprop)[i], las=2,
       col=ifelse(names(selprop)[i]%in%lasso_vars, yes="blue", no="grey"),
       col.axis=ifelse(names(selprop)[i]%in%lasso_vars, yes="blue", no="grey"))
}

#### STABILITY ENHANCED GRAPHICAL MODELS ####

# Running stability selection
t0=Sys.time()
out=CalibrateNetwork(data=proteins, K=100, tau=0.5, refine_calib_grid=FALSE, verbose=FALSE)
t1=Sys.time()
print(t1-t0)

CalibrationPlot(out)

# Adjacency matrix of the calibrated network
A=CalibratedAdjacency(out)

# Getting igraph object
mygraph=GetGraph(adjacency=A)

par(mar=rep(0,4))
set.seed(1)
plot(mygraph, layout=layout_with_fr(mygraph))

# highlight smoling related proteins from selprop
mynode_colour=rep("skyblue",ncol(proteins))
names(mynode_colour)=colnames(proteins)
smoking_related_proteins=names(selprop)[selprop>=hat_params[2]]
mynode_colour[smoking_related_proteins]="tomato"

mygraph=GetGraph(adjacency=A, node_color=mynode_colour)

par(mar=rep(0,4))
set.seed(1)
plot(mygraph, layout=layout_with_fr(mygraph))

# get the distances between the smoking-related proteins
print(mean(distances(mygraph)))
print(distances(mygraph)[smoking_related_proteins,smoking_related_proteins])

# Running stability selection under constraint
t0=Sys.time()
out=CalibrateNetwork(data=proteins, K=100, tau=0.5, PFER_thr=20,
                     refine_calib_grid=FALSE, verbose=FALSE)
t1=Sys.time()
print(t1-t0)

CalibrationPlot(out)

# Adjacency matrix of the calibrated network
A=CalibratedAdjacency(out)

# Getting igraph object
mygraph=GetGraph(adjacency=A)

par(mar=rep(0,4))
set.seed(1)
plot(mygraph, layout=layout_with_fr(mygraph))

print(sort(betweenness(mygraph), decreasing=TRUE))

#### MULTI-OMICs NETWORK ####

# Loading the data
proteins=readRDS("Data/Proteins_selected_denoised_re.rds")
transcripts=readRDS("Data/Transcripts_log_transformed.rds")
annot=readRDS("Data/Transcripts_annotation.rds")

# Combining measurements from the same individuals
ids=intersect(rownames(proteins), rownames(transcripts))
proteins=proteins[ids,]
transcripts=transcripts[ids,]
x=cbind(proteins, transcripts)

# Heatmap of correlation
pheatmap(cor(x), cluster_rows=FALSE, cluster_cols=FALSE, breaks=seq(-1,1,length.out=100), border=NA)

# Running stability selection under constraint
t0=Sys.time()
pk=c(ncol(proteins), ncol(transcripts))
Lambda=LambdaGridNetwork(data=x, pk=pk, lambda_path_refined_cardinal=20, PFER_thr=50)
out=CalibrateNetwork(data=x, pk=pk, Lambda=Lambda, lambda_other_blocks=apply(Lambda,2,min), 
                     K=100, tau=0.5, PFER_thr=50,
                     refine_calib_grid=FALSE, verbose=FALSE)
t1=Sys.time()
print(t1-t0)

block_names=c("Protein-protein", "Protein-transcript", "Transcript-transcript")
for (i in 1:3){
  print(block_names[i])
  CalibrationPlot(out, bi_dim = TRUE, block_id=i)
}

# Adjacency matrix of the calibrated network
A=CalibratedAdjacency(out)

# Using different colours for the proteins and transcripts
mynode_colours=c(rep("skyblue",ncol(proteins)),rep("forestgreen",ncol(transcripts)))

# Using gene symbols as labels (transcripts)
mynode_labels=c(colnames(proteins), annot[colnames(transcripts), "Symbol"])
names(mynode_labels)=colnames(x)

# Getting igraph object
mygraph=GetGraph(adjacency=A, node_color=mynode_colours, node_label=mynode_labels)
par(mar=rep(0,4))
set.seed(1)
plot(mygraph, layout=layout_with_fr(mygraph))

block_mat=GetBlockMatrix(pk=pk)
for (i in 1:3){
  print(block_names[i])
  print(sum(A[block_mat==i])/2)
}

#### COLOUR THE SMOKING RELATED TRANSCRIPTS ####

mynode_colours=c(rep("skyblue",ncol(proteins)),rep("forestgreen",ncol(transcripts)))
mynode_labels=c(colnames(proteins), annot[colnames(transcripts), "Symbol"])
names(mynode_labels)=mynode_labels
mynode_colours[which(mynode_labels %in% smoking_related_proteins)]="tomato"


mygraph=GetGraph(adjacency=A, node_color=mynode_colours, node_label=mynode_labels)
par(mar=rep(0,4))
set.seed(1)
plot(mygraph, layout=layout_with_fr(mygraph))
head(distances(mygraph))

