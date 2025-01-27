source('Auxiliar_functions.R')
library(kohonen)
library(caret)
library(aweSOM)
library(factoextra)
library(randomForestExplainer)
library(sf)
library(raster)
library(data.table)
library(GGally)
## Reading family data and pre-processing family data
family_data<-read.csv("family_data.csv",sep=";",header=T, row.names=1)
# Load family data and normalize by row totals
family_relative<-vegan::decostand(family_data,"total")
# Apply a fourth root transformation to reduce the effect of dominant taxa
family_sqrt4<-sqrt(sqrt(family_relative))

## Reading family data and pre-processing coordinates + depth data
coords_depth<-read.csv("Coords_Depth_data.csv",sep=";",row.names=1)
# Load coordinates and depth data and scale the features for consistency
cd_scaled<-scale(coords_depth)

###########################################################
#######                             #######################
#######   UNSUPERVISED LEARNING     #######################
#######                             #######################
###########################################################


## Prepare a list of two layers: environmental features (layer1) and family data (layer2)
training_list<-list(
  layer1=as.matrix(cd_scaled),
  layer2=as.matrix(family_sqrt4)
)

## Defining initial values for the codebook vectors: Setting the number of codebook vectors and initialize them randomly from the dataset
ncodes<-100
set.seed(42)
starters<-sample(1:nrow(cd_scaled), ncodes, replace = FALSE)
init<-lapply(training_list, function(x) x[starters, , drop = FALSE])


# Train a Self-Organizing Map (SOM) with two layers using hexagonal grid topology
set.seed(42)
som_model<-kohonen::supersom(
  data=training_list,
  init=init,
  radius=c(6.082763,0),
  alpha=c(0.05, 0.01),
  maxNA.fraction=c(0.001),
  dist.fcts=c("euclidean","euclidean"),
  grid=somgrid(x=10,y=10,topo="hexagonal",neighbourhood.fct="gaussian",toroidal=F),
  rlen=500,
  user.weights=c(0.05,0.95)
)

## SOM Evaluation metrics
# Assign sample names to SOM classification results
names(som_model$unit.classif)<-rownames(cd_scaled)
# Evaluate the SOM quality for layer1 (environmental features)
aweSOM::somQuality(som_model,as.matrix(cd_scaled))
# Temporarily switch to layer2 (family data) and re-evaluate SOM quality
som_model_temp<-som_model
som_model_temp$codes[[1]]<-som_model_temp$codes[[2]]
aweSOM::somQuality(som_model_temp,as.matrix(family_sqrt4))

## Hierarchical clustering
# Compute distances between SOM codebook vectors
dist_codes<-kohonen::object.distances(som_model,"codes")
# Perform hierarchical clustering with Ward's method into 5 clusters
hc_model<-factoextra::hcut(dist_codes,k=5,isdiss =T,hc_fun="hclust",hc_method = "ward.D2")
hcut_result<-hc_model$cluster
names(hcut_result)<-1:nrow(som_model$codes[[1]])

## Assigning samples to the SOM clusters: mapping each sample to a SOM cluster based on hierarchical clustering results
dfcut<-data.frame(neu=names(hcut_result),hcut_result)
list<-split(data.frame(id=names(som_model$unit.classif),neu=som_model$unit.classif),som_model$unit.classif)
res<-do.call(rbind,lapply(names(list),function(i){
  x<-list[[i]]
  x$hc<- dfcut[i,"hcut_result"]
  x
}))
newclass<-res$hc
names(newclass)<-rownames(res)
newclass<-newclass[rownames(som_model$data[[1]])]

## Ordering cluster levels by sample depth
fac<-cd_scaled[,"sample_depth",drop=F]
clusters<-newclass
newlevels<-names(sort(tapply(fac[,1],as.numeric(as.character(clusters)),mean)))
hc5<-factor(clusters,levels=rev(newlevels),labels=1:length(newlevels))
hcut_result<-factor(hcut_result,levels=rev(newlevels),labels=1:length(newlevels))

## Dendrogram
gg_dendrogram(hcut_result,hc_model, theme="theme_minimal")



## Indicators

set.seed(42)
indicators<-indicator_multipart(family_sqrt4, hc5,5,nperm=999)
indicator_taxa<-names(sort(filter_multipart_result(indicators,5)))
summary_indicators<-summary_multipatt(indicators)
#write.table(summary_indicators,'summary_indicators.csv',sep=";", col.names=NA)


###########################################################
#######                             #######################
#######     SUPERVISED LEARNING     #######################
#######                             #######################
###########################################################

# RF CLASSIFICATION - Microbial Associations

# Step 1: Split the data into training (80%) and test (20%) sets, stratifying by clusters
partition<-rep("test",length(hc5))
set.seed(42)
test_obs<-caret::createDataPartition(hc5,p=0.8)
partition[test_obs[[1]]]<-"training"
partition<-factor(partition)
environ_data<-read.csv("environ_data.csv",sep=";",header=T, row.names=1)
environ_datalist<-split(environ_data,partition)
hc5_list<-split(hc5,partition)
train_data<-environ_datalist$training
test_data<-environ_datalist$test
train_y<-hc5_list$training
test_y<-hc5_list$test

# Step 2: Train a Random Forest classification model with hyperparameter tuning and cross-validation
set.seed(42)
rf_model<-caret::train(
  x=train_data,
  y=train_y,
  method = "rf", ntree = 500L, replace = TRUE, nodesize = 1L,
  nPerm = 1L, norm.votes = TRUE, localImp = TRUE, keep.forest = TRUE,
  keep.inbag = TRUE,
  trControl = trainControl(method = "repeatedcv",
                           number = 5L,
                           repeats = 5L,
                           search = "grid",
                           p = 0.1,
                           initialWindow = NULL,
                           horizon = 1,
                           fixedWindow = TRUE,
                           skip = 0, verboseIter = FALSE,
                           returnData = TRUE,
                           returnResamp = "final",
                           savePredictions = "all"),
  tuneLength =5
)

## Feature importance
imp_rf_classication<-multipimp( rf_model,measures=c(c("mean_min_depth",
                                                      "accuracy_decrease",
                                                      "gini_decrease" ,
                                                      "no_of_nodes",
                                                      "times_a_root"),'p_value','no_of_trees'), mean_sample="top_trees")



# Step 3: Compute class proportions in the training data for stratified baseline
class_distribution <- table(train_y) / length(train_y)

# Step 4: Define stratified baseline prediction function
stratified_baseline <- function(test_data, class_distribution) {
  sample(names(class_distribution),
         size = nrow(test_data),
         replace = TRUE,
         prob = as.numeric(class_distribution)) # Use class probabilities
}

# Step 5: Generate stratified predictions based on the class distribution
set.seed(42) # Ensure reproducibility
stratified_predictions <- stratified_baseline(test_data, class_distribution)

# Step 6: Evaluate stratified baseline using confusion matrix
stratified_conf_matrix <- confusionMatrix(as.factor(stratified_predictions), test_y)

# Print stratified baseline performance
cat("Stratified Baseline Confusion Matrix:\n")
print(stratified_conf_matrix)

# Step 7: Obtain Random Forest predictions
rf_predictions <- predict(rf_model, test_data)

# Step 8: Evaluate Random Forest performance
rf_conf_matrix <- confusionMatrix(as.factor(rf_predictions), test_y)

# Print Random Forest performance
cat("Random Forest Confusion Matrix:\n")
print(rf_conf_matrix)

# Step 9: Perform McNemar test for error patterns between RF and stratified baseline
conf_matrix_rf_stratified <- table(
  RF = rf_predictions == test_y,
  Baseline = stratified_predictions == test_y
)

mcnemar_test_stratified <- mcnemar.test(conf_matrix_rf_stratified)
cat("McNemar Test p-value (RF vs Stratified Baseline):", mcnemar_test_stratified$p.value, "\n")

# Step 10: Bootstrap for accuracy difference (RF vs Stratified Baseline)
accuracy <- function(actual, predicted) {
  mean(actual == predicted)
}

set.seed(42)
bootstrap_diff_stratified <- replicate(1000, {
  sample_indices <- sample(seq_along(test_y), replace = TRUE)
  test_y_sample <- test_y[sample_indices]
  rf_sample_preds <- rf_predictions[sample_indices]
  stratified_sample_preds <- stratified_predictions[sample_indices]
  accuracy(test_y_sample, rf_sample_preds) - accuracy(test_y_sample, stratified_sample_preds)
})

# Step 11: Confidence interval for accuracy difference
ci_stratified <- quantile(bootstrap_diff_stratified, probs = c(0.025, 0.975))
cat("95% CI for Accuracy Difference (RF vs Stratified Baseline):", ci_stratified, "\n")
###################################################################################

# RF REGRESSION - CITOMETRY


# Load cytometry data and log-transform (log10) the response variables to normalize them
cito_data<-read.csv("Cito_data.csv",sep=";",row.names=1)
cito_log10<-log10(cito_data+1)
envi_cito<-environ_data[rownames(cito_log10),]


## SYNECHOCOCCUS
# Step 1: Split data into training (80%) and test (20%) sets for Synechococcus
synec<-cito_log10$SYNEC_CELL.ML
partition_synec<-rep("test",length(synec))
set.seed(42)
test_synec<-caret::createDataPartition(synec,p=0.8)
partition_synec[test_synec[[1]]]<-"training"
partition_synec<-factor(partition_synec)


envi_synec_datalist<-split(envi_cito,partition_synec)
train_x_synec<-envi_synec_datalist$training
test_x_synec<-envi_synec_datalist$test
synec_list<-split(synec,partition_synec)
train_y_synec<-synec_list$training
test_y_synec<-synec_list$test

# Step 2: Train a Random Forest regression model for Synechococcus
set.seed(42)
rf_synec<-caret::train(
  x=train_x_synec,
  y=train_y_synec,
  method = "rf", ntree = 500L, replace = TRUE, nodesize = 1L,
  nPerm = 1L, norm.votes = TRUE, localImp = TRUE, keep.forest = TRUE,
  keep.inbag = TRUE,
  trControl = trainControl(method = "repeatedcv",
                           number = 5L,
                           repeats = 5L,
                           search = "grid",
                           p = 0.1,
                           initialWindow = NULL,
                           horizon = 1,
                           fixedWindow = TRUE,
                           skip = 0, verboseIter = FALSE,
                           returnData = TRUE,
                           returnResamp = "final",
                           savePredictions = "all"),
  tuneLength =5
)

imp_rf_synec<-multipimp( rf_synec,measures=c(c('mean_min_depth', 'mse_increase', 'node_purity_increase', 'no_of_nodes', 'times_a_root'),'p_value','no_of_trees'), mean_sample="top_trees")
plot_feature_importance(imp_rf_synec)


# Step 3: Generate naive baseline predictions (mean of training data)
# Root mean square error function
rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2))
}

# Step 3: Generate naive baseline predictions (mean of training data)
naive_predictions_synec <- rep(mean(train_y_synec), length(test_y_synec))

# Step 4: Calculate RMSE for both Random Forest and naive baseline
rf_rmse_synec <- rmse(test_y_synec, predict(rf_synec, test_x_synec))
naive_rmse_synec <- rmse(test_y_synec, naive_predictions_synec)

# Step 5: Perform bootstrap to assess RMSE improvement for Random Forest over baseline
set.seed(42)
bootstrap_rmse_diff_synec <- replicate(1000, {
  sample_indices <- sample(seq_along(test_y_synec), replace = TRUE)
  test_y_sample <- test_y_synec[sample_indices]
  rf_sample_preds <- predict(rf_synec, test_x_synec[sample_indices, ])
  naive_sample_preds <- naive_predictions_synec[sample_indices]
  rmse(test_y_sample, rf_sample_preds) - rmse(test_y_sample, naive_sample_preds)
})

# Step 6: Calculate confidence interval for RMSE difference
ci_synec <- quantile(bootstrap_rmse_diff_synec, probs = c(0.025, 0.975))
cat("RF RMSE:", rf_rmse_synec, "\n")
cat("Naive RMSE:", naive_rmse_synec, "\n")
cat("95% CI for RMSE Difference:", ci_synec, "\n")

###################################################################################

#PROCHLOROCOCCUS
# Step 1: Split data into training (80%) and test (20%) sets for Prochlorococcus
proc<-cito_log10$PROC_CELL.ML
partition_proc<-rep("test",length(proc))
set.seed(42)
test_proc<-caret::createDataPartition(proc,p=0.8)
partition_proc[test_proc[[1]]]<-"training"
partition_proc<-factor(partition_proc)

envi_proc_datalist<-split(envi_cito,partition_proc)
train_x_proc<-envi_proc_datalist$training
test_x_proc<-envi_proc_datalist$test
proc_list<-split(proc,partition_proc)
train_y_proc<-proc_list$training
test_y_proc<-proc_list$test

# Step 2: Train a Random Forest regression model for Prochlorococcus
set.seed(42)
rf_proc<-caret::train(
  x=train_x_proc,
  y=train_y_proc,
  method = "rf", ntree = 500L, replace = TRUE, nodesize = 1L,
  nPerm = 1L, norm.votes = TRUE, localImp = TRUE, keep.forest = TRUE,
  keep.inbag = TRUE,
  trControl = trainControl(method = "repeatedcv",
                           number = 5L,
                           repeats = 5L,
                           search = "grid",
                           p = 0.1,
                           initialWindow = NULL,
                           horizon = 1,
                           fixedWindow = TRUE,
                           skip = 0, verboseIter = FALSE,
                           returnData = TRUE,
                           returnResamp = "final",
                           savePredictions = "all"),
  tuneLength =5
)
imp_rf_proc<-multipimp( rf_proc,measures=c(c('mean_min_depth', 'mse_increase', 'node_purity_increase', 'no_of_nodes', 'times_a_root'),'p_value','no_of_trees'), mean_sample="top_trees")
plot_feature_importance(imp_rf_proc)
# Step 3: Generate naive baseline predictions (mean of training data)
naive_predictions_proc <- rep(mean(train_y_proc), length(test_y_proc))

# Step 4: Calculate RMSE for both Random Forest and naive baselin
rf_rmse_proc <- rmse(test_y_proc, predict(rf_proc, test_x_proc))
naive_rmse_proc <- rmse(test_y_proc, naive_predictions_proc)

# Step 5: Perform bootstrap to assess RMSE improvement for Random Forest over baseline
set.seed(42)
bootstrap_rmse_diff_proc <- replicate(1000, {
  sample_indices <- sample(seq_along(test_y_proc), replace = TRUE)
  test_y_sample <- test_y_proc[sample_indices]
  rf_sample_preds <- predict(rf_proc, test_x_proc[sample_indices, ])
  naive_sample_preds <- naive_predictions_proc[sample_indices]
  rmse(test_y_sample, rf_sample_preds) - rmse(test_y_sample, naive_sample_preds)
})

# Step 6: Calculate confidence interval for RMSE difference
ci_proc <- quantile(bootstrap_rmse_diff_proc, probs = c(0.025, 0.975))
cat("RF RMSE:", rf_rmse_proc, "\n")
cat("Naive RMSE:", naive_rmse_proc, "\n")
cat("95% CI for RMSE Difference:", ci_proc, "\n")

###################################################################################

# HETEROTROPHIC BACTERIA
# Step 1: Split data into training (80%) and test (20%) sets for heterotrophic bacteria
bact_hetero<-cito_log10$BACT_HETERO_CELL.ML
partition_bact_hetero<-rep("test",length(bact_hetero))
set.seed(42)
test_bact_hetero<-caret::createDataPartition(bact_hetero,p=0.8)
partition_bact_hetero[test_bact_hetero[[1]]]<-"training"
partition_bact_hetero<-factor(partition_bact_hetero)
envi_bact_hetero_datalist<-split(envi_cito,partition_bact_hetero)
train_x_bact_hetero<-envi_bact_hetero_datalist$training
test_x_bact_hetero<-envi_bact_hetero_datalist$test
bact_hetero_list<-split(bact_hetero,partition_bact_hetero)
train_y_bact_hetero<-bact_hetero_list$training
test_y_bact_hetero<-bact_hetero_list$test
# Step 2: Train a Random Forest regression model for heterotrophic bacteria
set.seed(42)
rf_bact_hetero<-caret::train(
  x=train_x_bact_hetero,
  y=train_y_bact_hetero,
  method = "rf", ntree = 500L, replace = TRUE, nodesize = 1L,
  nPerm = 1L, norm.votes = TRUE, localImp = TRUE, keep.forest = TRUE,
  keep.inbag = TRUE,
  trControl = trainControl(method = "repeatedcv",
                           number = 5L,
                           repeats = 5L,
                           search = "grid",
                           p = 0.1,
                           initialWindow = NULL,
                           horizon = 1,
                           fixedWindow = TRUE,
                           skip = 0, verboseIter = FALSE,
                           returnData = TRUE,
                           returnResamp = "final",
                           savePredictions = "all"),
  tuneLength =5
)

imp_rf_bact_hetero<-multipimp( rf_bact_hetero,measures=c(c('mean_min_depth', 'mse_increase', 'node_purity_increase', 'no_of_nodes', 'times_a_root'),'p_value','no_of_trees'), mean_sample="top_trees")

plot_feature_importance(imp_rf_bact_hetero)

# Step 3: Generate naive baseline predictions (mean of training data)
naive_predictions_bact_hetero <- rep(mean(train_y_bact_hetero), length(test_y_bact_hetero))

# Step 4: Calculate RMSE for both Random Forest and naive baseline
rf_rmse_bact_hetero <- rmse(test_y_bact_hetero, predict(rf_bact_hetero, test_x_bact_hetero))
naive_rmse_bact_hetero <- rmse(test_y_bact_hetero, naive_predictions_bact_hetero)

# Step 5: Perform bootstrap to assess RMSE improvement for Random Forest over baseline
set.seed(42)
bootstrap_rmse_diff_bact_hetero <- replicate(1000, {
  sample_indices <- sample(seq_along(test_y_bact_hetero), replace = TRUE)
  test_y_sample <- test_y_bact_hetero[sample_indices]
  rf_sample_preds <- predict(rf_bact_hetero, test_x_bact_hetero[sample_indices, ])
  naive_sample_preds <- naive_predictions_bact_hetero[sample_indices]
  rmse(test_y_sample, rf_sample_preds) - rmse(test_y_sample, naive_sample_preds)
})

# Step 6: Calculate confidence interval for RMSE difference
ci_bact_hetero <- quantile(bootstrap_rmse_diff_bact_hetero, probs = c(0.025, 0.975))
cat("RF RMSE:", rf_rmse_bact_hetero, "\n")
cat("Naive RMSE:", naive_rmse_bact_hetero, "\n")
cat("95% CI for RMSE Difference:", ci_bact_hetero, "\n")

###################################################################################

# PROKARYOTES
Prokaryotes<-cito_log10$Prokaryotes_CELL.ML
# Step 1: Split data into training (80%) and test (20%) sets for prokaryotes
partition_Prokaryotes<-rep("test",length(Prokaryotes))
set.seed(42)
test_Prokaryotes<-caret::createDataPartition(Prokaryotes,p=0.8)
partition_Prokaryotes[test_Prokaryotes[[1]]]<-"training"
partition_Prokaryotes<-factor(partition_Prokaryotes)
envi_Prokaryotes_datalist<-split(envi_cito,partition_Prokaryotes)
train_x_Prokaryotes<-envi_Prokaryotes_datalist$training
test_x_Prokaryotes<-envi_Prokaryotes_datalist$test
Prokaryotes_list<-split(Prokaryotes,partition_Prokaryotes)
train_y_Prokaryotes<-Prokaryotes_list$training
test_y_Prokaryotes<-Prokaryotes_list$test

# Step 2: Train a Random Forest regression model for prokaryotes
set.seed(42)
rf_Prokaryotes<-caret::train(
  x=train_x_Prokaryotes,
  y=train_y_Prokaryotes,
  method = "rf", ntree = 500L, replace = TRUE, nodesize = 1L,
  nPerm = 1L, norm.votes = TRUE, localImp = TRUE, keep.forest = TRUE,
  keep.inbag = TRUE,
  trControl = trainControl(method = "repeatedcv",
                           number = 5L,
                           repeats = 5L,
                           search = "grid",
                           p = 0.1,
                           initialWindow = NULL,
                           horizon = 1,
                           fixedWindow = TRUE,
                           skip = 0, verboseIter = FALSE,
                           returnData = TRUE,
                           returnResamp = "final",
                           savePredictions = "all"),
  tuneLength =5
)
imp_rf_prok<-multipimp( rf_Prokaryotes,measures=c(c('mean_min_depth', 'mse_increase', 'node_purity_increase', 'no_of_nodes', 'times_a_root'),'p_value','no_of_trees'), mean_sample="top_trees")
plot_feature_importance(imp_rf_prok)

# Step 3: Generate naive baseline predictions (mean of training data)
naive_predictions_Prokaryotes <- rep(mean(train_y_Prokaryotes), length(test_y_Prokaryotes))

# Step 4: Calculate RMSE for both Random Forest and naive baseline
rf_rmse_Prokaryotes <- rmse(test_y_Prokaryotes, predict(rf_Prokaryotes, test_x_Prokaryotes))
naive_rmse_Prokaryotes <- rmse(test_y_Prokaryotes, naive_predictions_Prokaryotes)

# Step 5: Perform bootstrap to assess RMSE improvement for Random Forest over baseline
set.seed(42)
bootstrap_rmse_diff_Prokaryotes <- replicate(1000, {
  sample_indices <- sample(seq_along(test_y_Prokaryotes), replace = TRUE)
  test_y_sample <- test_y_Prokaryotes[sample_indices]
  rf_sample_preds <- predict(rf_Prokaryotes, test_x_Prokaryotes[sample_indices, ])
  naive_sample_preds <- naive_predictions_Prokaryotes[sample_indices]
  rmse(test_y_sample, rf_sample_preds) - rmse(test_y_sample, naive_sample_preds)
})

# Step 6: Calculate confidence interval for RMSE difference
ci_Prokaryotes <- quantile(bootstrap_rmse_diff_Prokaryotes, probs = c(0.025, 0.975))
cat("RF RMSE:", rf_rmse_Prokaryotes, "\n")
cat("Naive RMSE:", naive_rmse_Prokaryotes, "\n")
cat("95% CI for RMSE Difference:", ci_Prokaryotes, "\n")



## Influence Score


pred_synec <-  predict(rf_synec)
imp0<-caret::varImp(rf_synec, scale=T)$importance
inf<- cor(train_x_synec, pred_synec, use = "complete.obs")
influences<-imp0*inf
influences_synec<-data.frame(var=paste0("Synec_",rownames(influences)),value=influences[,1],group="Synechococcus")


imp0<-caret::varImp(rf_proc, scale=T)$importance
pred_proc <-  predict(rf_proc)
inf<- cor(train_x_proc, pred_proc, use = "complete.obs")
influences<-imp0*inf
influences_proc<-data.frame(var=paste0("Proc_",rownames(influences)),value=influences[,1],group="Prochlorococcus")


imp0<-caret::varImp(rf_bact_hetero, scale=T)$importance
pred_bact_hetero <-  predict(rf_bact_hetero)
inf<- cor(train_x_bact_hetero, pred_bact_hetero, use = "complete.obs")
influences<-imp0*inf
influences_bact<-data.frame(var=paste0("Bact_Hetero_",rownames(influences)),value=influences[,1],group="Bact_Hetero")

imp0<-caret::varImp(rf_Prokaryotes, scale=T)$importance
pred_Prokaryotes <-  predict(rf_Prokaryotes)
inf<- cor(train_x_Prokaryotes, pred_bact_hetero, use = "complete.obs")
influences<-imp0*inf
influences_prok<-data.frame(var=paste0("Prok_",rownames(influences)),value=influences[,1],group="Prokaryotes")


influences<-rbind(influences_bact,influences_synec,influences_proc,influences_prok)

{
  pdf("cito_influence_score.pdf", width=8, height = 8)
  ggplot(influences,aes(x=value,y=reorder(var,value), fill=value))+geom_bar(stat="identity",show.legend = F)+scale_y_discrete(
    breaks =  c(
      "Bact_Hetero_temp", "Bact_Hetero_sal", "Bact_Hetero_density", "Bact_Hetero_clo_a",
      "Bact_Hetero_spm", "Bact_Hetero_poc", "Bact_Hetero_pon", "Bact_Hetero_poc_pon",
      "Bact_Hetero_od_Lab", "Bact_Hetero_no2", "Bact_Hetero_no3", "Bact_Hetero_po4",
      "Bact_Hetero_siO4", "Bact_Hetero_nh4", "Bact_Hetero_ph_lab", "Bact_Hetero_cruise",
      "Synec_temp", "Synec_sal", "Synec_density", "Synec_clo_a",
      "Synec_spm", "Synec_poc", "Synec_pon", "Synec_poc_pon",
      "Synec_od_Lab", "Synec_no2", "Synec_no3", "Synec_po4",
      "Synec_siO4", "Synec_nh4", "Synec_ph_lab", "Synec_cruise",
      "Proc_temp", "Proc_sal", "Proc_density", "Proc_clo_a",
      "Proc_spm", "Proc_poc", "Proc_pon", "Proc_poc_pon",
      "Proc_od_Lab", "Proc_no2", "Proc_no3", "Proc_po4",
      "Proc_siO4", "Proc_nh4", "Proc_ph_lab", "Proc_cruise",
      "Prok_temp", "Prok_sal", "Prok_density", "Prok_clo_a",
      "Prok_spm", "Prok_poc", "Prok_pon", "Prok_poc_pon",
      "Prok_od_Lab", "Prok_no2", "Prok_no3", "Prok_po4",
      "Prok_siO4", "Prok_nh4", "Prok_ph_lab", "Prok_cruise"
    ),
    labels = c(
      expression(Temperature),              # Bact_Hetero_temp
      expression(Salinity),                 # Bact_Hetero_sal
      expression(Density),                  # Bact_Hetero_density
      expression(Chlorophyll~a),            # Bact_Hetero_clo_a
      expression(SPM),                      # Bact_Hetero_spm
      expression(POC),                      # Bact_Hetero_poc
      expression(PON),                      # Bact_Hetero_pon
      expression(POC:PON),                  # Bact_Hetero_poc_pon
      expression(OD),      # Bact_Hetero_od_Lab
      expression(NO[2]),                    # Bact_Hetero_no2
      expression(NO[3]),                    # Bact_Hetero_no3
      expression(PO[4]),                    # Bact_Hetero_po4
      expression(SiO[4]),                   # Bact_Hetero_siO4
      expression(NH[4]),                    # Bact_Hetero_nh4
      expression(pH),                 # Bact_Hetero_ph_lab
      expression(Cruise),              # Bact_Hetero_cruise
      expression(Temperature),              # Synec_temp
      expression(Salinity),                 # Synec_sal
      expression(Density),                  # Synec_density
      expression(Chlorophyll~a),            # Synec_clo_a
      expression(SPM),                      # Synec_spm
      expression(POC),                      # Synec_poc
      expression(PON),                      # Synec_pon
      expression(POC:PON),                  # Synec_poc_pon
      expression(OD),      # Synec_od_Lab
      expression(NO[2]),                    # Synec_no2
      expression(NO[3]),                    # Synec_no3
      expression(PO[4]),                    # Synec_po4
      expression(SiO[4]),                   # Synec_siO4
      expression(NH[4]),                    # Synec_nh4
      expression(pH),                 # Synec_ph_lab
      expression(Cruise),              # Synec_cruise
      expression(Temperature),              # Proc_temp
      expression(Salinity),                 # Proc_sal
      expression(Density),                  # Proc_density
      expression(Chlorophyll~a),            # Proc_clo_a
      expression(SPM),                      # Proc_spm
      expression(POC),                      # Proc_poc
      expression(PON),                      # Proc_pon
      expression(POC:PON),                  # Proc_poc_pon
      expression(OD),      # Proc_od_Lab
      expression(NO[2]),                    # Proc_no2
      expression(NO[3]),                    # Proc_no3
      expression(PO[4]),                    # Proc_po4
      expression(SiO[4]),                   # Proc_siO4
      expression(NH[4]),                    # Proc_nh4
      expression(pH),                 # Proc_ph_lab
      expression(Cruise),              # Proc_cruise
      expression(Temperature),              # Prok_temp
      expression(Salinity),                 # Prok_sal
      expression(Density),                  # Prok_density
      expression(Chlorophyll~a),            # Prok_clo_a
      expression(SPM),                      # Prok_spm
      expression(POC),                      # Prok_poc
      expression(PON),                      # Prok_pon
      expression(POC:PON),                  # Prok_poc_pon
      expression(OD),      # Prok_od_Lab
      expression(NO[2]),                    # Prok_no2
      expression(NO[3]),                    # Prok_no3
      expression(PO[4]),                    # Prok_po4
      expression(SiO[4]),                   # Prok_siO4
      expression(NH[4]),                    # Prok_nh4
      expression(pH),                 # Prok_ph_lab
      expression(Cruise)               # Prok_cruise
    )
  )+facet_wrap(~group,scales="free")+scale_fill_continuous()+theme_light()+theme(strip.text =  element_text(color = "black"))+ylab("")+xlab("Influence Score")
  graphics.off()
  shell('cito_influence_score.pdf')
}
