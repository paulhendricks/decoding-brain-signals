directory <- "E:/Brain_Competition_OnPrem"
functionfile <- paste(directory, "brain_competition_functions.R", sep="")
datafile <- paste(directory, "ecog_train_with_labels.csv", sep="")
outputtemplatefile <- paste(directory, "ecog_train_templates.csv", sep="")
model_rda_file <- paste(directory, "logitmodel.rda", sep="")

source(functionfile)
library(glmnet)
# Step 5.1: read the data into a dataframe dataset1 from a local file, defined in variable datafile, 
# assuming that you have downloaded the data from 
# http://az754797.vo.msecnd.net/competition/ecog/datasets/ecog_train_with_labels.csv

dataset1 <- read.csv(datafile, sep=",", header=T, stringsAsFactor=F)
summary(dataset1)

# Step 5.2, split the data into training and validation. For each patient, we take the first 150 stimulus presentation cycles into training,
# and the remaining 50 stimulus presentation cycles into validation.
# We name these two datasets as dataset_train, and dataset_valid
col_names <- colnames(dataset1)
training_portion <- 3/4 # Taking the first 75% stimulus presentation cycles as the training data, for each patient
unique_patients <- unique(dataset1[,1]) # get the list of unique patient ID
num_patients <- length(unique_patients) # get the number of unique patients
num_cols <- ncol(dataset1)
dataset_train <- NULL
dataset_valid <- NULL
for (i in 1:num_patients){
  data_i <- dataset1[dataset1[,1]==unique_patients[i],] #get the data of a patient
  nrows_i <- nrow(data_i)
  num_cols_i <- sum(data_i[1,] != -999999) # get the number of columns that are not valued -999999, which indicates that this column does not have signals. 
  signal_i <- data_i[,2:(num_cols_i-2)] # get the signals from valid channels for this patient.
  stim_i <- as.matrix(data_i[,(num_cols-1)]) # get the stimulus column of this patient, and convert it to matrix, which is required by function fh_get_events()
  
  ## create events vector. It returns two columns: col1: stimulus start time, col2: stimulus type 1 for house, 2 for face
  events_i=fh_get_events(stim_i); #get the events for this patient
  events_12_i= events_i[which(events_i[,3]!=0),] # only need the stimulus onset time
  events_12_i=events_12_i[,-2]
  num_stimulus_train <- floor(nrow(events_12_i)*training_portion) #number of events to be put in the training data
  train_last_time <- events_12_i[num_stimulus_train,1] + 399 #the last record in the training data should be the stimulus onset time of the last stimulus presentation cycle + 399 milliseconds
  dataset_train <- rbind(dataset_train, data_i[1:train_last_time,]) #add the training records of this patient to data.set
  dataset_valid <- rbind(dataset_valid, data_i[(train_last_time+1):nrows_i,])
}

dataset_train <- as.data.frame(dataset_train) #convert it to a data frame
dataset_valid <- as.data.frame(dataset_valid) #convert it to a data frame

# Step 5.3, create templates for each channel, stimulus type (house or face), and for each patient
ncols <- ncol(dataset_train) # get the number of columns of the training data. Training data has 67 columns. Col1: PatientID, Cols2-65, 64 channels, Col 66: Stimulus Type, Col 67: Stimulus Presentation Cycle ID
unique_patients <- unique(dataset_train[,1]) # Get the list of unique patients
num_patients <- length(unique_patients) # Get the number of unique patients
templates <- as.data.frame(matrix(NA, nrow=1201*num_patients, ncol=ncols-3)) # every patient will have 1201 rows, where the first row is the standard deviation of signals, which will be used to normalize the signal
# rows 2-601 will be the house template, and rows 602-1201 will be the face template. A template is defined as the average of the signal between 200 ms before the onset of stimulus,
# and 399 ms after the onset of stimulus. So totally 600 points for a house template, or a face template. 
PatientID_Col <- matrix(rep('',1201*num_patients), nrow=1201*num_patients,ncol=1) # A column of patientID, which is going to be cbind with templates

# Start building templates for each patient, channel, and stimulus class (house or face)
for (j in 1:num_patients){
  patient_id <- unique_patients[j] # Get the current patient ID
  PatientID_Col[((j-1)*1201+1):(j*1201),1] <- patient_id # Assign the same patient ID to the patientID column
  data_j <- dataset_train[dataset_train[,1]==patient_id,] # get the data of this specific patient
  ncols_j <- sum(data_j[1,] != -999999) # Determine how many valid columns this patient has (column -999999 means that this patient does not have that signal channel)
  signal_train <- as.matrix(data_j[,2:(ncols_j-2)]) # get the signal for this patient, excluding those -999999 channels
  signal_train <- apply(signal_train, 2, as.numeric) # convert the signal data to numeric, in case they might be treated as string features
  
  stim_train <- as.matrix(data_j[,ncols-1]) # get the column of stim for this patient
  events_train=fh_get_events(stim_train); # get the event matrix
  events_train= events_train[which(events_train[,3]!=0),] # Only keep the stimulus onset time in the trainign data
  events_train=events_train[,-2] # exclude the midway of events column
  train_t=c(1:nrow(stim_train)); # train_t is the row index of training data
  train_e=events_train; # make a copy of the stimulus onset time data 
  num_chans=ncol(signal_train); # get the number of channels this patient has.
  tlims=c(-200,399);# times to start and end erps, this is the time window we are going to add to the stimulus onset time later. The [stimulus onset time -200, stimulus onset time+399] is the time window we used to construct the templates
  erp_baseline=c(-200,49);# times to calcualte erp based upon(must be within tlims) # We take the [stimulus onset time -200, stimulus onset time+49] as the baseline for each stimulus presentation cycle, assuming that the brain has not reponded to the 
  # visual stimulus within 50 milliseconds after the stimulus onset
  train_chans_sd <- rep(0,num_chans) # initialite a variable to hold the standard deviation of signals for each channel. It will be used to normalize the signals.
  for (k in 1:num_chans){
    train_chans_sd[k] <- sd(signal_train[,k]) # get the standard deviation of each channel in the training data
    
    signal_train[,k] <- signal_train[,k]/sd(signal_train[,k]);# Normalize the scale of each signal by dividing by its own standard deviation
  }
  
  # This function generates the templates for each signal, and for house and face stimulus separately,
  # It is just the average of each signal in the 600 ms window (-200ms before stimulus, and 399 ms after stimulus for each stimulus type
  # over all stimulus presentation cycles in the training data.
  fh_sta = function(inputdata,events,fh_class,tlims) {
    cls_times= events[which(events[,2]==fh_class),1] # get the stimulus onset time of a specified class (this is the onset time of a stimulus, not the ending time of the previous ISI)
    sta=matrix(data=0, nrow=(tlims[2]-tlims[1]+1),ncol=ncol(inputdata))
    for (k in 1:length(cls_times)){
      sta=sta+inputdata[cls_times[k]+c(tlims[1]:tlims[2]),]; #accumulating the signals after realigning all stimulus presentation cycles along the stimulus onset time
    }
    sta=sta/k; # calculate the average
    return(sta) # output the average as the template
  } 
  
  #get sta templates
  
  sta_h=fh_sta(signal_train,train_e,1,tlims);# templates of house stimulus
  sta_f=fh_sta(signal_train,train_e,2,tlims);# templates of face stimulus
  
  # recenter stas w.r.t. baseline
  for (k in 1:num_chans) {
    sta_h[,k]=sta_h[,k]-mean(sta_h[c(erp_baseline[1]:erp_baseline[2])-tlims[1]+1,k]); #remove the baseline (the average between observation 1 to 250 of the template) from the template. These are the final templates.
    sta_f[,k]=sta_f[,k]-mean(sta_f[c(erp_baseline[1]:erp_baseline[2])-tlims[1]+1,k]);
  }
  
  train_chans_sd <- matrix(train_chans_sd,nrow = 1,ncol = num_chans)
  templates[1201*(j-1)+1,1:num_chans] <- train_chans_sd #indgest the calculated templates into the templates variable
  templates[(1201*(j-1)+2):(1201*(j-1)+601),1:num_chans] <- as.matrix(sta_h)
  templates[(1201*(j-1)+602):(1201*j),1:num_chans] <- as.matrix(sta_f)
  indx <- sapply(templates, is.factor)
  col_indx <- c(1:ncol(templates))[indx]
  templates[, col_indx[-1]] <- lapply(templates[,col_indx[-1]], function(x) as.numeric(as.character(x)))
}

col_names <- rep('',ncols-2)
col_names[1] <- 'PatientID'
for (i in 2:(ncols-2)){
  col_names[i] <- paste('Chanel ', (i-1), sep="")
}
templates <- as.matrix(templates)
templates <- data.frame(PatientID_Col, templates, stringsAsFactors=F)
colnames(templates) <- col_names #assign column names to the data frame before output
write.csv(templates, file = outputtemplatefile, row.names=FALSE) #write the template data to a local csv file
                                                                 #We will upload this file to AzureML workspace to build the predictive experiment

# Step 5.4. Project ECoG signals in both training and validation stimulus cycles to the templates to calculate the similarities, as features
# for model training and validation in next step
# Step 5.4.1 Project ECoG signals in training stimulus cycles to the templates
erp_train <- fh_project_2_templates(dataset_train, templates)
# Step 5.4.2 Project ECoG signals in validation stimulus cycles to the templates
erp_valid <- fh_project_2_templates(dataset_valid, templates)

# Step 5.5. Train a logistic regression model on training data, and valid it on validation data
col_names <- colnames(erp_train)
ncols <- length(col_names)
erp_train[,ncols-1] <- erp_train[,ncols-1] - 1 # convert label 1 and 2 to 0 and 1 for logistic regression
erp_valid[,ncols-1] <- erp_valid[,ncols-1] - 1 # convert label 1 and 2 to 0 and 1 for logistic regression
formula <- paste(col_names[2:(length(col_names)-2)], collapse="+")
formula <- paste("Stimulus_Type ~ ", formula, sep="")
glmnetmodel <- glmnet(x=as.matrix(erp_train[,2:(ncols-2)]), y=erp_train[,ncols-1], alpha=1, nlambda=1, lambda=0.01) #train a LASSO model, where lambda=0.01
summary(glmnetmodel)
#valid_pred <- predict(logitmodel, newdata = erp_valid, type="response")
valid_pred <- predict(glmnetmodel, newx = as.matrix(erp_valid[,2:(ncols-2)]), type="response")

valid_pred[valid_pred >= 0.5] <- 1
valid_pred[valid_pred < 0.5] <- 0

index <- erp_valid[,ncols-1] == valid_pred
print(paste("Validation accuracy = ", round(sum(index)/length(valid_pred)*100,4)), sep="")

# Step 6. Save the model object to a local rda file 
#save(logitmodel, file = model_rda_file) # model_rda_file is defined in the beginning of this R script file.
                                        # In this example, logitmodel.rda is the file name
save(glmnetmodel, file = model_rda_file)

