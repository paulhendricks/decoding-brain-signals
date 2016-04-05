library(abind)

# Define function fh_get_events() that will be widely used in the pipeline.
# It outputs events start time, midway from event start to the next event, and event type
# Event type has 3 values: 0 means it is the inter-stimulus interval start time, 
#                          1 means it is the house stimulus onset time
#                          2 means it is the face stimulus onset time
# It takes the stimulus column as the input, which records the stimulus type at each sampling point
fh_get_events = function(stim) {
	nrows <- nrow(stim)
	if (stim[nrows] != 0 & stim[nrows] != 101) {
		#if the last stimulus type is not 0 or 101, 
		#meaning that the stim column ends at the end of a house or a face stimulus presentation cycle. 
		#Padding a 0 after that, otherwise, the last stimulus onset event will not be recoganized. 
		stim <- rbind(stim, 0)
	}
	tmp = c(0, stim[1:((length(stim) - 1))])
	b <- which((stim - tmp) != 0) #Get the stimulus type changing time
	c <- floor(diff(b) / 2) #Get the half of the times between two consecutive events
	b <- b[1:length(b) - 1]
	d <- b + c # Get the midway between two consecutive events
	evs = matrix(data = NA, nrow = length(b), ncol = 3)
	evs[, 1] <- b
	evs[, 2] <- d
	evs[, 3] <- stim[d]
	evs <- evs[which(evs[, 3] != 0),] # If stimulus type is 0, it is not ISI, nor stimulus presentation time. Get rid of these events
	evs[which(evs[, 3] < 51), 3] <- 1 # Stimulus types 1 - 50 means it is house stimulus (stimulus class 1)
	evs[which(evs[, 3] == 101), 3] <- 0 # If original stimulus type is 101, meaning that it is ISI, relabel it as event 0
	evs[which(evs[, 3] > 50), 3] <- 2 # Stimulus types 51 - 100 means it is house stimulus (stimulus class 1)
	rm(b, c, d)
	return(evs)
	}


# Function that projects raw ECoG signals to templates.
# dataset1 is the raw ECoG signals, and dataset2 is the templates
fh_project_2_templates = function(dataset1, dataset2) {
  ncols1 <- ncol(dataset1)
  ncols2 <- ncol(dataset2)
  unique_patients <- unique(dataset1[,1]) # get the list of unique patients
  num_patients <- length(unique_patients) # get the number of unique patients
  train_data <- NULL
  Patient_IDs <- NULL
  Labels <- NULL
  Stimulus_ID <- NULL
  for (j in 1:num_patients){
    dataset2_j <- dataset2[dataset2[,1]==unique_patients[j],2:ncols2] # get the data for this specific patient
    dataset2_j <- apply(dataset2_j, 2, as.numeric)
    nonna_index <- !is.na(dataset2_j[1,])  # in the template data, if a channel does not have value for a patient, this column is missing. is.na() will be true on this column
    ncols2_j <- sum(nonna_index) #number of columns, not including the patientID column
    dataset1_j <- dataset1[dataset1[,1]==unique_patients[j],2:(ncols2_j+1)] # get the signals of this patient
    dataset1_j_stimulus_id <- dataset1[dataset1[,1]==unique_patients[j],ncols1] # get the stimulus presentation cycle ID of this patient
    train_chan_sd <- dataset2_j[1,1:ncols2_j] # for each patient, the first row in the template data is the standard deviation 
    sta_h <- dataset2_j[2:601,1:ncols2_j] # rows 2 to 601 are the house templates
    sta_f <- dataset2_j[602:1201,1:ncols2_j] # rows 602-1201 are the face templates
    signal_train <- as.matrix(dataset1_j)
    signal_train <- apply(signal_train, 2, as.numeric) # convert the signal to numeric
    stim_train <- as.matrix(dataset1[dataset1[,1]==unique_patients[j],ncols1-1]) # get the stimulus type column of this patient
    events_train <- fh_get_events(stim_train) # get the event onset time
    events_train <- events_train[which(events_train[,3]!=0),]
    events_train <- events_train[,-2]
    
    train_t=c(1:nrow(stim_train));
    train_e=events_train;
    num_stimulus <- nrow(train_e)
    Stimulus_ID_j <- matrix(0,nrow=num_stimulus, ncol=1)
    
    num_chans=ncols2_j;
    
    tlims=c(-200,399);# times to start and end erps
    erp_baseline=c(-200,49);# times to calcualte erp based upon(must be within tlims)
    
    for (k in 1:num_chans){
      signal_train[,k]=signal_train[,k]/train_chan_sd[k]; # Normalize the signals by dividing by the stand deviation in the training data
    }
    
    ## generate train data features
    # Originally, the author of the PLOS paper generates 4 training points between each stimulus for each stimulus presentation cycle. 
    # The purpose of doing so is to select features for machine learning models. 
    # That part is skipped in this sample training experiment. The feature selection part is left 
    # for participants to figure out
    
    f_template_train=matrix(0, nrow=nrow(train_e), ncol=num_chans);
    h_template_train=0*f_template_train;
    
    for (k in 1:nrow(train_e)){ # for each stimulus presentation cycle, calculate the similarity between the signal [stimulus onset time - 200, stimulus onset time + 399] and the templates
      Stimulus_ID_j[k,1] <- dataset1_j_stimulus_id[train_e[k,1]] # get the stimulus presentation cycle ID at the onset time of the stimulus
      # dot products to project the raw signal to templates, as similarity measurement. Very similar as calculating the cosine between two vectors
      for (chan in 1:num_chans){
        dt=signal_train[train_e[k,1]+c(tlims[1]:tlims[2]),chan];# select data for that channel
        dt=dt-mean(dt[c(erp_baseline[1]:erp_baseline[2])-tlims[1]+1]);#normalize by substracting the baseline mean of that signal
        
        f_template_train[k,chan]=sum(sta_f[,chan]*dt);# project to the template of face for that channel
        h_template_train[k,chan]=sum(sta_h[,chan]*dt);# project to the template of house for that channel
      }
    }
    
    traindata_f=f_template_train
    traindata_h=h_template_train
    trainlabels=train_e[,2]
    
    trainingdata=abind(traindata_f[which(trainlabels>0),],traindata_h[which(trainlabels>0),],along=2) #consolidate the face similarities and house similarities together into a single dataset
    training_label= as.matrix(trainlabels[which(trainlabels>0)]);
    nrows_j <- nrow(trainingdata)
    trainingdata_j <- matrix(NA, nrow=nrows_j, ncol=2*(ncols1-3))
    trainingdata_j[,1:ncols2_j] <- trainingdata[,1:ncols2_j]
    trainingdata_j[,(ncols1-3+1):(ncols1-3+ncols2_j)] <- trainingdata[,(ncols2_j+1):(2*ncols2_j)]
    train_data <- rbind(train_data,trainingdata_j)
    Patient_IDs <- rbind(Patient_IDs, matrix(unique_patients[j], nrow=nrows_j, ncol=1))
    Labels <- rbind(Labels, training_label)
    Stimulus_ID <- rbind(Stimulus_ID, Stimulus_ID_j)
  }
  nrows <- nrow(train_data)
  Patient_IDs <- as.matrix(Patient_IDs, nrow=nrows, ncol=1)
  train_data <- as.matrix(train_data)
  data.set <- cbind(Patient_IDs, train_data, Labels, Stimulus_ID)
  col_names <- rep('',(ncols1-3)*2+3)
  col_names[1] <- 'PatientID'
  for (i in 1:(ncols1-3)){
    col_names[i+1] <- paste('Chanel_',i,'_F', sep='') # the first 64 variables will be similarity features with face templates, and the next 64 variables will be with house templates
    col_names[i+1+ncols1-3] <- paste('Chanel_',i,'_H',sep='')
  }
  col_names[2*(ncols1-3)+2] <- 'Stimulus_Type'
  col_names[2*(ncols1-3)+3] <- 'Stimulus_ID'
  
  tmp <- apply(data.set[,2:ncol(data.set)], 2, as.numeric)
  tmp[is.na(tmp)] <- 0
  data.set <- data.frame(data.set[,1], tmp, stringsAsFactors = F)
  #data.set <- as.data.frame(data.set)
  colnames(data.set) <- col_names
  indx <- sapply(data.set, is.factor)
  col_indx <- c(1:length(col_names))[indx]
  data.set[, col_indx[-1]] <- lapply(data.set[,col_indx[-1]], function(x) as.numeric(as.character(x))) # convert numeric variables to numeric, not factors
  
  return(data.set)
}


