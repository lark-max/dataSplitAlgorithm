#' @title par.default
#' @description
#' The basic parameters defined by user
#'
#' @param include.inp Whether the input vector is included when calculating Euclidean distance between samples, which defaults to TRUE
#' @param seed Seed number, which defaults to 1000
#' @param sel.alg Algorithm selection, includes SOMPLEX, MDUPLEX, DUPLEX, SBSS.P, TC, which defaults to "MDUPLEX"
#' @param prop.Tr The proportion of data allocated to the training set, which defaults to 0.6
#' @param prop.Ts The proportion of data allocated to the test set, which defaults to 0.2
#' @param Train The name of the training set when the output file is required, which defaults to "Train.txt"
#' @param Test The name of the test set when the output file is required, which defaults to "Test.txt"
#' @param Validation The name of the validation set when the output file is required, which defaults to "Valid.txt"
#' @param loc.calib If using the TC algorithm, you can also specify the location of the splitting range, which defaults to c(0,0.6)
#' @param writeFile Whether to output the partition result to a txt file, which defaults to TRUE
#'
#' @return None
#'
par.default <- function(){
  list(
    include.inp = TRUE, #Whether the input vector is included when calculating Euclidean distance between samples
    seed = 1000, # Seed number
    sel.alg = "MDUPLEX", # Algorithm selection, default is MDUPLEX
    prop.Tr = 0.6, # The proportion of data allocated to the training set, which defaults to 0.6
    prop.Ts = 0.2, # The proportion of data allocated to the test set, which defaults to 0.2
    Train = "Train.txt", # The name of the training set when the output file is required, which defaults to "Train.txt"
    Test = "Test.txt", # The name of the test set when the output file is required, which defaults to "Test.txt"
    Validation = "Valid.txt", # The name of the validation set when the output file is required, which defaults to "Valid.txt"
    loc.calib = c(0,0.6), # If you use a simple continuous time series data partition, you can also specify the location of the partition segment
    writeFile = TRUE
  )
}


#' @title getAUC
#' @description
#' The similarity degree of distribution features between two data subsets are calculated and evaluated by auc index
#'
#' @param data1 A matrix or data.frame after data splitting
#' @param data2 A matrix or data.frame after data splitting
#'
#' @return a number, representing the average value of the auc
#'
#' @importFrom xgboost xgb.DMatrix xgboost
#' @importFrom pROC auc
#' @importFrom Matrix Matrix
#' @importFrom stats predict
#' @importFrom caret createFolds
#'
#' @export
#'
getAUC <- function(data1, data2){
  data <- as.data.frame(rbind(data1,data2))
  data$status <- c(rep(1,nrow(data1)),rep(0,nrow(data2)))
  data$Idex <- NULL

  splitNum <- 10
  auc_value <- c()

  # Set K folds
  folds <- createFolds(y=data$status, k=splitNum)

  # The classifier is trained for each fold and the auc value is calculated
  for(i in 1:splitNum){
    train <- data[-folds[[i]],]
    test <- data[folds[[i]],]

    # Data preprocessing
    trainData <- data.matrix(train[,-ncol(train)])
    trainData <- Matrix(trainData,sparse = T)
    dtrain <- xgb.DMatrix(data = trainData, label = train$status)
    testData <- data.matrix(test[,-ncol(test)])
    testData <- Matrix(testData,sparse = T)
    dtest <- xgb.DMatrix(data = testData, label = test$status)

    # xgboost classifier
    model <- xgboost(data = dtrain,verbosity = 0, max_depth=6, eta=0.5,nrounds=100,
                     objective='binary:logistic', eval_metric = 'auc')

    pre_xgb = round(predict(model,newdata = dtest))

    auc_value = c(auc_value, as.numeric(auc(test$status,pre_xgb,levels=c(1,0),direction=">")))
  }
  return(mean(auc_value))
}


getDist <- function(vect1, vect2) sqrt(sum((vect1 - vect2)^2))


#' @title selectData
#' @description
#' Built-in function: select data for splitting and decide whether to include input vectors based on
#' user-set parameter \code{include.inp}
#'
#' @param data Raw data set entered by the user in matrix or data.frame format
#' @param control The list of parameters given by the user
#'
#' @return Remove the subscript vector in the first column and return a matrix
#'
selectData <- function(data,control){
  if(control)
    return.data = as.matrix(data[,-1]) # Remove the index column
  else
    return.data = as.matrix(data[,ncol(data)]) # Only choose the Output column
  return(return.data)
}


#' @title standardise
#' @description
#' Build-in function: using for standardizing the raw data
#'
#' @param data Enter data in matrix or data.frame format
#'
#' @return A matrix with the internal data normalized by column
#'
standardise <- function(data){
  if(!is.matrix(data) && !is.data.frame(data)){
    stop(
      "[Error]:Invalid input data format, should be matrix or dataframe"
    )
  }
  if(ncol(data) < 2 || nrow(data) < 2)
    return(as.matrix(data))
  else
    return(as.matrix(apply(data,2,scale)))
}


#' @title DP.initialSample
#' @description
#' Built-in function: the initial sampling of DUPLEX algorithm
#'
#' @param split.info A list containing information such as data subscripts and classification subsets that need to be divided
#' @param choice The current sample set, a string, must have the same name as the subset in split.info
#'
#' @return A list to complete the current sampling
#'
DP.initialSample <- function(split.info, choice){
  if(length(split.info$index) < 2)
    return(split.info)

  a = 1; b = 2 # Set two points
  maxDistance = 0
  for(t in 1:(length(split.info$index)-1)){
    for(tt in (t+1):length(split.info$index)){
      cur.dist =  getDist(split.info$data[split.info$index[t],],
                          split.info$data[split.info$index[tt],])
      if(cur.dist > maxDistance){
        maxDistance = cur.dist
        a = t
        b = tt
      }
    }
  }
  eval(parse(text = paste("split.info$",choice,"=c(split.info$",choice,
                          ",split.info$index[c(a,b)])",sep="")))
  eval(parse(text = paste("split.info$index = split.info$index[c(-a,-b)]")))
  return(split.info)
}


#' @title DP.resample
#' @description
#' Built-in function: cyclic sampling of DUPLEX algorithm
#'
#' @param split.info A list containing information such as data subscripts and classification subsets that need to be divided
#' @param choice The current sample set, a string, must have the same name as the subset in split.info
#'
#' @return A list to complete the current sampling
#'
DP.resample <- function(split.info, choice){
  if(length(split.info$index) < 2){
    eval(parse(text = paste("split.info$",choice,"=c(split.info$",choice,
                            ",split.info$index)",sep="")))
    split.info$index <- c()
    return(split.info)
  }
  a = 1; b = 2
  distA = 0; distB = 0
  eval(parse(text = paste("curSetLen=length(split.info$",choice,")",sep="")))
  for(t in 1:length(split.info$index)){
    threshold.dist = 1.0e+6
    for(s in 1:curSetLen){
      eval(parse(text = paste(
        "cur.dist = getDist(split.info$data[split.info$index[t],],split.info$data[split.info$",
        choice,"[s],])",sep="")))
      if(cur.dist < threshold.dist)
        threshold.dist = cur.dist
    }

    if(threshold.dist >= distA){
      b = a
      distB = distA
      a = t
      distA = threshold.dist
    }else if(threshold.dist >= distB){
      b = t
      distB = threshold.dist
    }else
      next
  }
  eval(parse(text = paste("split.info$",choice,"=c(split.info$",choice,
                          ",split.info$index[c(a,b)])",sep="")))
  eval(parse(text = paste("split.info$index = split.info$index[c(-a,-b)]")))
  return(split.info)
}

#' @title checkFull
#' @description
#' Built-in function: check whether sampling of any two subsets has been completed,
#' and all remaining data is allocated to the last subset
#'
#' @param split.info A list containing information such as data subscripts and classification subsets that need to be divided
#' @param num.train The sampling set of the training set
#' @param num.test The sampling set of the test set
#' @param num.valid The sampling set of the validation set
#'
#' @return A list
#'
checkFull <- function(split.info, num.train, num.test, num.valid){
  res <- list(ini.info = split.info)
  if(!(length(split.info$trainKey) < num.train) & !(length(split.info$testKey) < num.test)){
    res$ini.info$validKey = c(res$ini.info$validKey, res$ini.info$index)
    res$ini.info$index = c()
    res$signal = TRUE
  }else if(!(length(split.info$trainKey) < num.train) & !(length(split.info$validKey) < num.valid)){
    res$ini.info$testKey = c(res$ini.info$testKey, res$ini.info$index)
    res$ini.info$index = c()
    res$signal = TRUE
  }else if(!(length(split.info$testKey) < num.test) & !(length(split.info$validKey) < num.valid)){
    res$ini.info$trainKey = c(res$ini.info$trainKey, res$ini.info$index)
    res$ini.info$index = c()
    res$signal = TRUE
  }else{
    res$signal = FALSE
  }
  return(res)
}


#' @title somCluster
#' @description
#' Built-in function: SOM cluster by using \code{[kohonen]{som}}
#'
#' @param data Given data set in matrix or data.frame format
#'
#' @return A List which contains data points clustered by each neuron in the SOM network
#' @importFrom kohonen som somgrid
#'
somCluster <- function(data){
  # Basic parameters for som cluster
  neuron.num = round(2*sqrt(nrow(data)))
  neuron.col = round(sqrt(neuron.num / 1.6))
  neuron.row = round(1.6 * neuron.col)
  neuron.info <- list(
    neuron.num = neuron.row * neuron.col,
    neuron.row = neuron.row,
    neuron.col = neuron.col
  )
  # som train, need loading the package: kohonen
  som.model <- som(data, grid = somgrid(neuron.info$neuron.row, neuron.info$neuron.col,
                                        "hexagonal"))
  # summary
  neuron.info$neuron.cluster <- list()
  for(i in 1:neuron.info$neuron.num){
    neuron.info$neuron.cluster[[i]] = which(som.model$unit.classif == i)
  }
  return(neuron.info)
}


#' @title getSampleNumofEachNeuron
#' @description
#' Built-in function: calculate the number of data segmentation for each neuron in the SOM grid
#'
#' @param som.info A List which contains data points clustered by each neuron in the SOM network
#' @param control The list of parameters given by the user
#'
#' @return A list containing the sampling number of three subsets in the SOM grid
#'
getSampleNumofEachNeuron <- function(som.info, control){
  sampleNum.eachNeuron = list()
  sampleNum.eachNeuron$Tr <- sampleNum.eachNeuron$Ts <- sampleNum.eachNeuron$Vd <- c()
  for(i in 1:som.info$neuron.num){
    sampleNum.eachNeuron$Tr[i] = round(length(som.info$neuron.cluster[[i]])* control$prop.Tr)
    sampleNum.eachNeuron$Ts[i] = round(length(som.info$neuron.cluster[[i]])* control$prop.Ts)
    sampleNum.eachNeuron$Vd[i] = length(som.info$neuron.cluster[[i]]) -
      sampleNum.eachNeuron$Tr[i] - sampleNum.eachNeuron$Ts[i]
  }
  return(sampleNum.eachNeuron)
}


#' @title TC algorithm
#' @description
#' Traditional data splitting algorithm: simply divide the raw data into two parts for calibration and evaluation
#'
#' @param data Raw data set entered by the user in matrix or data.frame format
#' @param control The list of parameters given by the user
#'
#' @return A list containing subsets of calibration and evaluation data
#'
TC <- function(data, control){
  len <- nrow(data)
  start.point = 1 + len * control$loc.calib[1]
  end.point = len * control$loc.calib[2]

  if(start.point<0 | start.point > len | end.point < 0 | end.point > len | end.point < start.point)
    stop(
      "[Error]:Parameter loc.calib is incorrect, should be set to [0,1]"
    )

  return(list(Calibration = data[start.point:end.point,],
              Evaluation = data[-(start.point:end.point),]))
}


#' @title SSsample
#' @description
#' Built-in function: core code of SS algorithm
#'
#' @param index A vector containing the index of the data
#' @param prop A number determining the sampling step
#'
#' @importFrom stats runif
#' @return A vector with sampled data index
#'
SSsample <- function(index, prop){
  sampleVec <- c()
  interval <- ceiling(1/prop)
  loc <- runif(1,0,32767) %% interval
  while(loc <= length(index)){
    k <- ceiling(loc)
    sampleVec <- c(sampleVec, index[k])
    loc = loc + interval
  }
  return(sampleVec)
}


#' @title remainUnsample
#' @description
#' Built-in function: the auxiliary function of SS algorithm is used to obtain the remaining unsampled data
#'
#' @param X A vector that needs to be sampled
#' @param Y A vector that have sampled
#'
#' @return A vector containing the remaining unsampled data
#'
remainUnsample <- function(X, Y){
  remainData <- c()
  ix <- iy <- 1
  while(ix <= length(X)){
    if(iy >= length(Y))
      break
    if(X[ix] == Y[iy])
      iy = iy + 1
    else
      remainData  = c(remainData, X[ix])
    ix = ix + 1
  }
  remainData = c(remainData, X[ix:length(X)])
  return(remainData)
}


#' @title SS algorithm
#' @description
#' The Systematic stratified(SS) data splitting method is a semideterministic method,
#' in which every kth sample from a random starting point is selected to form the training,
#' testing, and validation data sets. In implementing systematic sampling in this study,
#' the data are first ordered along the output variable dimension in increasing order.
#' Then, the sampling interval is determined based on the training and testing data proportions specified by the user.
#' Thereafter, a starting point is randomly selected and training samples are drawn first, followed by the testing samples.
#' Finally, unsampled data are allocated to the validation set.
#'
#' @references
#' Baxter, C. W., Stanley S. J., Zhang Q., and Smith D. W.(2000), Developing artificial neural network process models: A guide for drinking water utilities,
#' @references
#' paper presented at the 6th Environmental Engineering Society Specialty Conference of the CSCE, pp. 376–383, Canadian Society for Civil Engineering (CSCE),
#' @references
#' London, Ontario, Canada.
#'
#' @param data Raw data set entered by the user in matrix or data.frame format
#' @param control The list of parameters given by the user
#'
#' @return A list of three subsets
#'
SS <- function(data, control){
  # Basic parameters
  set.seed(control$seed)
  num.total = nrow(data)
  num.train = round(num.total*control$prop.Tr)
  num.test = round(num.total*control$prop.Ts)
  num.valid = num.total- (num.train + num.test)

  # deal with data
  select.data = selectData(data, control$include.inp)
  select.data.std = standardise(select.data)

  data.index <- seq(1,num.total,1)

  # Firstly get the output variable list
  outputVec <- select.data.std[,ncol(select.data.std)]

  # The data are first ordered along the output variable dimension in increasing order
  data.index <- data.index[order(outputVec),drop = FALSE]

  calibrateProp = control$prop.Tr + control$prop.Ts
  trainKey <- testKey <- validKey <- c()

  # The data is firstly divided into two parts,
  # which need to be sampled and those that do not need to be sampled
  sampleKey <- SSsample(data.index,calibrateProp)

  # The unsampled data are allocated to validating subset
  validKey <- remainUnsample(data.index, sampleKey)

  # Then split sample into systematic testing and training sets
  testKey <- SSsample(sampleKey, control$prop.Ts / calibrateProp)
  trainKey <- remainUnsample(sampleKey, testKey)

  print("SS sampling complete!")
  return(list(Train = data[trainKey,], Test = data[testKey,],
              Validation = data[validKey,]))
}


#' @title SBSS-P algorithm
#' @description
#' The self-organizing map based stratified sampling(SBSS) approach is a two-step data splitting method,
#' which has been implemented previously by (Bowden et al., 2002), (Kingston, 2006), and (May et al., 2010).
#' In the first step, multivariate stratified random sampling is performed to partition the data into K strata, where clustering is performed using a self-organizing map (SOM,Kohonen, 1998).
#' In the second step, uniform random intracluster sampling is applied to generate the data split.
#'
#' @references
#' May, R. J., Maier H. R., and Dandy G. C.(2010), Data splitting for artificial neural networks using SOM-based stratified sampling, Neural Netw, 23(2), 283-294.

#' @param data Raw data set entered by the user in matrix or data.frame format
#' @param control The list of parameters given by the user
#'
#' @return A list of three subsets
#'
SBSS.P <- function(data, control){
  # Basic parameters
  set.seed(control$seed)
  num.total = nrow(data)
  num.train = round(num.total*control$prop.Tr)
  num.test = round(num.total*control$prop.Ts)
  num.valid = num.total- (num.train + num.test)

  # deal with data
  select.data = selectData(data, control$include.inp)
  select.data.std = standardise(select.data)

  # SOM cluster
  som.info <- somCluster(select.data.std)

  # Determine the amount of sampled data for each neuron
  sampleNum.eachNeuron <- getSampleNumofEachNeuron(som.info, control)

  # Sampling
  trainSet <- testSet <- ValidSet <- c()
  print(paste("Total neuron:", som.info$neuron.num))

  # Create information list
  split.info <- list(
    trainKey = c(),
    testKey = c(),
    validKey = c()
  )

  # Sampling from each neuron
  for(i in 1:som.info$neuron.num){
    print(paste("sampling on neuron:",i))
    # Sampling for training dataset
    if(length(som.info$neuron.cluster[[i]]) > sampleNum.eachNeuron$Tr[i]){
      randomSample.index <- sample(c(1:length(som.info$neuron.cluster[[i]])), sampleNum.eachNeuron$Tr[i])
      split.info$trainKey <- c(split.info$trainKey, som.info$neuron.cluster[[i]][randomSample.index])
      som.info$neuron.cluster[[i]] = som.info$neuron.cluster[[i]][-randomSample.index]
    }
    # Sampling for test dataset
    if(length(som.info$neuron.cluster[[i]]) > sampleNum.eachNeuron$Ts[i]){
      randomSample.index <- sample(c(1:length(som.info$neuron.cluster[[i]])), sampleNum.eachNeuron$Ts[i])
      split.info$testKey <- c(split.info$testKey, som.info$neuron.cluster[[i]][randomSample.index])
      som.info$neuron.cluster[[i]] = som.info$neuron.cluster[[i]][-randomSample.index]
    }
    # The remain data points are all allocated to validation dataset
    split.info$validKey <- c(split.info$validKey, som.info$neuron.cluster[[i]])
    som.info$neuron.cluster[[i]] <- NA
  }
  print("SBSS-P sampling complete!")
  return(list(Train = data[split.info$trainKey,], Test = data[split.info$testKey,],
              Validation = data[split.info$validKey,]))

}


#' @title DUPLEX algorithm
#' @description
#' The DUPLEX data splitting method was developed by (Snee,1977) based on one of the earliest data splitting algorithms called CADEX or Kennard-Stone sampling (Kennard and Stone, 1969).
#' DUPLEX draws samples based on Euclidean distances. When applying DUPLEX, the two points which are farthest apart in terms of the Euclidean distance are assigned to the first data set.
#' The next pair of points that are farthest apart in the remaining list are assigned to the second data set. This process is repeated until both data sets are filled (Snee, 1977).
#' The original DUPLEX algorithm was used to divide data into two sets. (May et al., 2010) modified the original DUPLEX algorithm to sample data into three data sets in turn based on the proportions of the data sets specified by the user.
#' Thus, DUPLEX can be used to generate the training, testing, and validation data sets for ANN model development.
#'
#' @references
#' Snee, R. D. (1977), Validation of regression models: Methods and examples, Technometrics, 19(4), 415–428.
#'
#' @param data Raw data set entered by the user in matrix or data.frame format
#' @param control The list of parameters given by the user
#'
#' @return A list of three subsets
#'
DUPLEX <- function(data,control){
  # Basic parameters
  set.seed(control$seed)
  num.total = nrow(data)
  num.train = round(num.total*control$prop.Tr)
  num.test = round(num.total*control$prop.Ts)
  num.valid = num.total- (num.train + num.test)

  # deal with data
  select.data = selectData(data, control$include.inp)
  select.data.std = standardise(select.data)

  # Create information list
  split.info <- list(
    data = select.data.std,
    index = seq(1,num.total,1),
    trainKey = c(),
    testKey = c(),
    validKey = c()
  )

  # Step1: initial sampling
  print("Start the initial sampling...")
  if(num.train > 0)
    split.info = DP.initialSample(split.info,"trainKey")
  if(num.test > 0)
    split.info = DP.initialSample(split.info,"testKey")
  if(num.valid > 0)
    split.info = DP.initialSample(split.info,"validKey")
  print("Initial sampling successfully!")

  # Step2: sampling data through a cyclic sampling pool
  print("Start the loop sampling...")

  while(length(split.info$index) > 0){
    print(paste("Remaining unsampled data:",length(split.info$index)))
    if(length(split.info$trainKey) < num.train)
      split.info = DP.resample(split.info,"trainKey")
    if(length(split.info$testKey) < num.test)
      split.info = DP.resample(split.info,"testKey")
    if(length(split.info$validKey) < num.valid)
      split.info = DP.resample(split.info,"validKey")
    # Check full
    check.res = checkFull(split.info, num.train, num.test, num.valid)
    if(check.res$signal){ # The stop signal is TRUE
      split.info = check.res$ini.info
      print("Two of the datasets are full, and all remaining data is sampled to the other dataset")
      break;
    }
  }

  print("DUPLEX sampling complete!")
  return(list(Train = data[split.info$trainKey,], Test = data[split.info$testKey,],
              Validation = data[split.info$validKey,]))


}


#' @title Modified DUPLEX -- MDUPLEX algorithm
#' @description
#' DUPLEX have provided a mechanism for splitting the data in a deterministic manner,
#' a property that is attractive for practical applications. However,it can result in significant bias when used in the context of evaluating the performance of data-driven rainfall-runoff models where the size of the calibration subset is typically larger than that of the evaluation subset. Under such circumstances,
#' DUPLEX allocates a substantially larger number of normal (less extreme) data points to the larger (calibration) subset. This biases the calibration towards normal events and causes the model to have relatively poor performance on extreme events, resulting in pessimistic assessment of model evaluation performance
#' To overcome this problem, the authors of this package have developed an improved version of MDUPLEX algorithm, which can effectively solve the defects of traditional DUPLEX algorithm.
#' The details of the algorithm improvement are described in:
#'
#' @references
#' Chen, J., Zheng F., May R., Guo D., Gupta H., and Maier H. R.(2022), Improved data splitting methods for data-driven hydrological model development based on a large number of catchment samples, Journal of Hydrology, 613.
#' @references
#' Zheng, F., Chen J., MaierH. R., and Gupta H.(2022), Achieving Robust and Transferable Performance for Conservation‐Based Models of Dynamical Physical Systems, Water Resources Research, 58(5).
#'
#' @param data Raw data set entered by the user in matrix or data.frame format
#' @param control The list of parameters given by the user
#'
#' @return A list of three subsets
#'
MDUPLEX <- function(data,control){
  # Basic parameters
  set.seed(control$seed)
  num.total = nrow(data)
  num.train = round(num.total*control$prop.Tr)
  num.test = round(num.total*control$prop.Ts)
  num.valid = num.total- (num.train + num.test)

  # deal with data
  select.data = selectData(data, control$include.inp)
  select.data.std = standardise(select.data)

  # Parameters of the basic sampling pool
  poolSize = round(1/min(control$prop.Tr, control$prop.Ts,
                         1-(control$prop.Tr+control$prop.Ts)))
  samplingPool <- list(
    trainSize = round(poolSize * control$prop.Tr),
    testSize = round(poolSize * control$prop.Ts),
    validSize = round(poolSize * (1-control$prop.Tr-control$prop.Ts))
  )

  # Create information list
  split.info <- list(
    data = select.data.std,
    index = seq(1,num.total,1),
    trainKey = c(),
    testKey = c(),
    validKey = c()
  )

  # Step1: initial sampling
  print("Start the initial sampling...")
  if(num.train > 0)
    split.info = DP.initialSample(split.info,"trainKey")
  if(num.test > 0)
    split.info = DP.initialSample(split.info,"testKey")
  if(num.valid > 0)
    split.info = DP.initialSample(split.info,"validKey")
  print("Initial sampling successfully!")

  # Step2: sampling data through a cyclic sampling pool
  print("Start the loop sampling...")
  while (length(split.info$index)>0) {
    trainSize.cnt = samplingPool$trainSize
    testSize.cnt = samplingPool$testSize
    validSize.cnt = samplingPool$validSize
    print(paste("Remaining unsampled data:",length(split.info$index)))
    while(TRUE){
      stopSignal = TRUE
      if(trainSize.cnt!=0 & length(split.info$trainKey) < num.train){
        split.info = DP.resample(split.info,"trainKey")
        trainSize.cnt = trainSize.cnt - 1
        stopSignal = FALSE
      }
      if(testSize.cnt!=0 & length(split.info$testKey) < num.test){
        split.info = DP.resample(split.info,"testKey")
        testSize.cnt = testSize.cnt - 1
        stopSignal = FALSE
      }
      if(validSize.cnt!=0 & length(split.info$validKey) < num.valid){
        split.info = DP.resample(split.info,"validKey")
        validSize.cnt = validSize.cnt - 1
        stopSignal = FALSE
      }
      if(stopSignal)
        break
    }
  }
  print("MDUPLEX sampling complete!")
  return(list(Train = data[split.info$trainKey,], Test = data[split.info$testKey,],
              Validation = data[split.info$validKey,]))
}


#' @title SOMPLEX algorithm
#' @description
#' SOMPLEX is a multi-stage algorithm that combines SOM-based clustering with DUPLEX sampling of map units.
#' SOM clustering is used to identify distinct regions of the data space that can include both typical classes and examples of unique input–output cases.
#' However, overlooked by many data splitting applications is that the clusters may cover varying proportions of the database, data can be non-uniformly distributed within each partition,
#' and data in some partitions can be more widely spread than in others (May et al., 2010), due to which random sampling from within each cluster can be less than ideal.
#' Conversely, given any distribution of data, DUPLEX can generate a sample that uniformly covers the full range of the data (Snee, 1977).
#' By applying DUPLEX to each SOM-based partition, representative samples can be generated from each partition regardless of the distribution within the partition.
#'
#' @references
#' Chen, J., Zheng F., May R., Guo D., Gupta H., and Maier H. R.(2022), Improved data splitting methods for data-driven hydrological model development based on a large number of catchment samples, Journal of Hydrology, 613.
#'
#'
#' @param data Raw data set entered by the user in matrix or data.frame format
#' @param control The list of parameters given by the user
#'
#' @return A list of three subsets
#'
SOMPLEX <- function(data, control){
  # Basic parameters
  set.seed(control$seed)
  num.total = nrow(data)
  num.train = round(num.total*control$prop.Tr)
  num.test = round(num.total*control$prop.Ts)
  num.valid = num.total- (num.train + num.test)

  select.data = selectData(data, control$include.inp)
  select.data.std = standardise(select.data)

  # SOM cluster
  som.info <- somCluster(select.data.std)

  # Determine the amount of sampled data for each neuron
  sampleNum.eachNeuron <- getSampleNumofEachNeuron(som.info, control)

  # Sampling
  trainSet <- testSet <- ValidSet <- c()
  print(paste("Total neuron:", som.info$neuron.num))
  for(i in 1:som.info$neuron.num){
    # Create information list
    split.info <- list(
      data = select.data.std,
      index = som.info$neuron.cluster[[i]],
      trainKey = c(),
      testKey = c(),
      validKey = c()
    )
    print(paste("sampling on neuron:",i))
    if(sampleNum.eachNeuron$Tr[i] > 0)
      split.info = DP.initialSample(split.info,"trainKey")
    if(sampleNum.eachNeuron$Ts[i] > 0)
      split.info = DP.initialSample(split.info,"testKey")
    if(sampleNum.eachNeuron$Vd[i] > 0)
      split.info = DP.initialSample(split.info,"validKey")

    while(length(split.info$index) > 0){
      if(length(split.info$trainKey) < sampleNum.eachNeuron$Tr[i])
        split.info = DP.resample(split.info,"trainKey")
      if(length(split.info$testKey) < sampleNum.eachNeuron$Ts[i])
        split.info = DP.resample(split.info,"testKey")
      if(length(split.info$validKey) < sampleNum.eachNeuron$Vd[i])
        split.info = DP.resample(split.info,"validKey")
    }
    trainSet = c(trainSet, split.info$trainKey)
    testSet = c(testSet, split.info$testKey)
    ValidSet = c(ValidSet, split.info$validKey)
  }
  print("SOMPLEX sampling complete!")
  return(list(Train = data[trainSet,], Test = data[testSet,], Validation = data[ValidSet,]))
}


#' @title dataSplit
#' @description
#' The main program entry of data splitting algorithms, provides the interface of each algorithm
#' Data format requirements:
#' Refer to the data sample provided by the package:The first column should be the index, the second through Nth columns should be the input vector,
#' and the last column should be the output vector
#' Refer to \code{\link{par.default}} for user-defined parameters.
#'
#' @param data Raw data set entered by the user in matrix or data.frame format
#' @param control The list of parameters given by the user
#' @param ... Redundant parameter list
#'
#' @importFrom utils modifyList write.table
#' @return A list of three subsets
#' @export
#'
#' @author
#' Feifei Zheng \email{feifeizheng@zju.edu.cn}
#' @author
#' Junyi Chen \email{jun1chen@zju.edu.cn}
#'
#' @references
#' Chen, J., Zheng F., May R., Guo D., Gupta H., and Maier H. R.(2022).Improved data splitting methods for data-driven hydrological model development based on a large number of catchment samples, Journal of Hydrology, 613.
#' @references
#' Zheng, F., Chen J., Maier H. R., and Gupta H.(2022). Achieving Robust and Transferable Performance for Conservation‐Based Models of Dynamical Physical Systems, Water Resources Research, 58(5).
#' @references
#' Zheng, F., Chen, J., Ma, Y.,  Chen Q., Maier H. R., and Gupta H.(2023). A Robust Strategy to Account for Data Sampling Variability in the Development of Hydrological Models, Water Resources Research, 59(3).
#'
#' @examples
#' data("dataSplitAlgorithm_data_small")
#' result = dataSplit(dataSplitAlgorithm_data_small, control = list(sel.alg = "MDUPLEX",writeFile = FALSE))
#'
#' data("dataSplitAlgorithm_data_middle")
#' result = dataSplit(dataSplitAlgorithm_data_middle, control = list(sel.alg = "SBSS.P",writeFile = FALSE))
#'
#' data("dataSplitAlgorithm_data_large")
#' result = dataSplit(dataSplitAlgorithm_data_large, control = list(sel.alg = "SOMPLEX",writeFile = FALSE))
#'
dataSplit <- function(data,control = list(),...){
  # Check data format
  if(!is.matrix(data) & !is.data.frame(data)){
    stop(
      "[Error]:Invalid input data format!"
    )
  }
  # Check data integrity
  if(sum(is.na(data))!=0){
    stop(
      "[Error]:Missing values in the input data"
    )
  }
  # Check parameter list
  stopifnot(is.list(control))
  control <- modifyList(par.default(),control)
  isValid <- names(control) %in% names(par.default())
  if (any(!isValid)) {
    stop(
      "[Error]:Unrecognised options: ",
      toString(names(control)[!isValid])
    )
  }

  # Split data
  start.time = Sys.time()
  eval(parse(text = paste("obj=",control$sel.alg,"(data,control)",sep="")))
  obj$func.time = Sys.time() - start.time

  # Output data into txt files
  if(control$writeFile){
    if(control$sel.alg != "TC"){
      write.table(obj$Train,control$Train,row.names = F,col.names = T,sep='\t')
      write.table(obj$Test,control$Test,row.names = F,col.names = T,sep='\t')
      write.table(obj$Validation,control$Validation,row.names = F,col.names = T,sep='\t')
    }else{
      write.table(obj$Calibration,"Calibration.txt",row.names = F,col.names = T,sep='\t')
      write.table(obj$Evaluation,"Evaluation.txt",row.names = F,col.names = T,sep='\t')
    }
  }
  return(obj)
}
