#Copyright (C) 2018 Bing Zhu
# ===================================================
# SmoteENN: Smote+ENN
# ===================================================

SmoteENN<-
    function(x, y, percOver = 1400, k1 = 5, k2 = 3)
        # INPUTS
        #    x: A data frame of the predictors from training data
        #    y: A vector of response variable from training data
        #    percOver: Percent of new instance generated for each minority instance
        #    k1: Number of the nearest neighbors
        #    k2: Number of neighbours for ENN
        #  allowParallel: A logical number to control the parallel computing. If allowParallel = TRUE, the function is run using parallel techniques
    {
        source("packages/IRIC-master/R/Data level/SMOTE.R")
        newData <- SMOTE(x, y, percOver, k1)
        tgt <- length(newData)
        
        # remove option for parallel computing -- requires a library not available on hpc
        # indexENN  <- ENN(tgt, newData, k2,allowParallel)
        
        indexENN  <- ENN(tgt, newData, k2)  
        newDataRemoved <- newData[!indexENN, ]
        return(newDataRemoved)
    }


# ===================================================
#  ENN: using ENN rule to find the noisy instances
# ===================================================

ENN <-
    function(tgt, data, k)
    {
        # find column of the target 
        numRow  <- dim(data)[1]
        indexENN <- rep(FALSE, numRow)
        
        # transform the nominal data into  binary 
        source("packages/IRIC-master/R/Data level/Numeralize.R")
        dataTransformed <- Numeralize(data[, -tgt])
        classMode<-matrix(nrow=numRow)
        library("RANN")
        indexOrder <- nn2(dataTransformed, dataTransformed, k+1)$nn.idx
        
        # remove option for parallel computing
        # if  (allowParallel) {
        # 
        #     classMetrix <- matrix(data[indexOrder[,2:(k+1)], tgt], nrow = numRow)
        #     library("parallel")
        #     cl <- makeCluster(2)
        #     classTable   <- parApply (cl, classMetrix, 1, table)
        #     modeColumn   <- parLapply(cl, classTable, which.max)
        #     classMode    <- parSapply(cl, modeColumn, names)
        #     stopCluster(cl)
        #     indexENN[data[, tgt]!= classMode] <- TRUE
        # } else {
            
            for (i in 1:numRow)
            {
                classTable    <- table(data[indexOrder[i, ], tgt])
                classMode[i]  <- names(which.max(classTable))
            } 
        # }
        indexENN[data[, tgt]!= classMode] <- TRUE
        return(indexENN)
    }

