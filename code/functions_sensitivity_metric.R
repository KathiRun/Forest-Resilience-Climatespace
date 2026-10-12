
#weighted r2 function
r2ww <- function(x){
  SSe <- sum((x$w*x$resid)^2); #the residual sum of squares is weighted
  observed <- x$resid+x$fitted;
  SSt <- sum((x$w*observed-mean(x$w*observed))^2); #the total sum of squares is weighted
  value <- 1-SSe/SSt;
  return(value);
}

# colSds function
colSds<-function(x) apply(x, 2, function(x) sd(x, na.rm=T))

# detrend function
det.func<-function(x) {
  x<-t(x)
  x<- x-rowMeans(x, na.rm=T)	
  t(x)
}

# z score calc function
z.func<-function(x){
  x<-t(x)
  x<-x-rowMeans(x, na.rm=T)	
  x<-x/ apply(x, 1, function(x) sd(x, na.rm=T))
  t(x)
}

#vsi function
vsi <-function(ndvi.ts, temp.ts, rad.ts, wat.ts,init.y=2001,end.y=2020, ndviThresh=0.1, niter=50){
  #ndvi below 0.1 means dormant
  #here we use kNDVI
  
  
  # create the vector
  # Define variable names required in analysis
  years <- c(2001:2020)
  months <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug","sep", "oct", "nov", "dec")
  vars <- c("temp", "rad", "wat")
  
  if(any(!is.na(ndvi.ts)) && any(!is.na(temp.ts)) && any(!is.na(wat.ts)) && any(!is.na(rad.ts))){
    
    # Combine dataset	
    full <- list(ndvi.ts = ndvi.ts, temp.ts = temp.ts, rad.ts = rad.ts, wat.ts = wat.ts)
    
    # Replace no data values (-999) with NAs
    full <- lapply(full,function(x) replace(x,  which(x==-999), NA))
    
    # Build the year-month matrices	
    full <- lapply(full, function(x)matrix(x, nrow = 20, byrow = TRUE, dimnames = list(years, months)))
    
    # Calculate detrended time series
    full.det=lapply(full, function(x) return(resid(lm(data = data.frame(x=as.numeric(t(x)),t=1:length(as.numeric(x))),x~t,na.action = na.exclude))))
    
    full.det <- lapply(full.det, function(x)matrix(x, nrow = 20, byrow = TRUE, dimnames = list(years, months)))
    
    # Calculate deseasonalized and normalized matrices
    fullDet <- lapply(full.det, function(x)det.func(x))
    fullZ <- lapply(full.det, function(x)z.func(x))
    
    # Make t+1 variable	
    t <- as.vector(t(fullZ$ndvi.ts))
    t1 <- c(NA, t[-length(t)])
    fullZ$t1 <- matrix(t1, nrow=20, byrow = TRUE, dimnames=list(years, months))
    
    ######################################
    ####  specific for each period    ####
    ######################################
    
    years <- as.character(init.y:end.y)
    
    periodl <- length(years)
    
    # Combine dataset	
    full <- list(ndvi.ts = full$ndvi.ts[(init.y-2000):(end.y-2000),],
                 temp.ts = full$temp.ts[(init.y-2000):(end.y-2000),],
                 rad.ts = full$rad.ts[(init.y-2000):(end.y-2000),],
                 wat.ts = full$wat.ts[(init.y-2000):(end.y-2000),])
    
    fullDet <- list(ndvi.ts = fullDet$ndvi.ts[(init.y-2000):(end.y-2000),],
                    temp.ts = fullDet$temp.ts[(init.y-2000):(end.y-2000),],
                    rad.ts = fullDet$rad.ts[(init.y-2000):(end.y-2000),],
                    wat.ts = fullDet$wat.ts[(init.y-2000):(end.y-2000),])
    
    fullZ <- list(ndvi.ts = fullZ$ndvi.ts[(init.y-2000):(end.y-2000),],
                  temp.ts = fullZ$temp.ts[(init.y-2000):(end.y-2000),],
                  rad.ts = fullZ$rad.ts[(init.y-2000):(end.y-2000),],
                  wat.ts = fullZ$wat.ts[(init.y-2000):(end.y-2000),],
                  t1 = fullZ$t1[(init.y-2000):(end.y-2000),])
    
    # Calculate the climatological means +sds
    means <- lapply(full, function(x)colMeans(x, na.rm=T))
    sds <- lapply(full, function(x)colSds(x))#is that useful ?
    
    ##################
    #  masking step  #
    ##################
    
    # Find the months which should be masked due to dormancy (due to mean temperature being less than 0 and mean ndvi < 1000
    
    maskTemp <- which(means$temp.ts <= 0)
    maskNdvi <- which(means$ndvi.ts <= ndviThresh | is.na(means$ndvi.ts)) 
    
    maskMonths<-c(maskNdvi, maskTemp)
    maskMonths<- sort(unique(maskMonths)) # months to be excluded
    
    mask<-c(1,2,3,4,5,6,7,8,9,10,11,12)  
    if(length(maskMonths >0)) mask<- mask[-maskMonths] #months to be included
    
    ##### autocorrelation, instead of removing months, put NAs
    if(length(mask)<3){
      cortot <- NA
    }else{
      timesNA <- matrix(NA, nrow = periodl, ncol=12, byrow = TRUE, dimnames = list(years, months))
      timesNA[,mask] <- fullZ$ndvi.ts[,mask]
      reference <- as.vector(t(timesNA))
      timeslagNA <- matrix(NA, nrow = periodl, ncol=12, byrow = TRUE, dimnames = list(years, months))
      timeslagNA[,(mask%%12)+1] <- fullZ$t1[,(mask%%12)+1]
      lag <- as.vector(t(timeslagNA))
      
      if(sum(complete.cases(lag, reference))<35){
        cortot <- NA
      }else{
        e <- tryCatch(cor(reference, lag,  use = "complete.obs"), error = function(e) e = NA)
        if(is.na(e)){
          cortot <- NA
        }else{
          cortot <- cor(reference, lag,  use = "complete.obs")
        }
      }
    }
    
    # Mask the time series	
    full.mask <- lapply(full, function(x)x[,mask])
    fullDet.mask <- lapply(fullDet, function(x)x[,mask])	
    fullZ.mask <- lapply(fullZ, function(x)x[,mask])
    
    # Do an NA check- there are some pixels that flatline in ndvi or in wat. Want to identify these pixels and remove
    screeningMatrix <- matrix(unlist(fullZ.mask), ncol=5)
    NAcheck<- apply(screeningMatrix, 2, function(x) length(x[is.na(x)==T]))/nrow(screeningMatrix) # Gives you the proportion NAs in each column in the zscores.
    
    NAcheckResult <- length(which(NAcheck == 1)) # Could change the 1 to an argument- then can easily adjust this value in the settings to see the influence of these NAs. Currently removing if full of nas
    
    #################
    # P C A   R E G #
    #################
    
    # Make objects that are needed to sucessfully run the PCA regression and store the results	 
    pc1<-matrix(NA, nrow=periodl, ncol=12, dimnames=list(years, months))
    pc2<-matrix(NA, nrow=periodl, ncol=12, dimnames=list(years, months))
    pc3<-matrix(NA, nrow=periodl, ncol=12, dimnames=list(years, months))
    pc4<-matrix(NA, nrow=periodl, ncol=12, dimnames=list(years, months))
    
    PCAloadings<-array(NA, dim=c(4,4,12), dimnames=list(c("1","2", "3", "4"), c(vars, "t1"), months))
    rsquared<-rep(NA, 12); names(rsquared)<-months
    pVal<-rep(NA, 12); names(rsquared)<-months
    
    sig.var.inf.sum<-matrix(NA, nrow=4, ncol=12, dimnames=list(c(vars, "t1"), months))
    upCI.sum <- matrix(NA, nrow=4, ncol=12, dimnames=list(c(vars, "t1"), months))
    lowCI.sum <- matrix(NA, nrow=4, ncol=12, dimnames=list(c(vars, "t1"), months))
    
    if(NAcheckResult == 0) {
      
      # OK, run the PCA analysis for each month.
      
      PCmask<- seq(1, 12)
      PCmask[maskMonths]<-NA
      
      for(iii in 1:12){
        
        if(is.na(PCmask[iii])==FALSE){
          
          # select the relevant monthly data
          tempPCA <- fullZ$temp.ts[,iii] 
          radPCA <- fullZ$rad.ts[,iii]
          watPCA <- fullZ$wat.ts[,iii]
          t1PCA <- fullZ$t1[,iii]
          if(length(t1PCA[is.na(t1PCA)]) == periodl) t1PCA <- rep(0, periodl)				
          
          
          # combine these into a matrix
          varPCA<-cbind(tempPCA, radPCA, watPCA, t1PCA); colnames(varPCA)<-c(vars, "t1") 
          
          # In one pixel found ndvi Zscores NaNs because variance was zero in November. Replaced NaNs with 0
          varPCA[is.nan(varPCA)]<-0
          
          # remove rows with NAs in them (required for PCA)
          varPCA = na.omit(varPCA)					
          
          # run the PCA. Only doing a PCA if there are more than 4 complete rows
          if(nrow(varPCA) >4){
            
            pca.model<-princomp(varPCA) 
            
            # estimate the loadings of each of the variables on the 3 different PC axes
            PCAloadings[,,iii]<-rbind(pca.model$loadings[,1], pca.model$loadings[,2], pca.model$loadings[,3], pca.model$loadings[,4]) 
            PCAloadings1<-data.frame(rbind(pca.model$loadings[,1], pca.model$loadings[,2], pca.model$loadings[,3], pca.model$loadings[,4]))
            
            # identify which rows of PC tables should be filled
            pcTableMatch = match(row.names(pca.model$scores), row.names(pc1))
            
            # Find the scores for PC1, PC2, PC3, PC4
            pc1[pcTableMatch,iii]<-pca.model$scores[,1] 
            pc2[pcTableMatch,iii]<-pca.model$scores[,2]
            pc3[pcTableMatch,iii]<-pca.model$scores[,3]
            pc4[pcTableMatch,iii]<-pca.model$scores[,4]
            
            # Regress the PC scores against the ndvi data
            varTable <- data.frame(resp = fullZ$ndvi.ts[,iii], PC1 = pc1[,iii], PC2 = pc2[,iii], PC3 = pc3[,iii] , PC4 = pc4[,iii])
            varTable <- varTable[complete.cases(varTable), ]
            
            if(nrow(varTable)>=6){	
              #if the last axis is too low explanatory
              c <- tryCatch(rlm(resp ~ PC1 + PC2 + PC3 + PC4, data = varTable, maxit=niter), error = function(c) c = NA)
              f <- tryCatch(summary(rlm(resp ~ PC1 + PC2 + PC3 + PC4, data = varTable, maxit=niter)), error = function(f) f = NA)
              if(!is.list(c) || !is.list(f)){
                rm(c)
                #if the 3rd axis is too low explanatory
                d <- tryCatch(rlm(resp ~ PC1 + PC2 + PC3, data = varTable, maxit=niter), error = function(d) d = NA)
                g <- tryCatch(summary(rlm(resp ~ PC1 + PC2 + PC3, data = varTable, maxit=niter)), error = function(g) g = NA)
                if(!is.list(d) || !is.list(g)){
                  rm(d)
                  model <- rlm(resp ~ PC1 + PC2, data = varTable, maxit=niter)
                  
                  # extract the cofficients of the PC regression 
                  regCoeff <- model$coefficients[-1]
                  
                  # Find the rsquared and pvalue of the PC regression
                  
                  a <- tryCatch(r2ww(model), warning = function(w) w = NA)
                  if(is.na(a)){
                    rm(a)
                  } else {						
                    
                    rsquared[iii]<-r2ww(model)
                    if(rsquared[iii]>0.99){
                      pVal[iii] <- 0.00001
                    }else{
                      b <- tryCatch(f.robftest(model)$p.value, error = function(b) b = NA)
                      if(is.na(b)){
                        pVal[iii]<-NA
                      }else{
                        pVal[iii]<-f.robftest(model)$p.value
                      }
                    }
                    cInf <- confint.default(object = model, parm = c("PC1", "PC2", "PC3", "PC4"), level = 0.9)
                    
                    # if a pc axis is signicant at p = 0.1 it will tramsform these back to the original coefficients by mulitplying the loadings of that axis by the regression coefficients
                    sigCoef <- which(c(f.robftest(model, var="PC1")$p.value,
                                       f.robftest(model, var="PC2")$p.value,
                                       NA,
                                       NA) <0.1)
                    #end of if else where the 3th axis is weird
                    if(length(sigCoef)!=0){							
                      sig.var.inf<-PCAloadings1[sigCoef,]*regCoeff[sigCoef]
                      sig.var.inf.sum[,iii] <- apply(sig.var.inf, 2, sum)
                      # Find the upper and lower confidence in the same way
                      upCI <- PCAloadings1[sigCoef,]*cInf[sigCoef, 2]
                      upCI.sum[,iii] <- apply(upCI, 2, sum)
                      
                      lowCI <- PCAloadings1[sigCoef,]*cInf[sigCoef, 1]
                      lowCI.sum[,iii] <- apply(lowCI, 2, sum)
                      
                    } # closes the sig if statement	
                    
                  }# tryCatch if statement a
                }else{
                  
                  
                  model <- rlm(resp ~ PC1 + PC2 + PC3, data = varTable, maxit=niter)
                  
                  
                  # extract the cofficients of the PC regression 
                  regCoeff <- model$coefficients[-1]
                  
                  # Find the rsquared and pvalue of the PC regression
                  
                  a <- tryCatch(r2ww(model), warning = function(w) w = NA)
                  if(is.na(a)){
                    rm(a)
                  } else {						
                    
                    rsquared[iii]<-r2ww(model)
                    if(rsquared[iii]>0.99){
                      pVal[iii] <- 0.00001
                    }else{
                      b <- tryCatch(f.robftest(model)$p.value, error = function(b) b = NA)
                      if(is.na(b)){
                        pVal[iii]<-NA
                      }else{
                        pVal[iii]<-f.robftest(model)$p.value
                      }
                    }
                    cInf <- confint.default(object = model, parm = c("PC1", "PC2", "PC3", "PC4"), level = 0.9)
                    
                    # if a pc axis is signicant at p = 0.1 it will tramsform these back to the original coefficients by mulitplying the loadings of that axis by the regression coefficients
                    sigCoef <- which(c(f.robftest(model, var="PC1")$p.value,
                                       f.robftest(model, var="PC2")$p.value,
                                       f.robftest(model, var="PC3")$p.value,
                                       NA) <0.1)
                    #end of if else where the 4th axis is weird
                    if(length(sigCoef)!=0){							
                      sig.var.inf<-PCAloadings1[sigCoef,]*regCoeff[sigCoef]
                      sig.var.inf.sum[,iii] <- apply(sig.var.inf, 2, sum)
                      # Find the upper and lower confidence in the same way
                      upCI <- PCAloadings1[sigCoef,]*cInf[sigCoef, 2]
                      upCI.sum[,iii] <- apply(upCI, 2, sum)
                      
                      lowCI <- PCAloadings1[sigCoef,]*cInf[sigCoef, 1]
                      lowCI.sum[,iii] <- apply(lowCI, 2, sum)
                      
                    } # closes the sig if statement	
                    
                  }# tryCatch if statement a
                }
              }else{
                model <- rlm(resp ~ PC1 + PC2 + PC3 + PC4, data = varTable, maxit=niter)
                
                # extract the cofficients of the PC regression 
                regCoeff <- model$coefficients[-1]
                
                # Find the rsquared and pvalue of the PC regression
                
                a <- tryCatch(r2ww(model), warning = function(w) w = NA)
                if(is.na(a)){
                  rm(a)
                } else {						
                  
                  rsquared[iii]<-r2ww(model)
                  if(rsquared[iii]>0.99){
                    pVal[iii] <- 0.00001
                  }else{
                    b <- tryCatch(f.robftest(model)$p.value, error = function(b) b = NA)
                    if(is.na(b)){
                      pVal[iii]<-NA
                    }else{
                      pVal[iii]<-f.robftest(model)$p.value
                    }
                  }
                  cInf <- confint.default(object = model, parm = c("PC1", "PC2", "PC3", "PC4"), level = 0.9)
                  
                  # if a pc axis is signicant at p = 0.1 it will tramsform these back to the original coefficients by mulitplying the loadings of that axis by the regression coefficients
                  sigCoef <- which(c(f.robftest(model, var="PC1")$p.value,
                                     f.robftest(model, var="PC2")$p.value,
                                     f.robftest(model, var="PC3")$p.value,
                                     f.robftest(model, var="PC4")$p.value) <0.1)
                  #end of if else where the 4th axis is weird
                  if(length(sigCoef)!=0){							
                    sig.var.inf<-PCAloadings1[sigCoef,]*regCoeff[sigCoef]
                    sig.var.inf.sum[,iii] <- apply(sig.var.inf, 2, sum)
                    # Find the upper and lower confidence in the same way
                    upCI <- PCAloadings1[sigCoef,]*cInf[sigCoef, 2]
                    upCI.sum[,iii] <- apply(upCI, 2, sum)
                    
                    lowCI <- PCAloadings1[sigCoef,]*cInf[sigCoef, 1]
                    lowCI.sum[,iii] <- apply(lowCI, 2, sum)
                    
                  } # closes the sig if statement	
                  
                }# tryCatch if statement a
              }# tryCatch if statement c
            } # closes nrow masking statement
            
          }	# closes varPCA masking statement
          
        } # closes the NA masking if statement
      } #closes the monthly for loop for running the pcas/ regressions
      
      
    } # closes the NA check masking statement
    
    
    #############################################
    # # V A R I A N C E   B Y  SIG COEF #
    # #############################################
    monthCoefSel <- which(is.na(sig.var.inf.sum[1,])==FALSE)
    if(length(monthCoefSel) ==0){
      
      sdDetVecSig <- rep(NA, 4)
      meanVecSig <- rep(NA, 4)
      
    }else{
      
      detVecSig <- lapply(fullDet, function(x)c(t(x[,monthCoefSel])))
      fullVecSig <- lapply(full, function(x)c(t(x[,monthCoefSel])))
      
      sdDetVecSig<- unlist(lapply(detVecSig, function(x)sd(x, na.rm=T)))
      meanVecSig <- unlist(lapply(fullVecSig, function(x)mean(x, na.rm=T)))
      
    }
  }else{ # beginning if statement
    months <- months
    mask <- c()
    rsquared<-rep(NA, 12); names(rsquared)<-months
    pVal<-rep(NA, 12); names(rsquared)<-months
    
    sig.var.inf.sum<-matrix(NA, nrow=4, ncol=12, dimnames=list(c(vars, "t1"), months))
    upCI.sum <- matrix(NA, nrow=4, ncol=12, dimnames=list(c(vars, "t1"), months))
    lowCI.sum <- matrix(NA, nrow=4, ncol=12, dimnames=list(c(vars, "t1"), months))
    
    cortot <- NA
    sdDetVecSig <- rep(NA, 4)
    meanVecSig <- rep(NA, 4)
  }#end of the first if not everything is na
  # Store the results output. Processing of these takes place in the mapping script.
  
  results<-list(mask = mask, rsquared = rsquared, sig.var.inf.sum = sig.var.inf.sum, 
                pVal= pVal,  sdDetVecSig = sdDetVecSig, meanVecSig = meanVecSig, upCIsum = upCI.sum, lowCI.sum = lowCI.sum, autocor = cortot)
  
  results		
}



###below put in function file

makedf <- function(x, resultPath, var){
  resultMatrix <- x
  resultMatrix[is.na(resultMatrix)] <- -8888  # converts the NA values in the results to -8888 so that they can be distinguished
  write.table(resultMatrix, file = paste(resultPath, var, ".csv", sep="" ), sep= ",", col.names=T, row.names= FALSE)
}


cutAnomalies <- function(x){
  # Cut the anomolies			
  toPlot <- c(x)
  toPlot1<- toPlot[!is.na(toPlot)]
  toPlot <-toPlot1
  sortedLow<-sort(toPlot)
  lowCut<-sortedLow[length(sortedLow)*0.005]
  sortedHigh<-sort(toPlot, decreasing=T)
  highCut<-sortedHigh[length(sortedHigh)*0.005]
  x[which(x>highCut)]<-highCut
  x[which(x<lowCut)]<-lowCut
  x			
}


Proj <- function(resultMatrix, fullStand = TRUE, resultPath, raw = TRUE){
  
  # resultMatrix = results that need to be put in a df 
  # fullStand = Do you want to standardise between 0 and 100 against all the variables (TRUE), or for each variable individually (FALSE)
  # resultPath = path to save the results raster (as .csv file)
  # raw = Do you want to output the raw values when making the .csv matrix?
  
  print("Making projection to matrix")
  progBar <- txtProgressBar(min = 0, max = length(colnames(resultMatrix)), style = 3)
  
  if (raw == TRUE){
    anomCut <- apply(resultMatrix, 2, function(x)cutAnomalies(x)) # trim the anomalies
    varName <- colnames(anomCut)
    for(i in 1: length(varName)){
      makedf(anomCut[,i], resultPath, var= varName[i]) # Make the df of each variable
      setTxtProgressBar(progBar, i)	
    }
  } else {
    anomCut <- apply(resultMatrix, 2, function(x)cutAnomalies(x)) # trim the anomalies
    if(fullStand == FALSE) stand0_100 <- apply(anomCut, 2, function(x)range01(x)) else stand0_100 <- range01(anomCut) # standardise between 0 and 100
    varName <- colnames(stand0_100)
    for(i in 1: length(varName)) {
      makedf(stand0_100[,i], resultPath, var= varName[i]) # Make the df of each variable
      setTxtProgressBar(progBar, i)
    }
  }
}

compileCoefTiles <- function(results = resultsInput, absVal = TRUE){ # this should also be adapted  
  
  resultMatrix <- matrix(NA, nrow=1, ncol=4)
  
  results <- results

  vsi<-results
  rm(results)
  namesResults <- names(vsi[[1]])
  var <- "sig.var.inf.sum"
  dimension <- 4
  
  xToPlot <- matrix(NA, nrow=length(vsi), ncol= dimension) # Create your results matrix, in this instance will be n.land-pixels x 4.var
  varNum<-which(namesResults== var) # Identify correct item in list for each pixel
  
  for(i in 1:length(vsi)) {
    toPlot<-vsi[[i]] # extract a given pixel
    varCoef <- toPlot[[varNum]] # extract the coefficient matrix 
    
    if(absVal== TRUE){
      xToPlot[i,] <- apply(abs(varCoef), 1, function(x)mean(x, na.rm=TRUE)) 
    }else{
      xToPlot[i,] <- apply((varCoef), 1, function(x)mean(x, na.rm=TRUE)) # if don't want to use absolute coefficient values
    }		
    
  } # End of pixel forloop
  resultMatrix <- rbind(resultMatrix, xToPlot) # Combine all the results from the different variables
  
  resultMatrix <- resultMatrix[-1,] 
  colnames(resultMatrix) <- rownames(varCoef) 
  return(resultMatrix) # All you need to do now is line them up with the landsea mask (see next section)
}


compileAutocor <- function(results = resultsInput, absVal = TRUE){
  
  resultMatrix <- matrix(NA, nrow=1, ncol=1)
  
  vsi<-results
  rm(results)
  namesResults <- names(vsi[[1]])
  var <- "autocor"
  dimension <- 1
  
  xToPlot <- matrix(NA, nrow=length(vsi), ncol= dimension) # Create your results matrix, in this instance will be n.land-pixels x 4.var
  varNum<-which(namesResults== var) # Identify correct item in list for each pixel
  
  for(i in 1:length(vsi)) {
    toPlot<-vsi[[i]] # extract a given pixel
    varCoef <- toPlot[[varNum]] # extract the coefficient matrix 
    
    if(absVal== TRUE){
      xToPlot[i,] <- abs(varCoef)
    }else{
      xToPlot[i,] <- varCoef# if don't want to use absolute coefficient values
    }		
    
    
  } # End of pixel forloop
  resultMatrix <- rbind(resultMatrix, xToPlot) # Combine all the results from the different variables
  
  resultMatrix <- resultMatrix[-1,] 
  colnames(resultMatrix) <- rownames(varCoef) 
  return(resultMatrix) # All you need to do now is line them up with the landsea mask (see next section)
}


compileVarTilesAnon <- function(results = resultsInput, var = "meanVecSig"){
  
  resultMatrix <- matrix(NA, nrow=1, ncol=4)
  print(paste("Compiling the variable", var, sep=" ")) 
  progBar <- txtProgressBar(min =0, max = length(results), style = 3)
  
  namesResults <- names(results[[1]])
  dimension <- 4
  
  vsi <- results
  rm(results)
  xToPlot <- matrix(NA, nrow=length(vsi), ncol= dimension) # Create your results matrix, in this instance will be n.land-pixels x 3.var
  varNum<-which(namesResults== var) # Identify correct item in list for each pixel
  
  for(i in 1:length(vsi)) {
    toPlot<-vsi[[i]] # extract a given pixel
    varCoef <- toPlot[[varNum]] # extract the coefficient matrix 
    xToPlot[i,] <- varCoef
    setTxtProgressBar(progBar, i)
    
  } # End of pixel forloop
  resultMatrix <- rbind(resultMatrix, xToPlot) # Combine all the results from the different variables
  
  resultMatrix <- resultMatrix[-1,] 
  colnames(resultMatrix) <- c("ndvi", "temp", "rad", "wat") 
  # All you need to do now is line them up with the landsea mask (see next section)
  close(progBar)
  resultMatrix
}


anomCalcNew <- function(var = "temp", namefile = paste(mainDir, "/results/anom", sep= "")) {
  
  colSel <- match(var, colnames(mean_sig))
  
  x = mean_sig[,colSel]
  y = sd_sig[,colSel]

  # Compile the data to make the the plot
  require("nlme")
  
  # Residual result is the variable that you want map
  residualResult<-rep(NA, length(x))
  I <- which(is.na(x) == FALSE)
  
  nls2 <- nls(y~b0+b1*(x)+ b2*x^2,start=c(b0=0,b1=0, b2=0) )
  
  colList <- list(ndvi = rgb(100,100 ,100, 50, maxColorValue = 255), temp =  rgb(100, 0 ,0, 50, maxColorValue = 255),
                  wat = rgb(0, 0 , 100, 50, maxColorValue = 255), rad = rgb(0, 100 ,0, 50, maxColorValue = 255))
  
  smp <- if(length(I) < 1000){
    I
  }else{
    sample(I, 1000)
    }
  
  plot(x[smp], y[smp], col= eval(parse(text = paste("colList$", var, sep=""))), pch=16, xlab = "Mean", ylab = "Variance", main = var)
  
  xPred <- seq(floor(min(x, na.rm= TRUE)), ceiling(max(x, na.rm= TRUE)), length.out=100)
  yPred <- predict(nls2, newdata = data.frame(x= xPred))
  
  lines(xPred, yPred)
  
  residualResult[I] <-residuals(nls2)
  anomCut <-cutAnomalies(residualResult) # trim the anomalies
  
  #stand0_100 <- range01(anomCut)#dont change range to 0-100 because we compare the 2 periods
  stand0_100 <- anomCut
  makedf(x= stand0_100, resultPath= namefile, var)

  result <- list(xPred = xPred, yPred = yPred, x = x[smp], y = y[smp], varName = var) 
  result	
}