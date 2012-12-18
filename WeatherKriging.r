library(rgdal)
library(maptools)
library(gstat)
library(sp)
library(plyr)
library(SDMTools)
library(latticeExtra)

################
# Prediction   #
# Grid Creator #
################

makeGrid <- function(dir, res){
	#dir <- paste("C:\\Users\\dylan.mcleod\\Google Drive\\Study\\Summer Research Scholarship 2012\\Weather\\Shapefiles")
	setwd(dir)

	###import geographic border from shape file
	border<-readOGR("map1.shp", "map1")


	###Border bounds (in GPS degrees) stored in vals
	vals <- border@bbox

	###Height and width of grid
	deltaLong <- as.double((vals[1,2] - vals[1,1]))
	deltaLat <- as.double((vals[2,2] - vals[2,1]))

	###Grid reso (GPS degrees)
	gridRes <- res

	###Size of grid in cells
	gridSizeX <- deltaLong / gridRes
	gridSizeY <- deltaLat / gridRes

	###Grid Topology stored in grd
	grd <- GridTopology(vals[,1],c(gridRes,gridRes),c(gridSizeX,gridSizeY))

	###Coordinates and Spatial (bbox, proj4string) stored in SpatialPoints object
	pts <- SpatialPoints(coordinates(grd))

	###SpatialPoints, CRS stored in SpatialPointsDataFrame object pts1
	pts1 <- SpatialPointsDataFrame(as.data.frame(pts), data=as.data.frame(rep(1,nrow(as.data.frame(pts)))))

	###Overlay pts1 and border of grid
	Overlay=overlay(pts1,border)
	pts1$border=Overlay

	###Create grid coordinates vector
	nona<-na.exclude(as.data.frame(pts1))
	coordinates(nona)=~x+y
	gridded(nona) <- TRUE

	### Set coordinates system to GPS
	proj4string(nona)=CRS("+init=epsg:4326")

	###Write grid to ascii file
	writeAsciiGrid(nona,"prediction_grid.asc")

	###Show visualisation of grid and return it
	plot(nona)
	return(nona)

}

###################
# Load Prediction #
# Grid From File  #
###################

loadGrid <- function(dir){
	setwd(dir)
	nona <- as.data.frame(readAsciiGrid("prediction_grid.asc"))
	return(nona)
}


##############
# Load Data  #
##############

loadData <- function(dir){

	setwd(dir)

	###All files read into fileList dataframe list. Length of list stored in numFiles
	fileList <- list.files()
	allData <- lapply(fileList, read.table, header = TRUE, sep= "\t")
	return(allData)
}

##################
# Load Elevation #
# Data           #
##################

loadDem <- function(dir){

	setwd(dir)
	elevation<-read.asciigrid("elevation.asc")
	proj4string(elevation)=CRS("+init=epsg:4326")
	return(elevation)
	#image(elevation)
}


##############
# Date       #
# Randomiser #
##############


dateRandomiser <- function(){

	###Select random year and month between 2002 and 2011
	year <- round(runif(1, 2002, 2011))
	month <- round(runif(1, 1, 12))


	###Number of days established depending on month chosen
	if(month == 2)
	{
		numDays <- 28 
	} else if (month == 9 | 4 | 6 | 11)
	{
		numDays <- 30
	} else
	{
		numDays <- 31
	}

	###Select random day
	day <- round(runif(1, 1, numDays))
	return(c(year,month,day))
}

##################
# Get Sample Day #
##################


###Store subset of data for chosen day in data dataframe

getSample <- function(sampleDate, allData){

	data <- lapply(allData, subset, Year == sampleDate[1] & Month == sampleDate[2] & Day == sampleDate[3])
	data <- do.call(rbind, data)
	data <- subset(data, select = c(Lon, Lat, Maximum.temperature..Degree.C., Elevation))
	data <- na.omit(data)

	###coordinates set as Lon and Lat
	coordinates(data)=~Lon+Lat

	###coordinate system set to WGS84
	proj4string(data)=CRS("+init=epsg:4326")

	###Rename variable to MaxTemp
	names(data)[names(data)=="Maximum.temperature..Degree.C."] = "MaxTemp"

	return(data)

}

#################
# Test Splitter #
#################

testSplitter <- function(data,testFraction){
	valSets <- list()
	i<-sample(nrow(data),round(nrow(data)*testFraction))
	j<-c(1:length(data))
	j <- setdiff(j, i)
	test<-rbind(data[(i[1:length(i)]),])
	training<- rbind(data[(j[1:length(j)]),])
	valSets[[1]] <- test
	valSets[[2]] <- training
	return(valSets)
}

#################
# Date Iterator #
#################

dateIterator <- function(startDate){
	date <- startDate
	date <- as.Date(date)
	date <- date+1
	year <- as.numeric(format(date, format = "%Y"))
	month <- as.numeric(format(date, format = "%m"))
	day <- as.numeric(format(date, format = "%d"))
	return(c(year,month,day))
	
}

################################
# Ord Krig Variogram Modelling #
################################

ordKrigMod <- function(data, model){
	
	modChoice <- paste(model)

	###Plot semivariogram of lag distances
	samplevgm <- variogram(MaxTemp~1, data)
	#plot(samplevgm)

	###Variogram model generated. Sill set to variance, dist. Model selected based on parameter passed. Range estimated as 70 (Gau). Nugget set to 0.
	mod<-vgm(psill=var(data$MaxTemp),model=modChoice,range=70,nugget=0)
	# help(vgm)
	# print(mod)

	###Save model in dataframe fit 
	fit<-fit.variogram(samplevgm, model=mod, fit.ranges = TRUE)
	#plot(variogram(MaxTemp~1,data),fit,main="Model")
	return(fit)
}

#######################
# Ord Krig Validation #
#######################

ordKrigTester <- function(nona,valSets,varMod){
	
	varMod = model
	
	test <- valSets[[1]]
	training <- valSets[[2]]
	krig<-krige(MaxTemp~1,training,model=varMod,newdata=test)
	avgSqErr<-mean((test$MaxTemp-krig$var1.pred)^2)
	return(avgSqErr)
}

##############
# Ord Kriger #
##############

ordKriger <- function(data,varMod){
	cross1<-krige.cv(MaxTemp~1,data,model=fit)
	plot(data$MaxTemp,cross1$var1.pred,asp=1,xlab="Observed",ylab="Predicted")
	abline(0,1,col="red",cex=0.5)
}

#############
# Transform #
#############

trans <- function(dFrame)
{
	dFrame = transform(dFrame, Elevation = log(Elevation))
	coordinates(dFrame)=~Lon+Lat
	proj4string(dFrame)=CRS("+init=epsg:4326")
	return(dFrame)
}


##############
# Validation #
##############

validation <- function(type, nona, valSet, model, range = 70){
	#valSet <- sampleData
	#model <- "Lin"
	#type <- "cokrig"
	modChoice <- paste(model)
	typeChoice <- paste(type)
	rawResults <- matrix(0,length(valSet), 2)
	procResults <- vector(mode = "numeric", 4)

	if(type == "ordkrig"){
		for(i in 1:length(valSet)){
			test <- valSet[i,]
			training <- valSet[-i,]
			krig<-krige(MaxTemp~1,training,model=model,newdata=test)
			rawResults[i,1] <- test$MaxTemp
			rawResults[i,2] <- krig$var1.pred
		}
	} else if(type == "cokrig"){
		for(i in 1:length(valSet)){
			test <- valSet[i,]
			training <- valSet[-i,]
			test = trans(test)
			training = trans(training)
			gv <- gstat(id="MaxTemp",formula=MaxTemp~1,data=training)
			gv <- gstat(gv,id="Elevation",formula=Elevation~1,data=training)
			vario <- variogram(gv)
			gv <- gstat(gv,id=c("MaxTemp","Elevation"),model=vgm(psill=cov(training$MaxTemp,training$Elevation),model=modChoice,range=range, nugget=0))
			gv <- fit.lmc(vario,gv,model=vgm(psill=cov(training$MaxTemp,training$Elevation),model=modChoice,range=range,rangenugget=0))
			cokrig<-predict.gstat(gv,test)
			rawResults[i,1] <- test$MaxTemp
			rawResults[i,2] <- cokrig$MaxTemp.pred
		}
	}else{}
	# Mean residual (sqrt(sumofsquares))
	procResults[1] <- sqrt(mean((rawResults[,1]-rawResults[,2])^2))

	# Median absolute deviation
	procResults[2] <- mad(sqrt(((rawResults[,1]-rawResults[,2])^2)), center = median(sqrt((rawResults[,1]-rawResults[,2])^2)))
		
	# Standard deviation of residuals
	procResults[3] <- sd(sqrt((rawResults[,1]-rawResults[,2])^2))
		
	# Spearman's rank correlation coefficient
	procResults[4] <-as.numeric(cor.test(rawResults[,1],rawResults[,2], method = "spearman")$estimate)

	return(procResults)
}



#####################
# CoKrig Validation #
#####################

coKrigTester <- function(nona,valSets,model,range){

	#nona <- predGrid
	#model <- "Lin"
	modChoice <- paste(model)
	test <- valSets[[1]]
	training <- valSets[[2]]
	#over=overlay(elevation,test)
	#names(test)[names(test)=="elevation.asc"] <- "Elevation"
	test = trans(test)
	training = trans(training)

	###Create gstat object for training set
	gv<-gstat(id="MaxTemp",formula=MaxTemp~1,data=training)
	gv<-gstat(gv,id="Elevation",formula=Elevation~1,data=training)

	###Fit linear model of coregionalisation of training set
	gv<-gstat(gv,id=c("MaxTemp","Elevation"),model=vgm(psill=cov(training$MaxTemp,training$Elevation),model=modChoice,nugget=0))
	gv<-fit.lmc(variogram(gv),gv,model=vgm(psill=cov(training$MaxTemp,training$Elevation),model=modChoice,range=range,nugget=0))

	#plot(variogram(gv),gv$model)

	cokrig<-predict.gstat(gv,test)

	###Residual metrics
	#CVRSQR<-as.numeric(cor.test(test$MaxTemp,cokrig$MaxTemp.pred)$estimate)^2	#Pearson's Correlation Coefficient
	#CVRMSD<-sqrt(sum((test$MaxTemp-cokrig$MaxTemp.pred)^2)/length(test$MaxTemp))	#Root Mean Square Deviation

	###Plot observed versus residuals
	#plot(vario,g$model)

	###Return mean of square errors
	avgSqErr<-mean((test$MaxTemp-cokrig$MaxTemp.pred)^2)
	return(avgSqErr)

}

############
# CoKriger #
############

coKriger <- function(nona,data,elevation,model,range){

	#nona <- predGrid
	#model <- "Lin"
	#data <- sampleData
	#elevation<-elevData
	modChoice <- paste(model)

	over=overlay(elevation,nona)
	nona$elevation=over$elevation.asc

	nona <- transform(nona, elevation = log(elevation))
	coordinates(nona)=~x+y
	proj4string(nona)=CRS("+init=epsg:4326")

	data <- trans(data)
	

	###Create gstat object containing variate and covariate
	g<-gstat(id="MaxTemp",formula=MaxTemp~1,data=data)
	g<-gstat(g,id="Elevation",formula=Elevation~1,data=data)

	###Plot sample semivariogram of lag distances
	vario<-variogram(g)
	#plot(vario)

	###Fit linear model of coregionalisation
	g<-gstat(g,id=c("MaxTemp","Elevation"),model=vgm(psill=cov(data$MaxTemp,data$Elevation),model=modChoice,range=range,nugget=0))
	g<-fit.lmc(vario,g,model=vgm(psill=cov(data$MaxTemp,data$Elevation),model=modChoice,range=range,rangenugget=0))
	#plot(vario,g$model)
	k<-predict.gstat(g,nona)

	return(k)
}


###############
# Visualisers #
###############

ordKrigVisualise <- function(data,varMod,nona){
	map<-krige(MaxTemp~1, data, model=varMod, newdata=nona)
	spplot(map,"var1.pred",col.regions=colorRampPalette(c('dark blue','white','dark red')),main="Ord Krig Prediction Map (Temp)",scales=list(draw=T))
}

coKrigVisualise <- function(k){
	###Map predictions
	cuts <- seq(floor(min(k$MaxTemp.pred)), ceiling(max(k$MaxTemp.pred)), length.out=11)
	col.regions <- brewer.pal(n=11, "RdBu")
	col.regions <-rev(col.regions)
	print(spplot(k,"MaxTemp.pred",cuts=cuts,col.regions=col.regions,main="Co-Krig Prediction Map (Temp-Elevation)",colorkey=list(labels=list(at=cuts),at=cuts), pretty=TRUE, scales=list(draw=T)))

}



##########
## MAIN ##
##########

mapDir <- paste("C:\\Users\\dylan.mcleod\\Google Drive\\Study\\Summer Research Scholarship 2012\\Weather\\Shapefiles")
predGrid <- makeGrid(mapDir, 0.03)

dataDir <- paste("C:\\Users\\dylan.mcleod\\Google Drive\\Study\\Summer Research Scholarship 2012\\Weather\\Temperature")
allData <- loadData(dataDir)

covDir <- paste("C:\\Users\\dylan.mcleod\\Google Drive\\Study\\Summer Research Scholarship 2012\\Weather\\Elevation Data\\Cropped DEM")
elevData <- loadDem(covDir)
image(elevData, col = topo.colors(20))



iters <- 730
resultsKrig <- matrix(data = NA, nrow = iters, ncol = 4)
colnames(resultsKrig) = c("MeanResidual", "MedianAbsDev", "StdDevofResiduals", "Spearman")

resultsKrigGau <- matrix(data = NA, nrow = iters, ncol = 4)
colnames(resultsKrig) = c("MeanResidual", "MedianAbsDev", "StdDevofResiduals", "Spearman")

resultsCoKrig <- matrix(data = NA, nrow = iters, ncol = 4)
colnames(resultsCoKrig) = c("MeanResidual", "MedianAbsDev", "StdDevofResiduals", "Spearman")

resultsCoKrigGau <- matrix(data = NA, nrow = iters, ncol = 4)
colnames(resultsCoKrigGau) = c("MeanResidual", "MedianAbsDev", "StdDevofResiduals", "Spearman")

sampleDate <- dateIterator("2009-12-31")

for(i in 1:iters){
	sampleData <- getSample(sampleDate,allData)
	model <- ordKrigMod(sampleData, "Lin")
	resultsKrig[i,] <- validation("ordkrig", predGrid, sampleData, model)
	model <- ordKrigMod(sampleData, "Gau")
	
	setwd("C:\\Users\\dylan.mcleod\\Google Drive\\Study\\Summer Research Scholarship 2012\\Weather\\Output\\OrdKrig")
	name = (cbind(paste(sampleDate, collapse = "-"),".png"))
	png(filename= paste(name, collapse=""))
	ordKrigVisualise(sampleData,model,predGrid)
	dev.off()
	help(paste)

	resultsKrigGau[i,] <- validation("ordkrig", predGrid, sampleData, model, 80)
	resultsCoKrig[i,] <- validation("cokrig",predGrid, sampleData, "Lin")
	krigedGrid <- coKriger(predGrid,sampleData,elevData,"Lin",80)
	
	setwd("C:\\Users\\dylan.mcleod\\Google Drive\\Study\\Summer Research Scholarship 2012\\Weather\\Output\\CoKrig")
	png(filename= paste(name, collapse=""))
	coKrigVisualise(krigedGrid)
	dev.off()

	resultsCoKrigGau[i,] <- validation("cokrig", predGrid, sampleData, "Gau", 80)
	sampleDate <- dateIterator(paste(sampleDate, collapse = "-"))
}

resultsKrig
resultsKrigGau
resultsCoKrig
resultsCoKrigGau

meanResultsKrig <- c((mean(resultsKrig[,1])), (mean(resultsKrig[,2])), (mean(resultsKrig[,3])), (mean(resultsKrig[,4])))
meanResultsKrigGau <- c((mean(resultsKrigGau[,1])), (mean(resultsKrigGau[,2])), (mean(resultsKrigGau[,3])), (mean(resultsKrigGau[,4])))
meanResultsCoKrig <- c((mean(resultsCoKrig[,1])), (mean(resultsCoKrig[,2])), (mean(resultsCoKrig[,3])), (mean(resultsCoKrig[,4])))
meanResultsCoKrigGau <- c((mean(resultsCoKrigGau[,1])), (mean(resultsCoKrigGau[,2])), (mean(resultsCoKrigGau[,3])), (mean(resultsCoKrigGau[,4])))

#ordKrigVisualise(sampleData,model,predGrid)
#krigedGrid <- coKriger(predGrid,sampleData,elevData,"Lin",80)
#krigedGrid$MaxTemp.pred
#krigedGrid <- coKriger(predGrid,sampleData,elevData,"Gau",80)
#coKrigVisualise(krigedGrid)



#sampleDate <- dateRandomiser()
#sampleData <- getSample(sampleDate,allData)





