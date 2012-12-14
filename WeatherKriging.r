library(rgdal)
library(maptools)
library(gstat)
library(sp)
library(plyr)


################
# Prediction   #
# Grid Creator #
################

makeGrid <- function(dir, res){

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

################################
# Ord Krig Variogram Modelling #
################################

ordKrigMod <- function(data, model){
	
	modChoice <- paste(model)

	###Plot semivariogram of lag distances
	samplevgm <- variogram(MaxTemp~1, data)
	#plot(samplevgm)

	###Variogram model generated. Sill set to variance, dist. Model selected based on parameter passed. Range estimated as 60 (Gau). Nugget set to 0.
	mod<-vgm(psill=var(data$MaxTemp),model=modChoice,range=60,nugget=0)
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


#####################
# CoKrig Validation #
#####################

coKrigTester <- function(nona,valSets,model){

	#nona <- predGrid
	#model <- "Lin"
	modChoice <- paste(model)

	test <- valSets[[1]]
	training <- valSets[[2]]

	###Create gstat object for training set
	gv<-gstat(id="MaxTemp",formula=MaxTemp~1,data=training)
	gv<-gstat(gv,id="Elevation",formula=Elevation~1,data=training)

	###Fit linear model of coregionalisation of training set
	gv<-gstat(gv,id=c("MaxTemp","Elevation"),model=vgm(psill=cov(training$MaxTemp,training$Elevation),model=modChoice,nugget=0))
	gv<-fit.lmc(variogram(gv),gv,model=vgm(psill=cov(training$MaxTemp,training$Elevation),model=modChoice,nugget=0))

	#plot(variogram(gv),gv$model)

	cokrig<-predict.gstat(gv,test)
	help(predict.gstat)

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

coKriger <- function(nona,data,model){

	modChoice <- paste(model)
	
	###Create gstat object containing variate and covariate
	g<-gstat(id="MaxTemp",formula=MaxTemp~1,data=data)
	g<-gstat(g,id="Elevation",formula=Elevation~1,data=data)

	###Plot sample semivariogram of lag distances
	vario<-variogram(g)
	#plot(vario)

	###Fit linear model of coregionalisation
	g<-gstat(g,id=c("MaxTemp","Elevation"),model=vgm(psill=cov(data$MaxTemp,data$Elevation),model=modChoice,range=60,nugget=0))
	g<-fit.lmc(vario,g,model=vgm(psill=cov(data$MaxTemp,data$Elevation),model=modChoice,range=60,rangenugget=0))
	#plot(vario,g$model)
	k<-predict.gstat(g,nona)

	return(k)
}


###############
# Visualisers #
###############

ordKrigVisualise <- function(data,varMod,nona){
	map<-krige(MaxTemp~1, data, model=varMod, newdata=nona)
	spplot(map,"var1.pred",col.regions=heat.colors(20),main="Prediction Map",scales=list(draw=T))
}

coKrigVisualise <- function(k){
	###Map predictions
	spplot(k,"MaxTemp.pred",col.regions=heat.colors(20),main="Prediction Map",scales=list(draw=T))

	###Map errors
	#spplot(k,"MaxTemp.var",col.regions=terrain.colors(20),main="Error Map",scales=list(draw=T))
	#help(spplot)
}



##########
## MAIN ##
##########

mapDir <- paste("C:\\Users\\dylan.mcleod\\Google Drive\\Study\\Summer Research Scholarship 2012\\Weather\\Shapefiles")
predGrid <- makeGrid(mapDir, 0.03)

dataDir <- paste("C:\\Users\\dylan.mcleod\\Google Drive\\Study\\Summer Research Scholarship 2012\\Weather\\Temperature")
allData <- loadData(dataDir)

iters <- 100
sqErrKrig <- numeric(iters)

for(i in 1:iters){
	sampleDate <- dateRandomiser()
	sampleData <- getSample(sampleDate,allData)
	sampleData$Elevation <- NULL
	model <- ordKrigMod(sampleData, "Lin")
	#model <- ordKrigMod(sampleData, "Gau")
	valSets <- testSplitter(sampleData, 0.25)
	sqErrKrig[i] <- ordKrigTester(predGrid,valSets,model)
}
rootMeanSqrErr <- sqrt(mean(sqErrKrig, trim=.1))


ordKrigVisualise(sampleData,model,predGrid)



allData <- loadData(dataDir)
sqErrCoKrig <- numeric(iters)

for(i in 1:iters){
	sampleDate <- dateRandomiser()
	sampleData <- getSample(sampleDate,allData)
	valSets <- testSplitter(sampleData, 0.25)
	sqErrCoKrig[i] <- coKrigTester(predGrid,valSets,"Lin")
	#sqErrCoKrig <- coKrigTester(predGrid,valSets,"Gau")
}

rootMeanSqrErr <- sqrt(mean(sqErrCoKrig, trim=.1))


krigedGrid <- coKriger(predGrid,sampleData,"Lin")
#krigedGrid <- coKriger(predGrid,sampleData,"Gau")
coKrigVisualise(krigedGrid)