################################################################
## R code for **Browning et al., 2017. Breaks in MODIS time series 
## portend vegetation change – verification using long-term data 
## in an arid grassland ecosystem**
## Coding and commenting by Jonathan Maynard (USDA-ARS)
## Date: 03-13-2017
## Questions? jonathan.maynard@ars.usda.gov
###############################################################

###############################################################
## Required libraries
###############################################################
library(MODIS)
library(maptools)
library(rgdal)
library(raster)
library(maps)
library(doMPI)
library(foreach)
library(bfast)
library(strucchange)
library(Cairo)

#set up parrallel backend
cl <- startMPIcluster(count=24)
registerDoMPI(cl)
getDoParWorkers()
###############################################################
#Set repository for MODIS imagery
###############################################################
MODIS:::checkTools()
MODISoptions(localArcPath="~/MODIS/MODIS_ARC/",
             outDirPath="~/MODIS/MODIS_ARC/PROCESSES")
 

# #Vegetation Indicies (250 M)
#MODIS tile for USDA Jornada Experimental Range
getHdf(product="MOD13Q1",begin="2015321", tileH=9,tileV=5, collection='005') 
 

## Read in the Jornada boundary shapefile (specify source directory)
JRN <- readOGR(dsn="~/JRN_boundary",layer="JRN_boundary")
JRN.point <- readOGR(dsn="~/Jornada_1M_plot_centroid",layer="Jornada_1m_plot_centroid")



## Set up a raster "template" to use in rasterize()
##expand extent 20000 m on every side to prevent clipping across study area
ext <-  extent (JRN)+40000
xy <- abs(apply(as.matrix(bbox(ext)), 1, diff))  #this calculated the absolute distance in the x and y directions
r <- raster(ext, ncol=xy[1]/250, nrow=xy[2]/250) #this creates a raster box of the extended study area

## Rasterize the shapefile
JRN_raster <-rasterize(JRN, r) #we then create a raster of the JRN boundary within the extended study area
res(JRN_raster) <- 250 #assign resolution
plot(JRN_raster)

#find out which SDS numbers correspond to which bands by examining one downloaded image
#MOD13Q1
getSds(HdfName="MOD13Q1.A2000049.h08v04.005.2006269063040.hdf")

# $SDSnames
# [1] "250m 16 days NDVI"                       [2] "250m 16 days EVI"                      
# [3] "250m 16 days VI Quality"                 [4] "250m 16 days red reflectance"          
# [5] "250m 16 days NIR reflectance"            [6] "250m 16 days blue reflectance"         
# [7] "250m 16 days MIR reflectance"            [8] "250m 16 days view zenith angle"        
# [9] "250m 16 days sun zenith angle"           [10] "250m 16 days relative azimuth angle"   
# [11] "250m 16 days composite day of the year" [12]"250m 16 days pixel reliability" 

#run MODIS function runGdal to extract MODIS Data
#MOD13Q1
runGdal(job="JRN_Modis_VI", product="MOD13Q1", collection='005', extent=JRN_raster, begin="2000001", SDSstring="111111")

#set the path to get the correct data on the source directory: 
#path to 250M VI's
path <- paste(options("MODIS_outDirPath"),"/JRN_Modis_VI",sep="") 
ndvi <- preStack(path=path, pattern="*_NDVI.tif$")
JRN_NDVI_brick<- brick(stack(ndvi),  filename='/output_directory/JRN_NDVI_brick', overwrite=TRUE)


#Raster layers need to be in a vector list to parallelize accross a cluster.
stack2list <- function(stack){
  list<-vector("list", nlayers(stack))
  for(i in 1:nlayers(stack)) {
  list[[i]] <- stack[[i]]
  }
  return(list)
}

JRN_NDVI_list<- stack2list(stack(ndvi))
####################################################################################################################################

######################################################################################################
#MODIS QA Processing

## QA for MOD13 products
## MODIS13Q1 has 16 bit data so we want the fist 16 integers.
## MODIS QA data is parsed from right to left, and the individual bit groups (bit words) are read from left to right.
## Example
## as.integer(intToBits(7425)[16:1])
## produces 0 0 0 1 1 1 0 1 0 0 0 0 0 0 0 1 , the bit words or combinations of integers going from right to left mean different things
## for different MODIS products

###############MODIS13Q1 QA lookup table

QC_Data <- data.frame(Integer_Value = 0:65535,
                      Bit15 = NA,
                      Bit14 = NA,
                      Bit13 = NA,
                      Bit12 = NA,
                      Bit11 = NA,
                      Bit10 = NA,
                      Bit9 = NA,
                      Bit8 = NA,
                      Bit7 = NA,
                      Bit6 = NA,
                      Bit5 = NA,
                      Bit4 = NA,
                      Bit3 = NA,
                      Bit2 = NA,
                      Bit1 = NA,
                      Bit0 = NA,
                      QA_word1 = NA,
                      QA_word2 = NA,
                      QA_word3 = NA,
                      QA_word4 = NA,
                      QA_word5 = NA,
                      QA_word6 = NA,
                      QA_word7 = NA,
                      QA_word8 = NA,
                      QA_word9 = NA
                      
)

for(i in QC_Data$Integer_Value){
  AsInt <- as.integer(intToBits(i)[1:16])
  QC_Data[i+1,2:17]<- AsInt[16:1]
}

#VI Quality
QC_Data$QA_word1[QC_Data$Bit1 == 0 & QC_Data$Bit0==0] <- "VI GOOD"
QC_Data$QA_word1[QC_Data$Bit1 == 0 & QC_Data$Bit0==1] <- "VI Produced,Other Quality"
QC_Data$QA_word1[QC_Data$Bit1 == 1 & QC_Data$Bit0==0] <- "No Pixel,clouds"
QC_Data$QA_word1[QC_Data$Bit1 == 1 & QC_Data$Bit0==1] <- "No Pixel, Other QA"

#VI Usefullness
QC_Data$QA_word2[QC_Data$Bit5 == 0 & QC_Data$Bit4==0 & QC_Data$Bit3 == 0 & QC_Data$Bit2==0] <- "Highest Quality"
QC_Data$QA_word2[QC_Data$Bit5 == 0 & QC_Data$Bit4==0 & QC_Data$Bit3 == 0 & QC_Data$Bit2==1] <- "High Quality"
QC_Data$QA_word2[QC_Data$Bit5 == 0 & QC_Data$Bit4==0 & QC_Data$Bit3 == 1 & QC_Data$Bit2==0] <- "Good Quality"
QC_Data$QA_word2[QC_Data$Bit5 == 0 & QC_Data$Bit4==0 & QC_Data$Bit3 == 1 & QC_Data$Bit2==1] <- "Acceptable Quality"
QC_Data$QA_word2[QC_Data$Bit5 == 0 & QC_Data$Bit4==1 & QC_Data$Bit3 == 0 & QC_Data$Bit2==0] <- "Fair Quality"
QC_Data$QA_word2[QC_Data$Bit5 == 0 & QC_Data$Bit4==1 & QC_Data$Bit3 == 0 & QC_Data$Bit2==1] <- "Intermediate Quality"
QC_Data$QA_word2[QC_Data$Bit5 == 0 & QC_Data$Bit4==1 & QC_Data$Bit3 == 1 & QC_Data$Bit2==0] <- "Below Intermediate Quality"
QC_Data$QA_word2[QC_Data$Bit5 == 0 & QC_Data$Bit4==1 & QC_Data$Bit3 == 1 & QC_Data$Bit2==1] <- "Average Quality"
QC_Data$QA_word2[QC_Data$Bit5 == 1 & QC_Data$Bit4==0 & QC_Data$Bit3 == 0 & QC_Data$Bit2==0] <- "Below Average Quality"
QC_Data$QA_word2[QC_Data$Bit5 == 1 & QC_Data$Bit4==0 & QC_Data$Bit3 == 0 & QC_Data$Bit2==1] <- "Questionable Quality"
QC_Data$QA_word2[QC_Data$Bit5 == 1 & QC_Data$Bit4==0 & QC_Data$Bit3 == 1 & QC_Data$Bit2==0] <- "Above Marginal Quality"
QC_Data$QA_word2[QC_Data$Bit5 == 1 & QC_Data$Bit4==0 & QC_Data$Bit3 == 1 & QC_Data$Bit2==1] <- "Marginal Quality"
QC_Data$QA_word2[QC_Data$Bit5 == 1 & QC_Data$Bit4==1 & QC_Data$Bit3 == 0 & QC_Data$Bit2==0] <- "Low Quality"
QC_Data$QA_word2[QC_Data$Bit5 == 1 & QC_Data$Bit4==1 & QC_Data$Bit3 == 0 & QC_Data$Bit2==1] <- "No Atmopheric Correction"
QC_Data$QA_word2[QC_Data$Bit5 == 1 & QC_Data$Bit4==1 & QC_Data$Bit3 == 1 & QC_Data$Bit2==0] <- "Quality too low to be useful"
QC_Data$QA_word2[QC_Data$Bit5 == 1 & QC_Data$Bit4==1 & QC_Data$Bit3 == 1 & QC_Data$Bit2==1] <- "Not useful for other reasons"

#Aerosol quantity
QC_Data$QA_word3[QC_Data$Bit7 == 0 & QC_Data$Bit6==0] <- "Climatology"
QC_Data$QA_word3[QC_Data$Bit7 == 0 & QC_Data$Bit6==1] <- "Low"
QC_Data$QA_word3[QC_Data$Bit7 == 1 & QC_Data$Bit6==0] <- "Intermediate"
QC_Data$QA_word3[QC_Data$Bit7 == 1 & QC_Data$Bit6==1] <- "High"

#Adjacent Clouds
QC_Data$QA_word4[QC_Data$Bit8 == 0] <- "No"
QC_Data$QA_word4[QC_Data$Bit8 == 1] <- "Yes"

#Atomospheric BRDF correction
QC_Data$QA_word5[QC_Data$Bit9 == 0] <- "No"
QC_Data$QA_word5[QC_Data$Bit9 == 1] <- "Yes"

#Mixed Clouds
QC_Data$QA_word6[QC_Data$Bit10 == 0] <- "No"
QC_Data$QA_word6[QC_Data$Bit10 == 1] <- "Yes"

#Land/Water Mask
QC_Data$QA_word7[QC_Data$Bit13==0 & QC_Data$Bit12 == 0 & QC_Data$Bit11==0] <- "Shallow ocean"
QC_Data$QA_word7[QC_Data$Bit13==0 & QC_Data$Bit12 == 0 & QC_Data$Bit11==1] <- "Land"
QC_Data$QA_word7[QC_Data$Bit13==0 & QC_Data$Bit12 == 1 & QC_Data$Bit11==0] <- "Coastline/Shoreline"
QC_Data$QA_word7[QC_Data$Bit13==0 & QC_Data$Bit12 == 1 & QC_Data$Bit11==1] <- "Shallow inland water"
QC_Data$QA_word7[QC_Data$Bit13==1 & QC_Data$Bit12 == 0 & QC_Data$Bit11==0] <- "Ephemeral water"
QC_Data$QA_word7[QC_Data$Bit13==1 & QC_Data$Bit12 == 0 & QC_Data$Bit11==1] <- "Deep inland water"
QC_Data$QA_word7[QC_Data$Bit13==1 & QC_Data$Bit12 == 1 & QC_Data$Bit11==0] <- "Moderate or continental ocean"
QC_Data$QA_word7[QC_Data$Bit13==1 & QC_Data$Bit12 == 1 & QC_Data$Bit11==1] <- "Deep ocean"

#Possible snow/ice
QC_Data$QA_word8[QC_Data$Bit14 == 0] <- "No"
QC_Data$QA_word8[QC_Data$Bit14 == 1] <- "Yes"

#Possible shadow
QC_Data$QA_word9[QC_Data$Bit15 == 0] <- "No"
QC_Data$QA_word9[QC_Data$Bit15 == 1] <- "Yes"
##############################################################################################
#Run QA on MODIS 250m NDVI
#extract unique QA values from raster brick
JRN_QA_unique<- foreach(i=1:length(JRN_QA_list), .combine='c', .packages = c("raster", "rgdal")) %dopar% {   
  unique(JRN_QA_list[[i]])
}

#final list of unique QA integers from image brick
JRN_QA_unique<-unique(JRN_QA_unique)

#now subset QA_Data to extract only relevant QA integers
JRN_QC_Data <- subset(QC_Data, Integer_Value %in% JRN_QA_unique)

##Now subset based on the most important QA variables
#select all with VI quality of 00
VI.good <- subset(JRN_QC_Data, QA_word1=="VI GOOD")

#select subset of VI Usefullness 
VI.okay <- subset(subset(JRN_QC_Data, QA_word1=="VI Produced,Other Quality"), QA_word2=="Highest Quality" | QA_word2=="High Quality" | QA_word2=="Good Quality"
                  | QA_word2=="Fair Quality" | QA_word2=="Intermediate Quality" | QA_word2=="Below Intermediate Quality"| QA_word2=="Average Quality")

#combine good and okay integer values and then assign to a vector.
VI.QA <- rbind(VI.good, VI.okay)
QA.int<-VI.QA$Integer_Value

#Create JRN_QA raster list object to modify
JRN_QA <- JRN_QA_list

#First loop assigns all good pixels a value of 1 based on previous selection
for(i in 1:length(JRN_QA)){
  
  for(j in 1:length(QA.int)){
    JRN_QA[[i]][JRN_QA[[i]]==QA.int[j]] <-1
    
  }}

#Second loop assigns all values greater than 1 to NA
for(i in 1:length(JRN_QA)){
  JRN_QA[[i]][JRN_QA[[i]]>1] <-NA
}

## Final loop creates mask to mask out bad pixels
## This loop uses the foreach fucntion in a cluster configuration and thus requires 
## that the raster data be in a list format so that it can subset the data

#################Mask NDVI
JRN_NDVI_mask<-JRN_NDVI_list
JRN_NDVI_mask<- foreach(i=1:length(JRN_NDVI_mask), .packages = c("raster", "rgdal"))  %dopar% { 
  JRN_NDVI_mask[[i]]<-mask(JRN_NDVI_mask[[i]], (JRN_QA[[i]]))
}

#Final NDVI product
JRN_NDVI_mask <- brick(stack(JRN_NDVI_mask),  filename='/output_directory/JRN_NDVI_mask', overwrite=TRUE)

#Test it to see results
plot(JRN_NDVI_mask[[1]])
medianNDVI.mask <- calc(JRN_NDVI_mask, fun=function(x) median(x, na.rm = TRUE)) 
plot(medianNDVI.mask)
###################################################################################################
###################################################################################################
######################################################
####Now we process the MODIS data with the bfast package
## Processing raster bricks (satellite image time series of 16-day NDVI images) 

#modified code for function bfastts to account for leap year
bfastts.na <- function (data, dates, type = c("irregular", "16-day", "quarterly")) 
{
  yday365 <- function(x) {
    x <- as.POSIXlt(x)
    mdays <- c(31L, 28L, 31L, 30L, 31L, 30L, 31L, 31L, 30L, 
               31L, 30L, 31L)
    cumsum(c(0L, mdays))[1L + x$mon] + x$mday
  }
  if (type == "irregular") {
    zz <- zoo(data, 1900 + as.POSIXlt(dates)$year + (yday365(dates) - 
                                                       1)/365, frequency = 365)
  }
  if (type == "16-day") {
    z <- zoo(data, dates)
    yr <- as.numeric(format(time(z), "%Y"))
    jul <- as.numeric(format(time(z), "%j"))
    delta <- min(unlist(tapply(jul, yr, diff)))
    zz <- aggregate(z, yr + (jul - 1)/delta/23)
  }
  
  if (type == "quarterly") {
    zz <- zoo(data, dates)
  }
  #tso <- na.approx(as.ts(zz))
  tso <- na.approx(ts(coredata(zz),frequency=(365.25/16),start=start(zz))) ####modified part of bfastts code to account for leap year
  return(tso)
}

##############################################################################################################################
#Old code for bfast function, verions 1.4.4. This version uses the 'strucchange' package for breakpoint
# analysis whereas newer versions, e.g., version 1.5.7, do not and produce different a different result.
bfast_old <- function(Yt, h=0.15, season =c("dummy","harmonic","none"), max.iter = NULL, breaks = NULL, hpc = "none")
{
  season <- match.arg(season)
  ti <- time(Yt)
  f <- frequency(Yt)      # on cycle every f time points (seasonal cycle)
  if(class(Yt)!="ts")
    stop ("Not a time series object")
  ## return value
  output <- list()
  
  # Start the iterative procedure and for first iteration St=decompose result
  St <- stl(Yt, "periodic")$time.series[, "seasonal"]
  Tt <- 0
  
  # seasonal model setup
  if (season=="harmonic") {
    w <- 1/f # f = 23 when freq=23 :-)
    tl <- 1:length(Yt)
    co <- cos(2*pi*tl*w); si <- sin(2*pi*tl*w)
    co2 <- cos(2*pi*tl*w*2);si2 <- sin(2*pi*tl*w*2)
    co3 <- cos(2*pi*tl*w*3);si3 <- sin(2*pi*tl*w*3)
    smod <- Wt ~ co+si+co2+si2+co3+si3
  } else if (season=="dummy") {
    D <- seasonaldummy(Yt)
    D[rowSums(D)==0,] <- -1
    smod <- Wt ~ -1+D
  } else if (season=="none") {
    #print("No seasonal model will be fitted!")
  } else stop("Not a correct seasonal model is selected ('harmonic' or 'dummy') ")
  
  # number/timing of structural breaks in the trend/seasonal component
  Vt.bp <- 0
  Wt.bp <- 0 
  CheckTimeTt <- 1
  CheckTimeSt <- 1
  
  i <- 0
  while ( (!identical(CheckTimeTt,Vt.bp) | !identical(CheckTimeSt,Wt.bp)) & i < max.iter)
  {
    CheckTimeTt <- Vt.bp
    CheckTimeSt <- Wt.bp
    # TREND
    Vt <- Yt-St
    #         p.Vt <- sctest(efp(Vt ~ ti, h=h, type= "OLS-MOSUM"))
    #         if (p.Vt$p.value <=0.05) 
    #         {
    bp.Vt <- breakpoints(Vt ~ ti, h=h,breaks=breaks, hpc = hpc)
    nobp.Vt <- is.na(breakpoints (bp.Vt)[1])
    #         } 
    #         else 
    #         {
    #           nobp.Vt <- TRUE
    #           bp.Vt <- NA       
    #         }
    if (nobp.Vt)
    {
      fm0 <- lm(Vt ~  ti)
      Vt.bp <- 0      # no breaks times
      Tt <- ts(fitted(fm0))     # Data minus trend
      tsp(Tt) <- tsp(Yt)
      ci.Vt <- NA
    } 
    else
    {
      fm1 <- lm(Vt ~ breakfactor(bp.Vt)/ti)
      ci.Vt <- confint(bp.Vt, het.err = FALSE)
      Vt.bp <- ci.Vt$confint[,2]
      Tt <- ts(fitted(fm1))     # Data minus trend
      tsp(Tt) <- tsp(Yt)
    }
    
    # SEASONAL COMPONENT
    if (season=="none") {
      Wt <- 0
      St <- 0
      bp.Wt <- NA; ci.Wt <- NA; nobp.Wt<- TRUE
    } else
    {
      Wt <- Yt-Tt
      #        p.Wt <- sctest(efp(smod, h=h, type= "OLS-MOSUM"))      # preliminary test 
      #        if (p.Wt$p.value <=0.05) # OR statement 
      #        {
      bp.Wt <- breakpoints(smod, h=h,breaks=breaks, hpc = hpc) # Breakpoints in the seasonal component
      nobp.Wt <- is.na(breakpoints (bp.Wt)[1])
      #        } 
      #        else 
      #        {
      #            nobp.Wt <- TRUE
      #            bp.Wt <- NA       
      #        }
      if (nobp.Wt)
      {
        sm0 <- lm(smod)
        St <- ts(fitted(sm0))  #  The fitted seasonal component
        tsp(St) <- tsp(Yt)
        Wt.bp <- 0             # no seasonal breaks
        ci.Wt <- NA
      } 
      else
      {
        if(season=="dummy") sm1 <-lm(Wt ~ -1+D %in% breakfactor(bp.Wt))
        if(season=="harmonic") sm1 <- lm(Wt ~ (co+si+co2+si2+co3+si3) %in% breakfactor(bp.Wt)) 
        St <- ts(fitted(sm1))  #  The fitted seasonal component
        tsp(St) <- tsp(Yt)
        ci.Wt <- confint(bp.Wt, het.err = FALSE)
        Wt.bp <- ci.Wt$confint[,2] 
      }
    }
    i <- i+1
    output[[i]] <- list(Tt=Tt,St=St,Nt=Yt-Tt-St,
                        Vt=Vt, bp.Vt=bp.Vt, Vt.bp=Vt.bp, ci.Vt=ci.Vt,
                        Wt=Wt, bp.Wt=bp.Wt, Wt.bp=Wt.bp, ci.Wt=ci.Wt)
  }
  if (!nobp.Vt) # probably only works well for dummy model!
  {
    Vt.nrbp <- length(bp.Vt$breakpoints)
    co <- coef(fm1) # final fitted trend model
    Mag <- matrix(NA,Vt.nrbp,3)
    for (r in 1:Vt.nrbp) 
    {
      if (r==1) 
        y1 <- co[1]+co[r+Vt.nrbp+1]*ti[Vt.bp[r]]
      else 
        y1 <- co[1]+co[r]+co[r+Vt.nrbp+1]*ti[Vt.bp[r]]
      y2 <- (co[1]+co[r+1])+co[r+Vt.nrbp+2]*ti[Vt.bp[r]+1]
      Mag[r,1] <- y1
      Mag[r,2] <- y2
      Mag[r,3] <- y2-y1   
    }
    index <- which.max(abs(Mag[,3]))
    m.x <- rep(Vt.bp[index],2)
    m.y <- c(Mag[index,1],Mag[index,2]) #Magnitude position
    Magnitude <- Mag[index,3] # Magnitude of biggest change
    Time <- Vt.bp[index]
  } 
  else 
  {
    m.x <- NA; m.y <- NA
    Magnitude <- 0  # if we do not detect a break then the magnitude is zero
    Time <- NA # if we do not detect a break then we have no timing of the break
    Mag <- 0
  }
  return(structure(list(Yt=Yt,output=output,nobp=list(Vt=nobp.Vt,Wt=nobp.Wt),Magnitude=Magnitude,Mags=Mag,
                        Time=Time,jump=list(x=ti[m.x],y=m.y)),class="bfast"))  
}

###############################################################################################################
########################################################
####select sites based on following names: BASN CALI COLL EAST GRAV IBPE NORT RABB SAND SMAL SUMM TAYL TOBO WELL WEST
BASN <- JRN.point[JRN.point@data$SITE %in% c("BASN"), ]
CALI <- JRN.point[JRN.point@data$SITE %in% c("CALI"), ]
COLL <- JRN.point[JRN.point@data$SITE %in% c("COLL"), ]
EAST <- JRN.point[JRN.point@data$SITE %in% c("EAST"), ]
GRAV <- JRN.point[JRN.point@data$SITE %in% c("GRAV"), ]
IBPE <- JRN.point[JRN.point@data$SITE %in% c("IBPE"), ]
NORT <- JRN.point[JRN.point@data$SITE %in% c("NORT"), ]
RABB <- JRN.point[JRN.point@data$SITE %in% c("RABB"), ]
SAND <- JRN.point[JRN.point@data$SITE %in% c("SAND"), ]
SMAL <- JRN.point[JRN.point@data$SITE %in% c("SMAL"), ]
SUMM <- JRN.point[JRN.point@data$SITE %in% c("SUMM"), ]
TAYL <- JRN.point[JRN.point@data$SITE %in% c("TAYL"), ]
TOBO <- JRN.point[JRN.point@data$SITE %in% c("TOBO"), ]
WELL <- JRN.point[JRN.point@data$SITE %in% c("WELL"), ]
WEST <- JRN.point[JRN.point@data$SITE %in% c("WEST"), ]

#######extract plot points from raster brick: BASN CALI COLL EAST GRAV IBPE NORT RABB SAND SMAL SUMM TAYL TOBO WELL WEST

###############NDVI_mask
BASN.ndvi.points<-extract(JRN_NDVI_mask, BASN) #49 points from BASN plot
CALI.ndvi.points<-extract(JRN_NDVI_mask, CALI)
COLL.ndvi.points<-extract(JRN_NDVI_mask, COLL)
EAST.ndvi.points<-extract(JRN_NDVI_mask, EAST)
GRAV.ndvi.points<-extract(JRN_NDVI_mask, GRAV)
IBPE.ndvi.points<-extract(JRN_NDVI_mask, IBPE)
NORT.ndvi.points<-extract(JRN_NDVI_mask, NORT)
RABB.ndvi.points<-extract(JRN_NDVI_mask, RABB)
SAND.ndvi.points<-extract(JRN_NDVI_mask, SAND)
SMAL.ndvi.points<-extract(JRN_NDVI_mask, SMAL)
SUMM.ndvi.points<-extract(JRN_NDVI_mask, SUMM)
TAYL.ndvi.points<-extract(JRN_NDVI_mask, TAYL)
TOBO.ndvi.points<-extract(JRN_NDVI_mask, TOBO)
WELL.ndvi.points<-extract(JRN_NDVI_mask, WELL)
WEST.ndvi.points<-extract(JRN_NDVI_mask, WEST)

#mean across 49 points
BASN.ndvi.m <- colMeans(BASN.ndvi.points)
CALI.ndvi.m <- colMeans(CALI.ndvi.points)
COLL.ndvi.m <- colMeans(COLL.ndvi.points)
EAST.ndvi.m <- colMeans(EAST.ndvi.points)
GRAV.ndvi.m <- colMeans(GRAV.ndvi.points)
IBPE.ndvi.m <- colMeans(IBPE.ndvi.points)
NORT.ndvi.m <- colMeans(NORT.ndvi.points)
RABB.ndvi.m <- colMeans(RABB.ndvi.points)
SAND.ndvi.m <- colMeans(SAND.ndvi.points)
SMAL.ndvi.m <- colMeans(SMAL.ndvi.points)
SUMM.ndvi.m <- colMeans(SUMM.ndvi.points)
TAYL.ndvi.m <- colMeans(TAYL.ndvi.points)
TOBO.ndvi.m <- colMeans(TOBO.ndvi.points)
WELL.ndvi.m <- colMeans(WELL.ndvi.points)
WEST.ndvi.m <- colMeans(WEST.ndvi.points)
#########################################################

############################################################################################
### BASN
############################################################################################
## apply on MODIS composite point
ndvi.date <- (orgTime(ndvi,nDays=16,begin="2000049", end="2014177", pillow=0))$inputLayerDates  #, pillow=90  #original model run until end="2014177" but later cropped for figure
basn.ndvi <- bfastts.na(as.numeric(BASN.ndvi.m)/10000, ndvi.date, type = c("16-day")) 
dev.off()
plot(basn.ndvi)

#bfast function
#1. Create the time series data object from the MODIS imagery time series
basn.ndvi.fit<-stl(basn.ndvi, s.window=8)
colnames(basn.ndvi.fit$time.series)[2]<-"abrupt"
plot(basn.ndvi.fit)

basn.fit <- bfast_old(basn.ndvi, season="harmonic",h=0.15, max.iter=10)
plot(basn.fit) 

############################################################################################
### CALI
############################################################################################
## apply on one pixel for testing 
#ndvi.date <- (orgTime(ndvi,nDays=16))$inputLayerDates
cali.ndvi <- bfastts.na(as.numeric(CALI.ndvi.m)/10000, ndvi.date, type = c("16-day")) 
dev.off()
plot(cali.ndvi)

#bfast function
#1. Create the time series data object from the MODIS imagery time series
cali.ndvi.fit<-stl(cali.ndvi, s.window=8)
colnames(cali.ndvi.fit$time.series)[2]<-"abrupt"

#bfast model
cali.fit <- bfast_old(cali.ndvi, season="harmonic", max.iter=10) 

############################################################################################
### COLL
############################################################################################
## apply on one pixel for testing 
#ndvi.date <- (orgTime(ndvi,nDays=16))$inputLayerDates
coll.ndvi <- bfastts.na(as.numeric(COLL.ndvi.m)/10000, ndvi.date, type = c("16-day")) 
dev.off()
plot(coll.ndvi)

#bfast function
#1. Create the time series data object from the MODIS imagery time series
coll.ndvi.fit<-stl(coll.ndvi, s.window=8)
colnames(coll.ndvi.fit$time.series)[2]<-"abrupt"

#bfast model
coll.fit <- bfast_old(coll.ndvi, season="harmonic", max.iter=10) 

############################################################################################
### EAST
############################################################################################
## apply on one pixel for testing 
#ndvi.date <- (orgTime(ndvi,nDays=16))$inputLayerDates
east.ndvi <- bfastts.na(as.numeric(EAST.ndvi.m)/10000, ndvi.date, type = c("16-day")) 
dev.off()
plot(east.ndvi)

#bfast function
#1. Create the time series data object from the MODIS imagery time series
east.ndvi.fit<-stl(east.ndvi, s.window=18)
colnames(east.ndvi.fit$time.series)[2]<-"abrupt"

#bfast model
east.fit <- bfast_old(east.ndvi, season="harmonic",h=.17, max.iter=10) 

############################################################################################
### GRAV
############################################################################################
## apply on one pixel for testing 
#ndvi.date <- (orgTime(ndvi,nDays=16))$inputLayerDates
grav.ndvi <- bfastts.na(as.numeric(GRAV.ndvi.m)/10000, ndvi.date, type = c("16-day")) 
dev.off()
plot(grav.ndvi)

#bfast function
#1. Create the time series data object from the MODIS imagery time series
grav.ndvi.fit<-stl(grav.ndvi, s.window=8)
colnames(grav.ndvi.fit$time.series)[2]<-"abrupt"

#bfast model
grav.fit <- bfast_old(grav.ndvi, season="harmonic", max.iter=10) 

############################################################################################
### IBPE
############################################################################################
## apply on one pixel for testing 
#ndvi.date <- (orgTime(ndvi,nDays=16))$inputLayerDates
ibpe.ndvi <- bfastts.na(as.numeric(IBPE.ndvi.m)/10000, ndvi.date, type = c("16-day")) 
dev.off()
plot(ibpe.ndvi)

#bfast function
#1. Create the time series data object from the MODIS imagery time series
ibpe.ndvi.fit<-stl(ibpe.ndvi, s.window=8)
colnames(ibpe.ndvi.fit$time.series)[2]<-"abrupt"

#bfast model
ibpe.fit <- bfast_old(ibpe.ndvi, season="harmonic", max.iter=10) 

############################################################################################
### NORT
############################################################################################
## apply on one pixel for testing 
#ndvi.date <- (orgTime(ndvi,nDays=16))$inputLayerDates
nort.ndvi <- bfastts.na(as.numeric(NORT.ndvi.m)/10000, ndvi.date, type = c("16-day")) 
dev.off()
plot(nort.ndvi)

#bfast function
#1. Create the time series data object from the MODIS imagery time series
nort.ndvi.fit<-stl(nort.ndvi, s.window=8)
colnames(nort.ndvi.fit$time.series)[2]<-"abrupt"

#bfast model
nort.fit <- bfast_old(nort.ndvi, season="harmonic", max.iter=10) 

############################################################################################
### RABB
############################################################################################
## apply on one pixel for testing 
#ndvi.date <- (orgTime(ndvi,nDays=16))$inputLayerDates
rabb.ndvi <- bfastts.na(as.numeric(RABB.ndvi.m)/10000, ndvi.date, type = c("16-day")) 
dev.off()
plot(rabb.ndvi)

#bfast function
#1. Create the time series data object from the MODIS imagery time series
rabb.ndvi.fit<-stl(rabb.ndvi, s.window=8)
colnames(rabb.ndvi.fit$time.series)[2]<-"abrupt"

#bfast model
rabb.fit <- bfast_old(rabb.ndvi, season="harmonic", max.iter=10) 

############################################################################################
### SAND
############################################################################################
## apply on one pixel for testing 
#ndvi.date <- (orgTime(ndvi,nDays=16))$inputLayerDates
sand.ndvi <- bfastts.na(as.numeric(SAND.ndvi.m)/10000, ndvi.date, type = c("16-day")) 
dev.off()
plot(sand.ndvi)

#bfast function
#1. Create the time series data object from the MODIS imagery time series
sand.ndvi.fit<-stl(sand.ndvi, s.window=8)
colnames(sand.ndvi.fit$time.series)[2]<-"abrupt"

#bfast model
sand.fit <- bfast_old(sand.ndvi, season="harmonic", max.iter=10) 

############################################################################################
### SMAL
############################################################################################
## apply on one pixel for testing 
#ndvi.date <- (orgTime(ndvi,nDays=16))$inputLayerDates
smal.ndvi <- bfastts.na(as.numeric(SMAL.ndvi.m)/10000, ndvi.date, type = c("16-day")) 
dev.off()
plot(smal.ndvi)

#bfast function
#1. Create the time series data object from the MODIS imagery time series
smal.ndvi.fit<-stl(smal.ndvi, s.window=8)
colnames(smal.ndvi.fit$time.series)[2]<-"abrupt"

#bfast model
smal.fit <- bfast_old(smal.ndvi, season="harmonic", max.iter=10) 

############################################################################################
### SUMM
############################################################################################
## apply on one pixel for testing 
#ndvi.date <- (orgTime(ndvi,nDays=16))$inputLayerDates
summ.ndvi <- bfastts.na(as.numeric(SUMM.ndvi.m)/10000, ndvi.date, type = c("16-day")) 
dev.off()
plot(summ.ndvi)

#bfast function
#1. Create the time series data object from the MODIS imagery time series
summ.ndvi.fit<-stl(summ.ndvi, s.window=8)
colnames(summ.ndvi.fit$time.series)[2]<-"abrupt"

#bfast model
summ.fit <- bfast_old(summ.ndvi, season="harmonic", max.iter=10) 

############################################################################################
### TAYL
############################################################################################
## apply on one pixel for testing 
#ndvi.date <- (orgTime(ndvi,nDays=16))$inputLayerDates
tayl.ndvi <- bfastts.na(as.numeric(TAYL.ndvi.m)/10000, ndvi.date, type = c("16-day")) 
dev.off()
plot(tayl.ndvi)

#bfast function
#1. Create the time series data object from the MODIS imagery time series
tayl.ndvi.fit<-stl(tayl.ndvi, s.window=8)
colnames(tayl.ndvi.fit$time.series)[2]<-"abrupt"

#bfast model
tayl.fit <- bfast_old(tayl.ndvi, season="harmonic", max.iter=10) 

############################################################################################
### TOBO
############################################################################################
## apply on one pixel for testing 
#ndvi.date <- (orgTime(ndvi,nDays=16))$inputLayerDates
tobo.ndvi <- bfastts.na(as.numeric(TOBO.ndvi.m)/10000, ndvi.date, type = c("16-day")) 
dev.off()
plot(tobo.ndvi)

#bfast function
#1. Create the time series data object from the MODIS imagery time series
tobo.ndvi.fit<-stl(tobo.ndvi, s.window=8)
colnames(tobo.ndvi.fit$time.series)[2]<-"abrupt"

#bfast model
tobo.fit <- bfast_old(tobo.ndvi, season="harmonic", max.iter=10) 

############################################################################################
### WELL
############################################################################################
## apply on one pixel for testing 
#ndvi.date <- (orgTime(ndvi,nDays=16))$inputLayerDates
well.ndvi <- bfastts.na(as.numeric(WELL.ndvi.m)/10000, ndvi.date, type = c("16-day")) 
dev.off()
plot(well.ndvi)

#bfast function
#1. Create the time series data object from the MODIS imagery time series
well.ndvi.fit<-stl(well.ndvi, s.window=8)
colnames(well.ndvi.fit$time.series)[2]<-"abrupt"

#bfast model
well.fit <- bfast_old(well.ndvi, season="harmonic", max.iter=10) 

############################################################################################
### WEST
############################################################################################
## apply on one pixel for testing 
#ndvi.date <- (orgTime(ndvi,nDays=16))$inputLayerDates
west.ndvi <- bfastts.na(as.numeric(WEST.ndvi.m)/10000, ndvi.date, type = c("16-day")) 
dev.off()
plot(west.ndvi)

#bfast function
#1. Create the time series data object from the MODIS imagery time series
west.ndvi.fit<-stl(west.ndvi, s.window=8)
colnames(west.ndvi.fit$time.series)[2]<-"abrupt"

#bfast model
west.fit <- bfast_old(west.ndvi, season="harmonic", max.iter=10) 

####################################################################################################
#Function to extract Trend and Seasonal breaks
bfast.breaks <- function(bfast,date){
  
  trend <-list(list())
  if(is.na(bfast$output[[length(bfast$output)]]$bp.Vt$breakpoints)[1]==TRUE){
    trend <-NA
  } else {
    for(i in 1:length(bfast$output[[length(bfast$output)]]$bp.Vt$breakpoints)){
      trend[[i]] <- date[bfast$output[[length(bfast$output)]]$bp.Vt$breakpoints[i]]
      
    }
  }



  trend<-data.frame(as.character(as.Date(unlist(trend))))
  names(trend)<- c("Trend")
  season <-list(list())
  if(is.na(bfast$output[[length(bfast$output)]]$bp.Wt$breakpoints)[1]==TRUE){
    season <-NA
  }  else {
    for(i in 1:length(bfast$output[[length(bfast$output)]]$bp.Wt$breakpoints)){
      season[[i]] <- date[bfast$output[[length(bfast$output)]]$bp.Wt$breakpoints[i]]
      
    }
  }


  season<-data.frame(as.character(as.Date(unlist(season))))
  names(season)<- c("Season")
  breaks <- list(trend, season)
  names(breaks) <- c("Trend", "Season")
  return(breaks)
}


basn.breaks <- bfast.breaks(basn.fit, ndvi.date)
coll.breaks <- bfast.breaks(coll.fit, ndvi.date)
east.breaks <- bfast.breaks(east.fit, ndvi.date)
grav.breaks <- bfast.breaks(grav.fit, ndvi.date)
ibpe.breaks <- bfast.breaks(ibpe.fit, ndvi.date)
nort.breaks <- bfast.breaks(nort.fit, ndvi.date)
rabb.breaks <- bfast.breaks(rabb.fit, ndvi.date)
sand.breaks <- bfast.breaks(sand.fit, ndvi.date)
smal.breaks <- bfast.breaks(smal.fit, ndvi.date)
summ.breaks <- bfast.breaks(summ.fit, ndvi.date)
tayl.breaks <- bfast.breaks(tayl.fit, ndvi.date)
tobo.breaks <- bfast.breaks(tobo.fit, ndvi.date)
well.breaks <- bfast.breaks(well.fit, ndvi.date)
west.breaks <- bfast.breaks(west.fit, ndvi.date)

basn.breaks 
coll.breaks 
east.breaks 
grav.breaks 
ibpe.breaks 
nort.breaks
rabb.breaks
sand.breaks
smal.breaks
summ.breaks
tayl.breaks
tobo.breaks
well.breaks
west.breaks


###########################################################################
#function to extract CI interval dates
CI.interval.extract <- function(bfast,date, site){
      CI.trend <- data.frame(date[bfast$output[[length(bfast$output)]]$ci.Vt$confint[,1]],date[bfast$output[[length(bfast$output)]]$ci.Vt$confint[,2]],date[bfast$output[[length(bfast$output)]]$ci.Vt$confint[,3]])
      colnames(CI.trend) <- c("CI.start","Break", "CI.end")
      rownames(CI.trend) <- c(paste0(site,rep("T",nrow(CI.trend)),seq(1,nrow(CI.trend),1)))
      
      if(is.na(bfast$output[[length(bfast$output)]]$ci.Wt)[1]==TRUE){
      CI.final <- CI.trend
      } else {
      CI.season <- data.frame(date[bfast$output[[length(bfast$output)]]$ci.Wt$confint[,1]],date[bfast$output[[length(bfast$output)]]$ci.Wt$confint[,2]],date[bfast$output[[length(bfast$output)]]$ci.Wt$confint[,3]])
      colnames(CI.season) <- c("CI.start","Break", "CI.end")
      rownames(CI.season) <- c(paste0(site,rep("S",nrow(CI.season)),seq(1,nrow(CI.season),1)))
      CI.final <- rbind(CI.trend, CI.season)
      }
      return(CI.final)
  }


basn.CI.bounds <- CI.interval.extract( basn.fit,ndvi.date, "basn") 
cali.CI.bounds <- CI.interval.extract( cali.fit,ndvi.date, "cali")
coll.CI.bounds <- CI.interval.extract( coll.fit,ndvi.date, "coll") 
east.CI.bounds <- CI.interval.extract( east.fit,ndvi.date, "east") 
grav.CI.bounds <- CI.interval.extract( grav.fit,ndvi.date, "grav") 
ibpe.CI.bounds <- CI.interval.extract( ibpe.fit,ndvi.date, "ibpe") 
nort.CI.bounds <- CI.interval.extract( nort.fit,ndvi.date, "nort")
rabb.CI.bounds <- CI.interval.extract( rabb.fit,ndvi.date, "rabb")
sand.CI.bounds <- CI.interval.extract( sand.fit,ndvi.date, "sand")
smal.CI.bounds <- CI.interval.extract( smal.fit,ndvi.date, "smal")
summ.CI.bounds <- CI.interval.extract( summ.fit,ndvi.date, "summ")
tayl.CI.bounds <- CI.interval.extract( tayl.fit,ndvi.date, "tayl")
tobo.CI.bounds <- CI.interval.extract( tobo.fit,ndvi.date, "tobo")
well.CI.bounds <- CI.interval.extract( well.fit,ndvi.date, "well")
west.CI.bounds <- CI.interval.extract( west.fit,ndvi.date, "west")

CI.bounds <- rbind(basn.CI.bounds,  
                    coll.CI.bounds,  
                    cali.CI.bounds, 
                    east.CI.bounds,  
                    grav.CI.bounds,  
                    ibpe.CI.bounds,  
                    nort.CI.bounds, 
                    rabb.CI.bounds, 
                    sand.CI.bounds, 
                    smal.CI.bounds, 
                    summ.CI.bounds, 
                    tayl.CI.bounds, 
                    tobo.CI.bounds, 
                    well.CI.bounds, 
                    west.CI.bounds)

write.csv(CI.bounds, "JRN_CI_bounds.csv")

#############################################################################################################################################
###  plot.bfast function   ##################################################################################################################
# bfast uses the embeded function plot.bfast to plot results from timeseries decomposition. 
# To bring up code for plot.bfast: bfast:::plot.bfast
# I have modidied sections of this to allow plotting of NPP and precip data
#############################################################################################################################################

plot.bfast <- function (x, type = c("components", "all", "data", "seasonal", "trend", "noise"), sim = NULL, largest = FALSE, main, ANOVA = FALSE, npp, precip, ...) 
{
  type <- match.arg(type)
  realdata <- is.null(sim)
  Trend.bp <- !x$nobp$Vt
  if (type == "largest" & !Trend.bp) 
    stop("No trend breakpoints")
  title <- !missing(main)
  niter <- length(x$output)
  out <- x$output[[niter]]
  Tt <- out$Tt
  St <- out$St
  noise <- out$Nt
  npp <- npp
  precip <- precip
  if (type == "data") {
    if (!title) 
      main <- "Yt"
    plot(x$Yt, main = main, ...)
  }
  else if (type == "components") {
    ft <- cbind(seasonal = out$St, trend = out$Tt, remainder = out$Nt)
    tsp(ft) <- tsp(x$Yt)
    ft <- list(time.series = ft)
    if (!title) 
      main <- paste("no. iterations to estimate breakpoints:", 
                    niter)
    if (ANOVA == TRUE) {
      seasonal(ft, out, npp, precip, sim = sim, main = main, fit = x)
    }
    else {
      seasonal(ft, out, npp, precip,  sim = sim, main = main)
    }
  }
  else if (type == "noise") {
    require(forecast)
    if (!title) 
      main <- "Noise component"
    tsdisplay(noise, main = main)
  }
  else {
    if (type == "all") {
      idx <- 1:niter
      opar <- par(mfrow = c(2, niter))
    }
    else idx <- niter
    for (i in idx) {
      out <- x$output[[i]]
      if (type != "seasonal") {
        if (type == "trend" & !title) 
          main <- "Trend component"
        else if (!title) 
          main <- paste("Iteration ", i, ": Trend", sep = "")
        plot(out$Vt, main = main, ylab = "Vt", ...)
        lines(out$Tt, col = 4)
        if (Trend.bp) {
          lines(out$bp.Vt)
          lines(out$ci.Vt)
          legend("topright", paste("Time of BP(s)", paste(out$Vt.bp, 
                                                          collapse = ",")), col = 2)
        }
        if (!realdata) {
          lines(sim$time.series[, "abrupt"], col = 1, 
                lty = 2)
          legend("bottomleft", c("estimated", "simulated"), 
                 lty = c(1, 2), col = 1)
        }
        if (largest) {
          legend("bottomright", c("Magnitude of most sign change"), 
                 lty = c(1), col = 6)
          lines(x$jump, col = 6)
          points(x$jump, pch = 14, cex = 1, col = 6)
        }
      }
      if (type != "trend") {
        if (type == "seasonal" & !title) 
          main <- "Seasonal component"
        else if (!title) 
          main <- paste("Iteration ", i, ": Seasonal", 
                        sep = "")
        plot(out$Wt, main = main, ylab = "Wt", ...)
        lines(out$St, col = 2)
        Seas.bp <- !x$nobp$Wt
        if (Seas.bp) {
          lines(out$bp.Wt)
          lines(out$ci.Wt)
          legend("topright", paste("Time of BP(s)", paste(out$Wt.bp, 
                                                          collapse = ",")), col = 2)
        }
        if (!realdata) {
          lines(sim$time.series[, "seasonal"], col = 1, 
                lty = 2)
          legend("bottomleft", c("first run seasonality", 
                                 "first run estimated", "simulated"), lty = c(1, 
                                                                              1, 2), col = c(1, 2, 1))
        }
      }
    }
    if (type == "all") 
      par(opar)
  }
}





##################################################################################
#Plotting all components of the bfast output uses the seasonal function. I have modified this function to allow plotting of NPP and precip.
##Seasonal function


seasonal <- function (x, out, npp, precip, sim = NULL, labels = c(colnames(X), "Precipitation", "NPP"), set.pars = list(tck = -0.01, 
                                                                                                                        mar = c(0, 6, 0, 6), oma = c(6, 0, 4, 0)), main = NULL, range.bars = FALSE, 
                      col.range = "light gray", fit = NULL) 
{
  layout(matrix(c(1, 2, 3, 4,5,6), 6, 1, byrow = TRUE), heights = c(1, 
                                                                    0.75, 1, 0.75, 1, 1), TRUE) #modify trend plot height, change from 1.5 to 1
  sers <- x$time.series
  npp<-npp
  precip<-precip
  ncomp <- ncol(sers)
  data <- drop(sers %*% rep(1, ncomp))
  X <- cbind(data, sers)
  X<-ts(as.zoo(X)[1:291], start=2000.130, frequency=22.82812)
  colnames(X) <- c("Data", "Seasonal", "Trend", "Resid") #change plot x axis names from c("Yt", "St", "Tt", "et")
  nplot <- ncomp + 3
  if (range.bars) 
    mx <- min(apply(rx <- apply(X, 2, range), 2, diff))
  if (length(set.pars)) {
    oldpar <- do.call("par", as.list(names(set.pars)))
    on.exit(par(oldpar))
    do.call("par", set.pars)
  }
  if (is.null(fit) == F) {
    niter <- length(fit$output)
    out <- fit$output[[niter]]
    out_ANOVA <- array()
    out_breakdates <- array()
    if (out$Vt.bp[1] > 0) {
      breaks <- length(out$Vt.bp)
    }
    else {
      breaks <- 0
    }
    if (breaks > 0) {
      breakdates <- out$Vt.bp
      coefs <- coef(out$bp.Vt)
      sl <- coefs[, 2]
    }
    TS_anova <- fit$Yt - out$St
    dataframe <- data.frame(TIME = c(1:length(fit$Yt)), DATA = TS_anova)
    for (m in 1:(breaks + 1)) {
      startpoint <- if (m == 1) 
        1
      else breakdates[[m - 1]]
      endpoint <- if (m == (breaks + 1)) 
        length(fit$Yt)
      else breakdates[m] - 1
      df2 <- dataframe[startpoint:endpoint, ]
      model <- lm(DATA ~ TIME, data = df2)
      modelAnova <- anova(model)
      out_ANOVA[m] <- modelAnova$Pr[1]
      if (breaks == 0) {
        sl <- model$coefficients[2]
      }
    }
  }
  
  for (i in 1:nplot) {
    if (i == 6) {
      par(mar = c(0, 6, 0, 6))
    }
    if (i<5) {
      plot(X[, i], col = if (i == 1) 
        "black"
        else "red", ylim = if (i == 1 | i == 3) 
          range(X[, 1])
        else range(X[, i], sim$time.series[, i - 1], na.rm = TRUE), 
        type = if (i != 4) 
          "l"
        else "h", xlab = "", ylab = "", axes = FALSE) }
    
    if (i==5) {
      plot(precip, ann=FALSE, xaxt="n", yaxt="n", )
      
    }
    if (i==6) {
      plot(npp, type="b", ann=FALSE, yaxt="n", ylim=c(0,round(max(npp)+15))) 
      
    }    
    
    if (range.bars) {
      dx <- 1/64 * diff(ux <- par("usr")[1:2])
      y <- mean(rx[, i])
      rect(ux[2] - dx, y + mx/2, ux[2] - 0.4 * dx, y - 
             mx/2, col = col.range, xpd = TRUE)
    }
    if (i == 1 && !is.null(main)) {
      mtext(main, side = 3, font = 1.5, line = 1, cex = 1.1)
      if (!is.null(sim)) {
        lines(X[, i + 1] + X[, i + 2], col = "red", type = "l")
        legend("topleft", c("input", "estimated seasonal + trend "), 
               col = c("black", "red"), lty = 1, cex=0.7)
      }
    }
    if (i == 2) {
      lines(sim$time.series[, "seasonal"], col = "black")
      lines(out$bp.Wt)
      lines(out$ci.Wt)
    }
    if (i == 3) {
      lines(sim$time.series[, "abrupt"], col = "black")
      lines(out$bp.Vt)
      lines(out$ci.Vt)
      if (is.null(fit) == FALSE) {
        for (m in 1:(breaks + 1)) {
          x_coor <- out$bp.Wt$datatsp[[1]]
          if (m > 1) {
            x_coor <- x_coor + breakdates[[m - 1]]/frequency(fit$Yt)
          }
          y_range <- range(X[, 1])
          y_sl <- y_range[2] - (y_range[2] - y_range[1])/10
          y_Pr <- y_range[2] - (y_range[2] - y_range[1])/5
          beta <- formatC(sl[m], format = "f", digits = 3)
          text(x_coor, y_sl, bquote(beta == .(beta)), 
               pos = 4, cex=0.8)
          Pr <- formatC(out_ANOVA[m], format = "f", digits = 3)
          text(x_coor, y_Pr, bquote(p == .(Pr)), pos = 4, cex=0.8)
        }
      }
    }
    if (i == 4) {
      abline(h = 0)
      lines(sim$time.series[, "remainder"], col = "black")
    }
    if (i == 5) {
      lines(out$bp.Wt, col="blue", lty=2)    #this doesn't work
      lines(out$bp.Vt, col="green", lty=2)   #this doesn't work
    }
    if (i == 6) {
      lines(out$bp.Wt, col="blue", lty=2)   #this doesn't work
      lines(out$bp.Vt, col="green", lty=2)  #this doesn't work
    }
    
    
    box()
    right <- i%%2 == 0
    par(las=1)
    axis(2, labels = !right)
    axis(4, labels = right)
    #axis(1, labels = i == nplot)
    par(las=0)
    mtext(labels[i], side = 2, 3, cex=0.7)
  }
  mtext("Time", side = 1, line = 3, cex=0.8)
  invisible()
  if (is.null(fit) == FALSE) {
    return(data.frame(slope = sl, prob = out_ANOVA))
  }
  layout(matrix(1))
}

############################################################################################################################################
### Plotting NPP data by plant functional tyep ###############
############################################################################################################################################

jrn_func <- read.csv("G:/ARS_Static_Data/Jornada/JORNADA_NPP/biomass_funcgrp_site_2000_2012.csv", header=T, sep=',')

BASN.jrn_func <- jrn_func[jrn_func$SITE %in% c("BASN"), ]
CALI.jrn_func  <- jrn_func[jrn_func$SITE %in% c("CALI"), ]
COLL.jrn_func  <- jrn_func[jrn_func$SITE %in% c("COLL"), ]
EAST.jrn_func  <- jrn_func[jrn_func$SITE %in% c("EAST"), ]
GRAV.jrn_func  <- jrn_func[jrn_func$SITE %in% c("GRAV"), ]
IBPE.jrn_func  <- jrn_func[jrn_func$SITE %in% c("IBPE"), ]
NORT.jrn_func  <- jrn_func[jrn_func$SITE %in% c("NORT"), ]
RABB.jrn_func  <- jrn_func[jrn_func$SITE %in% c("RABB"), ]
SAND.jrn_func  <- jrn_func[jrn_func$SITE %in% c("SAND"), ]
SMAL.jrn_func  <- jrn_func[jrn_func$SITE %in% c("SMAL"), ]
SUMM.jrn_func  <- jrn_func[jrn_func$SITE %in% c("SUMM"), ]
TAYL.jrn_func  <- jrn_func[jrn_func$SITE %in% c("TAYL"), ]
TOBO.jrn_func  <- jrn_func[jrn_func$SITE %in% c("TOBO"), ]
WELL.jrn_func  <- jrn_func[jrn_func$SITE %in% c("WELL"), ]
WEST.jrn_func  <- jrn_func[jrn_func$SITE %in% c("WEST"), ]


######################################################################################################################################
#######################################################################################################################################
#####################################################################
############  Example plotting for BASN   ###############################################
#This code was used to create the NPP functional plots for the other 14 sites

#BASN Absolute Contribution
BASN.jrn_func$median_date <- as.Date(as.character(BASN.jrn_func$median_date))

ylim <- c(0, max(BASN.jrn_func$bio_annual + BASN.jrn_func$bio_pforb +  BASN.jrn_func$bio_pgrass +  BASN.jrn_func$bio_shrub))
xx <- c(BASN.jrn_func$median_date, rev(BASN.jrn_func$median_date))

basn.ab.ann <-BASN.jrn_func$bio_annual
basn.ab.pforb <- BASN.jrn_func$bio_pforb
basn.ab.grass <- BASN.jrn_func$bio_pgrass
basn.ab.shrub <- BASN.jrn_func$bio_shrub


######################################################order grass, shrub,  ann, pforb
CairoPDF("BASN_AbsCont_FunGrps.pdf",
         7, 5, bg = "transparent", pointsize=14 )
par(mar=c(6.5,4,3.5,4)+.1, las=1)
#order herb,  shrub, grass, ann, pforb
yy_grass <- c(rep(0, nrow(BASN.jrn_func)), rev(basn.ab.grass))
plot(x=BASN.jrn_func$median_date, y=basn.ab.grass, ylim=ylim, col='red', type='l', xaxt='n',
     ylab='', xlab='Date', main='BASN')
mtext(expression(paste('Mean Biomass (g ',m^2,')')), side=2, line=2.5, las=0)

polygon(xx, yy_grass, col="#00A600FF")

yy_shrub <- c( basn.ab.grass , rev(basn.ab.shrub) + rev(basn.ab.grass))
polygon(xx, yy_shrub, col="#E6E600FF")

yy_ann <- c(basn.ab.grass +basn.ab.shrub,  rev(basn.ab.grass) + rev(basn.ab.shrub) + rev(basn.ab.ann))
polygon(xx, yy_ann, col="#EBB25EFF")

yy_pforb <-  c( basn.ab.grass +basn.ab.shrub + basn.ab.ann,  rev(basn.ab.grass) + rev(basn.ab.shrub) + rev(basn.ab.ann) + rev(basn.ab.pforb))
polygon(xx, yy_pforb, col="#F0C9C0FF")

# x axis labels
labdates <- as.Date('1970-01-01') + min(BASN.jrn_func$median_date):max(BASN.jrn_func$median_date)
labdates <- labdates[ format(labdates, '%m') == '01']
labdates <- labdates[ format(labdates, '%d') == '01']
labdates <- unique(c(labdates, max(BASN.jrn_func$median_date)))

labnames <- format(labdates, '%Y')
axis(1, at=labdates, labels=labnames)

# # black lines first day of month
# for(a in labdates[ format(labdates, '%y') == '01']){
#   abline(v=a)
# }

# make the first 3 lower to avoid our legend "#00A600FF" "#63C600FF"  "#E6E600FF" "#EBB25EFF" "#F2F2F2FF"
cnt <- 0

# set alpha = 80 so it's relatively transparent
color <- rgb(190, 190, 190, alpha=80, maxColorValue=255)

# check if the first data point is a Sunday
BASN.jrn_func$mon <- format( BASN.jrn_func$median_date, '%m' )

# plot strip along sept biomass max
for( i in BASN.jrn_func[ BASN.jrn_func$mon == '09', ]$median_date ){
  rect(xleft=i-30, xright=i+30, ybottom=-1000, ytop=1.1*ylim[2], density=100, col=color)
  cnt <- cnt + 1
}
par(xpd=TRUE)
#order herb,  shrub, grass, ann, pforb
legend(x=BASN.jrn_func[1,]$median_date -100 , y=ylim[2]+5, c('Grass', 'Shrubs','Annuals', 'Forbs'), fill=c("#00A600FF",  "#E6E600FF", "#EBB25EFF" ,"#F0C9C0FF"), cex=0.67, bg = "white")
par(xpd=F)
dev.off()

######################################################################################################################################
#######################################################################################################################################


