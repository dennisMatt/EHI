#Code to tell R to only use the terra installation of PROJ (for performing projections)
plib<-Sys.getenv("PROJ_LIB")
prj<-system.file("proj",package="terra")[1]
Sys.setenv("PROJ_LIB"=prj)

#The following libraries are needed to run this script.

library(sf)
library(stars)
library(terra)
library(lwgeom)
library(raster)
library(gdistance)


###
# Linear interpolation function 'Lerp'
#   rescale value c from range a-b to range y-z
# This function was created to make sure that the sigmoid curve has a max of one
# It makes sure that values on the y axis in Fig.3 in the paper scale from 0-1


lerp <- function(c, a, b, y, z) {
  (c - a) * (z - y) / (b - a) + y
}



###
# Logistic function used to link matrix patch areas to distance decay of the edge effect
#   https://stackoverflow.com/questions/55725139/fit-sigmoid-function-s-shape-curve-to-data-using-python
# 
# Params:
#   `x` is the value for which you want to return the y
#   `L` is responsible for scaling the output range from [0,1] to [0,L]
#   `b` adds bias to the output and changes its range from [0,L] to [b,L+b]
#   `k` is responsible for scaling the input, which remains in (-inf,inf)
#   `x0` is the point in the middle of the Sigmoid, i.e. the point where Sigmoid should originally output the value 1/2 [since if x=x0, we get 1/(1+exp(0)) = 1/2].
# 
### Put simply, this function creates the sigmoid curve then we coerce this curve to scale 0-1 below

logistic <- function(x, L=1.0, k, x0) {
  L / (1.0 + exp(-k*(x-x0)))
}

###
# This makes a call to the logistic function and then scales the result to enforce a value of 1 for the 
#   stated maximum value ([max]), 0 for 0 and 0.5 where x = x0 (maximum val/2). This is achieved using a 
#   linear  interpolation function to scale either from 0-0.5 (where x<0.5) or 0.5-[max] (otherwise). This 
#   is necessary as scaling from 0-max can displace the case so that y!=0.5 where x=x0 if the adjustment 
#   required at x=0 and x=[max] are not equal.
# 
# Params:
#   `x` is the value for which you want to return the y
#   `fullEdgeEffectArea` is the x value that should be coerced to 1  

# 
###Finally, this function coerces the reults of the 'logistic' function uing the 'lerp' function above
coercedLogistic <- function(x, fullEdgeEffectArea) {
  
  # if outside of the parameter range, just return accordingly
  if (x < 0) stop("patch areas should always be positive!")
  if (x >= fullEdgeEffectArea) return(1)
  
  # work out order of magnitude and convert to scale factor
  #(this is arbitrary for now but changes according to fullEdgeEffectArea)
  
  k <- 1 / (10^(floor(log(fullEdgeEffectArea, 10))-1))
  
  # x0 is the midpoint of the scale
  x0 <- fullEdgeEffectArea*0.5
  
  # get the 'raw' logistic value for the current x
  vx <- logistic(x, k=k, x0=x0)
  
  # stretch it either up towards 1 or down towards 0 to coerce the scale to (x=0)=0, (x-x0)=0.5, (x=fullEdgeEffectArea)=1
  if (vx < 0.5) lerp(vx, logistic(0, k=k, x0=x0), 0.5, 0, 0.5) else lerp(vx, 0.5, logistic(fullEdgeEffectArea, k=k, x0=x0), 0.5, 1)
}


#load in the land-cover raster data using stars 

lcm<-read_stars("lcmSouthPennines.tif")
#convert to spatRaster
lcm<-rast(lcm)
#inspect
plot(lcm)

#load in polygons for the layer
lcmPoly<-st_read("lcmVector.shp")

#View first few cases
head(lcmPoly)

#create field to store land-cover classes
lcmPoly$hab<-NA

#name classes according to LCM legend
lcmPoly[lcmPoly$habNumber==1,]$hab="broadleaf"
lcmPoly[lcmPoly$habNumber==2,]$hab="conifer"
lcmPoly[lcmPoly$habNumber==3,]$hab="arable"
lcmPoly[lcmPoly$habNumber==4,]$hab="improved grass"
lcmPoly[lcmPoly$habNumber==6,]$hab="calcareous grass"
lcmPoly[lcmPoly$habNumber==7,]$hab="acid grass"
lcmPoly[lcmPoly$habNumber==9,]$hab="heather"
lcmPoly[lcmPoly$habNumber==10,]$hab="heather grass"
lcmPoly[lcmPoly$habNumber==11,]$hab="bog"
lcmPoly[lcmPoly$habNumber==12,]$hab="inland rock"
lcmPoly[lcmPoly$habNumber==14,]$hab="freshwater"
lcmPoly[lcmPoly$habNumber==20,]$hab="urban"
lcmPoly[lcmPoly$habNumber==21,]$hab="suburban"

##############################################Costs

#set costs according to those used in the paper (taken from Eycott et al. 2011)
lcmPoly$cost<-NA

lcmPoly[lcmPoly$habNumber==1,]$cost=0
lcmPoly[lcmPoly$habNumber==2,]$cost=3.03
lcmPoly[lcmPoly$habNumber==3,]$cost=10
lcmPoly[lcmPoly$habNumber==4,]$cost=10
lcmPoly[lcmPoly$habNumber==6,]$cost=4.35
lcmPoly[lcmPoly$habNumber==7,]$cost=4.35
lcmPoly[lcmPoly$habNumber==9,]$cost=2.22
lcmPoly[lcmPoly$habNumber==10,]$cost=2.22
lcmPoly[lcmPoly$habNumber==11,]$cost=2.44
lcmPoly[lcmPoly$habNumber==12,]$cost=5
lcmPoly[lcmPoly$habNumber==14,]$cost=10
lcmPoly[lcmPoly$habNumber==20,]$cost=5
lcmPoly[lcmPoly$habNumber==21,]$cost=5


##############################################Edge
#set edge effect values according to those used in the paper (taken from Eycott et al. 2011)
lcmPoly$edge<-NA

lcmPoly[lcmPoly$habNumber==1,]$edge=0
lcmPoly[lcmPoly$habNumber==2,]$edge=20.01
lcmPoly[lcmPoly$habNumber==3,]$edge=49.97
lcmPoly[lcmPoly$habNumber==4,]$edge=29.37
lcmPoly[lcmPoly$habNumber==6,]$edge=15.34
lcmPoly[lcmPoly$habNumber==7,]$edge=15.80
lcmPoly[lcmPoly$habNumber==9,]$edge=15.34
lcmPoly[lcmPoly$habNumber==10,]$edge=15.34
lcmPoly[lcmPoly$habNumber==11,]$edge=12.3
lcmPoly[lcmPoly$habNumber==12,]$edge=5
lcmPoly[lcmPoly$habNumber==14,]$edge=5
lcmPoly[lcmPoly$habNumber==20,]$edge=75.54
lcmPoly[lcmPoly$habNumber==21,]$edge=75.54


#get areas of all land-cover patches so we can apply our sigmoid function
lcmPoly$Area<-as.numeric(st_area(lcmPoly))



#create field for distance decay component (d) on the Y axis in Fig. 3. 
#We will set fullEdgeEffectArea to 100000 m2 (this means that any patch this size exerts a maximal edge effect at all distances up to the value in the "edge" field) 
  lcmPoly$log.i<-sapply(lcmPoly$Area,coercedLogistic, fullEdgeEffectArea=100000)
#create alpha from Equation 4 in paper (this sets the edge effect when multiplied by a distance)  
  lcmPoly$distEdge= -log(lcmPoly$log.i)/lcmPoly$edge
  #make sure that edge effects are always zero for woodland.
  lcmPoly[lcmPoly$habNumber==1,]$distEdge=0
  
  #create a separate object for the "habitat" class (woodland) 
  habPatch<-lcmPoly[lcmPoly$habNumber==1,]  
  
  #not ideal but let's keep things manageable by only focussing on woodland patches over 10 hectares
  habPatch<-habPatch[habPatch$Area>=100000,]
  
  #this function will create the edge surfaces for all land-use types and collect them in a list 
  funEdge<-function(x){
  
    
    matN.i<-lcmPoly[lcmPoly$edge==x,]#select all land-use polygons with the same edge value
    matN.iBuff<-st_buffer(matN.i,dist=x) # buffer these by the edge value
    
    if(unique(matN.i$edge)>0){# only continue if edge >0 (i.e. if NOT woodland)
    
    #rasterize the buffer o we can do some raster algebra (use the lcm object as a template)
    matN.iEdge<-rasterize(vect(matN.iBuff),lcm,field="distEdge") 
    
    #Now create another raster of the land-use to calculate distance away from all patches
    matN.iRas<-rasterize(vect(matN.i),lcm,field="habNumber")
    
    #create a distance raster
    buffDist<-gridDist(matN.iRas,target=NA)
    
    #now multiple the distance decay factor (aplha in Equation 4) by the distance to get the edge effect at each cell
    edge.i<-exp(-matN.iEdge*buffDist)
    
    #Let's tidy up by making sure all land-use (non-habitat) patches all have an edge effect value of 1
    #copy the land-use raster
    mask.i<-matN.iRas
    #set all values > 0 to 1
    mask.i[mask.i>0]<-1
    #set everything else to 0
    mask.i[is.na(mask.i)]<-0
    
    #add to the edge surface
    edge.i<-edge.i+mask.i
    
    #take a look
    plot(edge.i)
    #check progress
    print(x)
    
    return(edge.i)
    }
  }  
  
  #create object to hold edge values
  edge<-unique(lcmPoly$edge)
  
  #run the function
  edgeList<-lapply(edge,funEdge)
  
  #catch any layer that are null and remove
  rasL<-edgeList[!sapply(edgeList,is.null)]
  
  #now - use the sprc function to create a spatRaster collection object (do ?sprc to check documentation)
  edgeCol<-sprc(rasL)
  
  
  
  
  #now mosaic all edge raster surfaces together
  edgeFin<-mosaic(edgeCol,fun="sum") #add them all up
  edgeFin[is.na(edgeFin)]<-0 # any NAs should be core woodland so set edge to zero
  
  edgeFin[edgeFin>1]<-1 #makes no sense to have edge effect over 100% so set max to 1
  
  names(edgeFin)<-"edgeHab" # this is to know what to look for when we extract raster values later
  
  plot(edgeFin)#take a look
  
  #This is the EHI function (not slightly different arguments to the one linked from the paper as we already created edge surface before running)
  #Edge-weighted habitat function
  # 
  # Params:
  #   `matrix` = edge source polygons 
  #   `patches` = habitat patches 
  #   `specialism` = one of "interior", "edge" or "generalist"
  #   `maxDist` = maximum dispersal distance of the species being modeled
  #   `dispersalRate` = component setting rate of dispersal success
  
  
  EHI<-function(matrix, patches, specialism, dispersalRate, maxDist) {
    
   
    
    
    # EHI for interior species
    if (specialism == "interior") {
      habCells <- 1 - matrix
      
      # plot RH
      plot(habCells,axes=F,legend=T)
      extClump <- extract(habCells, patches, fun=sum, method="simple", bind=T)
      
      
      # extract raster values from within patches
      cellArea <- res(habCells)[1]^2 
      
      # create a vector for the new area
      extClump$areaMod <- extClump$edgeHab * cellArea
      
     
    }   # -- interior specialists
    
    #################  
    
    # Alternatively for edge species
    else if(specialism=="edge") {
      
      
      
      # for an edge specialist just leave as the edge habitat within the patch (i.e. inverse of above)
      habCells <- matrix
      
      # plot EH
      plot(habCells,axes=F,legend=T)
      
      # extract raster values from within patches
      extClump <- extract(matrix, patches, fun=sum, method="simple",bind=T)
      
      #get cell area
      cellArea <- res(habCells)[1]^2 
      head(extClump)
      # create a vector for the new area
      extClump$areaMod <- extClump$edgeHab * cellArea
      
   
    } #-- edge specialists
    
    #################  
    
    # for generalist just take habitat patch area
    else if (specialism=="generalist") {
      
      # for generalist all patch cells are one i.e. area not changed to just be the habitat patches  
      habCells <- st_rasterize(patches)
      
      # set patches to 1
      habCells<-rast(habCells[1,])
      habCells[habCells==0]<-1
      
      # plot EH
      plot(habCells,axes=F,legend=T)
      
      # extract raster values from within patches
      extClump <- extract(matrix, patches, fun=sum, method="simple", bind=T)
      
      
      #extract values
      extClump <- extract(matrix, patches, fun=sum, method="simple", bind=T)
      head(extClump)
      
      # extract raster values from within patches
      cellArea <- res(habCells)[1]^2 
      
      # create a vector for the new area
      extClump$areaMod <- extClump$edgeHab * cellArea
      
      
    }
    
    
    
    ###########################################Least Cost Path calculations##############################
    
    # convert matrix patches to cost raster (costs in col 5)
    costRast<-rasterize(lcmPoly,lcm,field="cost")
    
    
    #set habitat patches (NA in the data to 1)
    costRast[is.na(costRast)]<-1
    
    print("calculating transition layer, please wait")
    # gdistance takes raster package objects so convert 
    land_cost <- transition(raster(costRast), transitionFunction=function(x) 1 / mean(x), 8)
    print("transition layer done")
    # set destination points as centroids of the patches
    sites <- SpatialPoints(st_coordinates(st_point_on_surface(st_as_sf(extClump))))
    
    # init cost matrix and loop through each row
    costMat <- matrix(0, nrow=nrow(extClump), ncol=nrow(extClump))
    n.patch<-nrow(costMat)  
    for (i in 1:n.patch) {
      c <- gdistance::costDistance(land_cost, sites[i], sites) 
      costMat[i,] <- c
      print(i)
    }
    
    # this is a matrix of least cost distances between patches 
    
    #create basis for probability matrix
    distMat<-as.matrix(costMat, nrow=nrow(patches), nrow=nrow(patches)) 
    distMat <- apply(distMat, MARGIN=1, FUN=as.numeric) # make sure all elements are numeric here
    
    # set alpha which determines colonization probability of the species 
    alpha= -log(dispersalRate) / 10000
    
    # init empty matrix for adjacency matrix 
    A.prob <- matrix(0, nrow=nrow(distMat), ncol=ncol(distMat))
    
    # negative exponential of colonization kernel x distance get probability
    A.prob <- exp(-alpha * distMat) 
    
    # set diag to zero
    diag(A.prob) <- 1 
    
    # final matrix for connectivity graph
    A.prob <- as.matrix(A.prob)
    
    # final matrix for connectivity graph
    graph.Aprob <- graph_from_adjacency_matrix(A.prob, mode="undirected", weighted=T)
    
    ## calculate RHI
    
    # calculate all shortest paths between nodes
    pstar.mat <- shortest.paths(graph.Aprob, weights= -log(E(graph.Aprob)$weight))
    
    # back-transform to probabilities of connectedness
    pstar.mat <- exp(-pstar.mat)                                                   
    
    # get study area in m2
    AL <- ncell(matrix)*cellArea
    
    # get area vector 
    area <- extClump$areaMod
    
    # sum areas from the vector
    areaSum <- sum(area)
    
    # get product of all patch areas ij and multiply by probabilities above
    PCmat <- outer(area,area) * pstar.mat 
    
    # sum above
    pcMatSum <- sum(PCmat)
    
    # divide by total area of the study squared to get the PC metric  
    EHI <- sqrt(pcMatSum) / as.numeric(AL) 
    
      # compile results into list and return
    
    return(EHI)
    print("EHI done")
  } # -- EHI
  
  

  testInt<-EHI(matrix=edgeFin,patches=habPatch,dispersalRate = 0.05,maxDist = 10000,specialism = "interior")
  
  testInt  

  





