#Bee Landscape Population Modeling 

#Set working directory
setwd("~/Documents/GitHub/spatial-pop-modeling-test/Lonsdorf_Code_R/")

#Load libraries
library(tidyverse)
library(raster)
library(terra)
library(signal)

#Create loop function
run_species <- function(s){
  
  #Read in data layers (CDL floral resources, nesting, life history)
  reclass_floral <- read.csv("CA_bees_floral.csv") %>% na.omit()
  reclass_floral[,-1] <- reclass_floral[,-1] / 100
  
  reclass_nest <- read.csv("CA_bees_nesting.csv") %>% na.omit()
  reclass_nest[,-1] <- reclass_nest[,-1] / 100
  
  reclass_lifehist <- read.csv("CA_bees_forage.csv") %>% na.omit()
  
  #Create parameters
  cell_size <- 30 #30 x 30 m cell size
  spec_num <- nrow(reclass_lifehist) #7 species
  t_step <- 1 #time step
  max_nest <- reclass_lifehist$Nests[s] #max number of nests within a cell
  max_rep <- reclass_lifehist$Rep[s] #max number of queens produced per nest
  forage_dist <- reclass_lifehist$Disp[s]
  
  #Raster file path
  base <- "Tiff Files/cdl_Yolo_"
  raster_path <- paste0(base, (t_step + 2012), ".TIF") #raster file path
  
  #Read in raster layer and crop extent
  curr_LC <- rast(raster_path) #read in raster file
  curr_LC <- trim(curr_LC) #trim extra white space
  extent <- floor(ext(curr_LC)) #create extent to crop using
  curr_LC <- crop(curr_LC, extent) #crop raster
  
  nr <- nrow(curr_LC) #number of rows in raster
  nc <- ncol(curr_LC) #number of columns in raster
  
  #Initialize floral and nesting habitat 
  HF <- classify(curr_LC, reclass_floral[, c(1, s+1)]) #reclassify raster using reclass files 
  HN <- classify(curr_LC, reclass_nest[, c(1, s+1)]) #reclassify raster using reclass files
  HN_mat <- as.matrix(HN, wide = TRUE)
  
  #Create empty arrays 
  queen_array <- max_queen_array <- disperse <- array(0, dim = c(nr, nc, 20)) 
  queen_seed <- spec_pop <- pol_ls_score <- matrix(0, nrow = nr, ncol = nc)
  
  #Initialize queens
  set.seed(10) #for reproducibility
  queen_seed[611:620, 809:818] <- runif(n = 10 * 10) #generate random number in 10 x 10 grid in center of raster
  queen_array[,,1] <- ifelse(queen_seed > 0.8, 1, 0) #assign queen to cell if random number is greater than 0.8
  
  #Create foraging distance matrix
  radius <- 2*round(forage_dist/cell_size) #foraging radius in cells
  X <- outer(-radius:radius, rep(1, length(-radius:radius)), FUN = function(x,y) x)
  Y <- outer(-radius:radius, rep(1, length(-radius:radius)), FUN = function(x,y) y)
  
  dist_mat <- as.matrix(cell_size * sqrt(X^2 + Y^2)) #distance to center
  dist_mat_decay <- exp(-dist_mat/forage_dist) * (dist_mat <= 2 * forage_dist) #decay to make circle
  eff_forage <- sum(dist_mat_decay) #sum of circle distance matrix
  filter <- dist_mat_decay/eff_forage #normalized distance matrix (sums to 1), moving window
  
  for (t in 1:20) {
    print(paste0("Time step ", t)); flush.console()
    
    #Pull new raster file after 10 years of initializing 
    if(t > 10) {
      raster_path <- paste0(base, (t+2012-10),  ".TIF") #assign new raster path 
      curr_LC <- rast(raster_path) #current land cover raster
      curr_LC <- crop(curr_LC, extent) #crop to same extent
      
      #Initialize floral and nestin habitat
      HF <- classify(curr_LC, reclass_floral[, c(1, s+1)]) #reclassify raster using reclass files 
      HN <- classify(curr_LC, reclass_nest[, c(1, s+1)]) #reclassify raster using reclass files
      HN_mat <- as.matrix(HN, wide = TRUE)
    }
    
    #Moving window calculation
    forage <- focal(HF, filter, na.rm = TRUE) #convolution based on distance and available forage
    forage_mat <- as.matrix(forage, wide = TRUE)
    pol_ls_score <- HN * forage
    
    #Disperse queens
    if (t > 1) {
      max_queen_array[,,t] <- max_nest * HN_mat #calculate maximum number of queens per cell given landscape class
      queen_array[,,t] <- pmin(disperse[,,t-1], max_queen_array[,,t]) #take smallest value (can't exceed max)
      
      rep_temp <- matrix(rnorm(nr * nc), nrow = nr, ncol = nc) #random number generation to create rep variation
      range <- max(rep_temp) - min(rep_temp) #range of random variation
      standardized <- floor(max_rep * 2 * (rep_temp - min(rep_temp))/range) #standardize values
      rep <- queen_array[,,t] * (standardized * forage_mat) #generate rep
      rep_rast <- rast(rep) #convert to raster for focal calc
      disperse[,,t] <- as.matrix(focal(rep_rast, filter, na.rm = TRUE), wide = TRUE)
      
    } else {
      #Where t = 1 (initial conditions)
      rep_temp <- matrix(rnorm(nr * nc), nrow = nr, ncol = nc) #random number generation to create rep variation
      range <- max(rep_temp) - min(rep_temp) #range of random variation
      standardized <- floor(max_rep * 2 * (rep_temp - min(rep_temp))/range) #standardize values
      rep <- queen_array[,,t] * (standardized * forage_mat) #generate rep
      rep_rast <- rast(rep) #convert to raster for focal calc
      disperse[,,t] <- as.matrix(focal(rep_rast, filter, na.rm = TRUE), wide = TRUE)
    }
    
    spec_pop[,t] <- sum(queen_array[,,t], na.rm = TRUE) #population size of each species at a given time step
  }
  
  #Output
  return(list(queens = rast(curr_LC, nlyr = 20, vals = queen_array), 
              ls_score =  pol_ls_score))
  
}

#Run function
spec_1 <- run_species(1)
spec_1_q <- spec_1[[1]]
spec_1_ls <- spec_1[[2]]
plot(spec_1_q)

plet(spec_1_q[[20]])
