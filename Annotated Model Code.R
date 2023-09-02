### Code for "Cetacean-mediated nitrogen transport in the oceanic Gulf of Mexico" ###

## Citation: Woodstock, M.S.,J.J. Kiszka, M.R. Ramírez-León, T.T. Sutton, K. Fennel, B. Wang, Y. Zhang. (in review). Cetacean-mediated nitrogen transport in the oceanic Gulf of Mexico. Limnology and Oceanography

## Author (MS Woodstock) email: fishesofthedeep@gmail.com


#* Notes:
#* The new user will need to set the working directory (L 32) to where the "Cetacean Traits.xlsx file is.
#* The new user will have to set a results working directory (L75).
#* All other files will be created as the model runs.
#* The model will run in parallel (multiple iterations at once).
#* The model will export group-specific excretion rates, one group at a time. Each group file will have two .csv files (Consumption and Excretion Record.csv and Dive Table.csv). Consumption and Excretion Record.csv has all of the consumption and excretion information. Dive table is more of a diagnostic to make sure the group actually dove as it was supposed to.


## Clear all existing items in an environment
rm(list=(ls()))

### Load all of the libraries
  require(PBSmodelling)
  require(snowfall)
  require(parallel)
  require(snow)
  require(foreach)
  require(doSNOW)
  require(PBSadmb)
  require(tidyr)
	require(readxl)
  require(janitor)
  require(dplyr)


## Set global working directory SET YOUR OWN TO LOCATION WITH FILES
workdir <- "D:/Cetacean Mediated Nutrient Transport"
setwd(workdir)

## Import datasets with clean_names() function using the @janitor package
abundance <- clean_names(read_excel(paste(workdir,"/Input Data/Cetacean Traits.xlsx",sep=""),sheet = "Abundance"))
history <- clean_names(read_excel(paste(workdir,"/Input Data/Cetacean Traits.xlsx",sep=""),sheet = "Life History"))
fish_prox <- clean_names(read_excel(paste(workdir,"/Input Data/Cetacean Traits.xlsx",sep=""),sheet = "Fish Proximate"))
ceph_prox <- clean_names(read_excel(paste(workdir,"/Input Data/Cetacean Traits.xlsx",sep=""),sheet = "Cephalopod Proximate"))
crust_prox <- clean_names(read_excel(paste(workdir,"/Input Data/Cetacean Traits.xlsx",sep=""),sheet = "Crustacean Proximate"))
dive_table <- clean_names(read_excel(paste(workdir,"/Input Data/Cetacean Traits.xlsx",sep=""),sheet = "Dive Table"))


n_boot <-5000  #number of runs

#* Dataframe of input parameters. This is effectively useless when running in parallel because R will not retain information from previous iterations. However, it is kept here for non-parallel computing. The information can be gathered by analyzing the "groups" dataframes, which are output each iteration.
inputs <- data.frame(Simulation = rep(0,n_boot),Sperm_Whale_n = rep(0,n_boot),Rices_Whale_n = rep(0,n_boot),CBW_n = rep(0,n_boot),BBW_n = rep(0,n_boot),GBW_n = rep(0,n_boot),CBDO_n = rep(0,n_boot),PSD_n = rep(0,n_boot),STD_n = rep(0,n_boot),SPD_n = rep(0,n_boot),CD_n = rep(0,n_boot),FD_n = rep(0,n_boot),KW_n = rep(0,n_boot),FKW_n = rep(0,n_boot),PKW_n = rep(0,n_boot),DSW_n = rep(0,n_boot),PSW_n = rep(0,n_boot),MHW_n = rep(0,n_boot),RD_n = rep(0,n_boot),SFPW_n = rep(0,n_boot),fish_prot_prop = rep(0,n_boot),ceph_prot_prop = rep(0,n_boot),crust_prot_prop = rep(0,n_boot),n_pods = rep(0,n_boot))


n_ts <- (24*60)/1 #Number of minutes in a day divided by the length of the actual time step (e.g., X/1 = 1 minute time step);  1 minute makes sense because some surface intervals are very short. However, this model could be adjusted to run at different time steps


#* Protein concentration for fishes
pro_prop_fish_mean_init <- mean(fish_prox$mean_protein_concentration_percent_wet_weight)/100 #*Stickney and Torres (1989) %Wet Weight
pro_prop_fish_sd_init <- sd(fish_prox$mean_protein_concentration_percent_wet_weight)/100 #* Stickney and Torres (1989)

#* Protein concentration for cephalopods
pro_prop_ceph_mean_init <- mean(ceph_prox$mean_protein_concentration_percent_wet_weight)/100 #*Sinclair et al. 2015 %Wet Weight
pro_prop_ceph_sd_init <- sd(ceph_prox$mean_protein_concentration_percent_wet_weight)/100 #* Sinclair et al. 2015

#* Protein concentration for crustaceans
#** Some species had NAs for protein concentration, so this needs to be removed.
pro_prop_crust_mean_init <- mean(na.omit(crust_prox$mean_protein_concentration_percent_wet_weight))/100 #*Donnelly et al. 1993 %Wet Weight
pro_prop_crust_sd_init <- sd(na.omit(crust_prox$mean_protein_concentration_percent_wet_weight))/100 #* Donnelly et al. 1993


## Other paramters that are necessary. Congruent with Roman and McCarthy 2010
prot_n <- 0.17 #*% Nitrogen by weight for protein (Gaskin 1982)
prop_N_exc <- 0.8 #* Pretty standard assumption, but no great empirical estimate for this


#* Results output dataframe

## Set to wherever you want results to go
resdir <- paste(workdir,"/Results/Submit/",sep="") 


#Set up parallel computing
no_cores <- detectCores()  #determine number of cores on computer
cl<-makeCluster(no_cores)  #setup the size of your parallel cluster based on number of cores
registerDoSNOW(cl)         #register cluster to get it ready for parallel processing

#set up text progress bar....this is personal code just to have a progress bar letting your know % completion
pb <- txtProgressBar(max = n_boot, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# following is the main code for running in parallel
#ls is your parallel processing object, essentially the same as a for loop...need to include any packages that are required for use in the loop
#the code after setting up the parallel processing object is the loop you want to run on each core, this code will split the number of runs evenly across cores
# here the number of runs is set by nsim (so run the code following ls nsim times, and divide nsim/ncores)
#if need any info saved from each run, make sure that code explicitly saves that information to hardrive (e.g., as csv, etc.) otherwise values will not be saved into memory
ls=foreach(sim=1:n_boot,.options.snow = opts,.combine='rbind',.packages =c('PBSmodelling')) %dopar% {
  
	## Need to reload all libraries to be used in parallel
  require(readxl)
  require(janitor)
  require(dplyr)

	## Add simulation number to input dataframe
  inputs$Simulation[sim] <- sim
  
  ## Initial dataframe for species abundance; More of a placeholder for total species abundance
  general <- data.frame(Species = abundance$species,
                        Abundance_n = rep(0,n_distinct(abundance$species)))
  
  
  ## Apply abundance and biomass
  #* Our populations are assumed to all be within the GoM
  #* Function estimates a random value between minimum and maximum abundance (sensu Hayes et al. 2019)
  
  for (spec in 1:n_distinct(abundance$species)){
    while ((general$Abundance_n[spec] > abundance$abundance[spec]+(abundance$abundance[spec]*abundance$abundance_cv[spec])) || (general$Abundance_n[spec] < abundance$minimum_abundance[spec])){ #Need this loop so the bootstrapped value is not below the minimum abundance in the stock assessment
      general$Abundance_n[spec] <- round(rlnorm(1,meanlog = log(abundance$abundance[spec])-log(1+abundance$abundance_cv[spec]^2)/2,sdlog=sqrt(log(1+abundance$abundance_cv[spec]^2))),0)
    }
    inputs[sim,(spec+1)] <- general$Abundance_n[spec]
  }
  
	
  ## Sample pod sizes and distributions
  #* Dataframe containing information for each group. Data repeated for 5000 iterations as a placeholder since data are added to the dataframe at different points in the code, making @dplyrs add_row() function not applicable. 
  groups <- data.frame(Species = rep(abundance$species[1],5000),
                       Species_num = rep(0,5000),
                       Individuals_n = rep(0,5000),
                       Biomass_kg = rep(0,5000),
                       Latitude = rep(NA,5000),
                       Longitude = rep(NA,5000),
                       Bottom_Depth = rep(0,5000),
                       Meso_fish_Abun = rep(0,5000),
                       Meso_ceph_Abun = rep(0,5000),
                       Meso_crust_Abun = rep(0,5000))
  
  
  n_group <- 1 #* Counter for the number of groups in the model
  for (spec in 1:n_distinct(abundance$species)){
    spec_sum <- 0 #* Necessary to count species abundances being created and make sure we do not overshoot the estimate
    
    
    ## Fill groups dataframe until all groups for a species are made
    while (spec_sum < general$Abundance_n[spec]){ ## Will run until enough animals are assigned for the given species
      groups$Species[n_group] <- general$Species[spec] #Name of species
      groups$Species_num[n_group] <- spec ##Add a species ID
      ## Sample this until the model finds a suitable value for the group size
      while (groups$Individuals_n[n_group] < history$minimum_group_size[spec] || groups$Individuals_n[n_group] > history$maximum_group_size[spec]){
        
        ## Round to 0 (integer) and apply mean and SE values from Maze-Foley and Mullin 2006
        groups$Individuals_n[n_group]<- round(rnorm(1,mean = history$mean_group_size[spec],sd = history$se_group_size[spec]*sqrt(history$n_group_size[spec])),0)
        
      }
      
      #*Change final group size to match reserve to get to population size
      if ((groups$Individuals_n[n_group]+spec_sum) > general$Abundance_n[spec]){
        groups$Individuals_n[n_group] <- general$Abundance_n[spec] - spec_sum
      } 
      
      #Running counter so we know how many individuals of a species have been allocated
      spec_sum <- spec_sum + groups$Individuals_n[n_group]
      
      #Convert Abundance to biomass
      #* Assumes all individuals are the same weight
      groups$Biomass_kg[n_group] <- groups$Individuals_n[n_group] * history$average_individual_weight_kg[spec]
      
      #Continue the counter
      n_group <- n_group + 1
    }
  }
  
  ## Assign number of groups in model
  inputs$n_pods[sim] <- n_group
  
  groups <- groups[1:(n_group-1),] #Remove Excess Rows from the groups dataframe
  
  ## Resample Mesopelagic N proportion around mean 
  inputs$fish_prot_prop[sim] <- rnorm(1,pro_prop_fish_mean_init,pro_prop_fish_sd_init) #10% ± 3%
  inputs$ceph_prot_prop[sim] <- rnorm(1,pro_prop_ceph_mean_init,pro_prop_ceph_sd_init) #9% ± 2%
  inputs$crust_prot_prop[sim] <- rnorm(1,pro_prop_crust_mean_init,pro_prop_crust_sd_init) #9% ± 4%
  
  
  
  #### All of this is edited because the MaxENT results are not shareable
  
  # This normally doesn't exist and is for example
  for (group in 1:nrow(groups)){
  	groups$Latitude[group] <- 26
  	groups$Longitude[group] <- -90
  	groups$Bottom_Depth[group] <- 5000
  }
  
  
  ## The below code is normal with the .csv files available.
  
  ## Place each group in a designated latitude and longitude according to MaxENT model
  ## Load map for each group
  #n_group <- 1 #Count groups
  
  #for (spec in 1:n_distinct(abundance$species)){
  	
  	
  	### I rewrote this code so that the user does not need the habitat suitability maps ##
  	
  	
  	## A csv file per species that has spatially resolved habitat suitability vales that match the cetacean model grid
    #map <- read.csv(paste(workdir,"/Input Data/SDM Maps/SDM Excel/",history$map_file[spec],sep=""))
    
    ## Subset for one species to speed up code
    #sub <- groups %>% filter(Species_num == spec)
    
    ## Gathers random numbers (n = number of groups for the species) from 0 to 1
    #num <- runif(dim(sub)[1],0,map$Cumulative[dim(map)[1]])
    
    ## Assign latitude and longitude to groups weighted according to habitat suitability
    #for (a in 1:length(num)){
      
    #  if (num[a]>map$Cumulative[1]){
    #    sub2 <- map %>% filter(Cumulative <= num[a])
    #    groups$Latitude[n_group] <- sub2$y[dim(sub2)[1]]
    #    groups$Longitude[n_group] <- sub2$x[dim(sub2)[1]]
    #    groups$Bottom_Depth[n_group] <- sub2$mask_adj[dim(sub2)[1]]
    #    n_group <- n_group + 1
    #  } else {
    #    groups$Latitude[n_group] <- map$y[1]
    #    groups$Longitude[n_group] <- map$x[1]
    #    groups$Bottom_Depth[n_group] <- map$mask_adj[1]
    #    n_group <- n_group + 1
    #  }
    #}
  #}
  
  
  ############## Resume coding ##############
  
  
  ### Model is now initialized ###
  
  
  #### Model Run ----
  for (group in 1:dim(groups)[1]){ #Sample through the number of groups
    
    ## Dive table (1 per group) to be exported at the end
    dives <- data.frame(Species = rep("",n_ts),
                        Timestep = rep(0,n_ts),
                        Depth = rep(0,n_ts))
    
    ## Consumption Table (1 per group) to be exported at the end
    cons_exc <- data.frame(Species = character(),
                           Latitude = numeric(),
                           Longitude = numeric(),
                           Timestep = numeric(),
                           Consumed_kg_ind = numeric(),
                           Ration = numeric(),
                           Perc_BWGT = numeric(),
                           Depth = numeric(),
                           Z_bin = character(),
                           Consumed_N = numeric(),
                           Excreted_N = numeric())
    
    ### Determine species of group ###
    for (a in 1:dim(abundance)[1]){
      if (abundance$species[a]==groups$Species[group]){
        species_num <- a
      }
    }
    
    ## Set dives remaining for group at initial dive pattern
    deep_dives_remaining <- dive_table$deep_n_dives_per_day[species_num]
    shallow_dives_remaining <- dive_table$shallow_n_dives_per_day[species_num]
    
    ## Have a running counter to determine if animal should be resting ##
    deep_dive_capable <- TRUE
    shallow_dive_capable <- TRUE
    
    ## Counters for time spent in a dive/rest stage
    deep_dive_interval <- 0
    deep_surface_interval <- 0
    shallow_dive_interval <- 0
    shallow_surface_interval <- 0
    
    ## Set proportions of mesopelagics at night within the epipelagic zone
    
    #* All values represent the average and standard deviation from DEEPEND data
    groups$Meso_fish_Abun[group] <- rnorm(1,0.369862,0.032644)
    groups$Meso_ceph_Abun[group] <- rnorm(1,0.3999993,0.000813)
    groups$Meso_crust_Abun[group] <- rnorm(1,0.526189,0.052683)
    
    
    for (ts in 1:n_ts){ #Run through each minute of the day
      dives$Species[ts] <- groups$Species[group]
      dives$Timestep[ts]  <- ts
      
      ## Resample dive interval so that each interval is slightly different according to a coefficient of variation of 0.2
      deep_dive_int <- round(rnorm(1,dive_table$deep_dive_duration_min[species_num],dive_table$deep_dive_duration_min[species_num]*0.2),0)
      deep_surface_int <- round(rnorm(1,dive_table$deep_surface_interval_min[species_num],dive_table$deep_surface_interval_min[species_num]*0.2),0)
      
      ## Initialize the first shallow dive if necessary for the species
      if (!is.na(dive_table$shallow_day_mean_dive_depth[species_num])){
        shallow_dive_int <- round(mean(c(rnorm(1,dive_table$shallow_dive_duration_min[species_num],dive_table$shallow_dive_duration_min[species_num]*0.2))),0)
        shallow_surface_int <- round(rnorm(1,dive_table$shallow_surface_interval_min[species_num],dive_table$shallow_surface_interval_min[species_num]*0.2),0)
      }
      
      ## Initiate deep dive sequence if animal was at surface and can dive. Animal will also dive during the first time step
      if ((dives$Depth[ts-1]== 0 && deep_dive_capable==T) || ts == 1){
        
      	## Nightime diving sequence between 7 pm and 7 am
        if (ts/n_ts <= (420/1440) || ts/n_ts >= (1140/1440)){
        	
        	## Dive to depth in a normal distribution of dive depth and standard deviation.
          dives$Depth[ts] <- rnorm(1,dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_sd_dive_depth[species_num])
          
          ## 5% possibility of a deep dive (between normal upper limit and total max) ##            
          rnum <- runif(1,0,1)
          if (rnum <= 0.05){
            dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_max_deep_dive_depth[species_num])),dive_table$deep_night_sd_dive_depth[species_num])
            while (dives$Depth[ts] > dive_table$deep_night_max_deep_dive_depth[species_num]){
              #* Dive cannot be deeper than deepest recorded depth for the species, but with uncertainty in diving behaviors, it was possible and had to be constrained

              dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_max_deep_dive_depth[species_num])),dive_table$deep_night_sd_dive_depth[species_num])
            }
          }
          
          #Dive cannot be shallower than 0 meters but with a large enough dive depth uncertainty, it may be possible.
          while (dives$Depth[ts] < 0){
            dives$Depth[ts] <- rnorm(1,dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_sd_dive_depth[species_num])
          }
        } else {
          ## It is between 7 am and 7 pm, so the daytime dive sequence should be initiated
          
          dives$Depth[ts] <- rnorm(1,dive_table$deep_day_mean_dive_depth[species_num],dive_table$deep_day_sd_dive_depth[species_num])
          
          ## 5% possibility of a deep dive (between normal upper limit and total max) ##
          rnum <- runif(1,0,1)
          if (rnum <= 0.05){
            dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_day_mean_dive_depth[species_num],dive_table$deep_day_max_deep_dive_depth[species_num])),dive_table$deep_day_sd_dive_depth[species_num])
            
            while (dives$Depth[ts] > dive_table$deep_day_max_deep_dive_depth[species_num]){
              #* Dive cannot be deeper than deepest recorded depth
              dives$Depth[ts] <- rnorm(1,mean(c(dive_table$deep_day_mean_dive_depth[species_num],dive_table$deep_day_max_deep_dive_depth[species_num])),dive_table$deep_day_sd_dive_depth[species_num])
            }
          }
          
          #Dive cannot be shallower than 0 meters but with a large enough uncertainty, it is possible.
          while (dives$Depth[ts] < 0){
            dives$Depth[ts] <- rnorm(1,dive_table$deep_night_mean_dive_depth[species_num],dive_table$deep_night_sd_dive_depth[species_num])
          }
        }
        
      	## Animals deeper than bottom depth stay at maximum ocean depth in the location ##
      	if (dives$Depth[ts] > groups$Bottom_Depth[group]){
      		dives$Depth[ts] <- groups$Bottom_Depth[group]
      	}
      	
        ## Assure animal does not deep dive again until possible
        deep_dive_capable <- FALSE
        shallow_dive_capable <- FALSE
        cycle <- "deep"
        
        ### Initiate Foraging Sequence (All values in units kg) ###
        ## Gather Consumption rate for group that is between 1 and 4.5% bodyweight.
        ## The mean for all species is 3% and a 0.2 CV gives us between 1 and 4.5% bodyweight consumed
        #* Assumed a CV of 0.2
        ## Kg of biomass consumed in the dive
        cons <- rnorm(1,groups$Biomass_kg[group]*history$percent_bodyweight_consumed_percent[spec],groups$Biomass_kg[group]*history$percent_bodyweight_consumed_percent[spec]*0.2)
        
        # Calculate Ration (kg/day) by dividing consumption by number of dives animal will perform. The if/else statement is required to make the appropriate calculation for species that are going to make both deep dives and shallow dives
        if (is.na(dive_table$shallow_n_dives_per_day[species_num])){
        	rat <- cons/dive_table$deep_n_dives_per_day[species_num]
        } else {
        	rat <- cons/sum(dive_table$deep_n_dives_per_day[species_num],dive_table$shallow_n_dives_per_day[species_num])
        }
        
        ## Percent bodyweight consumed
        bwgt <- (rat/groups$Biomass_kg[group])*100

        
        ## Mesopelagic Nitrogen Consumed = Feces + Urine + Storage
        # Incorporate fish, cephalopod, and crustacean protein concentrations and then convert protein to N consumed. Equations are slightly different during night and day because the prey community changes. This means that the consumed and excreted N at night in the epipelagic zone is only considering the mesopelagic contribution. The consumption and excretion of epipelagic N during the day is still calculated here, but not included in analyses because there is no deep-sea contribution.
        
        if (ts/n_ts <= (420/1440) || ts/n_ts >= (1140/1440) && dives$Depth[ts] < 200){
          ## Nighttime feeding shallower than 200m
        	#* Only a proportion of the diet is mesopelagic
        	con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num]*groups$Meso_fish_Abun[group])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num]*groups$Meso_ceph_Abun[group])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num])*groups$Meso_crust_Abun[group])*prot_n  

        } else {
          ## Nighttime feeding below 200m and daytime feeding (i.e., no mesopelagic community adjuster)
          con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
          
        }
        
        #* Assumed 20% is stored
        exc <- con_n * prop_N_exc
        
        ## Assign depths ##
        
        depth <- dives$Depth[ts]
        
        ## Assign depth bins ##
        if (dives$Depth[ts] < 100){
          zbin <- "Upper Epipelagic"
        } else if (dives$Depth[ts]>=100 && dives$Depth[ts]<200){
          zbin  <- "Lower Epipelagic"
        } else if (dives$Depth[ts]>=200 && dives$Depth[ts]<600){
          zbin  <- "Upper Mesopelagic"
        } else if (dives$Depth[ts]>=600 && dives$Depth[ts]<1000){
          zbin  <- "Lower Mesopelagic"
        } else if (dives$Depth[ts]>=1000 && dives$Depth[ts] < groups$Bottom_Depth[group]){
          zbin  <- "Bathypelagic"
        }
        if (dives$Depth[ts] == groups$Bottom_Depth[group]){
          ## Determine if it is on the bottom
          zbin  <- "Benthopelagic"
        }
        
        ## Add to the consumption and excretion dataframe that is exported. Note that we are keeping track of the time (for day/night analyses), latitude and longitude (for spatial analyses), and depth bin (for vertical analyses). "Depth" here means the dive depth, which was used to identify the depth bin, but is not specificially used in analyses.
        cons_exc <- cons_exc %>% add_row(Species = groups$Species[group],Latitude = groups$Latitude[group],Longitude = groups$Longitude[group],Timestep = ts,Consumed_kg_ind = cons,Ration = rat,Perc_BWGT = bwgt,Depth = depth,Z_bin = zbin,Consumed_N = con_n,Excreted_N = exc)
        
      } else if (dives$Depth[ts-1] == 0 && deep_dive_capable == F){ ##Continuation of surface interval
        deep_surface_interval <- deep_surface_interval + 1
        dives$Depth[ts] <- 0
        
        
        ## End surface interval
        if (deep_surface_interval >= deep_surface_int){
          deep_surface_interval <- 0
          deep_dive_capable <- TRUE
        }
        
      } else if (dives$Depth[ts-1] > 0){ ## Continuation of dive interval
        deep_dive_interval <- deep_dive_interval + 1
        
        dives$Depth[ts] <- dives$Depth[ts-1]
        
        
        ## End dive interval
        if (deep_dive_interval >= deep_dive_int){
          deep_dive_interval <- 0
          dives$Depth[ts] <- 0
          cycle <- ""
        }
      }
      
      ## Initiate shallow dive sequence, if applicable ##
      ## All of this code is the same as the above code, but the dive sequence refers to shallow diving behavior, rather than deep. In practice, this is only used for the beaked whales, which exhibit two distinct diving behaviors. The words "deep" and "shallow" are a bit of a misnomer depending on your definition of "deep".
      
      if (!is.na(dive_table$shallow_day_mean_dive_depth[species_num])){
        if (shallow_dive_capable==TRUE) {
          ## Initiate shallow dive sequence ##
          
          ## It is between 7 pm and 7 am
          dives$Depth[ts] <- rnorm(1,dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_sd_dive_depth[species_num])
          
          ## 5% possibility of abnormally shallow dive (between normal upper limit and total max) ##
          rnum <- runif(1,0,1)
          if (rnum <= 0.05){
            dives$Depth[ts] <- rnorm(1,mean(c(dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_max_deep_dive_depth[species_num])),dive_table$shallow_night_sd_dive_depth[species_num])
            
            while (dives$Depth[ts] > dive_table$shallow_night_max_shallow_dive_depth[species_num]){
              #* Dive cannot be shallower than shallowest recorded depth
              dives$Depth[ts] <- rnorm(1,mean(c(dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_max_deep_dive_depth[species_num])),dive_table$shallow_night_sd_dive_depth[species_num])
            }
          }
          
          #Dive cannot be shallower than 0 meters but with a large enough uncertainty, it is possible.
          while (dives$Depth[ts] < 0){
            dives$Depth[ts] <- rnorm(1,dive_table$shallow_night_mean_dive_depth[species_num],dive_table$shallow_night_sd_dive_depth[species_num])
          }
          
          ## Assure animal does not shallow dive agin until possible
          shallow_dive_capable <- FALSE
          deep_dive_capable <- FALSE
          cycle <- "shallow"
          
          if (dives$Depth[ts] > groups$Bottom_Depth[group]){
            ## Animals deeper than bottom depth stay on bottom for benthopelagic coupling ##
            dives$Depth[ts] <- groups$Bottom_Depth[group]
          }
          
          ### Initiate Foraging Sequence (All values in units kg) ###
          ## Gather Consumption rate for group
          #* Assumed a CV of 0.2
          
          cons <- rnorm(1,groups$Biomass_kg[group]*history$percent_bodyweight_consumed_percent[spec],groups$Biomass_kg[group]*history$percent_bodyweight_consumed_percent[spec]*0.2)
          
          # Calculate Ration (kg/day). The if/else statement is required to make the appropriate calculation for species that are going to make both deep dives and shallow dives
          if (is.na(dive_table$shallow_n_dives_per_day[species_num])){
          	rat <- cons/dive_table$deep_n_dives_per_day[species_num]
          } else {
          	rat <- cons/sum(dive_table$deep_n_dives_per_day[species_num],dive_table$shallow_n_dives_per_day[species_num])
          }
          
          ## Percent bodyweight consumed
          bwgt <- (rat/groups$Biomass_kg[group])*100
          
          ## Mesopelagic Nitrogen Consumed = Feces + Urine + Storage
          # Incorporate fish, cephalopod, and crustacean protein concentrations and then convert protein to N consumed. Equations are slightly different during night and day because the prey community changes
          if (ts/n_ts <= (420/1440) || ts/n_ts >= (1140/1440) && dives$Depth[ts] < 200){
          	## Nighttime feeding shallower than 200m
          	#* Only a proportion of the diet is mesopelagic
          	con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num]*groups$Meso_fish_Abun[group])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num]*groups$Meso_ceph_Abun[group])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num])*groups$Meso_crust_Abun[group])*prot_n  
          	
          } else {
          	## Nighttime feeding below 200m and daytime feeding
          	con_n <- ((rat*inputs$fish_prot_prop[sim]*history$prop_fish_diet[species_num])+(rat*inputs$ceph_prot_prop[sim]*history$prop_ceph_diet[species_num])+(rat*inputs$crust_prot_prop[sim]*history$prop_crust_diet[species_num]))*prot_n
          	
          }  
          
          #* Assumed 20% is stored
          exc <- con_n * prop_N_exc
          
          
          ## Assign depths ##
          
          depth <- dives$Depth[ts]
          
          ## Assign depth bins ##
          if (dives$Depth[ts] < 100){
            zbin  <- "Upper Epipelagic"
          } else if (dives$Depth[ts]>=100 && dives$Depth[ts]<200){
            zbin  <- "Lower Epipelagic"
          } else if (dives$Depth[ts]>=200 && dives$Depth[ts]<600){
            zbin  <- "Upper Mesopelagic"
          } else if (dives$Depth[ts]>=600 && dives$Depth[ts]<1000){
            zbin  <- "Lower Mesopelagic"
          } else if (dives$Depth[ts]>=1000 && dives$Depth[ts] < groups$Bottom_Depth[group]){
            zbin  <- "Bathypelagic"
          } 
          
          if (dives$Depth[ts] == groups$Bottom_Depth[group]){
            zbin  <- "Benthopelagic"
          }
          
          
          cons_exc <- cons_exc %>% add_row(Species = groups$Species[group],Latitude = groups$Latitude[group],Longitude = groups$Longitude[group],Timestep = ts,Consumed_kg_ind = cons,Ration = rat,Perc_BWGT = bwgt,Depth = depth,Z_bin = zbin,Consumed_N = con_n,Excreted_N = exc)
          
          
        } else if (dives$Depth[ts-1] == 0 && shallow_dive_capable == F && ts != 1 && cycle !="deep"){
          ## Continuation of surface interval
          shallow_surface_interval <- shallow_surface_interval + 1
          dives$Depth[ts] <- 0
          
          ## End surface interval
          if (shallow_surface_interval >= shallow_surface_int){
            shallow_surface_interval <- 0
            shallow_dive_capable <- TRUE
          }
          
        } else if (dives$Depth[ts-1] > 0 && ts != 1 && cycle != "deep"){
          ## Continuation of dive interval
          shallow_dive_interval <- shallow_dive_interval + 1
          
          dives$Depth[ts] <- dives$Depth[ts-1]
          
          
          ## End dive interval
          if (shallow_dive_interval >= shallow_dive_int){
            shallow_dive_interval <- 0
            dives$Depth[ts] <- 0
            cycle <- ""
          }
        }
      }
      
      ## For diagnostic purposes to check if animal dove when needed to and vice versa.
      dives$Cap[ts] <- deep_dive_capable
      dives$Int[ts] <- deep_dive_interval
      
      ## End Time Step
    }
    setTxtProgressBar(pb,group)
    ## Create working directory if it doesn't exist
    if (!dir.exists(paste(resdir,"/Sim ",sim,sep=""))){
      dir.create(paste(resdir,"/Sim ",sim,sep=""))
    }
    if (!dir.exists(paste(resdir,"/Sim ", sim,"/Group ",group,sep=""))){
      dir.create(paste(resdir,"/Sim ", sim,"/Group ",group,sep=""))
    }
    setwd(paste(resdir,"/Sim ",sim,"/Group ",group,sep=""))
    #* Only need consumption and excretion values with observations
    cons_exc <- cons_exc %>% filter(Z_bin != "")
    write.csv(cons_exc,"Consumption and Excretion Record.csv")
    write.csv(dives,"Dive Table.csv")
    ## End Group
  }
  
  ## Export all data ##
  
  ## Export simulation-based dataframes ##
  setwd(paste(resdir,"/Sim ",sim,sep=""))
  
  ## The groups dataframe also contains maximum ocean depth in the location 
  write.csv(groups,"Groups Dataframe.csv")
  write.csv(inputs,"Input Directory.csv")
  
  ## End Iteration
  
  
}


# end parallel job
close(pb)
stopCluster(cl) #end the cluster for parallel processing
closeAllConnections() 

