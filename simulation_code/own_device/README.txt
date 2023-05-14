
This file contains all code used to re-generate the results of the simulation study on a personal device. Please follow the steps below to run the simulation 
on a personal device, we recommend only running the simulation on a personal device for a few iterations of scenarios which are not overly computationally intensive 
(Scenarios: 1,2,10,11) to test the code before sending to a high performance computer. 

Please note, our simulation study was conducted using high performance computing and we do not recommend running the full simulation on a personal device. 
Total computation time on a personal device is estimated to be approximately 1 calendar year, for the entire simulation.  

The contents of this folder are relevant mainly to reproduce the set up materials used in the simulation. 


Contents of this folder: 

1. data-generating-mechanism 

    dgm.Rmd     : R Markdown file containing all details for the data-generating mechanism 
    
    set.RData   : RData file containing a list with 18 objects. Each object specifies the data-generating parameters for 
                  a given simulation scenario (this file is generated in dgm.Rmd)
                
    sim_settings.csv  : A .csv file storing a data frame which specifies the simulation settings for each simulation scenario 
                       (this file is generated in dgm.Rmd)
                        

2. packages 

   In running our simulation on a high performance computer, we found four R packages necessary for 
   our simulation were unable to be compiled from their binaries on the high performance computer.  
   The source code for these packages (ebmc, IRIC, ROSE and simsalapar) was downloaded from the authors' 
   GitHub (IRIC) and CRAN (ebmc, ROSE, simsalapar) and sourced in the execute.R file for use in the simulation. 
   
   While these packages can be used on a personal device without downloading the source code, we wanted to maintain 
   identical implementation for the simulation on a personal device vs. the high performance computer. 


3. seeds 

   generate_seeds.R : a script file that generates all seeds necessary for the simulation study
   
   seedfarm.RData   : a .RData file containing a (2000 x 18) data frame which stores a unique seed for 
                      every iteration (2000) in all simulation scenarios (18).
    
    
4. execute.R 

   A script file used to run the simulation.  This file requires manual intervention to specify the 
   simulation scenario and iterations desired. 
   
   Please replace "X" with the desired simulation scenario (integer between 1 and 18), replace
   START with the desired iteration to begin with and STOP with the desired iteration to end with. 
   
   START and STOP may only take on integer values between 1 and 2000. 
   
   
5. functions.R 

   A script file that is sourced in execute.R containing all helper functions used in the simulation study. 
   
                      
6. sim_results 

   OUTPUT.txt  : In our simulation study, for each iteration four output files are generated.  Details regarding the 
		 content of each file and where it is stored are presented in OUTPUT.txt. 


   There are four folders here to store simulation results: iteration_info, per_iter_results, prediction, rep_checks
   

SET UP MATERIALS FOR SIMULATION STUDY: 

The following files are used in the simulation study, and if high performance computing is used to conduct the simulation, 
the user must ensure that the following four files exist in the folder simulation_code > send_to_hpc

If any changes are made to the following four files, please copy and paste the files to replace their equivalents in the
folder simulation_code > send_to_hpc 

	- execute.R
	- functions.R
	- set.RData
	- seedfarm.RData

 
TO RUN ON PERSONAL DEVICE: 

1. Run the Script execute.R, recall this requires manual specification of the simulation scenario and the desired iterations. 
   Also, it requires manual specification of a working directory (the path where execute.R lives on your computer). 
   Ensure all 2000 iterations are complete and that the output is stored appropriately (see OUTPUT.txt) in sim_results. 
2. Copy and paste the folder sim_results to the directory: results > simulation_results > raw_results > HERE 
3. Rename the file scX, where X is the integer representing the scenario number. 
4. Repeat for all 18 simulation scenarios. 

