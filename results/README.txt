

Results from the Simulation Study. 


Please find an outline of the contents of this folder introduced below: 

 
1. simulation_results:
   
    -  All simulation results are stored in this folder. 
    -  This folder contains two sub-directories: raw_results and recalibrated_results
    
       raw_results: stores all simulation output before re-calibration 
       recalibrated results: stores all simulation output after re-calibration
       
       
   simualtion_results > raw_results
      
    - contains 18 sub-directories, one for each simulation scenario
    - each sub-directories has four folders: 
    
       iteration_info   : one file per iteration (2000 total), storing warnings/errors for each simulation iteration. 
      
       per_iter_results : one file per iteration (2000 total), storing empirical performance metrics for each simulation iteration. 
      
       predictions      : one file per iteration (2000 total), storing the predicted probabilities and true class of the validation data 
                          for each simulation iteration.
      
       rep_check        : one file per iteration (2000 total), storing test values from simulated data so reproducibility may be assessed. 
    
       NOTE: the content of each of these four folders is detailed in file OUTPUT.txt: simulation_code > own_device > sim_results < OUTPUT.txt.    
        

   simulation_results > recalibrated_results
    
    - contains 18 subdirectories, one for each simulation scenario
    - each sub-directory has two folders: 
      
       per_iter_results : once results are processed it will store
        
                         one file per iteration (2000 total), storing empirical performance metrics for each simulation iteration. 
                         (after re-calibration)
      
       predictions      : once results are processed it will store
        
                         one file per iteration (2000 total), storing the predicted probabilities and true class of the validation data 
                         for each simulation iteration (after re-calibration)
                     
                     
                   
    
2. processed_results > plot_coords
    
    - contains two sub-directories: raw_results and re_calibrated results
    

      
   processed_results > plot_coords > raw_results
      
   - contains 18 sub-directories, one for each simulation scenario
   - once results are processed each sub-directory will contain:
      
     2000 files, each storing the plot coordinates of the flexible calibration curves for one iteration
     (without re-calibration) 
   
      
      
   processed_results > plot_coords > recalibrated_results 
      
   - contains 18 subdirectories, one for each simulation scenario 
   - once results are processed each sub-directory will contain:
      
     2000 files, each storing the plot coordinates of the flexible calibration curves for one iteration
     (after re-calibration) 

    
        
3. process_results.Rmd 

   Instructions and all code used to process the simulation results are included in this file. 

   TO REGENERATE THE RESULTS -- SEE INSTRUCTIONS IN THIS R MARKDOWN FILE.



4. process_results.R 
 
   Helper functions used in process_results.Rmd 



5. hpc_calibration_plots

   This folder contains all content necessary to generate calibration plots using a higher performance computer.  
   A README file is contained in this folder with further instructions. 







                   
                   