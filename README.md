**Class Imbalance Project:** 

This repository houses all necessary information to replicate the simulation study presented in our paper: 
The harms of imbalance corrections for calibration in machine learning: a simulation study. 

In this project, we investigated the effect of imbalance corrections on the predictive performance of machine learning models developed for dichotomous risk prediction. 
We considered 18 data-generating scenarios, each with 2000 iterations.  For each iteration we compared 30 prediction models.  Each prediction model was
developed using a unique combination of a class imbalance correction (used to pre-process data) and a machine learning algorithm (used to train the prediction model). 

Please consult the project manuscript for details, found here: manuscript > sim-manuscript > sim-manuscript.pdf.

Should you have any questions, please feel free to contact me (Alex Carriero) the owner of this GitHub repository. 

Email: alexjcarriero@gmail.com

NOTE: this repository does not contain all output from our simulation study (only the processed results that are necessary to generate our manuscript and shiny app are included). Instead, we demonstrate the file structure necessary to exactly replicate the simulation study by including all simulation output and processed results for one iteration from simulation scenario 1 as well as full instructions for how these folders may be populated with simulation results (more details below).  If you would like a link to download this repository populated with all from output and processed results from our simulation study, please contact me (Alex Carriero).

Thank you !

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    

This Repository Includes The Following: 

1.  LICENSE 

    All code for this project is published under a GNU general public license.  Please find details in the LICENSE file. 
 
    
    
2.  simulation_code 

     In this folder, you will find two sub-folders: own_device and send_to_hpc
     
     	own_device:    contains all files necessary to run the simulation on a personal computer 
     	send_to_hpc:   contains all files necessary to run the simulation on a high performance computer 
     
     README files are contained in each sub-folder for further orientation within each directory. 
     
     QUICK REFERENCES:  
     
     	- details regarding the data-generating mechanism are found in:  simulation_code > own_device > data-generating-mechanism > dgm.Rmd
     	- details regarding the output from the simulation study are found in:  simulation_code > own_device > sim_results > OUTPUT.txt
     

3.  results

     In this folder, you will find all output from the simulation and code used to process the results.  
     The folder contains a README with further details. 
   
   

4. manuscript 

     In this folder, you will find all necessary files to generate both our manuscript and supplementary materials. 
     This folder contains a README file with further details. 
     
     
    
5.  visualize_results 

    In this folder, you will find all files necessary to generate our Shiny App. 
    This folder contains a README with further details. 
   
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
    
To reproduce the entire project follow the procedure below: 

 PREFACE :  read the manuscript  (manuscript > sim-manuscript > sim-manuscript.pdf) for a detailed presentation of the simulation study

1.  Go to simulation_code > own_device 

     a. Open the README for details regarding running the simulation on a personal computer and a detailed overview of the contents of the folder.   
          If run on a personal device (not recommended) after step 1 please skip to step 3.
     
     	 NOTE: 
	 	
	 	- High performance computing is highly recommended, as the simulation is incredibly computationally intensive. 
		- The entire simulation run on a personal device would take approximately 1 calendar year to run. 
        
     b. Our simulation requires the file "set.RData".  To regenerate this file, run all code in the file: 
      	 simulation_code > own_device > data-generating-mechanism > dgm.Rmd
	 
     c. Our simulation requires the file "seedfarm.RData".  To regenerate this file, run all code in the R script: 
     	simulation_code > own_device > seeds > generate_seeds.R
     
     d. Paste the files from (b) and (c) to replace their equivalent files in the folder simulation_code > send_to_hpc
     
  
     
2.  Go to simulation_code > send_to_hpc 

    a. Open the README  for instructions regarding how to run the simulation on a high performance computer and how to store the results.  
    
        NOTE: 
        	
		- This requires manual specification of scenario number and iterations. 
		- The output of these simulation scenarios ranges heavily (from 947.9 MB to 33.71GB)
		- Depending on the HPC storage available to you, it is very likely that some simulation scenarios will need to be conducted in batches 
		   i.e., by running iterations 1-200, 201-400 .... 1801-2000 in independent batches, with results downloaded and cleared from the HPC before
		          subsequent batches are run. 
		          
        
     b. Run the simulation for all 18 simulation scenarios on a High Performance Computer and store results appropriately. 

     
     
3.  Go to results 

 	a. Open process_results.Rmd and follow the procedure outlined to process all results. 
	
	   NOTE: 
	   
	   	- generating calibration plots on a personal device is computationally intensive and can take up to a week for scenarios
	   	      3, 6, 9, 12, 15, and 18
		  
		- we recommend running Steps 4 and 6 on a high performance computer for the simulation scenarios mentioned above
		      
	
	b.  To generate the calibration plots on a high performance computer go to results > hpc_calibration_plots and open the README 
              for full instructions. 
              
              
	c. Ensure all steps in process_results.Rmd are complete. 
	
	

4.  Go to the folder  visualize_results and open the file app.Rmd.  Click "Run Document" to view a Shiny App displaying all results. 



5.  Finally go to the folder manuscript and open the README for details regarding the manuscript and 
    supplementary materials for this project. 
    
    
 
 ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
 
 Technical Requirements: 
 
The simulation was run on University Medical Center Utrecht's High Performance Computer and results were processed on a personal device.  Technical specifications for both of these devices, and
details regarding R package versions is available in Supplementary Materials section A. 

Supplementary Materials can be found in this repository: manuscript > supplementary_mat.pdf
 
 
 ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 
 Privacy / Ethics / Security:
 
 
Ethical Approval for the simulation study was granted by the Utrecht University Ethics committee and filed under number 22-1809.

All data are simulated, and not based on any empirical data.  Consequently, there are no privacy or security restrictions. 

 ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 
 
 THE END 
 
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
   
 
