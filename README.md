# Class Imbalance Project: 

This repository is a supplement to the paper: [The harms of class imbalance corrections for machine learning based prediction models: a simulation study](https://arxiv.org/abs/2404.19494). 

In this project, we investigated the effect of imbalance corrections on the predictive performance of various machine learning models. In this repository, you will find all code used in our project and full instructions for how to replicate the simulation study, case study, and the associated shiny app.  

Please note, due to size constraints, this repository does not contain the raw output from our simulation, only the processed results necessary for the creation of our shiny app. Instead, we include complete simulation output (and processed results) for one iteration of simulation scenario 1 to demonstrate the file structure necessary to replicate the simulation study.  This is accompanied by full instructions for how to populate the repository with all simulation results after replicating our simulation study (more details below).

If you would like a copy of this repository populated with all results, please contact me (Alex Carriero).

Email: a.j.carriero@umcutrecht.nl


## This Repository Includes The Following: 

| File / Folder                              | Contents                                                         |
| :----------------------------------------- | :--------------------------------------------------------------- |
| [`LICENSE`](./LICENSE)                     | The contents of this repository are made publicly available unter a `GNU General Public License v3.0`. <br> Please see this document for details.|
| [`simulation_code`](./simulation_code)     | In this folder, you will find two sub-directories: own_device and send_to_hpc. <br> <br> •  [`own_device`](./simulation_code/own_device):    contains all files necessary to run the simulation on a personal computer <br> <br> • [`send_to_hpc`](./simulation_code/send_to_hpc):  contains all files necessary to run the simulation on a high performance computer <br> <br> Each sub-directory contains a README with further details.|
| [`results`](./results)            |      This folder stores all output from the simulation and code used to process the results. |  
| [`case_study`](./case_study)            |      This folder contains the code for our case study based on MIMIC-III data. For privacy reasons the data are not included in our repository. |  
| [`visualize_results`](./visualize_results)               | In this folder, you will find all files necessary to generate our [`Shiny App`](https://alex-carriero.shinyapps.io/class_imbalance/). |
| [`supplementary_materials`](./supplementary_materials) | In this folder you will find the appendix and supplementary materials. |

   
 
    
## To reproduce the entire project please follow the procedure below:

1.  Go to [`simulation_code > own_device`](./simulation_code/own_device). In this step, we generate all files necessary for the simulation study. 

     a. Open the [`README`](./simulation_code/own_device/README.txt) for details regarding running the simulation on a personal computer and a detailed overview of the contents of the folder.  If run on a personal device (not recommended) after step 1 please skip to step 3. Please note: 
   
     - High performance computing is highly recommended, as the simulation is incredibly computationally intensive. 
     - The entire simulation run on a personal device would take approximately 1 calendar year to run. 
        
     b. Our simulation requires the file `set.RData`.  To regenerate this file, run all code in the R Markdown file: 
      	 [`dgm.Rmd`](./simulation_code/own_device/data-generating-mechanism/dgm.Rmd)
	 
     c. Our simulation requires the file `seedfarm.RData`.  To regenerate this file, run all code in the R script: 
     	[`generate_seeds.R`](./simulation_code/own_device/seeds/generate_seeds.R)
     
     d. Paste the files from (b) and (c) to replace their equivalent files in the folder simulation_code > send_to_hpc
  
     
2.  Go to [`simulation_code > send_to_hpc`](./simulation_code/send_to_hpc). In this step, we conduct the full simulation study using a High Performance Computer. 

    a. Open the [`README`](./simulation_code/send_to_hpc/README.txt) for instructions regarding how to run the simulation on a high performance computer and how to store the results. Please note: 
  
     - This requires manual specification of the simulation scenario and desired iterations. 
     - The output from the simulation scenarios ranges heavily (from 947.9 MB in scenario 1 to 33.71 GB in scenario 18)
     - Depending on the HPC storage available to you, it is very likely that some simulation scenarios will need to be conducted in batches 
       i.e., by running iterations 1-200, 201-400 .... 1801-2000 in independent batches, with results downloaded and cleared from the HPC before
	     subsequent batches are run. 
		          
     b. Run the simulation for all 18 simulation scenarios on a High Performance Computer and store results appropriately. 

     
     
3.  Go to [`results`](./results). In this step, all simulation results are processed. 

 	a. Open [`process_results.Rmd`](./results/process_results.Rmd) and follow the procedure outlined to process all results. Please note: 

	   - generating calibration plots on a personal device is computationally intensive and can take up to a week for scenarios 3, 6, 9, 12, 15, and 18 
	   - we recommend running Steps 4 and 6 on a high performance computer for the simulation scenarios mentioned above
		      
	
	b.  To generate the calibration plots on a high performance computer go to [`results > hpc_calibration_plots`](./results/hpc_calibration_plots) and open the [`README`](./results/hpc_calibration_plots/README.txt)
            for full instructions. 
              
              
	c. Ensure all steps in the file [`process_results.Rmd`](./results/process_results.Rmd) are complete. 
	
	

4.  Go to [`visualize_results`](./visualize_results). In this step, we generate our Shiny App and view the results in an interactive manner. 

	a. Open the R Markdown file [`app.Rmd`](./visualize_results/app.Rmd).  Click "Run Document" to view a Shiny App displaying all results. 

   
 
 ## Technical Requirements 
 
The simulation was run on University Medical Center Utrecht's High Performance Computer and results were processed on a personal device.  Technical specifications for both of these devices, are available in the Supplementary Materials section A.  All code for the simulation study and processing of results was run using the statistical software `R`. 


 ## Contact 
 
 The maintenance and public accessibility of this repository are managed by Alex Carriero. If you have any questions please feel free to contact me via email: a.j.carriero@umcutrecht.nl .
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
   
 
