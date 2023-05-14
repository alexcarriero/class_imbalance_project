# Class Imbalance Project: 

This repository houses all necessary information to replicate the simulation study presented in the paper: <br>
[`The harms of class imbalance corrections for calibration in machine learning: a simulation study`](./manuscript/sim-manuscript/sim-manuscript.pdf). 

In this project, we investigated the effect of imbalance corrections on the predictive performance of machine learning models developed for dichotomous risk prediction. Here in this repository, you will find all code used in our project and full instructions for how to replicate the simulation study, project manuscript and associated shiny app.  

Please note, due to size contstraints, this repository does not contain all output from our simulation, only the processed results necessary for the creation of our shiny app and project manuscript are present. Instead, we include complete simulation output (and processed results) for one iteration of simulation scenario to demsontrate the file structure necessary to replicate the simulation study.  This is accompanied by full instructions for how to populate the repository with all simulation results after replicating our simulation study (more details below).

To download this repository populated with all simulation output and processed results from our simulation study, please contact me (Alex Carriero) the owner and maintainer of this Github repository. 

Email: a.j.carriero@students.uu.nl


## This Repository Includes The Following: 

| File / Folder                              | Contents                                                         |
| :----------------------------------------- | :--------------------------------------------------------------- |
| [`LICENSE`](./LICENSE)                     | The contents of this repository are made publicly available unter a `GNU General Public License v3.0`. <br> Please see this document for details.|
| [`simulation_code`](./simulation_code)     | In this folder, you will find two sub-directories: own_device and send_to_hpc. <br> <br> •  [`own_device`](./simulation_code/own_device):    contains all files necessary to run the simulation on a personal computer <br> <br> • [`send_to_hpc`](./simulation_code/send_to_hpc):  contains all files necessary to run the simulation on a high performance computer <br> Each directory contains a README with further details. <br> <br> Highlights: <br> <br> •  details regarding the data-generating mechanism are found in the R Markdown file: [`dgm.Rmd`](./simulation_code/own_device/data-generating-mechanism/dgm.Rmd) <br> <br> •  details regarding the output from the simulation study are found in the .txt file: [`OUTPUT.txt`](./simulation_code/own_device/sim_results/OUTPUT.txt)|
| [`results`](./results)            |      This folder stores all output from the simulation and code used to process the results.  <br> The folder contains a README with further details. |
| [`manuscript`](./manuscript)     |  In this folder, you will find all necessary files to generate both our [`manuscript`](./manuscript/sim-manuscript/sim-manuscript.pdf) and [`supplementary materials`](./manuscript/supplementary_mat.pdf). <br> This folder contains a README file with further details. |           
| [`visualize_results`](./visualize_results)               | In this folder, you will find all files necessary to generate our [`Shiny App`](https://alex-carriero.shinyapps.io/class_imbalance/). <br> This folder contains a README with further details. |

   
 
    
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
  
     - This requires manual specification of scenario number and iterations. 
     - The output of these simulation scenarios ranges heavily (from 947.9 MB to 33.71GB)
     - Depending on the HPC storage available to you, it is very likely that some simulation scenarios will need to be conducted in batches 
       i.e., by running iterations 1-200, 201-400 .... 1801-2000 in independent batches, with results downloaded and cleared from the HPC before
	     subsequent batches are run. 
		          
     b. Run the simulation for all 18 simulation scenarios on a High Performance Computer and store results appropriately. 

     
     
3.  Go to [`results`](./results). In this step, all simulation results are processed and ready for the generation of the project manuscript and Shiny App. 

 	a. Open [`process_results.Rmd`](./results/process_results.Rmd) and follow the procedure outlined to process all results. Please note: 

	   - generating calibration plots on a personal device is computationally intensive and can take up to a week for scenarios 3, 6, 9, 12, 15, and 18 
	   - we recommend running Steps 4 and 6 on a high performance computer for the simulation scenarios mentioned above
		      
	
	b.  To generate the calibration plots on a high performance computer go to [`results > hpc_calibration_plots`](./results/hpc_calibration_plots) and open the [`README`](./results/hpc_calibration_plots/README.txt)
            for full instructions. 
              
              
	c. Ensure all steps in the file [`process_results.Rmd`](./results/process_results.Rmd) are complete. 
	
	

4.  Go to [`visualize_results`](./visualize_results). In this step, we generate our Shiny App and view the results in an interactive manner. 

	a. Open the R Markdown file [`app.Rmd`](./visualize_results/app.Rmd).  Click "Run Document" to view a Shiny App displaying all results. 



5.  Finally, go to [`manuscript`](./manuscript).  In this step, we generate all figures, tables, and supplementary materials as well as a manuscript in the style of the journal Statistics in Medicine. 

	a.  Open the [`README`](./manuscript/README.txt) for details. 
    
   
   
 
 ## Technical Requirements 
 
The simulation was run on University Medical Center Utrecht's High Performance Computer and results were processed on a personal device.  Technical specifications for both of these devices, are available in [`Supplementary Materials`](./manuscript/supplementary_mat.pdf) section A. 

All code for the simulation study and processing of results was conducted using the statistical programming language `R`. Session Info for the High Performance Computer used to run the simulation and the personal device used to process the results is available in [`Supplementary Materials`](./manuscript/supplementary_mat.pdf) section A. 
 
 
 
 ## Privacy / Ethics / Security
 
 
Ethical Approval for the simulation study was granted by the Utrecht University Ethics committee and filed under number 22-1809.

All data are simulated, and not based on any empirical data.  Consequently, there are no privacy or security restrictions. This repository is accessible to the public on GitHub under the license type of `GNU General Public License v3.0`.


 ## Contact 
 
 The maintenance and public accessibility of this repository are managed by Alex Carriero. If you have any questions please feel free to contact me via email: a.j.carriero@students.uu.nl .
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
   
 
