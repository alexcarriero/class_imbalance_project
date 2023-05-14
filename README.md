# Class Imbalance Project: 

This repository houses all necessary information to replicate the simulation study presented in the paper: <br>
[`The harms of imbalance corrections for calibration in machine learning: a simulation study`](./manuscript/sim-manuscript/sim-manuscript.pdf). 

In this project, we investigated the effect of imbalance corrections on the predictive performance of machine learning models developed for dichotomous risk prediction. Here in this repository, you will find all code used in our project and full instructions for how to replicate our simulation study.  Please note, this repository does not contain all output from our simulation, only the processed results necessary for the creation of our [`Shiny App`](https://alex-carriero.shinyapps.io/class_imbalance/) and [`project manuscript`](./manuscript/sim-manuscript/sim-manuscript.pdf) are present. Instead, we demonstrate the file structure necessary to exactly replicate the simulation study by including complete simulation output and processed results for one iteration from simulation scenario one.  This is accompanied by full instructions for how to populate this file structure with simulation results after replicating our simulation study (more details below).

If you would like a link to download this repository populated with all output and processed results from our simulation study, please contact me (Alex Carriero) the owner of this Github repository. 

Email: a.j.carriero@students.uu.nl


## This Repository Includes The Following: 

| File / Folder                              | Contents                                                         |
| :----------------------------------------- | :--------------------------------------------------------------- |
| [`LICENSE`](./LICENSE)                     | The contents of this repository are made publicly available unter a `GNU General Public License v3.0`.  Please see this document for details.|
| [`simulation_code`](./simulation_code)     | In this folder, you will find two sub-directories: own_device and send_to_hpc. <br> <br> •  own_device:    contains all files necessary to run the simulation on a personal computer <br> <br> • send_to_hpc:  contains all files necessary to run the simulation on a high performance computer <br> <br> <br> Highlights: <br> <br> •  details regarding the data-generating mechanism are found in:  <br> simulation_code > own_device > data-generating-mechanism > dgm.Rmd <br> <br> •  details regarding the output from the simulation study are found in: <br>  simulation_code > own_device > sim_results > OUTPUT.txt |
| [`results`](./results)            |      This folder stores all output from the simulation and code used to process the results.  The folder contains a README with further details. |
| [`manuscript`](./manuscript)     |  In this folder, you will find all necessary files to generate both our manuscript and supplementary materials. This folder contains a README file with further details. |           
| [`visualize_results`](./visualize_results)               | In this folder, you will find all files necessary to generate our Shiny App. This folder contains a README with further details. |

   
 
    
## To reproduce the entire project please follow the procedure below:

1.  Go to simulation_code > own_device 

     a. Open the README for details regarding running the simulation on a personal computer and a detailed overview of the contents of the folder.  If run on a personal device (not recommended) after step 1 please skip to step 3.
     
     NOTE:
     - High performance computing is highly recommended, as the simulation is incredibly computationally intensive. 
     - The entire simulation run on a personal device would take approximately 1 calendar year to run. 
        
     b. Our simulation requires the file "set.RData".  To regenerate this file, run all code in the R Markdown file: 
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
	
	

4.  Go to the folder  visualize_results and open the R Markdown file app.Rmd.  Click "Run Document" to view a Shiny App displaying all results. 



5.  Finally go to the folder manuscript and open the README for details regarding the manuscript and 
    supplementary materials for this project. 
    
   
 
 ## Technical Requirements 
 
The simulation was run on University Medical Center Utrecht's High Performance Computer and results were processed on a personal device.  Technical specifications for both of these devices, are available in Supplementary Materials section A.  All code for the simulation study and processing of results was conducted using the statistical programming language `R`. Session Info for the High Performance Computer used to run the simulation and the personal device used to process the results is available in Supplementary Materials section A. 

Supplementary Materials can be found in this repository: manuscript > supplementary_mat.pdf
 
 
 
 ## Privacy / Ethics / Security
 
 
Ethical Approval for the simulation study was granted by the Utrecht University Ethics committee and filed under number 22-1809.

All data are simulated, and not based on any empirical data.  Consequently, there are no privacy or security restrictions. This repository is accessible to the public on GitHub under the license type of `GNU General Public License v3.0`.


 ## Contact 
 
 The maintenance and public accessibility of the project repository are managed by Alex Carriero. If you have any questions please feel free to contact me via email: a.j.carriero@students.uu.nl .
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
   
 
