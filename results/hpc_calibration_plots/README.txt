

We used High Performance Computing to save the plot coordinates of all calibration curves in our simulation. 
This required uploading our simulation results to the HPC and subsequently, downloading the files containing 
the plot coordinates to a personal device.  

Please find the procedure below: 


1. Upload the entire folder "hpc_calibration_plots" to your user portal on the HPC.



2. We generated calibration plots for all 18 simulation scenarios with both the raw and re-calibrated results. 

   For each simulation scenario and specification of raw or re-calibrated results, it is necessary to 
   upload the corresponding simulation results to the HPC.
   
   For a given scenario X (where X is an integer from 1-18), select the raw or re-calibrated predictions 
   and upload them to the HPC by following steps (a) and (b). 


   a. From the simulation results: 
   
      data_archive > results > simulation results > SELECT A RAW OR RECALIBRATED > SELECT SCENARIO FOLDER (scX)
      
   
   b. From the folder scX, select the folder "predictions" and copy the contents into the appropriate path in the 
      hpc_calibration_plots folder. 
      
      For raw results: hpc_calibration_plots > sim_results > raw_results > PASTE HERE 
      
      For recalibrated results: hpc_calibration_plots > sim_results > recalibrated_results > PASTE HERE 
   
      NOTE: there will be 2000 files in the predictions folder, each corresponding to one simulation iteration. 

      

3. Open the file "shell_plots.sh"
   
    a. Select a scenario number (integer from 1-18) and replace "X" with that integer
    b. Select a range of iterations (integers from 1-2000) 
       -  replace START with the iteration you choose at the starting iteration 
       -  replace STOP  with the iteration you choose as the stopping iteration
       
      NOTE: iterations are numbered from 1-2000 and each iteration corresponds to a specific results file. 
      
      To conduct the full simulation it is insufficient to run iterations 1-1000 two times, 
      as this will repeat the same 1000 iterations. 
      
      All 2000 iterations must be covered for a given scenario i.e., by selecting 
      START = 1 and STOP  = 2000 

      or 
      
      START = 1 STOP  = 500 and subsequently, START = 501 STOP = 1000, etc. 


   
4. Open the file "execute.R" and replace "X" (line 17) with the integer you selected in step 2. This number represents the scenario number. 
   Then, you must specify if the calibration plot is for raw_results or recalibrated_results
   in the function save_plot_coords(..., recalibrated = T/F). Set the argument recalibrated equal to FALSE for the raw results and to 
   TRUE if the results are recalibrated. 
   



5. Call the file "shell_plots.sh" on the HPC to run the Rscript "execute.R"



6. Once completed, copy the folder "plot_coords" and paste it in the following directory: 

   data_archive > results > processed_results > plot_coords > SELECT RAW OR RECALIBRATED > PASTE HERE 
   
   Then, change the folder name to "scX", where X is an integer from 1-18 representing the simulation scenario selected. 
   
   

   
   