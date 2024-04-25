To Regenerate Simulation Results Using High Performance Computing: 

1. Upload the entire folder "send_to_hpc" to your user portal on the High Performance Computer (HPC). 


2. Open the file "shell.sh" 
  
    a. Select a scenario number (integer from 1-18) and replace "X" with that integer
    b. Select a range of iterations (integers from 1-2000) 
       -  replace START with the iteration you choose as the starting iteration 
       -  replace STOP  with the iteration you choose as the stopping iteration
       
      NOTE: iterations are numbered from 1-2000 and each iteration corresponds to a specific seed. 
      
      To conduct the 2000 iterations for a given scenario it is insufficient to run 
      iterations 1-1000 two times, as this will repeat the same 1000 iterations. 
      
      All 2000 iterations must be covered for a given scenario i.e., by selecting 
      START = 1 and STOP  = 2000 

      or 
       
      START = 1,    STOP = 500 and subsequently 
      START = 501,  STOP = 1000 
      ...
      START = 1501, STOP = 2000
      

3. Open the file "execute.R" and replace "X" (line 43) with the integer you selected in step 2. 
   This number represents the scenario number. 
   
   
4. Using the High Performance Computer call the .sh file.  This will run the R script execute.R on the HPC. 

   NOTE: ensure you have set the correct working directory (i.e., send_to_hpc). 
   
   
5. Once complete, copy the folder "sim_results" and paste in the directory: 

   results > simulation_results > raw_results > PASTE HERE 
   
   Then name the folder "scX", where X is an integer representing the scenario number. 
   


NOTE: all files used in the simulation are explained in the README found in simulation_code > own_device 


THE END

             