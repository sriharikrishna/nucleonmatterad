# nucleonmatterad
Algorithmic differentiation of either symmetric nuclean matter (SNM) or pure nucleon matter (PNM). 

## Setup

1. Obtain [Tapenade](http://www-sop.inria.fr/tropics/tapenade/downloading.html)

2. Change to tap branch
    `git checkout tap`
3. `cd src` 
4. Change the variable `TAPENADE_HOME` in MakefileTapf

## Execution
1. Use any of the following scripts. 
      Replace ID by an integer [1-60].
      Replace DELTA by 0.1 or 0.2.
      For the following scripts snm/DIFFSIZES.inc should declare    
      ```
             integer, parameter :: nbdirsmax=2
      ```
      ```
      ./script_dfo_snm.sh ID DELTA
      ./script_bfgs_snm.sh ID DELTA    
      ```
      For the following scripts pnm/DIFFSIZES.inc should declare    
      ```
             integer, parameter :: nbdirsmax=4      
      ```
      ```
      ./script_bfgs_pnm.sh ID DELTA        
      ./script_dfo_pnm.sh ID DELTA    
      ```
      For the following scripts snm/DIFFSIZES.inc should declare    
      ```
             integer, parameter :: nbdirsmax=7
      ```
      ```
      ./script_bfgs_snm_fullx.sh ID DELTA    
      ```

      For the following scripts snm/DIFFSIZES.inc should declare    
      ```
             integer, parameter :: nbdirsmax=7
      ```
      ```
       ./script_bfgs_pnm_fullx.sh ID DELTA	
      ```
	
    


