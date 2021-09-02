# Instructions

## Requirement: MATLAB, MCR 2018b 64 bit. 
MCR has been used for the purpose when the computers at different sites do not have MATLAB installed. In that case, install MCR on a computer that has MATLAB and create a standalone executable by writing **'mcc -e **.mlapp'** in the command prompt. This will create a standalone executable file. This executable file then can be used in any computer that has MCR istalled.  

## Instructions to run static tests: 
1. Download static localisation codes (App_Static_loc.mlapp and Staticloc.m).
2. If MATLAB has been installed, Open the MATLAB app (.mlapp file) and run the code. The UI requires the following inputs to perform the analysis:  
   * a) 'Parent Path' - The folder that contains KIM log files. Usually the Image folder contains 'Markerlocation_GA.txt' file/files which is/are needed for this analyis.
   * b) 'Coordinate file' - The patient coordinate file.  
      * Contains the [x y z] details of the marker positions in the patient (code assumes 3 markers per patient) followed by the [x y z] isocentre position in the patient
      * One line per marker, all positions in mm  
      * Final line is the isocentre in the same format as the markers  
      * The code requires this file to only contain numbers.  
   * c) 'Static shifts' - Applied couch shifts. 
   * d) Click on 'Compute Accuracy'. 
   * e) The code generate a file named 'Metrics.txt' with the information mean, standard deviation and percentiles in LR, SI and AP directions for all the data points and a figure with the KIM trace. 

## Instructions to run dynamic localisation test and treatment interruption test: 
1. Download dynamic localisation codes (App_Dynamic_loc.mlapp and Dynamicloc.m) and treatment interruption codes (App_treat_int.mlapp and TreatmentInt.m).
2. If MATLAB has been installed, Open the MATLAB app (.mlapp file) and run the code. The UI requires the following inputs to perform the analysis:  
   * a) 'Select KIM Traj folder' - The folder that contains KIM log files, these are the text file/s which include *MarkerLocation_GA* in the name. 
      * They can typically be found in the 'KIM-kV' folder
   * b) 'Select Robot Traj file' - The appropriate motion trace file (can be for the Robot or Hexamotion).  
   * c) 'Select coord file' - The patient coordinate file.  
      * Contains the [x y z] details of the marker positions in the patient (code assumes 3 markers per patient) followed by the [x y z] isocentre position in the patient
      * One line per marker, all positions in mm  
      * Final line is the isocentre in the same format as the markers  
      * The code requires this file to only contain numbers.    
   * d) 'Select param file' - Contains the parameters **a b c**; all three values are seconds  
      * The program shifts the timing of the KIM trace to match the timing of the robot trace between the times a and c with a step size of b until the SI positions match.  
      * Create a .txt file containing the parameters **a b c**, for example **-30 0.01 30**. In this example, the code will shift the trace from -30 sec to 30 sec in steps of 0.01 sec.  
   * e) Click on 'analyse'.
 
## Output
This will produce an output file with the information mean, standard deviation and percentiles in LR, SI and AP directions for all the data points in the KIM input trace and a figure with KIM trace and robot trace on top of each other. The name of the output file and the figure will be denoted by the KIM folder name stated in step 2a).  