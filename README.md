# SAS MACRO for GSD using TTE data that rely on PT assumption Version 1.0

We advise starting a new SAS session for each study run. We recommend at least 10,000 simulated samples.


For further information, contact:
Milind A. Phadnis, Department of Biostatistics and Data Science, University of Kansas Medical Center, 3901 Rainbow Blvd, Kansas City, KS 66103
Email: <mphadnis@kumc.edu> 
URL: <http://www.kumc.edu/school-of-medicine.html>

Code is available [here](github link)

## General Instructions

1) Parameters for the macro are entered at the bottom of the program.

- NumSimul        Number of simulated samples for the given sample size
- alpha           Type I error
- sides           1-sided or 2-sided test	
- lambda          Shape parameter of the Control Arm using GG distribution
- sigma           Scale parameter of the Control Arm using GG distribution
- med             User entered Median of the Control Arm using GG distribution
- evt_rate        Anticipated event rate for loss-to-follow-up (right censoring)
- seed            A random seed is chosen
- r               Allocation Ratio: (# in Treatment arm)/(# in Standard Arm)	
- Delta_PT_Ha     Under the alternative, PT is greater than 1
- a               Accrual time for the study
- a_type          Type of accrual pattern: "1" = Uniform, "2" = Truncated Exponential with parameter omega             
- a_omega         Parameter of "2" = truncated exponential distribution: >0 (convex) or <0 (concave); input will only be used with truncated exponential            
- t               Total time for the study = Accrual time + Follow-up time
- bind            Binding futility = 1; Non-binding futility=0
- num_look        Total number of looks (including the look at the end-of-study)
- look_points     1 = equally spaced looks, 2 = unequally spaced looks
- alpha_spend     Type of Alpha spending function: 1=Jennison-Turnbull, 2=Hwang-Shih-DeCani, 3=User defined spending
- rho             For both JT and HSD functions, rho=1 gives Pocock-type function
- beta            Type II error for the study
- beta_spend      Type of Beta spending function: 1=Jennison-Turnbull, 2=Hwang-Shih-DeCani, 3=User defined spending
- rho_f           For both JT and HSD functions, rho=1 gives Pocock-type function
- num_skip        Number of futility skips i.e. the first x looks at which futility is not tested
- maxiter         Maximum number of iterations. Change default value if algorithm does not converge
- convg           Convergence criteria. Change default value if algorithm does not converge
- direct          File path to store log file

2) The following datasets are needed to take advantage of the macro's user defined options.

- If unequal look points (look_points=2) are selected then UserDefTime will need to include the number of look points (numuser) followed by the time points. The timepoints need to be in cumulative, ascending order and the last time point must equal the total study time,t.
- If user defined alpha spending (alpha_spend=3) or beta spending (beta_spend=3) are selected then UserDefAlpha (type I error to be spent at each look) and UserDefBeta (type II error to be spent at each look) values must be split in cumulative, ascending order and the last entry must match the inputs for alpha and beta.

3) If the file is run through the SAS interface, the SAS log file is saved as a separate text file, Mydoc.log, under the user-defined filepath, while output tables are generated as an ods listing file in the output window. The first fifty errors will be printed at the bottom of the ods listing file.

4) If the file is run under a high computing cluster, the SAS log file is saved as a separate text file, Mydoc.log and the ods listing file is save as a separate .lst file.


## Instructions for utilizing the high performance computing cluster:

- Create an account through your institution to access the cluster.
- Download the following free software:

PuTTY:      <https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html>

FileZilla:  <https://filezilla-project.org/download.php?type=client>

Notepad++:  <https://notepad-plus-plus.org/downloads/>

- Login to the cluster through PuTTY and FileZilla. Check with your institution for instructions.
- Basic Linux commands:
1) cd work: load the working directory   
2) ls: list the files in the working directory
3) sbatch exponential.sh: submission script for a batch job
- Create a script file in Notepad++
1) Under Edit, EOL, select Unix Style
2) Select 'Save As' and save as a .sh file
3) Designate number of CPUs being utilized, RAM, and the SAS file that will be run. In this example, it is exponential.sas
- Submit the file 
1) Using FileZilla, move the exponential.sas file and your submission script (exponential.sh) to the working directory.
2) Type unix command `sbatch exponential.sh' to submit.
3) You will receive an email once the job is complete.



