# SAS MACRO for GSD using TTE data that rely on PT assumption Version 1.0

Code is available [here](https://github.com/thewan05/GSD_SAS_Macro/blob/main/gsd_sas_macro.sas?raw=true)

# General Instructions

- We advise starting a new SAS session for each study run.
- We recommend at least 10,000 simulated samples.
- Default values for some macro parameters have been provided in the text description.
- Error checks are in place to prevent impossible values from being used as macro parameters. If an impossible value is entered, the macro will stop executing and an error message is displayed with suggestions for correcting erroneous input values.
- Warning messages are also in place. If the number of simulations is < 10000, a warning message will display. If the first interim is conducted very early (relative to total study time) and with few events, the program may run into convergence problems.
- If any warnings or errors are detected, the SAS log file will display the first fifty messages.
- The macro will not run if simulations are <1000.


1) Parameters for the macro are entered at the bottom of the program.

- NumSimul:        Number of simulated samples for the given sample size
- alpha:           Type I error
- sides:           1-sided or 2-sided test	
- lambda:          Shape parameter of the Control Arm using GG distribution
- sigma:           Scale parameter of the Control Arm using GG distribution
- med:             User entered Median of the Control Arm using GG distribution
- evt_rate:        Anticipated event rate for loss-to-follow-up (right censoring)
- seed:            A random seed is chosen
- r:               Allocation Ratio: (# in Treatment arm)/(# in Standard Arm)	
- Delta_PT_Ha:     Under the alternative, PT is greater than 1
- a:               Accrual time for the study
- a_type:          Type of accrual pattern: "1" = Uniform, "2" = Truncated Exponential with parameter omega             
- a_omega:         Parameter of "2" = truncated exponential distribution: >0 (convex) or <0 (concave); input will only be used with truncated exponential            
- t:               Total time for the study = Accrual time + Follow-up time
- bind:            Binding futility = 1; Non-binding futility=0
- num_look:        Total number of looks (including the look at the end-of-study)
- look_points:     1 = equally spaced looks, 2 = unequally spaced looks
- alpha_spend:     Type of Alpha spending function: 1=Jennison-Turnbull, 2=Hwang-Shih-DeCani, 3=User defined spending
- rho:             For both JT and HSD functions, rho=1 gives Pocock-type function
- beta:            Type II error for the study
- beta_spend:      Type of Beta spending function: 1=Jennison-Turnbull, 2=Hwang-Shih-DeCani, 3=User defined spending
- rho_f:           For both JT and HSD functions, rho=1 gives Pocock-type function
- num_skip:        Number of futility skips i.e. the first x looks at which futility is not tested
- maxiter:         Maximum number of iterations. Change default value if algorithm does not converge
- convg:           Convergence criteria. Change default value if algorithm does not converge
- direct:          File path to store log file

2) The following datasets are needed to take advantage of the macro's user defined options.

- If unequal look points (look_points=2) are selected then UserDefTime will need to include the number of look points (numuser) followed by the time points. The timepoints need to be in cumulative, ascending order and the last time point must equal the total study time,t.
- If user defined alpha spending (alpha_spend=3) or beta spending (beta_spend=3) are selected then UserDefAlpha (type I error to be spent at each look) and UserDefBeta (type II error to be spent at each look) values must be split in cumulative, ascending order and the last entry must match the inputs for alpha and beta.

3) If the file is run through the SAS interface, the SAS log file is saved as a separate text file, Mydoc.log, under the user-defined filepath, while output tables are generated as an ods listing file in the output window. The first fifty errors will be printed at the bottom of the ods listing file.

4) If the file is run under a high computing cluster, the SAS log file is saved as a separate text file, Mydoc.log and the ods listing file is save as a separate .lst file.


# Instructions for utilizing the high performance computing cluster:

- Create an account through your institution to access the cluster.
- Download the following free software:

PuTTY:      <https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html>

FileZilla:  <https://filezilla-project.org/download.php?type=client>

Notepad++:  <https://notepad-plus-plus.org/downloads/>

- Login to the cluster through PuTTY and FileZilla. Check with your institution for instructions.
- Basic Linux commands:
  - cd work: load the working directory   
  - ls: list the files in the working directory
  - sbatch exponential.sh: submission script for a batch job
- Create a script file in Notepad++
  - Under Edit, EOL, select Unix Style
  - Select 'Save As' and save as a .sh file
  - Designate number of CPUs being utilized, RAM, and the SAS file that will be run. In this example, it is exponential.sas
- Submit the file 
  - Using FileZilla, move the exponential.sas file and your submission script (exponential.sh) to the working directory.
  - Type unix command `sbatch exponential.sh' to submit.
  - You will receive an email once the job is complete.

![exp](https://github.com/thewan05/GSD_SAS_Macro//blob/main/img/exponential.JPG?raw=true)

# Examples

These examples were first published under the [methodology paper](https://pubmed.ncbi.nlm.nih.gov/31571529/) (Phadnis and Mayo, 2020).
The examples are presented once more so the reader can easily reproduce them. There may be some minor variations,
depending on the seed used.

## These macro parameters are used to obtain the results for example one and figure one.

```
NumSimul=10000, alpha=0.025, sides=1, lambda=0.5, sigma=0.75, med=20, evt_rate=0.7,
seed=1729, r=1, Delta_PT_Ha=1.4, a=12, a_type=1, a_omega=1, t=60,
bind=1, num_look=3, look_points=2, alpha_spend=1, rho=1, beta=0.10, beta_spend=1,
rho_f=1, num_skip=0, maxiter=200, convg=1E3, direct=C:\Users\user1\Desktop
```

![example1](https://github.com/thewan05/GSD_SAS_Macro/blob/main/img/example1.PNG?raw=true)

![figure1](https://github.com/thewan05/GSD_SAS_Macro/blob/main/img/figure1.PNG?raw=true)




## These macro parameters are used to obtain the results for example two.

```
NumSimul=10000, alpha=0.025, sides=1, lambda=0.5 ,sigma=0.75,med=20, evt_rate=0.7,
seed=1729, r=1, Delta_PT_Ha=1.4, a=12,a_type=1, a_omega=1, t=60,
bind=1, num_look=3, look_points=2,alpha_spend=1, rho =3, beta=0.10, beta_spend=1,
rho_f=3,num_skip=0, maxiter=200, convg=1E3, direct=C:\Users\user1\DesktopUserDefTime dataset: 3 24 36 60
```

![example2](https://github.com/thewan05/GSD_SAS_Macro/blob/main/img/example2.PNG?raw=true)






## These macro parameters are used to obtain the results for example three and figure two.

```
NumSimul=10000, alpha =0.025, sides=1, lambda=0.5, sigma=0.75,med=20, evt_rate=0.7,
seed=1729, r=1, Delta_PT_Ha=1.4, a=12,a_type=1, a_omega=1, t=60,
bind=1, num_look=3, look_points=2,alpha_spend=3, rho=3, beta=0.10, beta_spend=3,
rho_f=3,num_skip=0, maxiter=200, convg=1E3,direct=C:\Users\user1\Desktop

UserDefTime dataset: 3 24 36 60
UserDefAlpha dataset: 0.0050 0.0125 0.0250
UserDefBeta dataset: 0.0100 0.0350 0.1000
```
![example3](https://github.com/thewan05/GSD_SAS_Macro/blob/main/img/example3.PNG?raw=true)

![figure2](https://github.com/thewan05/GSD_SAS_Macro/blob/main/img/figure2.PNG?raw=true)


## These macro parameters are used to obtain the results for example four.

```
NumSimul=10000, alpha=0.025, sides=1, lambda=1, sigma=1, med=1,evt_rate=1,
seed=1729, r =1, Delta_PT_Ha=1.75, a=1, a_type=1,a_omega=1, t=4,
bind=1, num_look=4, look_points=1, alpha_spend=2,rho=1, beta=0.20, beta_spend=2,
rho_f=1, num_skip=2,maxiter=200, convg=1E3, direct=C:\Users\user1\Desktop
```

![example4](https://github.com/thewan05/GSD_SAS_Macro/blob/main/img/example4.PNG?raw=true)




# Contributing to GSD SAS Macro

We gladly welcome suggestions! For further information, contact:
Milind A. Phadnis, Department of Biostatistics and Data Science, University of Kansas Medical Center, 3901 Rainbow Blvd, Kansas City, KS 66103
Email: <mphadnis@kumc.edu> 
URL: <http://www.kumc.edu/school-of-medicine.html>


# References

Phadnis MA, Mayo MS. Group sequential design for time-to-event data using the concept of proportional time. Stat Methods Med Res. 2020 Jul;29(7):1867-1890. doi: 10.1177/0962280219876313. Epub 2019 Oct 1. PMID: 31571529.

