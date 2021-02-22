* SAS Version 9.4;
* Copyright (c) 2002-2012 by SAS Institute Inc., Cary, NC, USA.;
* SAS (r) Proprietary Software 9.4 (TS1M2);
* Licensed to UNIVERSITY OF KANSAS - SFA T&R, Site 70081492;

* Copyright (C) 2020  Milind Phadnis

/* Topic: Group Sequential Designs

This program achieves the following:

A group sequential design using the concept of Proportional Time can be
used for time-to-event data where the assumptions of proportional hazards
or exponentially distributed lifetimes are not true. This method utilizes
the gamma ratio distribution to calculate the efficacy and safety boundaries
and be used for all distributions that are members of the generalized gamma
family applying an error spending approach.

Warning: At the start of this program, the log file, output window,
results window, and datasets in the working directory will be erased.

Recommended system requirements:

It takes 40 minutes to run 10,000 simulations on 24 GB RAM.

Runtime on a high performance computing cluster: 25 minutes, 2CPU, 8GB



/* Legal Notice

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; *without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

*/


/*Start new SAS session before you run this program*/;



dm 'log;clear;output;clear; ';	/*Clear the log file, results window, working directory*/;
dm 'odsresults; clear';;
PROC DATASETS LIB=work NOlist MEMTYPE=data kill; 
options nosymbolgen nomprint nomlogic;	/*system options set to default*/
proc printto; run;	/*ods listing will be redirected inside macro; set this to default again*/


%let _timer_start = %sysfunc(datetime());	/*Start timer*/

ods graphics off;
%macro main(NumSimul = , 				/* Number of simulated samples for the given sample size */
	alpha = ,							/* Type I error */
	sides = ,							/* 1-sided or 2-sided test */
	lambda = ,							/* Shape parameter of the Control Arm using GG distribution */
	sigma = ,							/* Scale parameter of the Control Arm using GG distribution */
	med = ,								/* User entered Median of the Control Arm using GG distribution */
	evt_rate = ,						/* Anticipated event rate for loss-to-follow-up (right censoring) */
	seed = ,							/* A random seed is chosen */
	r = ,								/* Allocation Ratio: (# in Treatment arm)/(# in Standard Arm) 	*/
	Delta_PT_Ha = ,						/* Under the alternative, PT is greater than 1 */
	a = ,								/* Accrual time for the study */
	a_type = ,							/* Type of accrual pattern: "1" = Uniform, "2" = Truncated Exponential with parameter omega */
	a_omega = ,							/* Parameter of "2" = truncated exponential distribution: >0 (convex) or <0 (concave); input will only be used with truncated exponential*/
	t = ,								/* Total time for the study = Accrual time + Follow-up time*/
	bind = ,								/* Binding futility = 1; Non-binding futility=0 */
	num_look = , 						/* Total number of looks (including the look at the end-of-study) */
	look_points = ,						/* 1 = equally spaced looks, 2 = unequally spaced looks */
	alpha_spend = ,						/* Type of Alpha spending function: 1=Jennison-Turnbull, 2=Hwang-Shih-DeCani */
	rho = ,								/* For both JT and HSD functions, rho=1 gives Pocock-type function */
	beta = ,							/* Type II error for the study */
	beta_spend = ,						/* Type of Beta spending function: 1=Jennison-Turnbull, 2=Hwang-Shih-DeCani */
	rho_f = ,							/* For both JT and HSD functions, rho=1 gives Pocock-type function */
	num_skip = , 						/* Number of futility skips i.e. the first x looks at which futility is not tested */
	maxiter= ,							/* Maximum number of iterations. Change default value if algorithm does not converge */
	convg= ,							/* Convergence criteria. Change default value if algorithm does not converge*/
	direct= );							/* File path to store log file */

	%let MyDist = gamma;						/* The member of the GG family. Can be changed if you know the distribution. Resulting runtime may be shorter.*/
	%let count = 0;								/* Variable specific to the Search-N algorithm */
	%let num = 0;								/* Variable specific to the Search-N algorithm */
	%let b_count = 0;							/* Variable specific to the Search-beta algorithm */
	%let b_num = 0;								/* Variable specific to the Search-beta algorithm */
	%let b_switch = 0;							/* Variable specific to the Search-beta algorithm */
	%let n1 = 128;								/* Initial sample size to start the Search-N algorithm */
	%let b1 = &beta ;							/* Initial value of type-II error used to the Search-beta algorithm */
	%let Delta_PT = 1;							/* Proportional Time (PT) ratio = {Theta1:new}/{Theta0:old} */


	/**********CALCULATE MU FROM SIGMA, LAMBDA, AND MEDIAN TIME**********/
	data gamma;
		do mu= -20 to 20 by 0.00001 until ((abs(x-0.5)<0.0001));  /* This steps reduces user-input burden */
			if &lambda = 0 then do;
				upperbound = (log(&med) - mu)/&sigma;
				x= cdf('Normal',upperbound);
				output;
				end;
			else do;
				kappa=1/(&lambda*&lambda);
				upperbound = exp(&lambda*((log(&med)-mu)/&sigma))*kappa;
				x=cdf('Gamma',upperbound,kappa);
				output;
				end;
		end;
		drop kappa upperbound x;
	run;

	data tmp;
		 set gamma nobs=nobs;
		 if _n_=nobs;
		 call symputx('mu',mu,'G'); /*create global macro variable with last observation of mu*/
	run;



	/**********MISSING VALUES**********/
	%if &NumSimul=. %then %let &NumSimul=10000;
	%if &alpha=. %then %let &alpha=0.05;
	%if &sides=. %then %let &sides=1;
	%if &r=. %then %let &r=1;
	%if &seed=. %then %let &seed=1726;
	%if &a_type=. %then %let &a_type=1;
	%if &alpha_spend=. %then %let &alpha_spend=1;
	%if &beta=. %then %let &beta=0.1;
	%if &beta_spend=. %then %let &beta_spend=1;
	%if &rho=. %then %let &rho=1;
	%if &rho_f=. %then %let &rho_f=1;

	


	/**********WARNINGS**********/
	%if &NumSimul < 10000 %then 
		%put WARNING: NumSimul should be at least 10000;


	/**********CREATE MACRO VARIABLE FOR ALPHA, BETA SPEND **********/
	%if &alpha_spend=3 %then %do;
		data UserDefAlpha2;
		set UserDefAlpha;
			do j=1 to &numuser;
				if j=i then
					call symputx(catt('alpha_look',j),alpha_look,'G');
			end;
		run;
	%end;

	%if &beta_spend=3 %then %do;
		data UserDefBeta2;
		set UserDefBeta;
			do j=1 to &numuser;
				if j=i then
					call symputx(catt('beta_look',j),beta_look,'G');
			end;
		run;
	%end;
	
	
	/**********TRANSFORM USER-DEFINED DATASET**********/ 
	%if &look_points=1 %then %do; 
		data UserDefTime;
			usertime1 = rand("Uniform");     /*Create a dummy UserDefTime dataset if look_points=1*/
			usertime2 = rand("Uniform");
		run;
	%end;


		

	%if &look_points=2 %then %do;
	*1) Error check;
		data UserDefTime;
			set UserDefTime;
				if usertime < lag(usertime) then do;
					put 'ERROR: User-defined times must be in ascending order.';
					ABORT;
				end;
				if usertime <0 then do;
					put 'ERROR: User-defined times must be greater than zero.';
					ABORT;
				end;
		run;

	*2) Transform datasets;
		proc transpose data=UserDefTime prefix=usertime out=UserDefTime;
		run;

		data UserDefTime;
			set UserDefTime;
			if _NAME_~='usertime' then delete;
		run;

	*3) More error Check;	
		data UserDefTime;
			set UserDefTime;
			if &t ~= usertime&num_look then do;
				put ERROR: "The user defined values of time do not add up to &t";
				ABORT;
			end;
		run;
	%end;


	%if &alpha_spend=3 %then %do;
		data UserDefAlpha;
			set UserDefAlpha;
			if alpha_look < lag(alpha_look) then do;
				put 'ERROR: User-defined alphas must be in ascending order.';
				ABORT;
			end;
			if alpha_look <0 then do;
				put 'ERROR: User-defined alphas must be greater than zero.';
				ABORT;
			end;
		run;

		proc transpose data=UserDefAlpha prefix=alpha_look out=UserDefAlpha;
		run;

		data UserDefAlpha;
			set UserDefAlpha;
			if _NAME_~='alpha_look' then delete;
		run;

		data UserDefAlpha;
			set UserDefAlpha;
			if &alpha ~= alpha_look&num_look then do;
				put ERROR: "The user defined values of alpha do not add up to &alpha";
				ABORT;
			end;
		run;
	%end;


	data UserDefBeta3;
		%do j=1 %to &num_look;		/*Create a dummy UserDefBeta dataset if beta spend doesn't equal 3*/
			b&j = rand("Uniform");
		%end;
		_NAME_=0;
		ERROR=0;
	run;

	

	%if &beta_spend=3 %then %do;
		data UserDefBeta;
			set UserDefBeta;
			if beta_look < lag(beta_look) then do;
				put 'ERROR: User-defined betas must be in ascending order.';
				ABORT;
			end;
			if beta_look <0 then do;
				put 'ERROR: User-defined betas must be greater than zero.';
				ABORT;
			end;
		run;

		proc transpose data=UserDefBeta prefix=beta_look out=UserDefBeta;
		run;

		data UserDefBeta;
			set UserDefBeta;
			if _NAME_~='beta_look' then delete;
		run;

		data UserDefBeta;
			set UserDefBeta;
			if &beta ~= beta_look&num_look then do;
				put ERROR: "The user defined values of beta do not add up to &beta";
				ABORT;
			end;
		run;


	/**********Dummy Dataset**********/
		data UserDefBeta3;
			set UserDefBeta(rename=(beta_look1-beta_look&numuser=b1-b&numuser));
		run;
	%end;


	

	/**********ERRORS**********/
	%if (%sysevalf(&mu = 20)) %then %do;
		%put ERROR: Mu is &mu. Rescale median time. For example, use days in place of hours or months in place of days.;
		%ABORT;
	%end;

	%if (%sysevalf(&mu = -20)) %then %do;
		%put ERROR: Mu is &mu. Rescale median time. For example, use days in place of hours or months in place of days.;
		%ABORT;
	%end;

	%if (%sysevalf(&NumSimul <= 0) or %sysevalf(&NumSimul <1000) or %datatyp(&NumSimul)=CHAR) %then %do; *comment out second part to check code for errors by using a small number of simulations;
		%put ERROR: The value "&NumSimul" for NumSimul is not allowed;
		%ABORT;
	%end;

	%if (%sysevalf(&alpha <= 0) or %sysevalf(&alpha >= 1) or %datatyp(&alpha)=CHAR) %then %do;
		%put ERROR: The value "&alpha" for ALPHA is not allowed. ALPHA should between 0 and 1;
		%ABORT;
	%end;

	%if (%sysevalf(&sides ne 1) and %sysevalf(&sides ne 2)) or %datatyp(&sides)=CHAR %then %do;
		%put ERROR: The value "&sides" for SIDES is not allowed. SIDES should be either 1 or 2;
		%ABORT;
	%end;

	%if %sysevalf(%sysfunc(abs(&lambda)) > 5) or %datatyp(&lambda)=CHAR  %then %do;
		%put ERROR: The value "&lambda" for LAMBDA is not allowed. LAMBDA should between -5 and 5;
		%ABORT;
	%end;

	%if %sysevalf(&sigma <=0) or %datatyp(&sigma)=CHAR %then %do;
		%put ERROR: The value "&sigma" for SIGMA is not allowed. SIGMA should be above 0;
		%ABORT;
	%end;

	%if (%sysevalf(%sysfunc(abs(&mu))>100) or %datatyp(&mu)=CHAR) %then %do;
		%put ERROR: The value "&mu" for MU is not allowed. MU should between -100 and 100;
		%ABORT;
	%end;

	%if (%sysevalf(&evt_rate <= 0) or %sysevalf(&evt_rate > 1) or %datatyp(&evt_rate)=CHAR) %then %do;
		%put ERROR: The value "&evt_rate" for EVENT RATE is not allowed. EVENT RATE should be between 0 and 1;
		%ABORT;
	%end;

	%if (%sysevalf(&seed <= 0) or %datatyp(&seed)=CHAR) %then %do;
		%put ERROR: The value "&seed" for Seed is not allowed. Seed should be greater than 0;
		%ABORT;
	%end;

	%if (%sysevalf(&r <= 0) or %sysevalf(&r >= 10) or %datatyp(&r)=CHAR) %then %do;
		%put ERROR: The value "&r" for ALLOCATION RATIO is not allowed. ALLOCATION RATIO should between 0 and 10;
		%ABORT;
	%end;

	%if %datatyp(&Delta_PT)=CHAR %then %do;
		%put ERROR: The value "&Delta_PT" for Seed is not allowed. Delta_PT should be numeric;
		%ABORT;
	%end;

	%if (%sysevalf(&Delta_PT_Ha < &Delta_PT) or %datatyp(&Delta_PT_Ha)=CHAR ) %then %do;
		%put ERROR: The value "&Delta_PT_Ha" for DELTA_PT UNDER THE ALTERNATE HYPOTHESIS is not allowed. Delta_PT_Ha should be >= Delta_PT;
		%ABORT;
	%end;

	%if (%sysevalf(&a <= 0) or %datatyp(&a)=CHAR ) %then %do;
		%put ERROR: The value "&a" for ACCRUAL TIME is not allowed. ACCRUAL TIME should be a positive number;
		%ABORT;
	%end;

	%if (%sysevalf(&a_type ne 1) and %sysevalf(&a_type ne 2)) or %datatyp(&a_type)=CHAR %then %do;
		%put ERROR: The value "&a_type" for ACCRUAL TYPE is not allowed. ACCRUAL TYPE should be either 1 or 2;
		%ABORT;
	%end;

	%if %datatyp(&a_omega)=CHAR %then %do;
		%put ERROR: The value "&a_omega" for Seed is not allowed. a_omega should be numeric;
		%ABORT;
	%end;

	%if (%sysevalf(&t <= &a) or %datatyp(&t)=CHAR ) %then %do;
		%put ERROR: The value "&t" for TOTAL TIME is not allowed. TOTAL TIME  should be a positive integer and should be greater than ACCRUAL TIME;
		%ABORT;
	%end;

	%if (%sysevalf(&num_look <= 0) or %sysfunc(int(&num_look)) ne &num_look or %datatyp(&num_look)=CHAR ) %then %do;
		%put ERROR: The value "&num_look" for NUMBER OF LOOKS is not allowed. NUMBER OF LOOKS should be at least 1 and be an integer;
		%ABORT;
	%end;

	%if (%sysevalf(&look_points ne 1) and %sysevalf(&look_points ne 2)) or %datatyp(&look_points)=CHAR %then %do;
		%put ERROR: The value "&look_points" for LOOK_POINTS is not allowed. LOOK_POINTS should be either 1 or 2;
		%ABORT;
	%end;

	%if (%sysevalf(&alpha_spend ne 1) and %sysevalf(&alpha_spend ne 2) and %sysevalf(&alpha_spend ne 3)) or %datatyp(&alpha_spend)=CHAR %then %do;
		%put ERROR: The value "&alpha_spend" for TYPE OF ALPHA SPENDING FUNCTION is not allowed. TYPE OF ALPHA SPENDING FUNCTION should be either 1 or 2 or 3;
		%ABORT;
	%end;

	%if %sysevalf(&alpha_spend  eq 1) and %sysevalf(&rho <= 0) or %sysevalf(&rho >= 10) or %datatyp(&rho)=CHAR %then %do;
		%put ERROR: If alpha_spend is "&alpha_spend" then the value "&rho" for RHO-EFFICACY PARAMETER should be greater than zero and less than ten;
		%ABORT;
	%end;

	%if %sysevalf(&alpha_spend  eq 2) and %sysevalf(%sysfunc(abs(&rho)) > 10) or %datatyp(&rho)=CHAR %then %do;
		%put ERROR: If alpha_spend is "&alpha_spend" then the value "&rho" for RHO-EFFICACY PARAMETER should be between -10 and 10;;
		%ABORT;
	%end;

	%if (%sysevalf(&beta <= 0) or %sysevalf(&beta >= 1) or %datatyp(&beta)=CHAR) %then %do;
		%put ERROR: The value "&beta" for BETA is not allowed. BETA should be between 0 and 1;
		%ABORT;
	%end;

	%if (%sysevalf(&beta_spend ne 1) and %sysevalf(&beta_spend ne 2) and %sysevalf(&beta_spend ne 3)) or %datatyp(&beta_spend)=CHAR %then %do;
		%put ERROR: The value "&beta_spend" for TYPE OF BETA SPENDING FUNCTION is not allowed. TYPE OF BETA SPENDING FUNCTION should be either 1 or 2 or 3;
		%ABORT;
	%end;


	%if %sysevalf(&beta_spend eq 1) and %sysevalf(&rho_f <= 0) or %sysevalf(&rho_f >= 10) or %datatyp(&rho_f)=CHAR %then %do;
		%put ERROR: If beta_spend is "&beta_spend" then the value "&rho_f" for RHO-EFFICACY PARAMETER should be greater than zero and less than ten;
		%ABORT;
	%end;

	%if %sysevalf(&beta_spend eq 2) and %sysevalf(%sysfunc(abs(&rho_f)) > 10) or %datatyp(&rho_f)=CHAR %then %do;
		%put ERROR: If beta_spend is "&beta_spend" then the value "&rho_f" for RHO-EFFICACY PARAMETER should be between -10 and 10;
		%ABORT;
	%end;

	%if (%sysevalf(&num_skip >= &num_look) or %sysevalf(&num_skip < 0) or %sysfunc(int(&num_skip)) ne &num_skip or %datatyp(&num_skip)=CHAR) %then %do;
		%put ERROR: The value "&num_skip" for NUMBER OF FUTILITY SKIPS is not allowed. NUMBER OF FUTILITY SKIPS should be a non-negative integer and be less than NUMBER OF LOOKS;
		%ABORT;
	%end;



	%macro PT_GSdesign(n);

	options nonotes symbolgen mprint mlogic;
	options missing = "";
	options linesize=100; *adjust ods listing parameters to determine how output tables are shown on your monitor;
	filename myfile "&direct\Mydoc.log"; /*Mydoc.log is stored under the specified file path*/
	proc printto log=myfile new;
	run;
	title "Group Sequential Test with PT assumption";
	ods listing close;
	

	data IGamma; 
		merge UserDefTime;
		k = 1/((&lambda)**2);
		mu2 = &mu + log(&Delta_PT); /*Formula for mu_star can be obtained from Phadnis et al (2020)*/
		mu3 = &mu + log(&Delta_PT_Ha);
		b_shape = abs(&lambda)/&sigma; 				
		Theta = exp(&mu)*(k**(-1/b_shape));
		Theta2 = exp(mu2)*(k**(-1/b_shape));
		Theta3 = exp(mu3)*(k**(-1/b_shape));
		factor_num = &evt_rate**(1/k);
		factor_denom= 1 - (&evt_rate**(1/k));
		factor = (factor_num/factor_denom)**(1/b_shape);
		Cens_theta = Theta*factor;
		Cens_theta2 = Theta2*factor;
		Cens_theta3 = Theta3*factor;

		/* Readers can refer to the original paper by Phadnis et al (2020) Group Sequential design for time-to-event
		data using the concept of PT, published under Statistical Methods in Medical Research OR
		Relia wiki article on the generalized gamma distribution for the pdf and alternative paramaterization:
		http://reliawiki.org/index.php/The_Generalized_Gamma_Distribution*/

		do SimulID = 1 to &NumSimul;
	
			Array LOOK[*] Look1 - Look&Num_Look ;
			Array L_EVT[*] L_evt1 - L_Evt&Num_Look ;
			Array Y[*] Y1 - Y&Num_Look ;
			Array UserTime[&Num_Look] UserTime1-UserTime&Num_Look; 
			
			
			/* This do loop calculates survival times (yi) under the null hypothesis for the standard arm*/
			do j = 1 to &n by 1;

				unif = ranuni(j);										
				arm = 'Std';
				

				/*Readers can refer to of Relia Wiki article for the GG Reliable life formula
				This formula is used to simulate event times from the 3-parameter GG*/

				if abs(&lambda) > 0.1 then 
				do;
					q = quantile('GAMMA',unif,k);
					time_init = exp(&mu + (&sigma/&lambda)*log(q/k));
					time_cens = Cens_theta*((ranexp(&seed+j))**(1/b_shape));				
					time = min(time_init, time_cens);
					if time = time_init then event1 = 1;
					else event1 = 0;
				end;

				if abs(&lambda) <= 0.1 then 
				do;
					q = quantile('normal',unif);
					time_init =  exp(&mu + &sigma*q);
					time_cens = Cens_theta*((ranexp(&seed+j))**(1/b_shape));				
					time = min(time_init, time_cens);
					if time = time_init then event1 = 1;
					else event1 = 0;
				end;

				/*The user-defined accural time follows either a uniform or a truncated exponential with paramter omega*/

				if &a_type = 1 then a_time = &a*ranuni(j) ;
				else a_time = -log(1-ranuni(j)*(1-exp(-&a_omega*&a)))/&a_omega ;
		

				do m=1 to &Num_Look;
				
					if &look_points=1 then Look[m] = m*&t/&num_look;
					else if &look_points=2 then Look[m] = UserTime[m];

					/*This step is for subject who have been recruited amd event occurs within the lookpoint period*/

					if ((a_time < Look[m])  and (time < Look[m] - a_time)) then 
					do;
						L_EVT[m] = 1;
						Y[m] = time;
						if event1 = 0 then L_EVT[m] = 0;
					end;

					/*This step is for subject who have been recruited but event does not occur within lookpoint (censored)*/

					if ((a_time < Look[m])  and (time >= Look[m] - a_time)) then 
					do;
						L_EVT[m] = 0;
						Y[m] = Look[m] - a_time;
						if event1 = 0 then L_EVT[m] = 0;
					end;

					/*This step is for subject who have not been recruited by the lookpoint, */

					if (a_time > Look[m])  then 
					do;
						L_EVT[m] = .;
						Y[m] = .;
					end;
			
				end;
				output;
			end;

			/*This DO loop calculates survival times (yi) under the null hypothesis for the treatment arm*/

			do j = &n+1 to ceil(&n*(1+&r)) by 1;															

				unif2 = ranuni(j);
				arm = 'New';

				if abs(&lambda) > 0.1 then 
				do;
					q = quantile('GAMMA',unif2,k);
					time_init = exp(mu2 + (&sigma/&lambda)*log(q/k));
					time_cens = Cens_theta2*((ranexp(&seed+j))**(1/b_shape));				
					time = min(time_init, time_cens);
					if time = time_init then event2 = 1;
					else event2 = 0;
				end;

				if abs(&lambda) <= 0.1 then 
				do;
					q = quantile('normal',unif2);
					time_init =  exp(mu2 + &sigma*q);
					time_cens = Cens_theta2*((ranexp(&seed+j))**(1/b_shape));				
					time = min(time_init, time_cens);
					if time = time_init then event2 = 1;
					else event2 = 0;
				end;

				if &a_type = 1 then a_time = &a*ranuni(j) ;
				else a_time = -log(1-ranuni(j)*(1-exp(-&a_omega*&a)))/&a_omega ;

				do i=1 to &Num_Look;
	
					if &look_points=1 then Look[i] = i*&t/&num_look;
					else if &look_points=2 then Look[i] = UserTime[i];

					if ((a_time < Look[i])  and (time < Look[i] - a_time)) then 
					do;
						L_EVT[i] = 1;
						Y[i] = time;
						if event2 = 0 then L_EVT[i] = 0;
					end;
					if ((a_time < Look[i])  and (time >= Look[i] - a_time)) then 
					do;
						L_EVT[i] = 0;
						Y[i] = Look[i] - a_time;
						if event2 = 0 then L_EVT[i] = 0;
					end;
					if (a_time > Look[i])  then 
					do;
						L_EVT[i] = .;
						Y[i] = .;
					end;
			
				end;

			output;
			end;
		

			/*This DO loop calculates survival times (yi) under the alternative hypothesis*/
			do j = &n*(1+&r)+1 to ceil(&n*(1+&r*2)) by 1;													

				unif3 = ranuni(j);
				arm = 'HA';

				if abs(&lambda) > 0.1 then 
				do;
					q = quantile('GAMMA',unif3,k);
					time_init = exp(mu3 + (&sigma/&lambda)*log(q/k));
					time_cens = Cens_theta3*((ranexp(&seed+j))**(1/b_shape));				
					time = min(time_init, time_cens);
					if time = time_init then event3 = 1;
					else event3 = 0;
				end;

				if abs(&lambda) <= 0.1 then 
				do;
					q = quantile('normal',unif3);
					time_init =  exp(mu3 + &sigma*q);
					time_cens = Cens_theta3*((ranexp(&seed+j))**(1/b_shape));				
					time = min(time_init, time_cens);
					if time = time_init then event3 = 1;
					else event3 = 0;
				end;

				if &a_type = 1 then a_time = &a*ranuni(j) ;
				else a_time = -log(1-ranuni(j)*(1-exp(-&a_omega*&a)))/&a_omega ;

				do i=1 to &Num_Look;
	
					if &look_points=1 then Look[i] = i*&t/&num_look;
					else if &look_points=2 then Look[i] = UserTime[i];

					if ((a_time < Look[i])  and (time < Look[i] - a_time)) then 
					do;
						L_EVT[i] = 1;
						Y[i] = time;
						if event3 = 0 then L_EVT[i] = 0;
					end;
					if ((a_time < Look[i])  and (time >= Look[i] - a_time)) then 
					do;
						L_EVT[i] = 0;
						Y[i] = Look[i] - a_time;
						if event3 = 0 then L_EVT[i] = 0;
					end;
					if (a_time > Look[i])  then 
					do;
						L_EVT[i] = .;
						Y[i] = .;
					end;
			
				end;

				output;
			end;

		end;
	run;


	/* This dataset contains survival times for the control arm and new arm under the null hypothesis*/
	data IGamma0 (where = (arm NE 'HA'));
	set IGamma ;
	run;

	
	/* This dataset contains the survival times for the control arm under the null, and survival times under the alternative*/
	data IGamma1 (where = (arm NE 'New'));
	set IGamma ;
	run;

	

	%macro GrpSeq;
	ods listing close;

	%do s = 1 %to &num_look;

		/* This step produces a frequency table by treatment arm and events under each look point*/
		ods output CrossTabFreqs = Freq_est&s;
		proc freq data = IGamma;
			by SimulID ;
			tables arm*L_EVT&s;
		run;

		/* This steps isolates survival times under the control arm*/
		data Count_evt_std&s (where = (L_EVT&s = 1 and arm = 'Std'));
		set Freq_est&s ;
		by SimulID;
			n1_interim = frequency;
			keep SimulID L_EVT&s arm n1_interim;
		run;

		/* This step isolates treatment times under the null*/
		data Count_evt_new&s (where = (L_EVT&s = 1 and arm = 'New'));
		set Freq_est&s ;
		by SimulID;
			n2_interim = frequency;
			keep SimulID L_EVT&s arm n2_interim;
		run;

		
		/*This step isolates treatment times under the alternative*/
		data Count_evt_HA&s (where = (L_EVT&s = 1 and arm = 'HA'));
		set Freq_est&s ;
		by SimulID;
			n3_interim = frequency;
			keep SimulID L_EVT&s arm n3_interim;
		run;

		
		/*Calculations for the JT and HSD for different values of rho are undertaken*/
		data count_evt&s (keep = SimulID LookNumber n1_interim n2_interim n3_interim Prop_0 Prop_1 Alpha_Look&s Beta_Look&s);
		merge Count_evt_std&s Count_evt_new&s Count_evt_HA&s;
		by SimulID;
			LookNumber = &s;

			/* This is the fraction of information observed at the jth interim analysis for the alpha spending function
			used in calculating the efficiacy boundary. Refer Phadnis et al (2020)*/
			Prop_0 = &s/&num_look;

			if LookNumber <= &num_skip then Prop_1 = 0 ;
			/* prop_1 is fraction of information at jth interim analysis for beta spending function
			used for calculating the futility boundary. Refer Phadnis et al (2020)*/
			else Prop_1 = &s/&num_look;

			if &alpha_spend = 1 then
			do;
				Alpha_Look&s = (&alpha/&sides)*(Prop_0**&rho);
			end;
			if &alpha_spend = 2 then
			do;
				if &rho = 0 then Alpha_Look&s = (&alpha/&sides)*Prop_0 ;
				else Alpha_Look&s = (&alpha/&sides)*((1-exp(-&rho*Prop_0))/(1-exp(-&rho)));
			end;
			%if &alpha_spend = 3 %then
			%do;
				Alpha_Look&s = &&alpha_look&s;			
			%end;
			
	
			if &beta_spend = 1 then										
			do;
				Beta_Look&s = (&beta)*(Prop_1**&rho_f);
			end;
			if &beta_spend = 2 then
			do;
				if &rho_f = 0 then Beta_Look&s = (&Beta)*Prop_1 ;
				else Beta_Look&s = (&beta)*((1-exp(-&rho_f*Prop_1))/(1-exp(-&rho_f)));
			end;
			%if &beta_spend = 3 %then
			%do;
				Beta_Look&s = &&beta_look&s;			
			%end;			
		run;


		proc sort data = count_evt&s;
		by SimulID;
		run;
	
		ods output Summary = Evt_interim&s;
		proc means data = count_evt&s mean;
			var n1_interim n2_interim n3_interim;
		run;

		/*The mean of events under each arm is calculated for each look point*/
		data Evt_interim&s (keep = Evt_n1_H0 Evt_n2_H0 Evt_n2_HA);
		set Evt_interim&s;
			Evt_n1_H0 = n1_interim_mean;
			Evt_n2_H0 = n2_interim_mean;
			Evt_n2_HA = n3_interim_mean;
		run;
	%end;

	data Alpha_all (keep = SimulID Alpha_Look1-Alpha_Look&num_look );
	merge count_evt1-count_evt&num_look;
	by SimulID;
	run;

	/* This steps calculates the incremental increase in alpha from one lookpoint to another*/
	data Alpha_all2;
	set Alpha_all;
		Array ALPHA_L[&num_look] Alpha_Look1 - Alpha_Look&Num_Look ;
		Array ALPHA_T[*] ALPHA_T1 - ALPHA_T&Num_Look ;
		ALPHA_T[1] = Alpha_Look1;
		do q = 2 to &num_look;
			ALPHA_T[q] = ALPHA_L[q] - ALPHA_L[q-1];
		end;
	run;

	
	/* This step calculates the incremental increase in beta from one lookpoint to another*/
	data Beta_all (keep = SimulID Beta_Look1-Beta_Look&num_look Beta_Look1-Beta_Look&num_look);
	merge count_evt1-count_evt&num_look;
	by SimulID;
	run;

	data Beta_all2;
	set Beta_all;
		Array BETA_L[&num_look] Beta_Look1 - Beta_Look&Num_Look ;
		Array BETA_T[*] BETA_T1 - BETA_T&Num_Look ;
		BETA_T[1] = Beta_Look1;
		do q = 2 to &num_look;
			BETA_T[q] = BETA_L[q] - BETA_L[q-1];
		end;
	run;


	/* Efficacy and futility boundaries are calculated. Efficiacy boundaries are obtained by calculating the change in proportional time
		under the null hypothesis. Values are exponentiated for each study arm. Futility boundaries are calculated by looking at the change
		in proportional times under the alternative hypothesis. Values are also exponentiated.*/


	%do s = 1 %to &Num_Look ;

		%let u = %eval(&s-1) ;

		/*  If w=0 then we are looking at efficacy boundaries. Only looking at alpha spending.
			If w=1 then we are looking at futility boundaries. Only looking at beta spending.
			Use proc lifereg to obtain parameter estimates by arm for each simulation ID.
			ParamTable will save parameter estimates for efficacy at each lookpoint. */

		/* Use LIFEREG to find MLEs of GG distribution parameters */
		%do w = 0 %to 1;
			ods output ParameterEstimates = ParamTable&w&s;
			proc lifereg data = IGamma&w ;										
			by SimulID;
				class  arm ;
				model Y&s*L_evt&s(0)= arm/dist = &MyDist maxiter=&maxiter CONVG=&convg ; /* MAXITER may come into play when dealing with small sample sizes */		
			run;																				


			/* Only want instaces where parameter = arm and DF=1.
			This step will save Mu1_hat-Mu0_hat into a table*/
			data Table_Mu&w&s ;													
			set ParamTable&w&s (where =(Parameter='arm' and DF=1));
				Est_Q = exp(Estimate);
				est_se_mu = StdErr;
			run;
			
			/* Merge tables. Keep only alpha spending and beta spending relevant to each look point.
			Both done at the same time using indicator function. */
			data Calc&w&s ;
			merge count_evt&s Table_Mu&w&s Alpha_all2 Beta_all2;
			run;

			data Calc&w&s(keep = SimulID LookNumber Prop_0 Prop_1 Alpha_Look&s ALPHA_T&s Est_Q 
								F_prob&w Beta_Look&s Beta_T&s);
			set Calc&w&s;
				F_prob&w = ((1-ALPHA_T&s)**(1-&w))*((1-BETA_T&s)**&w);										
			run;

			/*Sort so we obtain cut points that are below the F_prob*/
			proc sort data = Calc&w&s;
			by Est_Q;
			run;

			proc sort data = Table_Mu&w&s ;
			by Est_Q;
			run;		
		%end;

	
		%if %eval(&s) = 1 %then 
		%do;
			%do w = 0 %to 1;
				data NewCalc&w&s; /*If s=1, then u=0 and there is no cutting*/
				set Calc&w&s;
				run;
			%end;
		%end;

		%if %eval(&s) > 1 %then 
		%do;
			%do w = 0 %to 1;
				data Cut&w&u ;
				set NewCalc&w&u nobs = total&w&u;
					if (&w=1 and Beta_T&u = 0) then 
					do;
						Est_0 = 0;
						call symput ('Trim_Q',Est_0);
					end;

					else
					do;
						if _n_ = ((total&w&u - ceil(&NumSimul*Alpha_T&u))**(1-&w))*((ceil(&NumSimul*max(0.00001,Beta_T&u)))**&w); 
						call symput('Trim_Q', Est_Q);
					end;
				run;



			/*
			Under efficacy and futility, calculate the data to be cut at each lookpoint.
			Lagged NewCalc is set since NewCalc is created in the previous step.
			If futility and beta spending=0 then there is no futility boundary at the look point.
			Does not apply for efficiacy case.
			*/

				data Temp&w&u (keep=SimulID trim&w&u);
				set Calc&w&u;
					if ((-1)**&w)*(Est_Q) >= (1-2*&w)*(&Trim_Q) then trim&w&u = . ;
					else trim&w&u = 1; *Trimwu =1 are all marked for cutting.;
				run;

				proc sort data = Temp&w&u;
				by SimulID;
				run;

				proc sort data = Calc&w&s;
				by SimulID;
				run;

			%end;


			%do w=0 %to 0; /*binding futility*/

				%if &bind=0 %then
				%do;
					data NewCalc&w&s;
					merge Calc&w&s Temp01-Temp0&u;
					by SimulID; 
						if cmiss (of _all_) then delete;
					run;
				%end;

				%if &bind=1 %then
				%do;
					data NewCalc&w&s;
						merge Calc&w&s Temp01-Temp0&u Temp11-Temp1&u;
						by SimulID; 
							if cmiss (of _all_) then delete;
					run;
				%end;
			%end;

			%do w=1 %to 1; /*non-binding futility*/

				data NewCalc&w&s;
				merge Calc&w&s Temp01-Temp0&u Temp11-Temp1&u;
				by SimulID; 
					if cmiss (of _all_) then delete;
				run;

			%end;
		
			%do w = 0 %to 1;
				proc sort data = NewCalc&w&s;
				by Est_Q;
				run;

				proc sort data = Calc&w&s; 
				by Est_Q;
				run;
			%end;
		
		%end;





		/* Creating nice-looking output with prefix MyOutput*/
		%do w = 0 %to 1;
			data MyOutput&w&s (keep= LookNumber Prop_&w Est_Q Alpha_Spent Cumul_Alpha_Spent F_prob&w f_crit_look&w 
									 Q_crit_look&w check&w Beta_Spent Cumul_Beta_Spent);
			set NewCalc&w&s nobs=total_obs&w&s;

			if (&w=1 and Beta_T&s = 0) then
			do;
				Est_Q = 0;
				Q_crit_look1 = 0;
				f_crit_look1 = 0;
				check1 = 0;
				if _n_ = 1;
			end;
			else
			do;
				check&w = ((total_obs&w&s - ceil(&NumSimul*Alpha_T&s))**(1-&w))*((ceil(&NumSimul*max(0.00001,Beta_T&s)))**&w);  
				if check&w <= _N_ <= check&w;

				Q_crit_look&w = Est_Q;
				f_crit_look&w = round(Q_crit_look&w**(abs(&lambda)/&sigma),0.001);
			end;

			Alpha_Spent = ALPHA_T&s;
			Cumul_Alpha_Spent = Alpha_Look&s;

			Beta_Spent = BETA_T&s;
			Cumul_Beta_Spent = Beta_Look&s;

			run;
		%end;

		data MyOutput&s;
		merge MyOutput0&s MyOutput1&s Evt_interim&s;
		run;

		data StopEFF&s;
		set MyOutput&s;
			call symput('Eff_cut', Q_crit_look0);
		run;

		data StopFUT&s;
		set MyOutput&s;
			call symput('Fut_cut', Q_crit_look1);
		run;

		%do w = 0 %to 1;

			data Table_Mu&w&s;
			set Table_Mu&w&s;
				mark = _n_/&NumSimul ;
			run;

			data Upper&w&s (keep = mark_up&w);
			set Table_Mu&w&s;		
				if ((Est_Q > &Eff_Cut) and (lag(Est_Q)<= &Eff_Cut));
				mark_up&w = mark;
			run;

			data Lower&w&s (keep = mark_down&w);
			set Table_Mu&w&s;
					if ((Est_Q >= &Fut_Cut) and (lag(Est_Q)< &Fut_Cut));
					mark_down&w = mark - 1/&NumSimul;
			run;
		%end;

		data MyOutput&s;
		merge MyOutput&s Upper0&s Lower0&s Upper1&s Lower1&s;
			Cumul_Stop_Prob_H0 = (1 - mark_up0)+ mark_down0 ;
			Cumul_Stop_Prob_HA = (1 - mark_up1)+ mark_down1 ;
		run;

		data LT&s (keep = LookTimes);
		set IGamma;
			if _n_ = 1;
			if &look_points = 1 then LookTimes = Look&s;
			if &look_points = 2 then LookTimes = UserTime&s;
		run;

	%end;

	data LTpoints;
	set LT1-LT&num_look;
	run;

	proc sort data=IGamma1;
	by arm;
	run;

	ods output Summary = Sum_Table;
	proc means data = IGamma1 sum;
	by arm;
		var Y1-Y&num_look;
	run;

	data Cum_Time_h0 (where = (arm='Std')); 
	set Sum_Table;	
	run;

	proc transpose data = Cum_Time_h0 out=Time_H0;
	run;

	data Cum_Time_ha (where = (arm='HA')); 
	set Sum_Table;	
	run;

	proc transpose data = Cum_Time_ha out=Time_HA;
	run;

	data Time_H0 (keep = Cum_Sub_Time_H0);
	set Time_H0;
		Cum_Sub_Time_H0 = COL1/&NumSimul;
	run;
	data Time_HA (keep = Cum_Sub_Time_HA);
	set Time_HA;
		Cum_Sub_Time_HA = COL1/&NumSimul;
	run;

	data FinalOutput;
	set MyOutput1-MyOutput&num_look ;
		df1_look = 2*Evt_n1_H0/((&lambda)**2);
		df2_look = 2*Evt_n2_H0/((&lambda)**2);
		df3_look = 2*Evt_n2_HA/((&lambda)**2);
		Pval_crit0 = round(1-cdf('F',f_crit_look0,df1_look,df2_look),0.0001);
		Pval_crit1 = round(1-cdf('F',f_crit_look1,df1_look,df3_look),0.0001);
		F_prob_cum = 1-Cumul_Alpha_Spent;

		if LookNumber = &num_look then 
		do;
			Q_fixed = (quantile('F',F_prob_cum,df1_look,df2_look))**(&sigma/abs(&lambda));
			SampleSize = &n;
		end;
		else 
		do;
			Q_fixed = .;
			SampleSize = .;
		end;
	
	run;

	data FinalOutput;
	merge FinalOutput Time_H0 Time_HA LTpoints;
	
		lag_0 = lag(Cumul_Stop_Prob_H0);
		if LookNumber = 1 then Stop_Prob_H0 = Cumul_Stop_Prob_H0;
		else Stop_Prob_H0 = Cumul_Stop_Prob_H0 - lag_0;

		lag_1 = lag(Cumul_Stop_Prob_HA);
		if LookNumber = 1 then Stop_Prob_HA = Cumul_Stop_Prob_HA;
		else Stop_Prob_HA = Cumul_Stop_Prob_HA - lag_1;

		if f_crit_look1 = 0 then
		do;
			f_crit_look1 = .;
			Q_crit_look1 = .;
			Pval_crit1 = .;
		end;
	run;


	ods listing ;

	data FinalOutput;
	set FinalOutput;

		Avg_n_H0 = round((Evt_n1_H0 + Evt_n2_H0)*Stop_Prob_H0/2,0.01);
		Avg_n1_Ha = round(Evt_n1_H0*Stop_Prob_Ha,0.01);
		Avg_n2_Ha = round(Evt_n2_Ha*Stop_Prob_Ha,0.01);

		lag_1 = lag(Cumul_Stop_Prob_HA);
		if LookNumber = 1 then AvgEvtNum_fixedN = round(((Evt_n1_H0 + Evt_n2_Ha)/2)*(min(1,Stop_Prob_Ha)),0.001) ;
		if LookNumber = &num_look then AvgEvtNum_fixedN = round(((Evt_n1_H0 + Evt_n2_Ha)/2)*(1-min(1,lag_1)),0.001)  ;
		if LookNumber > 1 and LookNumber < &num_look then AvgEvtNum_fixedN = round(((Evt_n1_H0 + Evt_n2_Ha)/2)*(min(1,Cumul_Stop_Prob_HA)-min(1,lag_1)),0.001)  ;

		Evt_n1_H0 = round(Evt_n1_H0,0.01);
		Evt_n2_H0 = round(Evt_n2_H0,0.01);
		Evt_n1_Ha = round(Evt_n1_Ha,0.01);
		Evt_n2_Ha = round(Evt_n2_Ha,0.01);

		Alpha_Spent = round(Alpha_Spent,0.00001);
		Cumul_Alpha_Spent = round(Cumul_Alpha_Spent,0.00001);
		Q_crit_look0 = round(Q_crit_look0,0.001);
		Stop_Prob_H0 = round(Stop_Prob_H0,0.0001);
		Cumul_Stop_Prob_H0 = round(Cumul_Stop_Prob_H0,0.0001);
		Q_fixed = round(Q_fixed,0.001);

		Beta_Spent = round(Beta_Spent,0.00001);
		Cumul_Beta_Spent = round(Cumul_Beta_Spent,0.00001);
		Q_crit_look1 = round(Q_crit_look1,0.001);
		Stop_Prob_Ha = round(Stop_Prob_Ha,0.0001);
		Cumul_Stop_Prob_Ha = round(Cumul_Stop_Prob_Ha,0.0001);

		Cum_Sub_Time_H0 = round(Cum_Sub_Time_H0,0.0001);
		Cum_Sub_Time_Ha = round(Cum_Sub_Time_Ha,0.0001);

	run;

	data LastQ0 (keep = Q_crit_look0);					
	set FinalOutput;
		if _n_ = &num_look;
		call symputx('last_e', Q_crit_look0,'G');
	run;

	data LastQ1 (keep = Q_crit_look1);					
	set FinalOutput;
	if _n_ = &num_look;
	call symputx('last_f', Q_crit_look1,'G');
	run;

	title 'Average Number of Events under H0';
	proc means data = FinalOutput sum;
		var Avg_n_H0;
	run;

	title 'Average Number of Events N1 under HA';
	proc means data = FinalOutput sum;
		var Avg_n1_Ha;
	run;

	title 'Average Number of Events N2 under HA';
	proc means data = FinalOutput sum;
		var Avg_n2_Ha;
	run;

	title 'Average Number of Events (both arms combined) for Fixed N and Effect Size Delta';
	proc means data = FinalOutput sum;
		var AvgEvtNum_fixedN;
	run;

	title 'Group Sequential Design using Proportional Time';
	proc print data = FinalOutput;
		var LookNumber LookTimes Evt_n1_H0 Evt_n2_H0 Evt_n2_HA Alpha_Spent Cumul_Alpha_Spent Q_crit_look0 Check0 Stop_Prob_H0 
			Cumul_Stop_Prob_H0 Q_fixed Beta_Spent Cumul_Beta_Spent Q_crit_look1 Check1 Stop_Prob_Ha Cumul_Stop_Prob_Ha Cum_Sub_Time_H0 
			Cum_Sub_Time_HA SampleSize;
	run;


	%mend;

	dm 'log;
	clear; 
	';

	%GrpSeq;

	data final_N;
	set FinalOutput;
		if _n_ = &num_look;
		call symputx('N_final',put(SampleSize,8.),'G');
	run;

%mend;

%PT_GSdesign(&n1);



************************************************************************************************************************************
/************************ Doubling-Halving algorithm starts here: *****************************************************************/;

/* Find the smallest samples size n such that Q0 and Q1 under each look point are equal*/
%macro Search_N ;
	%let switch = 0;

	%do v = 1 %to 50;
	
		%if %sysevalf(&last_e > &last_f) %then
		%do;
			%if &v = %eval(&count + 1)  %then 
			%do;
				%let n2 = %eval(2*&n1) ;
			%end;
			%else %do;
				%let n2 = %sysevalf(%sysevalf((&n1 + &m2)/2),ceil);
			%end;
			%PT_GSdesign(&n2);
			%if %sysevalf((%sysfunc(abs(%eval(&n2) - %eval(&n1)))) < 2) and (&last_e <= &last_f) %then 
			%do;
				%return;
			%end;
			%else %do;
				%let m1 = &n1;
				%let n1 = &n2;
				%let count = %eval(&count + 1);
			%end;		
		%end;

		%if %sysevalf(&last_e <= &last_f) and %sysevalf(&v ne %eval(&count)) %then
		%do;
			%if (&v = %eval(&num + 1)) and %sysevalf(%eval(&switch) = 0) %then 
			%do;
				%let n2 = %sysevalf((&n1/2),ceil) ;
			%end;
			%else %do;
				%let n2 = %sysevalf(%sysevalf((&n1 + &m1)/2),ceil);
			%end;
			%PT_GSdesign(&n2);
			%if %sysevalf((%sysfunc(abs(%eval(&n2) - %eval(&n1)))) < 2) and (&last_e <= &last_f) %then 
			%do;
				%return;
			%end;
			%if %sysevalf((%sysfunc(abs(%eval(&n2) - %eval(&n1)))) >= 2) and (&last_e <= &last_f) %then
			%do;
				%let m2 = &n2;
				%let n1 = &n2;
				%let num = %eval(&num + 1);
			%end;
			%if (&last_e > &last_f) %then
			%do;
				%let m2 = &n1;
				%let n1 = &n2;
				%let switch = 1;
				%let num = %eval(&num + 1);
			%end;

		%end;

	%end;

%mend;

%Search_N;


/***************************************** Doubling-Halving algorithm ends here****************************************************/
************************************************************************************************************************************


************************************************************************************************************************************
/**************************************** The macro to adjust Type II (beta) error starts here - 09/02/2018 ************************************/;

/* Our intention is to get Q0 = Q1 under each lookpoint. Keeping n fixed, decrease the type II error until a four-digit decimal value of
beta is found that yields Q0m = Q1m. This completes the search for n and beta.*/

%macro BetaAdjust(b);
ods listing close;

	%do s = 1 %to &num_Look ;
		data Beta&s (keep=Beta_Look&s Beta_T&s);
		set Calc0&s;
		run;

		%do w = 0 %to 1;
			data Calc&w&s (drop = Beta_Look&s Beta_T&s);
			set Calc&w&s;
			run;
		%end;
	%end;

	data Beta_Only;
	merge Beta1-Beta&num_look;
	run;



	data Beta_all_new(drop=x1-x&num_look);
	merge Beta_Only UserDefBeta3;								
		
		%do j=1 %to &num_look;									
			retain X&j;
			if not missing(b&j) then X&j=b&j;
			b&j = X&j;
			if missing(X&j) then delete;
		%end;

		Array new_prop[&num_look] new_prop1 - new_prop&num_look;
		Array BETA_L[&num_look] Beta_Look1 - Beta_Look&Num_Look ;
		Array BETA_T[*] BETA_T1 - BETA_T&Num_Look ;
		Array B[&num_look]	B1 - B&num_look;						

		do i = 1 to &num_look;

			if i <= &num_skip then new_prop[i] = 0 ;
			else new_prop[i] = i/&num_look;

			if &beta_spend = 1 then										
			do;
				Beta_L[i] = (&b)*(new_prop[i]**&rho_f);
			end;
			if &beta_spend = 2 then
			do;
				if &rho_f = 0 then Beta_L[i] = (&b)*new_prop[i] ;
				else Beta_L[i] = (&b)*((1-exp(-&rho_f*new_prop[i]))/(1-exp(-&rho_f)));
			end;
			if &beta_spend = 3 then 					
			do;
				Beta_L[i] = (&b)*(B[i]/&beta);
			end;
		end;

		BETA_T[1] = Beta_Look1;
		do q = 2 to &num_look;
			BETA_T[q] = BETA_L[q] - BETA_L[q-1];		
		end;

	drop ERROR _NAME_ b1 - b&num_look;					

	run;


	%do s = 1 %to &Num_Look ;
		%do w = 0 %to 1;
			data Calc&w&s;
			merge Calc&w&s Beta_all_new ;
			run;
		%end;
	%end;

	%do s = 1 %to &Num_Look ;

		%let u = %eval(&s-1) ;
	
		%if %eval(&s) = 1 %then 
		%do;
			%do w = 0 %to 1;
				data NewCalc&w&s;
				set Calc&w&s;
				run;
			%end;
		%end;

		%if %eval(&s) > 1 %then 
		%do;
			%do w = 0 %to 1;
				data Cut&w&u ;
				set NewCalc&w&u nobs = total&w&u;
					if (&w=1 and Beta_T&u = 0) then 
					do;
						Est_0 = 0;
						call symput ('Trim_Q',Est_0);
					end;
					else
					do;
						if _n_ = ((total&w&u - ceil(&NumSimul*Alpha_T&u))**(1-&w))*((ceil(&NumSimul*max(0.00001,Beta_T&u)))**&w);
						call symput('Trim_Q', Est_Q);
					end;
				run;

				data Temp&w&u (keep=SimulID trim&w&u);
				set Calc&w&u;
					if ((-1)**&w)*(Est_Q) >= (1-2*&w)*(&Trim_Q) then trim&w&u = . ;
					else trim&w&u = 1;
				run;

				proc sort data = Temp&w&u;
				by SimulID;
				run;

				proc sort data = Calc&w&s;
				by SimulID;
				run;

			%end;
		
			%do w=0 %to 0; /*binding futility*/

				%if &bind=0 %then
				%do;
					data NewCalc&w&s;
					merge Calc&w&s Temp01-Temp0&u;
					by SimulID; 
						if cmiss (of _all_) then delete;
					run;
				%end;

				%if &bind=1 %then
				%do;
					data NewCalc&w&s;
						merge Calc&w&s Temp01-Temp0&u Temp11-Temp1&u;
						by SimulID; 
							if cmiss (of _all_) then delete;
					run;
				%end;
			%end;

			%do w=1 %to 1; /*non-binding futility*/

				data NewCalc&w&s;
				merge Calc&w&s Temp01-Temp0&u Temp11-Temp1&u;
				by SimulID; 
					if cmiss (of _all_) then delete;
				run;

			%end;
		
			%do w = 0 %to 1;
				proc sort data = NewCalc&w&s;
				by Est_Q;
				run;

				proc sort data = Calc&w&s; 	
				by Est_Q;
				run;
			%end;
		
		%end;

		%do w = 0 %to 1;
			data MyOutput&w&s (keep= LookNumber Prop_&w Est_Q Alpha_Spent Cumul_Alpha_Spent F_prob&w f_crit_look&w 
									 Q_crit_look&w check&w Beta_Spent Cumul_Beta_Spent); 
			set NewCalc&w&s nobs=total_obs&w&s;

				if (&w=1 and Beta_T&s = 0) then
				do;
					Est_Q = 0;
					Q_crit_look1 = 0;
					f_crit_look1 = 0;
					check1 = 0;
					if _n_ = 1;
					end;
				else
				do;
					check&w = ((total_obs&w&s - ceil(&NumSimul*Alpha_T&s))**(1-&w))*((ceil(&NumSimul*max(0.00001,Beta_T&s)))**&w);  
					if check&w <= _N_ <= check&w;

					Q_crit_look&w = Est_Q;
					f_crit_look&w = round(Q_crit_look&w**(abs(&lambda)/&sigma),0.001);
				end;

				Alpha_Spent = ALPHA_T&s;
				Cumul_Alpha_Spent = Alpha_Look&s;

				Beta_Spent = BETA_T&s;
				Cumul_Beta_Spent = Beta_Look&s;

			run;
		%end;

		data MyOutput&s;
		merge MyOutput0&s MyOutput1&s Evt_interim&s;
		run;

		data StopEFF&s;
		set MyOutput&s;
			call symput('Eff_cut', Q_crit_look0);
		run;

		data StopFUT&s;
		set MyOutput&s;
			call symput('Fut_cut', Q_crit_look1);
		run;

		%do w = 0 %to 1;

			data Table_Mu&w&s;
			set Table_Mu&w&s;
				mark = _n_/&NumSimul ;
			run;

			data Upper&w&s (keep = mark_up&w);
			set Table_Mu&w&s;		
				if ((Est_Q > &Eff_Cut) and (lag(Est_Q)<= &Eff_Cut));
				mark_up&w = mark;
			run;

			data Lower&w&s (keep = mark_down&w);
			set Table_Mu&w&s;
				if ((Est_Q >= &Fut_Cut) and (lag(Est_Q)< &Fut_Cut));
				mark_down&w = mark - 1/&NumSimul;
			run;
		%end;

		data MyOutput&s;
		merge MyOutput&s Upper0&s Lower0&s Upper1&s Lower1&s;
			Cumul_Stop_Prob_H0 = (1 - mark_up0)+ mark_down0 ;
			Cumul_Stop_Prob_HA = (1 - mark_up1)+ mark_down1 ;
		run;

		data LT&s (keep = LookTimes);
		set IGamma;
			if _n_ = 1;
			if &look_points = 1 then LookTimes = Look&s;
			if &look_points = 2 then LookTimes = UserTime&s;
		run;

	%end;

	data LTpoints;
	set LT1-LT&num_look;
	run;

	proc sort data=IGamma1;
	by arm;
	run;

	ods output Summary = Sum_Table;
	proc means data = IGamma1 sum;
		by arm;
		var Y1-Y&num_look;
	run;

	data Cum_Time_h0 (where = (arm='Std')); 
	set Sum_Table;	
	run;

	proc transpose data = Cum_Time_h0 out = Time_H0;
	run;

	data Cum_Time_ha (where = (arm='HA')); 
	set Sum_Table;	
	run;

	proc transpose data = Cum_Time_ha out=Time_HA;
	run;

	data Time_H0 (keep = Cum_Sub_Time_H0);
	set Time_H0;
		Cum_Sub_Time_H0 = COL1/&NumSimul;
	run;

	data Time_HA (keep = Cum_Sub_Time_HA);
	set Time_HA;
		Cum_Sub_Time_HA = COL1/&NumSimul;
	run;

	data FinalOutput;
	set MyOutput1-MyOutput&num_look ;
		df1_look = 2*Evt_n1_H0/((&lambda)**2);
		df2_look = 2*Evt_n2_H0/((&lambda)**2);
		df3_look = 2*Evt_n2_HA/((&lambda)**2);
		Pval_crit0 = round(1-cdf('F',f_crit_look0,df1_look,df2_look),0.0001);
		Pval_crit1 = round(1-cdf('F',f_crit_look1,df1_look,df3_look),0.0001);
		F_prob_cum = 1-Cumul_Alpha_Spent;

		if LookNumber = &num_look then 
		do;
			Q_fixed = (quantile('F',F_prob_cum,df1_look,df2_look))**(&sigma/abs(&lambda));
			SampleSize = &N_final;
		end;
		else do;
			Q_fixed = .;
			SampleSize = .;
		end;
	
	run;

	data FinalOutput;
	merge FinalOutput Time_H0 Time_HA LTpoints;
		lag_0 = lag(Cumul_Stop_Prob_H0);
		if LookNumber = 1 then Stop_Prob_H0 = Cumul_Stop_Prob_H0;
		else Stop_Prob_H0 = Cumul_Stop_Prob_H0 - lag_0;

		lag_1 = lag(Cumul_Stop_Prob_HA);
		if LookNumber = 1 then Stop_Prob_HA = Cumul_Stop_Prob_HA;
		else Stop_Prob_HA = Cumul_Stop_Prob_HA - lag_1;

		if f_crit_look1 = 0 then
		do;
			f_crit_look1 = .;
			Q_crit_look1 = .;
			Pval_crit1 = .;
		end;
	run;


	ods listing ;

	data FinalOutput;
	set FinalOutput;

		Avg_n_H0 = round((Evt_n1_H0 + Evt_n2_H0)*Stop_Prob_H0/2,0.01);
		Avg_n1_Ha = round(Evt_n1_H0*Stop_Prob_Ha,0.01);
		Avg_n2_Ha = round(Evt_n2_Ha*Stop_Prob_Ha,0.01);

		Evt_n1_H0 = round(Evt_n1_H0,0.01);
		Evt_n2_H0 = round(Evt_n2_H0,0.01);
		Evt_n1_Ha = round(Evt_n1_Ha,0.01);
		Evt_n2_Ha = round(Evt_n2_Ha,0.01);

		Alpha_Spent = round(Alpha_Spent,0.00001);
		Cumul_Alpha_Spent = round(Cumul_Alpha_Spent,0.00001);
		Q_crit_look0 = round(Q_crit_look0,0.001);
		Stop_Prob_H0 = round(Stop_Prob_H0,0.0001);
		Cumul_Stop_Prob_H0 = round(Cumul_Stop_Prob_H0,0.0001);
		Q_fixed = round(Q_fixed,0.001);

		Beta_Spent = round(Beta_Spent,0.00001);
		Cumul_Beta_Spent = round(Cumul_Beta_Spent,0.00001);
		Q_crit_look1 = round(Q_crit_look1,0.001);
		Stop_Prob_Ha = round(Stop_Prob_Ha,0.0001);
		Cumul_Stop_Prob_Ha = round(Cumul_Stop_Prob_Ha,0.0001);

		Cum_Sub_Time_H0 = round(Cum_Sub_Time_H0,0.0001);
		Cum_Sub_Time_Ha = round(Cum_Sub_Time_Ha,0.0001);

	run;

	data LastQ0 (keep = Q_crit_look0);					
	set FinalOutput;
		if _n_ = &num_look;
		call symput('last_e', Q_crit_look0);
	run;

	data LastQ1 (keep = Q_crit_look1);					
	set FinalOutput;
		if _n_ = &num_look;
		call symput('last_f', Q_crit_look1);
	run;


	title 'Average Number of Events under H0';
	proc means data = FinalOutput sum;
		var Avg_n_H0;
	run;

	title 'Average Number of Events N1 under HA';
	proc means data = FinalOutput sum;
		var Avg_n1_Ha;
	run;

	title 'Average Number of Events N2 under HA';
	proc means data = FinalOutput sum;
		var Avg_n2_Ha;
	run;

	title 'Group Sequential Design using Proportional Time';
	proc print data = FinalOutput;
		var LookNumber LookTimes Evt_n1_H0 Evt_n2_H0 Evt_n2_HA Alpha_Spent Cumul_Alpha_Spent Q_crit_look0 Check0 Stop_Prob_H0 
			Cumul_Stop_Prob_H0 Q_fixed Beta_Spent Cumul_Beta_Spent Q_crit_look1 Check1 Stop_Prob_Ha Cumul_Stop_Prob_Ha Cum_Sub_Time_H0 
			Cum_Sub_Time_HA SampleSize;
	run;

%mend;

dm 'log;
clear; 
';


***************************************** The macro to adjust Type II (beta) error ends here*************************************
/*********************************************************************************************************************************/



************************************************************************************************************************************
/************************ Beta-searching algorithm starts here *********************************************************************/;

%let b_count = 0;
%let b_num = 0;
%let b1 = &beta ;
%let b_switch = 0;


%macro Search_Beta ;

	%do v = 1 %to 50;

		%if %sysevalf(&last_e = &last_f) and %sysevalf(&b_count=0) %then
		%do;
			%return;
		%end;

		%if %sysevalf(&last_e > &last_f) %then
		%do;
			%if &v = %eval(&b_count + 1) %then 
			%do;
				%let b2 = %sysevalf(&b1 + 0.005) ;
			%end;
			%else %do;
				%let b2 = %sysfunc(round(%sysevalf((&b1 + &m2)/2),0.0001));
			%end;
			%BetaAdjust(&b2);
			%if %sysevalf((%sysfunc(abs(%sysevalf(&b2) - %sysevalf(&b1)))) < 0.0002) and (&last_e < &last_f) %then 
			%do;
				%return;
			%end;
			%else %do;
				%let m1 = &b1;
				%let b1 = &b2;
				%let b_count = %eval(&b_count + 1);	
			%end;	
		%end;

		%if %sysevalf(&last_e <= &last_f) %then
		%do;
			%if &v = %eval(&b_num + 1) and %sysevalf(%eval(&b_switch) = 0) %then 
			%do;
				%let b2 = %sysevalf(&b1 - 0.005) ;
			%end;
			%else %do;
				%let b2 = %sysfunc(round(%sysevalf((&b1 + &m1)/2),0.0001));
			%end;
			%BetaAdjust(&b2);
			%if %sysevalf((%sysfunc(abs(%sysevalf(&b2) - %sysevalf(&b1)))) < 0.0002) and (&last_e <= &last_f) %then 
			%do;
				%return;
			%end;
			%if %sysevalf((%sysfunc(abs(%sysevalf(&b2) - %sysevalf(&b1)))) >= 0.0002) and (&last_e <= &last_f) %then
			%do;
				%let m2 = &b2;
				%let b1 = &b2;
				%let b_num = %eval(&b_num + 1);
			%end;
			%if (&last_e > &last_f) %then
			%do;
				%let m2 = &b1;
				%let b1 = &b2;
				%let b_switch = 1;
				%let b_num = %eval(&b_num + 1);
			%end;
		%end;

	%end;

%mend;

%Search_Beta;

/************************ Beta-searching algorithm ends here: Added on 09-02-2018 *********************************************/
********************************************************************************************************************************/;


/**************************************************************************************************************************************
*********************************** The macro for graphs starts here *****************************************************************/

%macro PowerGraph;

%do s = 1 %to &Num_Look;

data PlotOutput&s;
set MyOutput&s ;
	df1_look = 2*Evt_n1_H0/((&lambda)**2);
	df2_look = 2*Evt_n2_H0/((&lambda)**2);
	df3_look = 2*Evt_n2_HA/((&lambda)**2);		
	F_prob_cum = 1-Cumul_Stop_Prob_H0;									
	f_crit_cum = quantile('F',F_prob_cum,df1_look,df3_look);
	do g = 1 to 3000 by 1; 
		increment = g/1000;
		Effect_Size = &Delta_PT+increment;	
		output;
		Ha_prob = f_crit_cum/(Effect_Size**(abs(&lambda)/&sigma));
		Cumul_Stop_Prob_Ha = 1-cdf('F',Ha_prob,df1_look,df3_look);
		output;
	end;
run;

%end;


data FinalOutput;
	set FinalOutput;
	if LookNumber=1 then do;
		maxQ=max(Q_crit_look0, Q_crit_look1, Q_fixed)+0.1; 
			call symputx('maxQ', maxQ, 'G');
		minQ=min(Q_crit_look0, Q_crit_look1, Q_fixed)-0.1; 
			call symputx('minQ', minQ, 'G');
	end;
run;


title "Boundary Plot for Test Statistic vs Look Times";
proc sgplot data=FinalOutput;
	series x = LookNumber y = Q_crit_look0 /markers;
	series x = LookNumber y = Q_crit_look1 /markers;
	series x = LookNumber y = Q_fixed /markers;
	XAXIS LABEL = 'Look Times' GRID VALUES = (1 TO &Num_Look BY 1) ;
	YAXIS LABEL = 'Test Statistic Boundary Value' GRID VALUES = (&minQ TO &maxQ BY 0.1) ; 
run;


%mend;
%PowerGraph;

/*********************************** The macro for graphs ends here *****************************************************************/
/**************************************************************************************************************************************/

data runtime; *Stop timer;
	 run_time = datetime() - &_timer_start;
run;

proc print data=runtime;
	var run_time;
	format run_time time11.2;
run;


proc printto; run; *set this to default again to check sas log file for errors and warnings;


data warning; 
   infile "&direct\Mydoc.log" truncover; 
   input error_and_warning $2000.; 
   if index(_infile_, 'ERROR:') = 1 then output; 
   if index(_infile_, 'WARNING:') = 1 then output; 
   run;

proc print data=warning (obs=50); *print just the first 50 observations;
run;


%mend;


data UserDefTime; *when look_points=2;
	input numuser @;
	call symputx("numuser",numuser,'G');
	do i=1 to numuser;
		input usertime @;
		output;
	end;
	*first value of datalines must be num_look (the number of looks);
	datalines;
	3 24 36 60
	;
run;


data UserDefAlpha; *when alpha_look=3;
	do i=1 to &numuser;
		input alpha_look @;
		output;
	end;
	datalines;
	0.0050 0.0125 0.0250
	;
run;


data UserDefBeta; *when beta_look=3;
	do i=1 to &numuser;
		input beta_look @;
		output;
	end;
	datalines;
	0.0100 0.0350 0.1000
	;
run;




*%main(NumSimul =10000 ,alpha =0.025 ,sides =1,lambda =0.5 ,sigma =0.75 ,med =20 ,evt_rate =0.7 ,seed =1729,r =1,Delta_PT_Ha =1.4 ,a =12 ,a_type =1 ,a_omega =1 ,t =60 , bind=1, num_look =3 ,look_points =2 ,alpha_spend =1,rho =1 ,beta =0.10, beta_spend =1 ,rho_f =1 ,num_skip =0, maxiter=200, convg=1E3, direct=C:\Users\jss);

%main(NumSimul =10000 ,alpha =0.025 ,sides =1,lambda =0.5 ,sigma =0.75 ,med =20 ,evt_rate =0.7 ,seed =1729,r =1,Delta_PT_Ha =1.4 ,a =12 ,a_type =1 ,a_omega =1 ,t =60 , bind=1, num_look =3 ,look_points =2 ,alpha_spend =1,rho =3 ,beta =0.10, beta_spend =1 ,rho_f =3 ,num_skip =0, maxiter=200, convg=1E3, direct=C:\Users\jss);

*%main(NumSimul =10000 ,alpha =0.025 ,sides =1,lambda =0.5 ,sigma =0.75 ,med =20 ,evt_rate =0.7 ,seed =1729,r =1,Delta_PT_Ha =1.4 ,a =12 ,a_type =1 ,a_omega =1 ,t =60 , bind=1, num_look =3 ,look_points =2 ,alpha_spend =3,rho =3 ,beta =0.10, beta_spend =3 ,rho_f =3 ,num_skip =0, maxiter=200, convg=1E3, direct=C:\Users\jss);


/**********END OF CODE**********/
