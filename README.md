<img width="468" height="126" alt="image" src="https://github.com/user-attachments/assets/57de6f44-0d3d-434b-aa51-4de23f77e72a" /># The long-term epidemiological impacts of doxycycline post-exposure prophylaxis and vaccination against multidrug-resistant Neisseria gonorrhoeae: a mathematical modelling study

Zihao Wang, Dariya Nikitin, George T.B. Young, Matan Yechezkel, Liang En Wee, Martin T.W. Chio, Lin Geng, Rayner Kay Jin Tan, Yi Wang, David N. Fisman, Joseph A. Lewnard, Lilith K. Whittles, Jue Tao Lim </br>

Requires: RStan (version 2.32.7), R (version 4.5.0), deSolve package (version 1.40) </br>

### File description
* run_mcmc_amr_fixedinitialstate_burntin_UK.R and amr_fixedinitialstate_burntin.stan contain R scripts and Stan code, respectively, used to calibrate the strain-specific gonorrhoea transmission model to data (annual gonorrhoea diagnoses, tests, symptomatic diagnoses, asymptomatic diagnoses, percentage ceftriaxone-resistant, and percentage tetracycline-resistant among MSM) in England
* run_amr_6years+covid.R includes R script to forward-simulate strain-specific gonorrhoea transmission dynamics under baseline conditions (without doxy-PEP and vaccination) and various intervention strategies (i.e., doxy-PEP standalone, vaccination standalone, and dual interventions) for England
* run_amr_6years+covid_failure.R includes R script to forward-simulate strain-specific gonorrhoea transmission dynamics with an adjusted ceftriaxone treatment failure rate
* run_heatmap_uptake_failure.R and run_heatmap_uptakes.R include R scripts for sensitivity analysis on intervention uptake rates and ceftriaxone treatment failure rate

### Objective
* To project the long-term population-level AMR dynamics of N. gonorrhoeae over a 15-year horizon within the epidemiological context of England
* To quantify both the long-term programmatic efficiency and the evolutionary consequences of doxy-PEP and vaccination by jointly evaluating infections averted per enrolment and strain-specific resistance dynamics

### Methodology
* Please refer to Supplementary Information for full implementation details
  
![alt text](https://github.com/killingbear999/amr_gonorrhoea/blob/main/amr_gonorrhea.png)

### Results
* Although doxy-PEP initially reduces gonorrhoea incidence, its long-term programmatic efficiency (0.072 [95% CrI: 0.018 – 0.23] averted infections per enrolment over the 15-year horizon) is substantially lower than vaccination (0.73 [95% CrI: 0.17 – 1.63] overall infections per enrolment) because the selective expansion of tetracycline-resistant strains progressively erodes its epidemiological benefit.
* Vaccination substantially improves long-term outcomes by reducing transmission without introducing additional antimicrobial selection pressure, and combined intervention strategies (0.83 [95% CrI: 0.22 – 1.71] averted infections per enrolment) mitigate resistance expansion more effectively than doxy-PEP alone. 
* The long-term emergence of ceftriaxone-resistant and dual-resistant strains is governed primarily by background ceftriaxone treatment failure (510 [95% CrI: 38 –258600] excess dual-resistant infections across the 15-year horizon (146.18% [95% CrI: 31.30 – 920.86%] cumulative increase) when ceftriaxone treatment failure is 20%) rather than intervention uptake, highlighting preservation of frontline treatment efficacy as the critical determinant of sustainable AMR control.
