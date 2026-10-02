# The long-term epidemiological impacts of doxycycline post-exposure prophylaxis and vaccination against multidrug-resistant Neisseria gonorrhoeae: a mathematical modelling study

Zihao Wang, Dariya Nikitin, George T.B. Young, Matan Yechezkel, Liang En Wee, Martin T.W. Chio, Lin Geng, Rayner Kay Jin Tan, Yi Wang, David N. Fisman, Joseph A. Lewnard, Lilith K. Whittles, Jue Tao Lim </br>

Requires: RStan (version 2.32.7), R (version 4.5.0), deSolve package (version 1.40) </br>

### File description
* run_mcmc_amr_fixedinitialstate_burntin_UK.R and amr_fixedinitialstate_burntin.stan contain R scripts and Stan code, respectively, used to calibrate the strain-specific gonorrhoea transmission model to data (annual gonorrhoea diagnoses, tests, symptomatic diagnoses, asymptomatic diagnoses, percentage ceftriaxone-resistant, and percentage tetracycline-resistant among MSM) in England
* run_amr_6years+covid.R includes R script to forward-simulate strain-specific gonorrhoea transmission dynamics under baseline conditions (without doxy-PEP and vaccination) and various intervention strategies (i.e., doxy-PEP standalone, vaccination standalone, and dual interventions) for England
* run_amr_6years+covid_failure.R includes R script to forward-simulate strain-specific gonorrhoea transmission dynamics with an adjusted ceftriaxone treatment failure rate
* run_heatmap_uptake_failure.R and run_heatmap_uptakes.R include R scripts for sensitivity analyses on intervention uptake rates and ceftriaxone treatment failure rate, respectively
* run_sensitivity_doxypep_efficacy.R and run_sensitivity_efficacy.R include R scripts for sensitivity analyses on doxy-PEP efficacy and vaccine efficacy, respectively

### One sentence summary
Combining vaccination with doxy-PEP can improve gonorrhoea control while mitigating, but not necessarily eliminating, doxy-PEP-associated selection for antimicrobial resistance.

### Abstract
Doxycycline post-exposure prophylaxis (doxy-PEP) reduces bacterial sexually transmitted infections but may select for antimicrobial-resistant Neisseria gonorrhoeae. We developed a strain-specific transmission model calibrated to surveillance data from men who have sex with men in England to evaluate population-level implementation of doxy-PEP and/or 4CMenB vaccination over 2027–2041. The model represented susceptible, tetracycline-resistant (Tet-R), ceftriaxone-resistant, and dual-resistant strains. Doxy-PEP initially reduced gonorrhoea incidence but selected for Tet-R strains, yielding only a 3.38% (95% credible interval [CrI]: 0.54 – 17.93%) cumulative reduction in infections over 15 years. At the primary assumption of 40% vaccine efficacy, vaccination alone and the combined intervention produced larger and more sustained reductions in transmission; the combined intervention achieved a 44.10% (95% CrI: 4.18 – 86.62%) cumulative reduction, with annual incidence in 2041 of 94,300 (95% CrI: 0 – 518,300) versus 203,600 (95% CrI: 44,800 – 573,400) without intervention. However, cumulative Tet-R infections under the combined intervention remained higher than under no intervention. Vaccine-efficacy sensitivity analysis showed that 30% was the lowest tested efficacy eliminating the median rebound in total incidence, whereas 50% was required for the median cumulative Tet-R burden to fall below the no-intervention counterfactual, indicating that lower protection may prevent incidence rebound while greater efficacy may be needed to offset doxy-PEP-associated Tet-R expansion. When ceftriaxone treatment failure was increased to 20%, doxy-PEP alone generated 510 excess dual-resistant infections (95% CrI: 38 – 258,600) over 15 years. These findings indicate that vaccination can enhance gonorrhoea control during doxy-PEP implementation, but long-term benefit depends on vaccine efficacy and preservation of ceftriaxone effectiveness. 

### Methodology
* Please refer to Supplementary Information for full implementation details
  
![alt text](https://github.com/killingbear999/amr_gonorrhoea/blob/main/amr_gonorrhea.png)
