## README

---

Data cleaning and analysis scripts for work on estimating the COI for MIP Pf samples
in the DRC and surrounding countries. 

Analysis files are in `analysis` directory. This directory also includes the `data`
used in the analysis and any `data` outputs, as well as the `plots` produced during 
investigations and the final figures placed in `manuscript_figures`. 

---

Overview of file roles:

* `00_meta_cleaning.R` - Initial data cleaning to add lat/longs to samples that were
collected from DRC DHS clusters as well as those also from clusters but not included
in DHS due to age of individuals. 
* `01_data_input.R` - Reading in MIP data, assigning locations and creating allele frequency
tables for the samples.
* `02_mccoil_submission.R` - Using output of 01_data_input.R to submit set of jobs to 
cluster to estimate the COI using McCOILR 
* `03_mccoil_retrieval.R` - Grab COI REAL McCOIL outputs and save to file. 
* `04_mccoil_sensitivity.R` - Investigate role of heterozygote thresholds and impact of
different location allocation to confirm little effect on whether samples are monoclonal
or not. Save outputs of sample COI estimates and which samples are monoclonals. 
* `05_coi_vs_prev.R` - Use COIs to plot and investigate COI differences between countries
and for patterns in COI against malaria prevalence in the DRC. 

---

Additional work not in final figures is contained in the Rmd `99_report.Rmd`
