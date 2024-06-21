File and folder structure of reproducibility materials for: Bayesian Methods for Extended Batch Marking

Reproducibility Files
├── 01 - Case Study 1
|   ├── Cowen2017_data.csv             <--- Data used in Case Study 1
|   ├── MLE_compute_time.R             <--- Reproduce Table 2 and Table 3 (c) and (d)
|   ├── Marked Likelihood
|   |   ├── Lmarked_Lj_phi_p.R         <--- Reproduce Table 2 and Table 3 (a) row 1
|   |   ├── Lmarked_Lj_phi_p_t.R       <--- Reproduce Table 2 and Table 3 (a) row 2
|   |   ├── Lmarked_Lj_phi_t_p.R       <--- Reproduce Table 2 and Table 3 (a) row 3
|   |   └── Lmarked_Lj_phi_t_p_t.R     <--- Reproduce Table 2 and Table 3 (a) row 4
|   └── Combined Likelihood
|       ├── Lcombined_Lj_phi_p.R       <--- Reproduce Table 2 and Table 3 (b) row 1
|       ├── Lcombined_Lj_phi_p_t.R     <--- Reproduce Table 2 and Table 3 (b) row 2
|       ├── Lcombined_Lj_phi_t_p.R     <--- Reproduce Table 2 and Table 3 (b) row 3
|       └── Lcombined_Lj_phi_t_p_t.R   <--- Reproduce Table 2 and Table 3 (b) row 4, and Figure 1 and Figure 2
└── 02 - Case Study 2
    ├── Data
    |   └── SandLanceData.xls          <--- Data used in Case Study 2, summarized in Table 4 and Table 5
    └── fit model.R                    <--- Reproduce results of Case Study 2: Table 6, Table 7

Data Descriptions
├── Cowen2017_data.csv
|
| Weatherloach data table.
|
| See Huggins et al. (2010) for additional details on the weatherloach data table.
| First row contains the captured unmarked at each sampling occasion. First column starting row 2 contains the 
| total released at each sampling occasion. Remaining entries are the matrix of recapture counts for each 
| group (row number-1) and each sampling occasion (column number).
└──
├── SandLanceData.xls
|   ├── Sheet 1 (SandLanceData)
|   |   └── Columns: Date, DeltaDays, Colour, Sample Occasion (g), Number Sampled (Sg), Number Marked (Rg), Number Unmarked (ut), Temperature, Salinity
|   ├── Sheet 2 (Batch Mark Data Array)
|   |   └── Recapture Matrix: Groups (g) \ Recapture Occassions (t)
|   └── Sheet 3 (Batch Mark Covariates)
|       └── Columns: Batch (g), Sampling Occasion (t), Length, Mass
|
| Pacific sand lance data tables. 
└──
