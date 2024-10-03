# County-Level Associations between Drinking Water PFAS Contamination and COVID-19 Mortality in the United States

This repository holds replication code for the following study:  

Liddie, J.M., Bind, M.A., Karra, M., Sunderland, E.M. County-Level Associations between Drinking Water PFAS Contamination and COVID-19 Mortality in the United States. Journal of Exposure Science & Environmental Epidemiology (2024) DOI: 10.1038/s41370-024-00723-5 doi: 10.1038/s41370-024-00723-5

Files are numbered according to the order in which they are run for the main analyses, sensitivity analyses (including propensity score matching), and table/figure generation. Replication datasets and additional descriptions of the data are available on the Harvard Dataverse here[https://doi.org/10.7910/DVN/PN0RI5].

# Contents
## Modelling  
- **1_Modelling.R**: primary models (including counties in the statewide sampling dataset and the Third Unregulated Contaminant Monitoring Rule) for the association between PFAS drinking water contamination and cumulative COVID-19 mortality 
- **2a_Modelling_matching.R**: propensity score matched analyses (unadjusted and adjusted) as robustness/sensitivity analyses
- **2b_Sensitivity_analyses.R**: all other sensitivity analyses related to modeling choices, outcome measures, exposure measures, and confounding adjustments

## Tables and figures  
- **3_Final tables and figures.RMD**: all primary and supplemental figures and tables (excluding love plots)

# Authors  

- [Jahred Liddie](https://scholar.harvard.edu/jmliddie), Department of Environmental Health, Harvard T.H. Chan School of Public Health
- Marie-Abèle Bind, Biostatistics Center, Massachusetts General Hospital
- Mahesh Karra, Frederick S. Pardee School of Global Studies, Boston University
- [Elsie M. Sunderland](https://bgc.seas.harvard.edu/), Department of Environmental Health, Harvard T.H. Chan School of Public Health; Harvard John A. Paulson School of Engineering and Applied Sciences, Harvard University

# Contact

Jahred Liddie
Harvard T.H. Chan School of Public Health
Department of Environmental Health
jliddie@g.harvard.edu
