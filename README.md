# FWER_Knockoff_filter

This repository contains all codes and data for reproducing simulation experiments and data analysis of our Biometrics paper "Summary Statistics Knockoffs Inference with Family-wise Error Rate Control" (http://).

The folder "Functions" contains three scripts that contain functions that will be used in all simulations, including "ESGWAS.R", "GhostKnockoff_hello.R" and "KnockoffScreen_AL_vMar16.R".

The folder "Section 3.2" contains scripts that perform comparisons of computational efficiency in Section 3.2. Specifically, "Computational time (Our method).R" is implemented to count the computational time of our method and so are "Computational time (Ren et al. 2023).R" and "Computational time (Trivial approach).R" analogously. "Summary.R" is implemented to summarize all computation time returned.

The folder "Section 4.1" contains scripts that perform simulation in Section 4.1. Specifically, "AR(1) 4.1.R" performs experiments under the AR(1) correlation structure and "Compound Symmetry 4.1.R" performs experiments under the compound symmetry correlation structure

The folder "Section 4.2" contains scripts that perform simulation in Section 4.2. Specifically, "AR(1) Sample Size 4.2.R" and "AR(1) Number of nonnull features 4.2" perform experiments under the AR(1) correlation structure and "Compound Symmetry Sample Size 4.2.R" and "Compound Symmetry Number of nonnull features 4.2" perform experiments under the compound symmetry correlation structure


The folder "Section 4.3" contains scripts that perform simulation in Section 4.3. Specifically, "AR(1) 4.3.R" performs experiments under the AR(1) correlation structure and "Compound Symmetry 4.3.R" performs experiments under the compound symmetry correlation structure


The folder "Section 5" contains scripts that perform analysis in Section 5.
Specifically, 
1. File "AD_Zscores_candidate.csv" contains all information of the 7,963 variants we analyze.
2. Files "Cor_1.csv"~"Cor_22.csv" contains all information of the correlation among the 7963 variants, separated by chromosomes.
3. "Ghost.R" and "Deran.R" perform the analysis using our method and the procedure of Ren et al. (2023).
4. "Plot_Ghost.R" and "Plot_Deran.R" visualize the analysis results of our method and the procedure of Ren et al. (2023).
5. "Reference_panel.csv", "Candidate_variants_info_other.csv" contain the information required for drawing plots.
