# MetaBinQ #

## Program
* metabinQ.py - main script
* tnCounter.c - c module for tetranuclotide counting
* setup.py - script for seting up tnCounter module

## Directories:

### scripts - other scripts used during the project 
    - analysis.ipynb - script for analysis of the prediction results
    - details.ipynb - script for quick visulization of bins
    - model_contamination.py - script producing contaminatd bins
    - genomes.txt - list of genomes used to generate viral data set
    
### results - results for some data sets

    - Metabat - Bins from Metabat repo binned with Metabat
        - metabat_high.csv - real values of contamination for high complexity data set
        - metabat_med.csv - real values of contamination for medium complexity data set
        - metabat_low.csv - real values of contamination for low complexity data set
        - metabat_high_results.csv - predicted contamination for high complexity data set
        - metabat_med_results.csv - predicted contamination for medium complexity data set
        
    - Maxbin - Bins from Metabat repo binned with Maxbin
        - maxbin_high.csv - real values of contamination for high complexity data set
        - maxbin_med.csv - real values of contamination for medium complexity data set
        - maxbin_low.csv - real values of contamination for low complexity data set
        - maxbin_high_results.csv - predicted contamination for high complexity data set
        - maxbin_med_results.csv - predicted contamination for medium complexity data set

    - mouse - bins from cami mouse gut data set, binned using Metabat
        - cami_mouse.csv - real values of contamination
        - new_mouse_gut_res.csv - predicted contamination

    - viral - simulated bins using bacterial genomes contaminated with viruses
        - vscont.csv - real values of contamination
        - new_bac_vir.csv - predicted contamination
        
