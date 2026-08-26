# ASTEP 400 data

The extracted photometry is in:

    betapic_astep_2017.csv
    betapic_astep_2018.csv

The field of Beta Pictoris along with the reference stars are labelled in:

    betapic_field.pdf

The files contain the following columns:

    BJD: Barycentric Julian Date
    FCAL0: Beta Pic flux divided by its daily median 
    FCAL1: Beta Pic flux with optimal calibration
    FCAL2: Beta Pic flux with a simple calibration using reference star #3 to calculate the daily median and reference star #2 for variations on shorter timescales
    FCAL3: Beta Pic flux with a simple calibration using reference star #5 to calculate the daily median and reference star #2 for variations on shorter timescales
    SUNELEV: Sun elevation (in degrees)
    AIRMASS: Beta Pic airmass
    SKY: Level of the sky (in ADU)
    FLAG: 0 if the stars were not correctly detected or the guiding was poor. 1 otherwise. 

The column FCAL1 is used for this paper. It is calculated as follows:

For BJD < 57970, we use the reference star #3 to calculate the daily median, and for BJD > 57970 the reference star #5 is used to calculate the daily median, due to telescope adjustments made around that epoch.

Between BJD > 57891 & BJD < 57907, we use reference star #3 to calculate the daily median with a scaling in the flux of 0.985, to account for ice stuck on the window of the camera box, partially obscuring the light and affecting stars differentially. 

