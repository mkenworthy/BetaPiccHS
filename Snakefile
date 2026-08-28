rule combine_astep:
    input:
        "src/data/astep/betapic_astep_2017.csv",
        "src/data/astep/betapic_astep_2018.csv"
    output:
        "src/data/astep_all.fits"
    conda:
        "environment.yml"
    script:
        "src/scripts/02_Astep_write_out_data.py"

rule combine_brite:
    input:
        "src/data/brite/betaPic_2015-16-17-18-BHr.dat",
        "src/data/brite/betaPic_2019_BHr.dat",
        "src/data/brite/betaPic_2021-BHr-all.dat" 
    output:
        "src/data/brite/brite_all_R.fits"
    conda:
        "environment.yml"        
    script:
        "src/scripts/03_Combine_brite_data.py"
        
rule bin_BRITE:
    input:
        "src/data/brite/brite_all_R.fits"  
    output:
        "src/data/binned_BRITE.dat"
    conda:
        "environment.yml"
    script:
        "src/scripts/04_Bin_data_BRITE.py"

rule bin_BRING:
    input:
        "src/data/bring/Reduced_betaPic.fits"  
    output:
        "src/data/binned_BRING.dat"
    conda:
        "environment.yml"
    script:
        "src/scripts/04_Bin_data_BRING.py"

rule bin_ASTEP:
    input:
        "src/data/astep_all.fits"  
    output:
        "src/data/binned_ASTEP.dat"
    conda:
        "environment.yml"
    script:
        "src/scripts/04_Bin_data_ASTEP.py"