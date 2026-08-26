rule my_script1:
    input:
        "src/scripts/02_Astep_write_out_data.py"
    output:
        "src/data/astep/astep_all.fits"
    shell:
        "python {input}"

rule my_script2:
    input:
        "src/scripts/03_Combine_brite_data.py"  
    output:
        "src/data/brite/brite_all_R.fits"
    shell:
        "python {input}"
        
rule my_script3:
    input:
        "src/scripts/04_Bin_data.py"  
    output:
        "src/data/binned_BRITE.dat"
        "src/data/binned_BRING.dat"
        "src/data/binned_ASTEP.dat"
    shell:
        "python {input}"
