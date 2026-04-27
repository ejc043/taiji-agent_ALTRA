# import necessary libraries
import sys as sys
import os as os
import numpy as np
import pandas as pd

### name key variables
# import passed variables, static variables
input_fp = sys.argv[1]
if not os.path.exists(input_fp):
    print("%s doesn't exist!" %input_fp)
    sys.exit([1])

taiji_config_formulafile_fp = "/stg3/data1/eunice/Taiji/Code/Taiji_uniqueGen/Taiji_uniqueGen/Taiji/taiji_config_formulafile.yml"
taiji_arrayscript_fp = "/stg3/data1/eunice/Taiji/Code/Taiji_uniqueGen/Taiji_uniqueGen/Taiji/Dependencies/Taiji_UniGen_array.sh"

# name variables
master_dir = input_fp[0:input_fp.rfind("/")]
masterinput_dir = master_dir + "/Input"
masteroutput_dir = master_dir + "/Output"
partialoutput_dir = masteroutput_dir + "/Partial"


### create necessary directories
# import input_fp as dataframe
input_xls = pd.ExcelFile(input_fp)
metadata_df = pd.read_excel(input_xls, "active_metadata")
taijiinput_df = pd.read_excel(input_xls, "Active")


# create master and partial directories
try:
    os.mkdir(masterinput_dir)
except OSError:
    print ("Creation of the directory %s failed" % masterinput_dir)
else:
    print ("Successfully created the directory %s " % masterinput_dir)

try:
    os.mkdir(masteroutput_dir)
except OSError:
    print ("Creation of the directory %s failed" % masteroutput_dir)
else:
    print ("Successfully created the directory %s " % masteroutput_dir)

try:
    os.mkdir(partialoutput_dir)
except OSError:
    print ("Creation of the directory %s failed" % partialoutput_dir)
else:
    print ("Successfully created the directory %s " % partialoutput_dir)


# create input and output directories for each sample
#metadata_noNA_df = metadata_df[metadata_df.vcf_F1 != "NA"]
metadata_noNA_df = metadata_df
samples_list = metadata_df["Submitter_ID"].to_list()

for sample in samples_list:
    sampleinput_dir = masterinput_dir + "/" + sample + "_input"
    sampleoutput_dir = partialoutput_dir + "/" + sample + "_output"
    try:
        os.mkdir(sampleinput_dir)
    except OSError:
        print ("Creation of the directory %s failed" % sampleinput_dir)
    else:
        print ("Successfully created the directory %s " % sampleinput_dir)
    try:
        os.mkdir(sampleoutput_dir)
    except OSError:
        print ("Creation of the directory %s failed" % sampleoutput_dir)
    else:
        print ("Successfully created the directory %s " % sampleoutput_dir)


### create unique config and input files for each sample
# create config files
for sample in samples_list:
    sampleconfig_fp = masterinput_dir + "/" + sample + "_input" + "/" + sample + "_config.yml"
    sampleinputfile_fp = masterinput_dir + "/" + sample + "_input" + "/" + sample + "_input.tsv"
    sampleoutput_dir = partialoutput_dir + "/" + sample + "_output"
    samplegenome_fp = metadata_df.loc[metadata_df.Submitter_ID == sample].vcf_Location.values[0]
    with open(taiji_config_formulafile_fp, "r") as f:
        sample1 = f.read().replace("[insert_input_filepath_here]", str(sampleinputfile_fp), 1)
        sample2 = sample1.replace("[insert_output_directory_here]", str(sampleoutput_dir), 1)
        sample3 = sample2.replace("[insert_genome_filepath_here]", str(samplegenome_fp), 1)
    with open(sampleconfig_fp, "w") as f:
        f.write(sample3)


# create input files
for sample in samples_list:
    #sampleinput_df = taijiinput_df[taijiinput_df["group"].str.contains(sample)]
    sampleinput_df = taijiinput_df[taijiinput_df["group"]==(sample)]
    sampleinputfile_fp = masterinput_dir + "/" + sample + "_input" + "/" + sample + "_input.tsv"
    sampleinput_df.to_csv(sampleinputfile_fp, sep = '\t', index = False)


# create file of all config files to be Taiji'ed
taiji_config_list_fp = master_dir + "/" + "taiji_config_files.txt"
filecount = 0
with open(taiji_config_list_fp, "w") as f:
    for root, dirs, files in os.walk(masterinput_dir):
        for file in files:
            if file.endswith(".yml"):
                filepath = os.path.join(root, file)
                print(filepath, file = f)
                filecount += 1


# pass sbatch array command to slurm

if filecount > 15:
    paralelljobs = 15
else:
    paralelljobs = filecount
arraycommand = "1-" + str(filecount) + "%" + str(paralelljobs)
print(arraycommand)
os.system("sbatch -a %s %s %s" % (arraycommand, taiji_arrayscript_fp, taiji_config_list_fp))
