import multiprocessing
import os
import subprocess
import pandas as pd
import statistics

# samples and dirs
data = pd.read_excel("Data/SRA_IDs.xlsx")
data["path"] = "Results/" + data["organism"] + "/" + data["tpntorg"]

check_dirs = data.path.unique()

for dir in check_dirs:
    checker=os.path.isdir(dir)
    if not checker:
        os.makedirs(dir)
        print("Created dir: ", dir)

#Main function
def my_pipeline(x):
    organism = data.iloc[x, 0]
    sra_id = data.iloc[x, 1]
    layout = data.iloc[x, 2]
    path = "./Results/"+data.iloc[x, 0]+"/"+data.iloc[x, 3]+"/"+sra_id

# Calling SRAtoolkit
    print ("Generating fastq for: " + sra_id)
    fasterq_dump = "~/sratoolkit.3.0.1-ubuntu64/bin/fasterq-dump -O ./Temp " + sra_id
    print ("The command used was: " + fasterq_dump)
    subprocess.call(fasterq_dump, shell=True)

#Computing mean and sd for single-end data
    if layout == "SINGLE":
        n_row = 0
        lens = []
        with open("./Temp/" + sra_id + ".fastq") as fastq:
            for line in fastq:
                n_row += 1
                if n_row % 2 == 0 and n_row % 4 != 0:
                    lens.append(len(line.strip("N\n")))
            length = str(statistics.mean(lens))
            sd = str(statistics.stdev(lens))

#Calling kallisto
    if layout == "SINGLE":
        reads = "./Temp/" + sra_id + ".fastq"
        kallisto_params="-i ./Indexes/" + organism + "_index.idx --single -l " + length + " -s " + sd
    elif layout == "PAIRED":
        reads = "./Temp/" + sra_id + "_1.fastq ./Temp/" + sra_id + "_2.fastq"
        kallisto_params="-i ./Indexes/" + organism + "_index.idx"

    kallisto = "kallisto quant " + kallisto_params + " -o " + path + " " + reads
    print ("The command used was: " + kallisto)
    subprocess.call(kallisto, shell=True)
    subprocess.call("rm " + reads, shell=True)

n_cores = multiprocessing.cpu_count()

# Run in parallel
with multiprocessing.Pool(n_cores-1) as pool:
    pool.map(my_pipeline, range(0, len(data)))
    