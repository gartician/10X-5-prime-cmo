import re
import os
import glob
import pandas as pd
from collections import defaultdict

configfile: "config/config.yaml"

# establish experiment, config file, and cmo relations with a DF --------------
exp_ids = list(config["multi"].keys())
multi_df = pd.DataFrame({
    "Exp": exp_ids,
    "Config": [ config["multi"][i]["config"] for i in exp_ids ],
    "CMO": [ config["multi"][i]["cmo"] for i in exp_ids ]
})
multi_df.index = exp_ids

print("\nCellRanger Multi Configuration")
print(multi_df)

experiment_sample = pd.read_csv(config["experiment_sample"])

print("\nExpected Samples from Each Experiment")
print(experiment_sample)

exp_sample = defaultdict(list)
for i in range(0, experiment_sample.shape[0]):
    exp = experiment_sample.iloc[i, ]["experiment"]
    smp = experiment_sample.iloc[i, ]["sample"]
    exp_sample[exp].append(smp)

# define constants
project_id = config["project_id"]
experiments_list = experiment_sample["experiment"]
samples_list = experiment_sample["sample"]

wildcard_constraints:
    experiment="[a-zA-Z0-9-_]+",
    sample="[a-zA-Z0-9-_]+"

def get_reads(wildcards):

    """
    Given a metrics_summary.csv output from CellRanger, determine number of GEX reads.
    """

    # import file
    df = pd.read_csv(wildcards)

    # retrieve number of GEX reads
    nreads = df[ (df["Library Type"] == "Gene Expression") & (df["Metric Name"] == "Number of reads")]

    if len(pd.unique(nreads["Metric Value"])) != 1:
        os.stop("Discordant number of GEX reads in a sequencing run")

    nread = nreads["Metric Value"].str.replace(",", "").astype(int).unique()[0]
    nread = str(round(nread, -6)) # round to the nearest million, then turn into a string integer

    return(nread)

def get_cells(wildcards):

    """
    Given a metrics_summary.csv output from CellRanger, determine number of expected cells to feed into CellRanger count
    """

    infile = f"data/multi/{wildcards.experiment}/outs/per_sample_outs/{wildcards.sample}/metrics_summary.csv"
    df = pd.read_csv(infile)
    num_cells = int(df[df["Metric Name"] == "Cells"]["Metric Value"].str.replace(",", "").astype(int)[0])
    return(num_cells)

def get_fastqs(wildcards):

    """
    Determine which set of BAM files to use after running bamtofastq in order to feed to
    a second round of CellRanger count
    """

    headers = f"data/headers/{wildcards.experiment}/{wildcards.sample}.txt"

    # determine which bamtofastq folder holds your GEX FASTQs via the gex_code
    with open(headers) as fh:
        for line in fh:
            if "Gene Expression" in line:
                gex_code = re.search('"library_id":.', line).group()[-1]
                break

    # determine the folders holding your fastqs, then extract a seq code
    folders_with_fastqs = glob.glob(f"data/bamtofastq/{wildcards.experiment}/{wildcards.sample}/*", recursive=True)
    
    # use this adjustment when bamtofastq has not finished running yet - good for printing dry runs and such
    if len(folders_with_fastqs) == 0:
        f = "TBD"
    elif len(folders_with_fastqs) >= 1:

        # identify a flow cell code
        seq_code = list(set([ i.split("_")[-1] for i in folders_with_fastqs  ]))[0]

        # output the folder containing GEX reads
        f = f"data/bamtofastq/{wildcards.experiment}/{wildcards.sample}/{wildcards.experiment}_{gex_code}_1_{seq_code}"

    return(f)

rule all:
    input:
        # initial demultiplex with CR multi
        expand([
            "data/multi/{experiment}/outs/per_sample_outs/{sample}/metrics_summary.csv",
            "data/multi/{experiment}/outs/per_sample_outs/{sample}/count/sample_alignments.bam"], zip,
            experiment = experiment_sample["experiment"],
            sample = experiment_sample["sample"]
        ),

        # bamtofastq
        expand("data/bamtofastq/{experiment}/{sample}", zip,
            experiment = experiments_list,
            sample = samples_list
        ),

        # BAM headers in preparation for CR count
        expand("data/headers/{experiment}/{sample}.txt", zip,
            experiment = experiments_list,
            sample = samples_list
        ),

        # CR count GEX reads
        expand("data/count/{experiment}/{sample}/outs/molecule_info.h5", zip,
            experiment = experiments_list,
            sample = samples_list
        )

rule cellranger_multi:
    output:
        expand("data/multi/{{experiment}}/outs/per_sample_outs/{sample}/metrics_summary.csv", sample = samples_list),
        expand("data/multi/{{experiment}}/outs/per_sample_outs/{sample}/count/sample_alignments.bam", sample = samples_list)
    params:
        cellranger = config["cellranger"],
        reference = config["reference"],
        config = lambda wildcards: config["multi"][wildcards.experiment]["config"]
    shell:
        "{params.cellranger} multi " 
        "--id {wildcards.experiment} "
        "--csv {params.config} "
        "--localmem 64 "
        "--localcores 8; "
        "rm -rf data/multi/{wildcards.experiment}; "
        "mv {wildcards.experiment} data/multi"

rule bamtofastq:
    input:
        bam = "data/multi/{experiment}/outs/per_sample_outs/{sample}/count/sample_alignments.bam",
        metrics_summary = "data/multi/{experiment}/outs/per_sample_outs/{sample}/metrics_summary.csv"
    output:
        directory("data/bamtofastq/{experiment}/{sample}")
    params:
        bamtofastq = config["bamtofastq"],
        num_reads = lambda wildcards, input: get_reads(input.metrics_summary)
    shell:
        "{params.bamtofastq} --reads-per-fastq {params.num_reads} {input.bam} {output}"

rule samtools_header:
    input:
        "data/multi/{experiment}/outs/per_sample_outs/{sample}/count/sample_alignments.bam"
    output:
        "data/headers/{experiment}/{sample}.txt"
    params:
        samtools = config["samtools"]
    shell:
        "{params.samtools} view -H {input} > {output}"

rule cellranger_count:
    input:
        bam = rules.bamtofastq.output,
        headers = "data/headers/{experiment}/{sample}.txt",
        metrics_summary = "data/multi/{experiment}/outs/per_sample_outs/{sample}/metrics_summary.csv"
    output:
        "data/count/{experiment}/{sample}/outs/molecule_info.h5"
    params:
        cellranger = config["cellranger"],
        reference = config["reference"],
        fastqs = get_fastqs,
        chemistry = config['chemistry'],
        num_cells = get_cells
    threads: 8
    shell:
        "{params.cellranger} count "
        "--id {wildcards.sample} " 
        "--transcriptome {params.reference} "
        "--fastqs {params.fastqs} "
        "--create-bam true "
        "--force-cells {params.num_cells} "
        "--check-library-compatibility false "
        "--disable-ui "
        "--nosecondary "
        "--chemistry {params.chemistry} "
        "--localmem 64 "
        "--localcores 8; "
        "rm -rf data/count/{wildcards.experiment}/{wildcards.sample}; "
        "mv {wildcards.sample} data/count/{wildcards.experiment} "