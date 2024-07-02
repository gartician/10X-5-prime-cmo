# How to demultiplex 10X 5' scRNA-Seq data with CMOs

This SnakeMake pipeline processes a special scenario of single-cell RNA-Seq (scRNA-Seq) data. In this situation, single cells are tagged with Cell Multiplexing Oligomers (CMO) on the cell surface followed by 5' gene expression characterization with a 10X chromium machine. 

While the combination of 3' scRNA-Seq with CMO is officially supported by CellRanger, 5' scRNA-Seq with CMO requires a workaround to correctly demultiplex cells. The source of this process can be found [here](https://www.10xgenomics.com/analysis-guides/demultiplexing-and-analyzing-5%E2%80%99-immune-profiling-libraries-pooled-with-hashtags). 

# Quick Start

The inputs of this pipeline includes the following:

* A filled out `config/config.yaml` file
* A set of GEX and CMO FASTQ files that is yet to be demultiplexed (see below example)

Example input FASTQs in one experiment/sequencing run
```
path/to/RunID_LD_43_GEX_Library_S1_L001_I1_001.fastq.gz
path/to/RunID_LD_43_GEX_Library_S1_L001_I2_001.fastq.gz
path/to/RunID_LD_43_GEX_Library_S1_L001_R1_001.fastq.gz
path/to/RunID_LD_43_GEX_Library_S1_L001_R2_001.fastq.gz
path/to/RunID_LD_44_Feature_Barcode_Library_S2_L001_I1_001.fastq.gz
path/to/RunID_LD_44_Feature_Barcode_Library_S2_L001_I2_001.fastq.gz
path/to/RunID_LD_44_Feature_Barcode_Library_S2_L001_R1_001.fastq.gz
path/to/RunID_LD_44_Feature_Barcode_Library_S2_L001_R2_001.fastq.gz
```

After filling out the inputs and CellRanger configuration files, the pipeline will perform the following: 

1. Align the 5' GEX and CMO FASTQs to a reference genome.
2. Automatically retrieve the number of reads and expected cells per sample.
3. Detect the GEX BAM file per sample, and convert that back to FASTQs.
4. Re-align the GEX FASTQs using `cellranger count` from step 3 using the parameters obtained from step 2. Note that the alignment program `cellranger multi` can be used if there are additional omics data available in your experiment such as feature barcodes, TCR, BCR, etc. Please see Limitations below.

The main output of this pipeline is one set of CellRanger counts results per sample, which can be fed into Seurat for further processing.

# Limitations

* This pipeline does not support using feature barcodes, but modifications/pull requests are encouraged.
* Only one transcriptome reference and one CellRanger executable version is allowed. Again, pull requests are encouraged. 
* This pipeline will not run CellRanger aggr, but the pipeline can be extended to include that step. 

Only one transcriptome reference and one cellranger executable is allowed.

# Running the Pipeline

Please run this pipeline interactively (preferably on a head node) with a SLURM profile.