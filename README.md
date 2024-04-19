# tprK pipeline
This pipeline was designed to take Illumina and PacBio files straight off the sequencer to a final comparison table of all the different variable regions with their relative frequencies, as well as various pretty plots along the way.

>[!WARNING]
Fastqs much be gzipped for both Illumina and PacBio runs
>
>Fastqs for Illumina should be paired and trimmed before being run through the pipeline

## Table of Contents
* [Installation](#Setup)
* [Formating Metadata Table](#Metadata)
* [WorkFlow](#WorkFlow)
    * [Example Commands](#Commands)
* [Example](#Example)
    * [Illumina](#Illumina)
    * [Pacbio](#Pacbio)
* [Output](#Output)


## Installation

1. Install [nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation).
   - Make sure you move nextflow to a directory in your PATH variable.
2. Install [docker](https://docs.docker.com/get-docker/). The first time running this program will take a while, as the docker image will take some time to build, but this is a one time thing!
3. Setup [nextflow tower](https://seqera.io/)

> [!WARNING]
This pipeline was written with DSL1 and must be run with an older version of nextflow before your nexflow run command use:
>
>NXF_VER=22.10.4

## Metadata

A metadatafile is passed in to give sample names and locations.

> [!NOTE]
> For Sample name column remove ".fastq.gz" from the end of the sample name and "Ill_" or "PB_" from the start of the filename.

>For Illumina Runs add "Ill_" before the sample name for the fastq for the file and in the metadata

|SampleName | Illumina | PacBio | BioSample | Replicate |
| ---       |  ---     | ---    | ---       | ---       | 
|Example_1	|Ill_Example_1.fastq.gz	|	
|Example_2	|Ill_Example_2.fastq.gz	|
|Example_3	|Ill_Example_3.fastq.gz	|		
|Example_4	|Ill_Example_4.fastq.gz	|
|Example_5	|Ill_Example_5.fastq.gz	|	

>For Pacbio Runs add "PB_" before the sample name for the fastq for the file and in the metadata

>[!WARNING]
>For Pacbio runs a newline must be added to the end of the metadata file or you will receive an incomplete final line error for the createPacBioTree step.

|SampleName | Illumina | PacBio | BioSample | Replicate |
| ---       |  ---     | ---    | ---       | ---       | 
|Example_1	||PB_Example_1.fastq.gz	|	
|Example_2	||PB_Example_2.fastq.gz	|
|Example_3	||PB_Example_3.fastq.gz	|		
|Example_4	||PB_Example_4.fastq.gz	|
|Example_5	||PB_Example_5.fastq.gz	|
||||	

>[!IMPORTANT]
>If running the fastqs in a specific directory the directory must be specified in the --INPUT option and in the metadata file. It is recommended to run the pipeline in the directory with your fastqs and pass "./" for --INPUT.

|SampleName | Illumina | PacBio | BioSample | Replicate |
| ---       |  ---     | ---    | ---       | ---       | 
|Example_1	|Example_Directory/Ill_Example_1.fastq.gz	|	
|Example_2	|Example_Directory/Ill_Example_2.fastq.gz	|
|Example_3	|Example_Directory/Ill_Example_3.fastq.gz	|		
|Example_4	|Example_Directory/Ill_Example_4.fastq.gz	|
|Example_5	|Example_Directory/Ill_Example_5.fastq.gz	|

>For running with technical replicates name any replicates with the same biosample name in the biosample column but differentiate with A B for what replicate the sample is

|SampleName | Illumina | PacBio | BioSample | Replicate |
| ---       |  ---     | ---    | ---       | ---       | 
|R501B1_R1_S5_merged|	Ill_R501B1_R1_S5_merged.fastq.gz||		KO_Time1_R501|	A
|R501B1_R2_S13_merged|	Ill_R501B1_R2_S13_merged.fastq.gz||		KO_Time1_R501|	B
|R501B2_R1_S21_merged|	Ill_R501B2_R1_S21_merged.fastq.gz||		KO_Time2_R501	|A
|R501B2_R2_S29_merged|	Ill_R501B2_R2_S29_merged.fastq.gz||		KO_Time2_R501|	B

## Options

| Command  | Description |
| ---      | ---         | 
| `--INPUT`  | Input folder where fastqs are located. For current  directory, `./` can be used.
| `--OUTDIR` | Output folder where .bams and consensus fastas will be piped into.
| `--METADATA` | Path to metadata file with specific format. 
| `--PACBIO` | Specify that there are only PacBio files to be read.
| `--ILLUMINA` | Specify that there are only Illumina files to be read.
| `--LARGE` | Splitting of visualizations will be done, heatmap and variable site comparison will produce a separate file for each variable region rather than combining them. Recommended for all runs as output is cleaner.
| `--REFERENCE` | Specify Illumina sample name (not file), to compare others to for dot-line plots. Can be used in tandem with --LARGE.
| `--TECHNICAL_REPLICATE` | Generates additional statistics if using technical replicates. 
| `--RF_FILTER` | Specify relative frequency filter. Default is 0.00001.
| `--COUNT_FILTER` | Specify count filter. Default is 0.
| `--ILLUMINA_FILTER` | Specify whether PacBio reads should be filtered to only include files supported by Illumina reads that reach the cutoff.
| `-resume`  | nextflow will pick up where it left off if the previous command was interrupted for some reason.
| `-with-docker ubuntu:18.04` | Runs command with Ubuntu docker.
| `-with-trace` | Outputs a trace.txt that shows which processes end up in which work/ folders. 
|`-profile`|`standard`: For less computationally intensive systems run locally, not reccommended
||`standardMORE`: For running on the more computationally strong local systems, used for lab i9 imacs 
||`laptop`: For running on very low power computational systems, very slow
||`Cloud`: For running on the cloud adds more computational power for memory intensive steps, recommended
|`-c`|Add you nextflow config file to access cloud
|`-with-tower`|Monitor your run with nextflow tower 

This pipeline can be run with both Illumina and PacBio samples at the same time skip --ILLUMINA and --PACBIO option and fill out the metadata table with both PacBio and Illumina samples, each sample must have a PacBio and an Illumina sample

### Commands

Illumina run in cloud:

```
NXF_VER=22.10.4 nextflow run main.nf  \
    --METADATA Metadata_Example_Illumina.csv \
    --ILLUMINA \
    --INPUT ./ \
    --OUTDIR Example_Illumina/ \
    -with-docker ubuntu:18.04 \
    -profile Cloud \
    -c your_nextflow_aws.config \
    -with-tower
```

Pacbio run in cloud:

```
NXF_VER=22.10.4  nextflow run main.nf \
    --METADATA Metadata_Example_PACBIO.csv \
    --PACBIO \
    --INPUT Example_PacBio_fastq/ \
    --OUTDIR Example_PACBIO_Out/ \
    -with-docker ubuntu:18.04 \
    -profile Cloud \
    -c your_nextflow_aws.config \
    -with-tower 
```

Illumina run local standard cpu usage:

```
NXF_VER=22.10.4 nextflow run main.nf  \
    --METADATA Metadata_Example_Illumina.csv \
    --ILLUMINA \
    --INPUT ./ \
    --OUTDIR Example_Illumina/ \
    -with-docker ubuntu:18.04 \
    -profile standard
```

Pacbio run local higher cpu usage:

```
NXF_VER=22.10.4  nextflow run main.nf \
    --METADATA Metadata_Example_PACBIO.csv \
    --PACBIO \
    --INPUT Example_PacBio_fastq/ \
    --OUTDIR Example_PACBIO_Out/ \
    -with-docker ubuntu:18.04 \
    -profile standardMORE
```

## WorkFlow

Default workflow for Illumina and Pacbio runs

```mermaid
flowchart TD;
    subgraph  
    visualizeAllData;
    end
    subgraph Illumina;
    createAllAssignments-->createFrequencyTables_Illumina;
    createFrequencyTables_Illumina-->summaryStats_Illumina;
    summaryStats_Illumina-->filterReads;
    classDef bar stroke:#0f0
    filterReads-->createFrequencyPlots_Illumina;
    createFrequencyPlots_Illumina-->createVariableRegionComparisons;
    createVariableRegionComparisons-->subsetReads;
    subsetReads-->visualizeAllData:::bar;
    end;
    subgraph PacBio;
    denoisePacBioFiles-.->createFrequencyTables_PacBio;
    createFrequencyTables_PacBio-.->summaryStats_PacBio;
    summaryStats_PacBio-.->createFrequencyPlots_PacBio;
    createFrequencyPlots_PacBio-.->createPacBioTree;
    createPacBioTree-.->visualizeAllData;
    end
```

## Example

Install [sratoolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit)

On Mac you can use brew:
```bash
brew install sratoolkit
```

### Illumina

Make an _Example_Illumina_fastq_ folder

```bash
mkdir Example_Illumina_fastq
```

Enter Example folder 
```bash
cd Example_Illumina_fastq 
```

Download samples from SRA and place them in the _Example_Illumina_fastq_ folder.

```bash
fasterq-dump --split-files SRR24317964
fasterq-dump --split-files SRR24317966
fasterq-dump --split-files SRR24317967	
fasterq-dump --split-files SRR24317968	
fasterq-dump --split-files SRR24317969
```

Cat R1 and R2 samples together

```bash
for i in *_1.fastq
do
base=$(basename $i _1.fastq)
cat ${base}_1.fastq ${base}_2.fastq > ${base}.fastq
done
```

Remove unpaired reads
```bash
rm *_1.fastq
rm *_2.fastq
```

gzip reads
```bash
gzip *.fastq
```

Add "Ill_: to the front of all read names
```bash
for file in *; do mv "$file" "Ill_$file"; done
```

Exit folder 
```bash
cd ..
```

Run example workflow or use one from above
```
NXF_VER=22.10.4 nextflow run greninger-lab/TPRK-Pipeline -r main --METADATA Metadata_Example_Illumina.csv --ILLUMINA --INPUT Example_Illumina_fastq/ --OUTDIR Example_Illumina/ -with-docker ubuntu:18.04 -profile standard
```

### Pacbio

Make an _Example_PacBio_fastq_ folder

```bash
mkdir Example_PacBio_fastq
```

Enter Example folder 
```bash
cd Example_PacBio_fastq
```

Download samples from SRA and place them in the _Example_PacBio_fastq_ folder.

```bash
fasterq-dump SRR10294254		
fasterq-dump SRR10294238
```

gzip reads
```bash
gzip *.fastq
```

Add "PB_: to the front of all read names
```bash
for file in *; do mv "$file" "PB_$file"; done
```

Exit folder 
```bash
cd ..
```

Run example workflow or use one from above
```
NXF_VER=22.10.4  nextflow run greninger-lab/TPRK-Pipeline -r main --METADATA Metadata_Example_PACBIO.csv --PACBIO --INPUT Example_PacBio_fastq/ --OUTDIR Example_PACBIO_Out/ -with-docker ubuntu:18.04 -profile standard
```

## Output

Illumina Output
```
Example_Illumina_Out
├── Figures
│   ├── Relative_Frequency_Plots
│   │   ├── Ill_SRR24317964_RelativeFreqPlot.html                               # Bar chart comparing variability between variable sites for each sample 
│   │   ├── Ill_SRR24317964_filtered_RelativeFreqPlot.html                      # Bar chart comparing variability between variable sites for each sample with filtering
│   │   ├── Ill_SRR24317966_RelativeFreqPlot.html
│   │   ├── Ill_SRR24317966_filtered_RelativeFreqPlot.html
│   ├── Variable_Region_Comparisons
│   │   ├── SRR24317964_vs_SRR24317966_VariableRegions_DotLine_Filtered.pdf     # Plot comparing variability of reads between samples 
│   │   └── Variable_region_compare.RData                                       # R data to be loaded into R of variable region comparisons
│   ├── all_heatmap.html                                                        # Heatmap of all variable region variability for all samples, will be separated by --large option
│   └── all_variable_regions.html                                               # Barchart of all variable region variability for all samples, will be separated by --large option
├── Tables
│   ├── Ill_SRR24317964_final_data.csv                                          # separated table of counts as nucleotide for each sample
│   ├── Ill_SRR24317964_overcount_final_AA_data.csv                             # separated table of counts as amino acid for each sample
│   ├── Ill_SRR24317964_summary_statistics.csv                                  # Total read counts per variable region per sample
│   ├── Ill_SRR24317966_final_data.csv
│   ├── Ill_SRR24317966_overcount_final_AA_data.csv
│   ├── Ill_SRR24317966_summary_statistics.csv
│   ├── all_assignments.csv                                                     # dataframe of reads uncounted per variable site for all samples
│   ├── all_summary_stats.csv                                                   # counts for all samples for all variable regions
│   ├── allreads.csv
│   ├── allreads_abs.csv                                                        # all reads csv with only counts
│   ├── allreads_filt_abs.csv                                                   # all reads csv with only counts filtered
│   ├── allreads_filt_rel.csv                                                   # all reads csv with only frequency filtered
│   ├── allreads_rel.csv                                                        # all reads csv with only frequency
│   └── compare_illumina_df.csv
├── allreads.csv                                                                # Counts and Frequencies for each variable site 
├── allreads_filtered.csv                                                       # Filters dataframe if RF or Count filter option is applied
└── allreads_filtered_heatmap.csv                                               # CSV used to generate heatmap with less filtering (with default filtering values allreads csvs will be identical)
```

Pacbio Output

```
Example_PACBIO_OUT
├── Figures
│   ├── Relative_Frequency_Plots
│   │   ├── PB_SRR10294238_RelativeFreqPlot.html                                # Bar chart comparing variability between variable sites
│   │   ├── PB_SRR10294238_filtered_RelativeFreqPlot.html                       # Bar chart comparing variability between variable sites for each sample
│   │   ├── PB_SRR10294254_RelativeFreqPlot.html
│   │   └── PB_SRR10294254_filtered_RelativeFreqPlot.html
│   └── Tree
│       ├── Isolates_aa_filt_fullORFs.aln.aln.tree.nwk                          # PacBio newick tree file
│       ├── Isolates_aa_filt_fullORFs.fasta                                     # Fasta converted to AA filtered on RF cuttoff, default is identical to non filtered                                    
│       ├── Isolates_aa_fullORFs.fasta                                          # Fasta converted to AA filtered on RF cuttoff
│       ├── PacBio_Tree_Filtered.RData                                          # Pacbio tree R data
│       ├── PacBio_Tree_Filtered.pdf                                            # Pacbio tree
│       └── Table_allAAfilt_fullORFs.tsv                                        # Table of counts per sequence
├── Tables
│   └── Frequency_Tables
│       ├── PB_SRR10294238.noprimers.filtered.RAD.nolines.fix_final_data.csv                # Counts and relative frequency of sequences for sample as nucleic acids
│       ├── PB_SRR10294238.noprimers.filtered.RAD.nolines.fix_overcount_final_AA_data.csv   # Counts and relative frequency of sequences for sample as AA
│       ├── PB_SRR10294238.noprimers.filtered.RAD.nolines.fix_summary_statistics.csv        # Counts per variable region
│       ├── PB_SRR10294254.noprimers.filtered.RAD.nolines.fix_final_data.csv
│       ├── PB_SRR10294254.noprimers.filtered.RAD.nolines.fix_overcount_final_AA_data.csv
│       ├── PB_SRR10294254.noprimers.filtered.RAD.nolines.fix_summary_statistics.csv
│       ├── all_assignments.csv # All amino acids and nucleic acid sequences for samples
│       └── compare_pacbio_df.csv # all reads csv with counts and frequency
└── denoised_fastas
    ├── PB_SRR10294238.noprimers.filtered.RAD.nolines.fix.fasta                             # Denoised fasta for each sample
    └── PB_SRR10294254.noprimers.filtered.RAD.nolines.fix.fasta
```
