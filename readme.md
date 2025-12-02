This repository contains code used to create figures for NanoRAPID paper.

Data however is not uploaded and need to be generated and pasted into the correct folders inside data directory. 

Most of the code is inside the "scripts" folder - see there for relevant code.

If you have any questions, please contact renzh6@mail2.sysu.edu.cn

---


# NanoRAPID V1.0

NanoRAPID is a deep learning framework designed for direct RNA structure profiling using nanopore direct RNA sequencing. By extracting structural signatures directly from raw electrical signals, NanoRAPID identifies probe-modified sites at single-molecule resolution without requiring external training standards. The model incorporates iterative noise-reduction and self-distillation strategies, enabling reliable detection of subtle structural features even under high signal variability. NanoRAPID generalizes well across diverse RNA backgrounds and can be readily extended to multiple types of chemical probing. Compared with traditional statistical approaches, it provides higher sensitivity to weak or context-dependent structural signals and offers a scalable solution for transcriptome-wide RNA structure analysis.

---

#### System Requirements

* **Operating System:** Linux
* **Python Version:** ≥ 3.8

#### Dependencies:

* dorado, guppy_basecaller
* nanopolish
* minimap2, samtools, cutadapt, seqkit2

> Ensure all dependencies are installed and available in the system PATH.

#### Required Files:

* **Reference genome sequence (.fa)**
  The reference FASTA file used for read mapping and transcript annotation.
* **Annotation file (.gtf)**
  Gene and transcript annotation file in GTF format, required for isoform assignment and expression quantification.

---

## Analysis Workflow

##### Step 1. DRS-basecall-mapping

```
./basecaller.sh /path/to/my_project gpu

# 
./mapping.sh  /path/to/my_project  /path/to/reference

#
eventalign.sh  /path/to/my_project  /path/to/reference
```

* **`./basecaller.sh /path/to/my_project 0`**
  * **Purpose:** This step takes the raw electrical signal data (often FAST5 files) and converts it into DNA/RNA sequences (FASTQ files). This process is called ​**basecalling**​.
  * **Arguments:**
    * `/path/to/my_project`: The directory containing the raw sequencing data.
    * `gpu`: Specifies the **GPU ID** to use for accelerated basecalling.
* **`./mapping.sh /path/to/my_project /path/to/reference`**
  * **Purpose:** This step aligns the newly generated basecalled reads to a ​**reference genome/transcriptome**​. This produces an alignment file, typically in BAM format.
  * **Arguments:**
    * `/path/to/my_project`: The directory containing the basecalled reads.
    * `/path/to/reference`: The path to the reference sequence file (e.g., FASTA file).
* **`eventalign.sh /path/to/my_project /path/to/reference`**
  * **Purpose:** In ONT analysis, **event alignment** is crucial. It aligns the raw electrical current **events** (from the FAST5 files) directly to the mapped reference coordinates. This is often a prerequisite for analyzing modifications (like RNA modifications) which are detected as subtle changes in the electrical current.
  * **Arguments:**
    * `/path/to/my_project`: The directory containing the raw data and/or alignment files.
    * `/path/to/reference`: The path to the reference sequence file.

##### Step 2.  Extract Features

##### Step 3.  NanoRAPID Predict

```
python -u NanoRAPID-predict.py -d {motif} -m {motif_model_50.pth.tar} -g {gpu} -o ${motif_predict}
```

-d {motif} Specify the RNA motif (sequence context) to analyze
-m {motif\_model\_50.pth.tar} Path to the pre-trained deep learning model for modification prediction
-g {gpu}  Specify the GPU ID to use for accelerated computation (e.g., 0, 1)
-o \${motif\_predict} Define the output file or directory path for the prediction results

---

This software was independently developed by Ze-Hui Ren and integrates multiple open-source tools into a unified analysis pipeline.

All third-party software copyrights remain the property of their original authors.

This software is intended for research and educational purposes only and may not be used for commercial purposes without authorization.

**Developer:** Ze-Hui Ren






