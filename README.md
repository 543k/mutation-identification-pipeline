# Goal: Identify SNPs causing colored mold in each patient's ears, and generate a report.
# Author: Vivek Mathew
# Email: vivekmathew@brandeis.edu
# Date: 11/25/2024


## Python Pipeline Script ###

### Description ###
A single Python script named pipeline.py that generates trimmed fastq files for each patient, aligns each fastq file to the reference FASTA file (dgorgon_reference.fa), and then generates sorted BAM files for variant discovery.
A final report is written for the variant discovery (report.txt).

Input:

A txt file containing clinical data (named harrington_clinical_data.txt).
A FASTA file called dgorgon_reference.fa containing a reference sequence for alignment.
A FASTQ file called hawkins_pooled_sequences.fastq containing pooled sequence data of all patients.

Output:

FASTQ and sorted BAM files associated with each patients, found in the fastqs and bams directories respectively.

How to Run:

To execute the script, navigate to the directory containing the script file (ensuring the above input files are in the same directory) and run:

chmod +X pipeline.py
python3 pipeline.py
