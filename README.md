## SPECTRA 

#### Peptide and Protein identification from LC-MS/MS data


## What is Proteomics?
Proteomics involves the applications of technologies for the identiﬁcation and quantiﬁcation of overall proteins present content of a cell, tissue or an organism. 
It supplements the other “omics” technologies such as genomic and transcriptomics to expound the identity of proteins of an organism, and to cognize the structure and functions of a particular protein. Proteomics-based technologies are utilized in various capacities for different research settings such as detection of various diagnostic markers, candidates for vaccine production, understanding pathogenicity mechanisms, alteration of expression patterns in response to different signals and interpretation of functional protein pathways in different diseases. 

## What is Mass spectrometry?
Mass spectrometry is a popular method to analyse bio-molecules by measuring the intact mass-to-charge ratios of their in-situ generated ionised forms or the mass-to-charge ratios of in-situ-generated fragments of these ions. The resulting mass spectra are used for a variety of purposes, among which is the identification, characterization, and absolute or relative quantification of the analysed molecules. 
Mass spectrometry with LC–MS-MS is the central among current proteomics.

## What does SPECTRA do?

* Preview the relevant spectrum info of the high throughput LC-MS/MS data (.mzML) 
* Implement De Novo sequencing to determine most probable peptides from one scan of MS2 data
* Match the resulting peptides with known proteins from reliable database

## Workflow
* Extract MS2 spectrum file; 
* M_Z/ intensity/ precursor info from one MS2 scan (here you should specify your candidate scan number, e.g. 0 or 3);
* MS2 de-isotoped spectrum plotting
* Peptide identification:  De novo sequencing
* API query: Match peptides to candidate proteins --> API request

  UniProt/ PIR (Protein Information Resource)
* Construct frontend GUI (SPECTRA) to show results
* Construct Docker images and containers

## Data

**1. reduced_data.mzML** 
* Raw LC-MS/MS data file
* Reduced version from one large raw data file 
* Can be used to extract MS2 info
* Organism: human serum

**2. ms2.mzML**
* MS2 data file
* Can also be generated by: **Parse_data_MS2** module
* Contains important info of precursor M/Z and charge, fragment M/Z and intensity
* **Used for peptide identification and GUI uploading**


## Project_spectra_temp
**1. project_spectra package**

* Parse_data_MS2.py: MS2 file generation and Parse data (M/Z, intensity, precursor M/Z and charge)

* Spec_peptide.py: peptide identification (The output is a list contains possible peptides, peak intensity related to each AA and sum of relative intensity as SCORE.)

* identifier.py: candidate protein prediction (The output is a dictionary: Accession number as key, Protein full name and Organism in a list as value)

* cli.py: CLI functions

**2. test folder**
* test functions


**3. The setup script for package**
* setup.py

## Frontend
**1. template**
* html files for GUI

**2. static**
* CSS and image folder for html

**3.script**
* run.py

## Dockerfile
* Docker script for Container

## README.md
* Description about SPECTRA project and related info

## Parameter Setting in peptide identification

* peaks searching window: 57-187
* C-terminal AA identification tolerance error: 1 Da
* peptide to precursor matching tolerance error: 1 Da

## CLI 

**1. Ms2_cli**: 

python cli.py get_ms2 /path/to/data/reduced_data.mzML /path/to/data/ms2.mzML

**2. check_scan_size**

python cli.py check_scan_size /path/to/data/ms2.mzML

**3. parse_cli**:

python cli.py parse /path/to/data/ms2.mzML 0

**4. get_plot_spectrum**:

python cli.py plot /path/to/data/ms2.mzML 0 --verbose

**5. identify_C_terminal**:

python cli.py check_C /path/to/data/ms2.mzML 0

**6. identify_N_terminal**:

python cli.py check_N /path/to/data/ms2.mzML 0

**7. get_peptide**:

python cli.py peptide_list /path/to/data/ms2.mzML 0 --verbose

**8. get_protein**:

python cli.py protein /path/to/data/ms2.mzML 0

## GUI
**Input**: 
* In the web interface, you should upload a **MS2 mzML** file (which can also be generated using our function),
and one specific scan number to submit your query.

**Output**: 
* Precursor info
* A table of number of peaks in this scan, max/min M/Z peaks, peak with highest intensity.
* A plot of de-isotoped spectrum.
* A table of candidate peptides derived from this scan spectrum and its sum intensity/ score; lists of candidate peptide and its corresponding peak intensity
* A table of all candidate proteins.

## Docker
$ docker build . -t project_4:latest

$ docker run --name project_test -p 5005:5005 -d project_4:latest

## Useful Link about De Novo Sequencing
* **Wikipedia**

   [De novo peptide sequencing](https://en.wikipedia.org/wiki/De_novo_peptide_sequencing)
* **Youtube - Matthew Padula**

   [De novo peptide sequencing from mass spectrometry data](https://www.youtube.com/watch?v=zZsI8GM45gc)

## Library and API

* [pyopenms](https://pyopenms.readthedocs.io/en/latest/): Parse data

* [UniProt](https://www.ebi.ac.uk/proteins/api/proteins/): UniProt API for protein identification

* [PIR](https://research.bioinformatics.udel.edu/peptidematchapi2/match_get?peptides=): PIR peptide match API

## Comments
* For programming lab exercise, please **do not** take it as commercial use case.
* We recommend you to use our data for testing related functions and CLI, since some .mzML data from Web do not have information of precursor charge. 
* Unfortunately not all scans can be identified by our algorithm, you can try scan 0 or 3 for testing.


