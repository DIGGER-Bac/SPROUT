# INTRODUCTION

SPROUT (Small RNA Profiling, Reconstruction, and Optimization Utility Tool) contains four programs: Class I, Class II, sRNA\_create and sRNA\_visualize.

With this tool, class I and II small RNAs (sRNAs) can be analyzed with regard to motifs within their scaffold sequences (function sRNAlyze in Class\_I and Class\_II).
Based on this analysis, new synthetic sRNA scaffolds can be generated based on the consensus motifs (function sRNAgen in Class\_I and Class\_II).
Additionally, synthetic sRNAs can be generated from scratch and the consensus motifs can be modified as needed by the user (sRNA\_create).
Another functionality is the visualization of the analyzed/created sRNAs (sRNA\_visualize).



Class I and II are very similar to operate, but differ in the Hfq binding motifs they recognize.
Class I: UA-rich sequence, gap, stem-loop, poly-U sequence
Class II: ARN motif, gap, stem-loop, poly-U sequence

sRNA\_create allows to completely customize an sRNA with the help of a text-based code.
sRNA\_visualize shows the 2D-structure of the input sequences and the positional entropy of the structure.



You can find the exact instructions how to use the programs below.



# REQUIREMENTS

* Podman or Docker Desktop installed
* folder SPROUT\_Linux downloaded



# FIRST TIME START

### in Podman:

* load the image within the container application:
quay.io/biocontainers/viennarna:2.6.4--py38pl5321h5cf8b27\_0

### in Docker:

* load the image using the command line input:
docker pull quay.io/biocontainers/viennarna:2.6.4--py38pl5321h5cf8b27\_0



# INSTRUCTIONS

## CREATE CONTAINER

* while creating the container using the desktop interface in Docker or Podman,
mount the host path "\~/SPROUT\_Linux" (the "\~" represents the complete path on your computer where you stored the folder) on the container path "/home"
* after starting the container, openpyxl has to be installed. You can do this with the command: "pip install openpyxl"



## CLASS I

### ANALYZE EXISTING SEQUENCES (sRNAlyze)

* save the sRNA sequences you want to analyze in .fasta format in "\~/SPROUT\_Linux/Class\_I/sRNAlyze/Data"
(a database with sRNAs from enterobacteria can be found in "\~/SPROUT\_Linux/database\_sequences")
* inside the container terminal, navigate to "/home/Class\_I/sRNAlyze"
* start analyzing with the command "python3 sRNA\_analyze\_I.py"
(sRNAlyze compares the saved sRNA sequences and finds common motifs)
* the information of the consensus motifs is stored in "\~/SPROUT\_Linux/Class\_I/sRNAlyze/Output/consensus\_conclusion.txt"



### GENERATE ONE SET OF NEW sRNA SCAFFOLD SEQUENCES (sRNAgen)

* the key data of the identified motifs (from sRNAlyze) are directly written to the config file of sRNAgen ("\~/SPROUT\_Linux/Class\_I/sRNAgen/config\_single.yml")
* you can adjust the config file if needed ("\~/SPROUT\_Linux/Class\_I/sRNAgen/config\_single.yml")
* inside the container terminal, navigate to "/home/Class\_I/sRNAgen"
* start generating with the command "python3 sRNA\_generate\_single\_I.py"
* the generated sequences can be found in the Output folder in two formats: .fasta and .xlsx



### GENERATE MULTIPLE SETS OF NEW sRNA SCAFFOLD SEQUENCES (sRNAgen)

* if you want to generate more than one set of sRNA scaffolds, you can copy the consensus output lines from different "\~/SPROUT\_Linux/Class\_I/sRNAgen/config\_single.yml" files and paste them into "\~/SPROUT\_Linux/Class\_I/sRNAgen/config\_multiple.yml"
* after that, navigate to "/home/Class\_I/sRNAgen"
* start generating with the command "python3 sRNA\_generate\_multiple\_I.py"
* the generated sequences can be found in the Output folder in two formats: .fasta and .xlsx



## SHORT INSTRUCTIONS CLASS I (example for one set of new scaffolds)

### command line input inside the container:

* cd home/Class\_I/sRNAlyze
* python3 sRNA\_analyze\_I.py
* cd ..
* cd sRNAgen
* python3 sRNA\_generate\_single\_I.py



## CLASS II

### ANALYZE EXISTING SEQUENCES (sRNAlyze)

* save the sRNA sequences you want to analyze in .fasta format in "\~/SPROUT\_Linux/Class\_II/sRNAlyze/Data"
(a database with sRNAs from enterobacteria can be found in "\~/SPROUT\_Linux/database\_sequences")
* inside the container terminal, navigate to "/home/Class\_II/sRNAlyze"
* start analyzing with the command "python3 sRNA\_analyze\_II.py"
(sRNAlyze compares the saved sRNA sequences and finds common motifs)
* the information of the consensus motifs is stored in "\~/SPROUT\_Linux/Class\_II/sRNAlyze/Output/consensus\_conclusion.txt"



### GENERATE ONE SET OF NEW sRNA SCAFFOLD SEQUENCES (sRNAgen)

* the key data of the identified motifs (from sRNAlyze) are directly written to the config file of sRNAgen ("\~/SPROUT\_Linux/Class\_II/sRNAgen/config\_single.yml")
* you can adjust the config file if needed ("\~/SPROUT\_Linux/Class\_II/sRNAgen/config\_single.yml")
* inside the container terminal, navigate to "/home/Class\_II/sRNAgen"
* start generating with the command "python3 sRNA\_generate\_single\_II.py"
* the generated sequences can be found in the Output folder in two formats: .fasta and .xlsx



### GENERATE MULTIPLE SETS OF NEW sRNA SCAFFOLD SEQUENCES (sRNAgen)

* if you want to generate more than one set of sRNA scaffolds, you can copy the consensus output lines from different "\~/SPROUT\_Linux/Class\_II/sRNAgen/config\_single.yml" files and paste them into "\~/SPROUT\_Linux/Class\_II/sRNAgen/config\_multiple.yml"
* after that, navigate to "/home/Class\_II/sRNAgen"
* start generating with the command "python3 sRNA\_generate\_multiple\_II.py"
* the generated sequences can be found in the Output folder in two formats: .fasta and .xlsx



## SHORT INSTRUCTIONS CLASS II (example for one set of new scaffolds)

### command line input inside the container:

* cd home/Class\_II/sRNAlyze
* python3 sRNA\_analyze\_II.py
* cd ..
* cd sRNAgen
* python3 sRNA\_generate\_single\_II.py



## sRNA\_create

* with the file "\~/SPROUT\_Linux/sRNA\_create/blueprint.txt", you can design and combine sRNA motifs (the exact procedure is described in this file)
* after that, navigate to "/home/sRNA\_create"
* start creating with the command "python3 sRNA\_create.py"
* the generated sequences can be found in the Output folder in two formats: .fasta and .xlsx



## sRNA\_visualize

* save the sRNA sequences you want to visualize in .fasta format in "\~/SPROUT\_Linux/sRNA\_visualize/Input"
* after that, navigate to "/home/sRNA\_visualize"
* start visualizing with the command "python3 sRNA\_visualize.py"
* the visualization can be found in the Output folder as a .xlsx file



## Overview of the folder structure

[Structure.txt](https://github.com/DIGGER-Bac/MOISten/files/14463396/Structure.txt)

