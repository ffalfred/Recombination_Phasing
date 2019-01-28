# Recombination_Phasing
Simple tool for phasing human PB-oocyte trios, as a programming assignment for the candidate informatician at Hoffmann Lab.
It is the pilot stage of the phasing algorithm for phasing  human oocytes from the SNP array data. This is the phasing by a common reference. The algorithm is described in Ottolini et al. 2016 and Ottolini et al. 2015. It reaches the step 45 of the Protocol described in Ottolini et al. 2016, previous to remove the common crossover. The program follows all the rules described in the Protocolo of Ottolini et al. 2016, until the step 45, except for the ones regarding the choosing of a reference.

## Files 
* **phasing_class.py**: Class containing the algorithms
* **phasing_trio.py**: File that uses the algorithms in "phasing_class.py" to produce the output
* **input_data_trios.txt**:  File containing example data to run a test. Contains maternal genotype and genotypes of four trios from the same individual for chromosome 1. The data is tab-delimited. 


## Getting started
Follow the instructions below to use the program.

### Input
The input data contains maternal genotype and genotypes of trios from the same individual. The data is tab-delimited. Each row is an SNP, and contains information about it like the chromosome, the position, the name of the SNP, the maternal genotype of the SNP, the genotype of the SNP for each cell.

### Output
Simple BED file per cell with chr, start coordinate, end coordinate and the phase. The name of the file is: "Trio[number of the trio]Cell[type of cell].bed"

### Prerequisits
The code has been build and uses Python 3.6.5, and the following Python packages:
* math
* numpy
* pandas
* matplotlib
* operator
* argparse

### Clone or downlowad
Clone or download this repo to your local machine using https://github.com/ffalfred/Recombination_Phasing

### Usage
The file that have to be launched is the phasing_trio.py. There are two options/arguments required by the program:
* -i/--input : file containing the input data in format .txt tablimited.
* -m/--min_snps : Amount of SNPs required to define an haploblock

#### Example/test
This is an example of the usage of the program with the example data, and requiring 20 consecutive SNPs to define an haploblock:
```
python phasing_trio.py -i input_data_trios.txt -m 20
```

## Author
Alfred Ferrer Florensa ; f.f.alfred@gmail.com 



