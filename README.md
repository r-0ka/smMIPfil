# smMIPfil

What this script does:

(1) Reads in mutation positions provided as comma separated text file ([mutation position]).

(2) Extracts all the reads covering the position from BAM file.

(3) Extracts UMI (unique molecular identifier) sequences and the base at the position for each read.

(4) Count the number of reads containing the same nucleotide for each UMI. 
    If there are less than 5 reads containing a UMI, the UMI is excluded.  

(5) For the nucleotide supported by the majority of reads, calculates the fraction of 
    the nucleotide-supporting-reads in all the reads for the UMI.  

    For example, if a UMI contained 5 A's, 0 C's, 250 G's and 45 T's, the script calculates 
    fraction for G because most reads were sequenced as G at the mutation position.  

    The fraction = read count for G / read count for this UMI
                 = 250 / (5 + 0 + 250 + 45)
                 ~ 0.83

(6) When the fraction was above 0.7, then the UMI is kept with the nucleotide.
    ALl the UMI's with lower fraction are discarded.  

    With the example above, the UMI would be kept as G-supporting.  

(7) At the step (6), the number of supported nucleotides is counted.

(8) UMI VAF for each mutation position is calculated. 



# Dependency

samtools



# Before running (script requirement)

This is a stand-alone perl script.  Except that this is dependent on the samtools, no installation required.  

Please modify the working directory in the script ($wk_dir).
The script assumes there is a tmp directory under the working directory and actually this is the only place this working directory is used. 



# Before running (input requirement)

The UMI sequences need to be included in the RX tag header in the BAM file. 


# Usage:

	perl smMIP_UMI_filtering.pl [mutation position] [input bam] [sample name] > [output]

[mutation position]:

  Comma separated text file with "chrom,position,refBase,mutBase" 

[input bam]:

  SORTED BAM file.
  Before the alignment, the smMIP fastq has to be clipped so that the 
  UMI's could be found back.

[sample name]:

  Unique name for each sample runs to prevent the intermediate files to be overwritten.
  Does not affect on the output. 



# Output

Tab delimited text file with the following columns: 
  
* chr

  Chromosome number (provided by the input file)
  
* pos

  Chromosome coordinate (provided by the input file)
  
* ref

  Reference base (provided by the input file)
  
* alt

  Mutated variant base (provided by the input file)
  
* totalUMI

  Total number of UMI covering the position
  
* A

  Total number of UMI supporting A at the position
 
* T

  Total number of UMI supporting T at the position
  
* G

  Total number of UMI supporting G at the position
  
* C

  Total number of UMI supporting C at the position
  
* VAF

  Allele frequency for alt-base calculated from the number of UMI's
