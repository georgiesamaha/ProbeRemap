# ProbeRemap 
This workflow can be used to map SNP probes from either an Illumina genotyping array manifest or PLINK .bim file to a supplied reference assembly using BLAST, filter results for best hit, and output updates genomic positions. This workflow is a python implementation of Willet & Haase, 2015 (DOI: 10.1111/age.12169). It accepts a series of markers with known genomic positions, collects flanking sequences up and downstream of the target SNP, BLASTs them against a provided reference assembly, and outputs updated marker positions, as well as unmapped and multi-mapping sequences. The output updated positions can be used to update a cohort array dataset in PLINK. 

These scripts have been used in Samaha et al. 2020 (DOI: 10.1038/s41598-020-76166-3) to remap array markers to the felCat9 reference assembly.  

## Workflow

<img src="https://github.com/georgiesamaha/ProbeRemap/blob/main/probe_remap_workflow.png" width="80%" height="80%"> 

## Software
- python3
- blast+/2.9.0

## Set up 
### Clone repository 

`git clone https://github.com/georgiesamaha/ProbeRemap.git`

### Prepare input
#### Manifest 
- Expected input is either a [PLINK .bim file](https://www.cog-genomics.org/plink/1.9/formats#bim) or an Illumina Infinium array manifest.   
- These scripts have been tested on canine and feline Infinium arrays. Manifests are available for download on [Illumina's website](https://sapac.support.illumina.com/array/downloads.html). 
- Other inputs such as BED format files can be converted to a relevant input file format using a conversion method [described here](https://sapac.support.illumina.com/bulletins/2016/05/basespace-sequence-hub-how-to-convert-a-custom-bed-file-to-a-manifest-file-for-enrichment-analysis.html).  

#### Reference 
- Reference file in fasta format for target genome to remap probe sequences 

## Quickstart guide
1. Run `python create_probes.py <manifest.csv> --illumina`

This will take the manifest file, create probes using flanking sequence for each marker and output two files `probes.seq` and `probes.fasta` in a directory called `./Output`. A FASTA file will be created using the SNV identifiers and flanking genomic sequence upstream and downstream of each SNV in top orientation are collected. The longest flanking sequence will be retained and output in FASTA format in `probes.fasta`. 

2. Run `bash probeblast.sh <ref.fasta>` 

This will take the probe sequence files output to `./Output` and the provided reference sequence and use them to create a custom BLAST database using BLAST+ to make the reference FASTA sequence searchable. A nucleotide BLAST search will then be performed and the results will be reformatted and output to `./Output/<reference.fasta>.BLASTn_Output`. 

3. Run `filter_blast.py ./Output/<reference.fasta>.BLASTn_Output`

This will filter BLAST results for all probe sequences to collect the single best hit for each sequence. Hits that failed to return any result, hits that failed to reach the base position immediately adjacent to the SNV and those that had an E-value greater than 1e−05 will rejected. You can edit Evalue within the script. A summary log will be printed to the screen and files will be output for unmapped markers `./Output/<reference>.failedBlast.txt`, multi-mapped markers `./Output/<reference>.multipleBLAST.txt`, and successfully updated positions `./Output/<reference>.newpos.txt`.   

## Acknowledgements 
Written by Mitchell O'Brien and Georgie Samaha 
