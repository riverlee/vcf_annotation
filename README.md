---
output:
  pdf_document: default
  html_document: default
---
# A simple tools to annotate VCF files

## Requirement and implementation to each points

For this challenge, you are asked to prototype a variant annotation tool. We will provide you with
a VCF file, and you will create a small software program to annotate each variant in the file.

Each variant must be annotated with the following pieces of information:

  1. Type of variation (substitution, insertion, CNV, etc.) and their effect (missense, silent,
intergenic, etc.). If there are multiple effects, annotate with the most deleterious
possibility. 
  2. Depth of sequence coverage at the site of variation.
  3. Number of reads supporting the variant.
  4. Percentage of reads supporting the variant versus those supporting reference reads.
  5. Allele frequency of variant from ExAC API (API documentation is available here:
http://exac.hms.harvard.edu/).
  6. Any additional annotations that you feel might be relevant.
  
## Usage

```
usage: vcf_annotation.py [-h] --infile INFILE --outfile OUTFILE [--bulk BULK]
                         [--vtype_tag VTYPE_TAG] [--depth_tag DEPTH_TAG]
                         [--altdepth_tag ALTDEPTH_TAG] [--max_try MAX_TRY]
                         [--sleep SLEEP]

Annotation VCF file format

optional arguments:
  -h, --help            show this help message and exit
  --infile INFILE       Input vcf file
  --outfile OUTFILE     Output file
  --bulk BULK           ExAC bulk query size(default 1000)
  --vtype_tag VTYPE_TAG
                        Variant type TAG in INFO(default TYPE)
  --depth_tag DEPTH_TAG
                        Variant depth TAG in INFO(default DP)
  --altdepth_tag ALTDEPTH_TAG
                        Variant alt allele depth TAG in INFO(default AO)
  --max_try MAX_TRY     Maximum query for each bulk(default 3)
  --sleep SLEEP         sleep time before next query(default 5s)
```

## Example

```
cd test
python ../vcf_annotation.py --infile ../Challenge_data_\(1\).vcf --outfile result.tsv 2>&1|tee run.log
```


## Details to each requirment

Note: **For variants in the VCF which are multiallelic, we will break it down into multiple variant records. For those variant records, the ALT_DEPTH+REF_DEPTH is not equal to DEPTH in the output file**


  1. Type of variation is implemented by taking from the **TYPE in the INFO column**. Their effect is implemented by taking from the **ExAC annotation (vep_annotations->Consequence)**. 
  2. Depth of sequence coverage at the site of variation is implemented by taking from the **DP in the INFO column**.
  3. Number of reads supporting the variant is implemented by taking from the **AO in the INFO column**.
  4. Percentage of reads supporting the variant versus those supporting reference reads is implemented as **ALT_DEPTH/REF_DEPTH**.
  5. Allele frequency of variant from ExAC API - **Use bulk query http://exac.hms.harvard.edu//rest/bulk/variant/variant**.
  6. Any additional annotations that you feel might be relevant - **In the output file, last column is effect;gene;acid_change;protein_change;SIFT;PolyPhen;BIOTYPE**.
  
## Things can be improved

- Use try/exception on the depth(DP) and alt depth (AO) when converting them into integer(int)
- Current load all the vcf records into a list. When vcf is large(millions records), we can implemete it in a way that read a certain number of variants and  process them, then read the next certaion number of variants,process them until read all variants.
- Type of variation could be implemeetd as a function instead of reading it from INFO
- Variants effection could be implemeted local as a function by taking an gene annotation file. (E.g, GTF format)
- Nice to add testing code, either unittest or pytest
- 