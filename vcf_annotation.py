#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 15:07:53 2020

@author: jiang
"""

"""
         Tempus Bioinformatics Technical Challenge
         
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

Please upload all relevant code (written in whatever language you like) to a public repo on
GitHub. Be sure to include the annotated variants in a csv/tsv file. Please also upload a
README with the link to the repo via dropbox. Note that work will be assessed based on quality
of code and documentation more-so than the annotation.

"""
import argparse
import logging
import sys
from collections import namedtuple
import json
import requests
import time

## For log
logging.basicConfig(level=10,
        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
        datefmt='%a, %d %b %Y %H:%M:%S',
        stream=sys.stdout,
        filemode="w"
        )

def read_vcf(infile,type_tag="TYPE",depth_tag="DP",alt_depth_tag="AO"):
    """
    read a vcf file and return a list of namedtuple
    When variants is multiallelic variant, break the variants into multiple records (ALT column has ",")
    For example:
    chr1  1000 . A C,G 
    
    Parameters
    ----------
    infile : string
        Input vcf file name

    Returns
    -------
    variants : list
        A list of variants in namedtuple object(CHROM,POS,ID,REF,ALT,VTYPE,EFFECT,DEPTH,ALT_DEPTH,REF_DEPTH,ALT2REF,ALLELE_FREQ,ANNO)
        
    CHROM        chromosome from vcf
    POS          position from vcf
    ID           id from vcf, if value is "" or ".", will be updated based on the information get from "http://exac.hms.harvard.edu//rest/bulk/variant/variant
    REF          reference allele
    ALT          if the current recode is multiallelic variant, will break down to non multiallelic variant
    VTYPE        get from the INFO column based on 'TYPE' tag
    EFFECT       information from http://exac.hms.harvard.edu//rest/bulk/variant/variant (Consequence)
    DEPTH        get from the INFO column based on 'DP' tag
    ALT_DEPTH    get from the INFO column based on 'AO' tag,
    REF_DEPTH    calculated by DEPTH-all ALT_DEPTH, if the variants in the vcf is multiallelic variant, in the output DEPTH will not equal to ALT_DEPTH+REF_DEPTH
    ALT2REF      ALT_DEPTH/REF_DEPTH
    ALLELE_FREQ  information from http://exac.hms.harvard.edu//rest/bulk/variant/variant  (allele_freq)
    ANNO         information from http://exac.hms.harvard.edu//rest/bulk/variant/variant  ('Consequence','SYMBOL','HGVSc','HGVSp','SIFT','PolyPhen','BIOTYPE')
        
    """
    ## namedtuple for each of the variants' records
    VariantRecord=namedtuple('VariantRecord',"CHROM,POS,ID,REF,ALT,VTYPE,EFFECT,DEPTH,ALT_DEPTH,REF_DEPTH,ALT2REF,ALLELE_FREQ,ANNO")

    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	normal	vaf5
    variants=[]  ## return value
    with open(infile,"r") as f:
        for line in f:
            ## Only variant record
            if not line.startswith("#"):
                lst=line.rstrip().split("\t")  ## convert to list
                chrom,pos,ID,ref,alt,info=lst[0],lst[1],lst[2],lst[3],lst[4],lst[7]
                
                if "," in alt:
                ## Need to break down current line into multiple variant records
                    alts = alt.split(",")
                    depth=int(get_value_from_INFO(info,depth_tag))
                    vtypes=[]      ## store variant type for each record
                    alt_depths=[]  ## store alternative allele depth for each record
                    for i in range(len(alts)):
                        vtype=get_value_from_INFO(info,type_tag,i)               ## variant type, either snp,del,ins,mmp or complex
                        alt_depth=int(get_value_from_INFO(info,alt_depth_tag,i)) ## alt allele depth
                        vtypes.append(vtype)
                        alt_depths.append(alt_depth)
                    
                    ## Now get the ref_depth
                    ref_depth=depth-sum(alt_depths)
                    ## Loop again and append to variants
                    for i in range(len(alts)):
                        alt2ref=None
                        if ref_depth>0:
                            alt2ref=alt_depths[i]/ref_depth
                        vv = [chrom,pos,ID,ref,alts[i],vtypes[i],None,depth,alt_depths[i],ref_depth,alt2ref,None,None]
                        variants.append(VariantRecord._make(vv))
                        
                else:
                    depth=int(get_value_from_INFO(info,depth_tag)) ## 
                    vtype=get_value_from_INFO(info,type_tag)  ## variant type, either snp,del,ins,mmp or complex
                    alt_depth=int(get_value_from_INFO(info,alt_depth_tag)) ## alt allele depth
                    ref_depth=depth-alt_depth
                    alt2ref=None
                    if ref_depth>0:
                        alt2ref=alt_depth/ref_depth
                        VariantRecord=namedtuple('VariantRecord',"CHROM,POS,ID,REF,ALT,VTYPE,EFFECT,DEPTH,ALT_DEPTH,REF_DEPTH,ALT2REF,ALLELE_FREQ,ANNO")

                    vv = [chrom,pos,ID,ref,alt,vtype,None,depth,alt_depth,ref_depth,alt2ref,None,None]
                    variants.append(VariantRecord._make(vv))
    return variants
   
def get_value_from_INFO(info,TAG="DP",idx=None):
    """
    Get value from the INFO based on TAG, default to get  Depth of sequence coverage at the site of variation
    
    DP for Depth of sequence coverage at the site of variation
    AO for Number of reads supporting the variant.

    example of INFO looks like 
    AB=0;ABP=0;AC=0;AF=0;AN=6;AO=95;CIGAR=1X;DP=4124;DPB=4124;DPRA=0.999031;EPP=9.61615;EPPR=316.776;GTI=0;LEN=1;MEANALT=1;MQM=59.7579;MQMR=65.2274;NS=2;NUMALT=1;ODDS=591.29;PAIRED=0.989474;PAIREDR=0.966741;PAO=0;PQA=0;PQR=0;PRO=0;QA=3774;QR=160284;RO=4029;RPL=51;RPP=4.13032;RPPR=101.278;RPR=44;RUN=1;SAF=40;SAP=8.15326;SAR=55;SRF=1663;SRP=269.369;SRR=2366;TYPE=snp'

    Parameters
    ----------
    info : string
        info from INFO column of VCF.
    TAG : string, optional
        TAG in the INFO. The default is "DP".
    idx: int
        if not None, we need to idxth value of the tag which is seprated by ","

    Returns
    -------
    d : string
        values of the TAG.

    """
    d = None    
    for pair in  info.split(";"):
        k,v=pair.split("=")
        if k == TAG:
            if idx is None:
                d = v
            else:
                d =v.split(",")[idx]
            break
    return d
    
def fix_chr(chrom, with_chr):
    """
    Fix the chrom, whether want to return with 'chr1' or just "1"

    Parameters
    ----------
    chrom : string
        Original chromosome name.
    with_chr : bool
        Whether want to have 'chr' or not.

    Returns
    -------
    string
        

    """
    if with_chr:
        return chrom if chrom.startswith('chr') else 'chr' + chrom
    else:
        return chrom if not chrom.startswith('chr') else chrom[3:]
        return chrom if not chrom.startswith('chr') else chrom[3:]


def info_from_ExAC(lst,url="http://exac.hms.harvard.edu//rest/bulk/variant/variant"):
    """
    Use ExAC bulk query API get get variants information.

    Parameters
    ----------
    lst : list
        a list of variants.
    url : string, optional
        url address for ExAC bulk query API. The default is "http://exac.hms.harvard.edu//rest/bulk/variant/variant".

    Raises
    ------
    Exception
        query failed.

    Returns
    -------
    r : dict
        return value from API.
    """
    ## Make json array for the api post
    chr_pos_ref_alt=[]
    for v in lst:
        chr=fix_chr(v.CHROM,False)
        chr_pos_ref_alt.append(chr+"-"+v.POS+"-"+v.REF+"-"+v.ALT) 
    
    jsondata=json.dumps(chr_pos_ref_alt)
    response = requests.post(url, data=jsondata)
    r = response.json()
    if "errors" in r:
        raise Exception(str(r["errors"]))
    return r


def divide_list_to_chunks(lst,size):
    ## Loop with step of size
    for i in range(0,len(lst),size):
        yield lst[i:i+size]
        
def add_annotation_to_variants(chunk,anno):
    """
    Add annotation from ExAC to the current variants, which will be the allele_freq and ANNO

    Parameters
    ----------
    chunk : list
        a list of variants.
    anno : dict
        annotation from ExAC, 'chr-pos-ref-alt' as the key.

    Returns
    -------
    list
        a list of variant with annotatin added from ExAC.

    """
    
    if anno is None:
        return chunk
    vep_wants=['Consequence','SYMBOL','HGVSc','HGVSp','SIFT','PolyPhen','BIOTYPE']
    rsid='rsid'
    r=[]
    for v in chunk:
        chrom=fix_chr(v.CHROM,False)
        k = chrom+"-"+v.POS+"-"+v.REF+"-"+v.ALT
        if k in anno:
            ## Get allele frequency etc information from http://exac.hms.harvard.edu//rest/bulk/variant/variant
            if 'allele_freq' in anno[k]:
                v = v._replace(ALLELE_FREQ=anno[k]['allele_freq'])           
            ## Update ID
            if (v.ID=="." or v.ID=="") and rsid in anno[k]:
                v = v._replace(ID=anno[k][rsid])           
            ## VEP annotation
            if 'vep_annotations' in anno[k]:
                veps = anno[k]['vep_annotations']
                ## Get most frequence Consequence(variants type)
                Consequences=[]
                annotations=[]
                for vep in veps:
                    tmp=[] ## For full annotation
                    for want in vep_wants:
                        if want in vep:
                            tmp.append(vep[want])
                        else:
                            tmp.append("")
                    annotations.append(";".join(tmp))
                   
                    ## Get Consequences
                    if 'Consequence' in vep:
                        Consequences.append(vep['Consequence'])
                   
                ## Now update the v
                if Consequences:
                    v=v._replace(EFFECT=Consequences[0])
                if annotations:
                    v=v._replace(ANNO="|".join(annotations))
        r.append(v)
    
    return r
               
            
               
            
           
def RUN(infile,outfile,size=1000,type_tag="TYPE",depth_tag="DP",alt_depth_tag="AO",maxtry=3,sleeptime=5):
    ## Read variants from vcf file
    variants = read_vcf(infile,type_tag,depth_tag,alt_depth_tag)
    logging.info("Loading {} variants".format(len(variants)))
    updated_variants=[]
    ## Break down to a list of list, the maxsize of the sub-list is 100(size), 
    ## This is used to do bulk annotation from http://exac.hms.harvard.edu//rest/bulk/variant/variant
    variants_chunks=list(divide_list_to_chunks(variants,size))
    n_chunks=len(variants_chunks)
    
    for i in range(n_chunks):
        chunk = variants_chunks[i]
        logging.info("ExAC annotates {}/{} chunks with {} variants".format(i+1,n_chunks,len(chunk)))
        tried=0
        anno=None
        ## Will try maxtry time to query
        while tried<maxtry:
            tried+=1
            time.sleep(sleeptime)
            try:
                anno=info_from_ExAC(chunk)
                break
            except Exception as e:
                logging.error(e)
                logging.info("Try agian - {} times".format(tried))
         
        if anno is None:
            for v in chunk:
                updated_variants.append(v)
        else:
            for v in add_annotation_to_variants(chunk,anno):
                updated_variants.append(v)
    
    ## Now write out
    logging.info("Write out to '{}'".format(outfile))
    with open(outfile,"w") as f:
        f.write("CHROM\tPOS\tID\tREF\tALT\tTYPE\tEFFECT\tDEPTH\tALT_DEPTH\tREF_DEPTH\tALT2REF\tALLELE_FREQ\tANNOTATION\n")
        for v in updated_variants:
            conv = lambda i : i or ''
            vv = [str(conv(i)) for i in list(v)]
            outline="\t".join(vv)
            f.write(outline+"\n")
   
    logging.info("Finished")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Annotation VCF file format', )
    parser.add_argument("--infile",dest="infile",help="Input vcf file",required=True)
    parser.add_argument("--outfile",dest="outfile",help="Output file ",required=True)
    parser.add_argument("--bulk",dest="bulk",help="ExAC bulk query size(default 1000)",default=1000)
    parser.add_argument("--vtype_tag",dest="vtype_tag",help="Variant type TAG in INFO(default TYPE)",default="TYPE")
    parser.add_argument("--depth_tag",dest="depth_tag",help="Variant depth TAG in INFO(default DP)",default="DP")
    parser.add_argument("--altdepth_tag",dest="altdepth_tag",help="Variant alt allele depth TAG in INFO(default AO)",default="AO")
    parser.add_argument("--max_try",dest="max_try",help="Maximum query for each bulk(default 3)",default=3)
    parser.add_argument("--sleep",dest="sleep",help="sleep time before next query(default 5s)",default=5)


    args = parser.parse_args()
    RUN(args.infile,args.outfile,args.bulk,args.vtype_tag,args.depth_tag,args.altdepth_tag,args.max_try,args.sleep)
