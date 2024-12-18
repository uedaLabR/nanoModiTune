import numpy as np
import pandas as pd


def parseResult(infostr):

    strand =  False
    if "STRAND=True" in infostr:
        strand = True
    infos = infostr.split(",")
    genename = None
    ID = None
    af = 0
    for i in infos:
        if "gene_name=" in i:
            genename = i.replace("gene_name=","")
        if "ID=" in i:
            ID = i.replace("ID=","")
        if "AF=" in i:
            af = float(i.replace("AF=",""))

    return strand,af,ID,genename


def parse_transcript_line(line):
    # Extract basic fields
    chromosome = line[0]
    start = line[3]
    end = line[4]
    strand = line[6]
    attributes_str = line[8]

    # Parse the attributes field
    attributes = {}
    # Split the attributes string into key-value pairs
    attrs = [attr.strip() for attr in attributes_str.strip().split(';') if attr.strip()]
    for attr in attrs:
        # Split each attribute into key and value
        if ' ' in attr:
            key, value = attr.split(' ', 1)
            # Remove quotes from the value
            value = value.strip('"')
            attributes[key] = value

    # Extract specific attributes
    reference_id = attributes.get('reference_id')
    ref_gene_name = attributes.get('ref_gene_name')
    TPM = attributes.get('TPM')

    # Return the extracted information
    return (reference_id,chromosome,start,end,strand,ref_gene_name,TPM)

import csv
def getTPMdict(file_path):

    infodict = {}
    with open(file_path, 'r', encoding='utf-8') as file:

        for line in file:
            # Split the line into fields
            fields = line.strip().split('\t')
            print(fields)
            if len(fields) > 2:
                if fields[2] == 'transcript':
                    result = parse_transcript_line(fields)
                    gene = result[5]
                    infodict[gene] = result

    return infodict

def toDf(annodict,tpmdict):

    data = []
    for key in annodict:

        if key is None or len(key) < 2:
            continue

        sl = []
        v = annodict[key]
        counter, aflist = v
        ID, chrom, txstart, txend,strand,genename,tpm = tpmdict[key]
        sl.append(genename)
        sl.append(tpm)
        sl.append(chrom)
        sl.append(strand)
        sl.append(txstart)
        sl.append(txend)
        sl.extend(counter)
        data.append(sl)

    df = pd.DataFrame(data, columns=['gene', 'TPM', 'chrom', 'strand', 'start', 'end','m6A','Y','Inosine','m5C'])
    return df


def output(writer,annodict,tpmdict):

    summaryDf  = toDf(annodict,tpmdict)
    summaryDf = summaryDf.sort_values(by='m6A', ascending=False)
    summaryDf.to_excel(writer, sheet_name='summary',  index_label='Row Label')



def summary(vcf,gtf,out):

    tpmdict = getTPMdict(gtf)
    writer = pd.ExcelWriter(out, engine='openpyxl')
    mod_kinds = ["a","17802","17596","m"]
    annodict = {}

    with open(vcf, 'r') as read_file:
        cnt = 0
        for line in read_file:

            print(line)
            if cnt ==0:
                cnt+=1
                continue

            fields = line.strip().split('\t')
            if fields[6] == "FAIL" or fields[6] == "False":
                continue
            modkind = fields[4]
            modidx = mod_kinds.index(modkind)
            strand,af,id,genename = parseResult(fields[7])
            print(strand,id,genename)
            if genename in tpmdict:

                if genename in annodict:

                    counter, aflocallist = annodict[genename]
                    annodict[genename] = counter, aflocallist
                    aflocallist.append((af, modidx))

                else:

                    aflocallist = []
                    counter = np.zeros(4)
                    aflocallist.append((af, modidx))
                    annodict[genename] = counter, aflocallist

                counter[modidx] += 1

    #output
    output(writer,annodict,tpmdict)
    writer.close()

def run():


    vcf = "/mnt/share/ueda/RNA004/hek293/result_filter.vcf"

    gtf = "/mnt/ssdnas07/nanozero/rna/nanomoditune_v01/HEK293T_DR13/HEK293T_DR13/HEK293T_DR13.gtf"
    genetsv = "/mnt/share/ueda/Docker/source/gencode_v44.tsv"
    out = "/mnt/share/ueda/RNA004/hek293/hek293_summary_all.xlsx"
    summary(vcf,gtf,out)

run()