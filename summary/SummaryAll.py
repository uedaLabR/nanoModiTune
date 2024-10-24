import numpy as np



def parseResult(infostr):

    asr =  False
    if "ASR=True" in infostr:
        asr = True


    annostr = ""
    print(infostr)
    if "Anno=(" in infostr:
        annostr = infostr[infostr.index("Anno=(")+6:infostr.index(")")]
        print(annostr)
    #
    strs = infostr.split(",")
    #Anno = ('CDS', 'ENST00000521466.5', 'RPS14', -873)
    maints = ""
    af = 0
    dp = 0
    for astr in strs:
        # if "Anno=" in annostr:
        #
        #     annostr = annostr.replace("Anno = (","")
        #     print(annostr)
        #     if "None" in annostr:
        #         return None
        #     annostr = annostr.replace(")","")
        if "MainTs=" in astr:

            maints = astr.replace("MainTs=","")
        if "AF=" in astr:

            afs = astr.replace("AF=","")
            af = float(afs)

        if "MOD_Count=" in astr:

            hmc = astr.replace("MOD_Count=","")
            hmc = int(hmc)

        if "DP=" in astr:

            dps = astr.replace("DP=","")
            dp = int(dps)

    print("annostr",annostr)
    annos = annostr.split(",")
    if len(annos) < 2:
        return None
    flg = annos[0].replace("'", "")
    id = annos[1].replace("'", "")
    gene = annos[2].replace("'", "")
    seq = annos[4].replace("'", "")

    return id,gene, flg, asr, seq, maints, af,dp, hmc



def create_histogram(data, bin_size):

    # print(data)
    num_bins = int(1 / bin_size)
    # Create the bins and initialize the histogram
    histogram = {f"{bin_size * i:.1f}-{bin_size * (i + 1):.1f}": 0 for i in range(num_bins)}

    # Bin the data
    for number in data:
        if number > 0.99:
            number = 0.99
        bin_index = int(number / bin_size)
        bin_key = f"{bin_size * bin_index:.1f}-{bin_size * (bin_index + 1):.1f}"
        histogram[bin_key] += 1

    # print(histogram)
    return histogram

def sep(aflist):

    aflist1=[]
    aflist2 = []
    aflist3 = []

    for af,flg in aflist:

        if flg ==0:
            aflist1.append(af)
        elif flg ==1:
            aflist2.append(af)
        elif flg ==2:
            aflist3.append(af)

    return aflist1,aflist2,aflist3

import statistics
def toDf(annodict,tpmdict,idtogenesymbol):

    data = []
    datagenomic = []

    cnt = 0
    tpmwritten = set()
    counterlen = 0
    for key in annodict:
        print(cnt)
        v = annodict[key]
        counter,aflist = v
        sl =[]
        gene = key.split(":")[0]
        id = key.split(":")[1]


        if id =="ge":
            id = "genomic"


        sl.append(id)
        sl.append(gene)

        # ID\tTPM\tChrom\tSTRAND\ttxStart\ttxEnd\tGENELEN"
        tpm,chrom,strand,txstart,txend,genelen = 0,"","",0,0,0
        regtranscript = False
        if id in tpmdict:
            tpm,chrom,strand,txstart,txend,genelen = tpmdict[id]
            tpmwritten.add(id)
            regtranscript = True

        sl.append(tpm)
        sl.append(chrom)
        sl.append(strand)
        sl.append(txstart)
        sl.append(txend)
        sl.append(genelen)
        sl.append(sum(counter))
        aflist1,aflist2,aflist3 = sep(aflist)
        sl.append(len(aflist1))
        sl.append(statistics.median(aflist1) if aflist1 else 0)
        sl.append(max(aflist1) if aflist1 else 0)
        sl.append(sum(aflist1) if aflist1 else 0)
        sl.append(len(aflist2))
        sl.append(statistics.median(aflist2) if aflist2 else 0)
        sl.append(max(aflist2) if aflist2 else 0)
        sl.append(sum(aflist2) if aflist2 else 0)
        sl.append(len(aflist3))
        sl.append(statistics.median(aflist3)  if aflist3 else 0)
        sl.append(max(aflist3) if aflist3 else 0)
        sl.append(sum(aflist3) if aflist3 else 0)

        counterlen = len(counter)
        sl.extend(counter)
        if regtranscript:
            data.append(sl)
        else:

            sl2 = []
            sl2.append(id)
            sl2.append(gene)
            sl2.append(sum(counter))
            sl2.extend(counter)
            datagenomic.append(sl2)

        cnt+=1

    #write other transcript
    idkeys = tpmdict.keys()
    for id in idkeys:

        print(cnt)
        if id == "ge":
            continue
        if id not in tpmwritten:

            tpm,chrom,strand, txstart, txend, genelen = tpmdict[id]
            if float(tpm) == 0:
                continue
            sl = []
            sl.append(id)
            gene =""
            if id in idtogenesymbol:
                gene = idtogenesymbol[id]
            sl.append(gene)
            sl.append(tpm)
            sl.append(chrom)
            sl.append(strand)
            sl.append(txstart)
            sl.append(txend)
            sl.append(genelen)

            zero_list = [0] * (counterlen+13)
            sl.extend(zero_list)
            data.append(sl)
            cnt += 1

    # order = ['CDS', 'INTRON', '3pUTR', '5pUTR', 'upstream', 'downstream', 'noncoding', 'novel_ncRNA', 'genomic',
    #          'lookuperror']
    df = pd.DataFrame(data, columns=['ID','gene','TPM','Chrom','strand','txStart','txEnd','GENE_LEN',
                                     'TotalModPos',
                                     'm6A_counts','m6A_AFmedian','m6A_AFmax','m6A_AFsum',
                                     'Y_counts','Y_AFmedian','Y_AFmax','Y_AFsum',
                                     'm5C_counts', 'm5C_AFmedian', 'm5C_AFmax', 'm5C_AFsum',

                                     'CDS_m6A', 'INTRON_m6A', '3pUTR_m6A', '5pUTR_m6A',
                                     'upstream_m6A', 'downstream_m6A','noncoding_m6A', 'novel_transcript_m6A',
                                     'inter_genic_m6A','CDS_Y', 'INTRON_Y', '3pUTR_Y', '5pUTR_Y', 'upstream_Y',
                                     'downstream_Y','noncoding_Y', 'novel_transcript_Y', 'inter_genic_Y',
                                     'CDS_m5C', 'INTRON_m5C', '3pUTR_m5C', '5pUTR_m5C', 'upstream_m5C',
                                     'downstream_m5C','noncoding_m5C', 'novel_transcript_m5C', 'inter_genic_m5C'])
    dfgenomic = pd.DataFrame(datagenomic, columns=['ID', 'neighbor_gene', 'TotalModPos', 'CDS_m6A', 'INTRON_m6A', '3pUTR_m6A', '5pUTR_m6A',
                                     'upstream_m6A', 'downstream_m6A', 'noncoding_m6A', 'novel_transcript_m6A',
                                     'inter_genic_m6A', 'CDS_Y', 'INTRON_Y', '3pUTR_Y', '5pUTR_Y', 'upstream_Y',
                                     'downstream_Y', 'noncoding_Y', 'novel_transcript_Y', 'inter_genic_Y', 'CDS_m5C',
                                     'INTRON_m5C', '3pUTR_m5C', '5pUTR_m5C', 'upstream_m5C', 'downstream_m5C',
                                     'noncoding_m5C', 'novel_transcript_m5C', 'inter_genic_m5C'])
    return df,dfgenomic

def toDf2(rdict):

    data = []
    for key in rdict:
        v = rdict[key]
        sl =[]
        sl.append(key)
        sl.append(v)
        data.append(sl)

    df = pd.DataFrame(data, columns=['AF_range', 'count'])
    return df

def toSummaryDf(array):

    row1 = array[:9]
    row2 = array[9:9*2]
    row3 = array[9*2:]

    data = []
    data.append(row1)
    data.append(row2)
    data.append(row3)
    df = pd.DataFrame(data, columns=['CDS', 'INTRON', '3pUTR', '5pUTR', 'upstream', 'downstream','noncoding', 'novel_transcript', 'inter_genic'])
    df.index =  [ 'm6A','Y','m5C']
    return df

def output(writer,annodict,af_histogram1,af_histogram2,af_histogram3,tpmdict,allcounter,idtogenesymbol):
#def output(writer,annodict,annodictAsr,af_histogram,seqlist,tpmdict,allcounter,allcounterSpr):

    print("finish 0")
    summaryDf = toSummaryDf(allcounter)
    print("finish 1")
    annoDf,genomicDf = toDf(annodict,tpmdict,idtogenesymbol)
    print("finish 2")
    dfhist1 = toDf2(af_histogram1)
    print("finish 3")
    dfhist2 = toDf2(af_histogram2)
    print("finish 4")
    dfhist3 = toDf2(af_histogram3)
    print("finish 5")

    summaryDf.to_excel(writer, sheet_name='summary',  index_label='Row Label')
    annoDf.to_excel(writer, sheet_name='All', index=False)
    # genomicDf.to_excel(writer, sheet_name='nongene_region', index=False)
    dfhist1.to_excel(writer, sheet_name='AF_m6A', index=False)
    dfhist2.to_excel(writer, sheet_name='AF_Y', index=False)
    dfhist3.to_excel(writer, sheet_name='AF_m5C', index=False)

import csv
def getTPMdict(file_path):
    infodict = {}
    with open(file_path, 'r', encoding='utf-8') as file:
        tsv_reader = csv.reader(file, delimiter='\t')
        for row in tsv_reader:
            id = row[0]
            if "#" in id:
                continue
            # ID\tTPM\tChrom\tSTRAND\ttxStart\ttxEnd\tGENELEN"
            infodict[id] = (row[1],row[2],row[3],row[4],row[5],row[6])
    return infodict

def getidtogenesymbol(genetsv):

    d = {}
    with open(genetsv, 'r') as file:
        tsv_reader = csv.reader(file, delimiter='\t')
        ididx = 0
        gsidx = 16
        isFirst = True
        target_title = "hg38.kgXref.geneSymbol"
        for row in tsv_reader:

            if isFirst:
                try:
                    gsidx = row.index(target_title)
                except ValueError:
                    pass
                isFirst = False
                continue

            id = row[ididx]
            genesymbol = row[gsidx]
            d[id] = genesymbol
    print("d",d)
    return d

order = ['CDS', 'INTRON', '3pUTR', '5pUTR', 'upstream', 'downstream','noncoding','novel_ncRNA','genomic']
import pandas as pd
def summary(vcf,tpm,genetsv,out):

    tpmdict = getTPMdict(tpm)
    idtogenesymbol = getidtogenesymbol(genetsv)
    writer = pd.ExcelWriter(out, engine='openpyxl')
    annodict = {}
    allcounter = np.zeros(27)
    aflist1 =[]
    aflist2 = []
    aflist3 = []
    mod_kinds = ["m6A","Y","m5C"]

    with open(vcf, 'r') as read_file:
        cnt = 0
        for line in read_file:

            if cnt ==0:
                cnt+=1
                continue

            fields = line.strip().split('\t')
            if fields[6] == "FAIL" or fields[6] == "False":
                continue

            modkind = fields[4]
            modidx = mod_kinds.index(modkind)
            print("modidx",modidx,line)
            #
            parsedresult = parseResult(fields[7])

            if parsedresult is not None:
                id,gene, flg, asr, seq, maints, af,dp, hmc =  parsedresult

                if id is None:
                    continue
                keygene = gene+":"+maints
                idx = order.index(flg)
                idx = (modidx *9)+ idx
                if modidx == 0:
                    aflist1.append(af)
                elif modidx == 1:
                    aflist2.append(af)
                elif modidx == 2:
                    aflist3.append(af)

                if keygene in annodict:

                    counter,aflocallist = annodict[keygene]
                    annodict[keygene] = counter, aflocallist
                    aflocallist.append((af,modidx))

                else:

                    aflocallist = []
                    counter = np.zeros(27)
                    aflocallist.append((af,modidx))
                    annodict[keygene] = counter,aflocallist

                counter[idx] += 1
                allcounter[idx] +=1



    bin_size = 0.1  # Bin size of 0.1
    af_histogram1 = create_histogram(aflist1, bin_size)
    af_histogram2 = create_histogram(aflist2, bin_size)
    af_histogram3 = create_histogram(aflist3, bin_size)

    #output
    output(writer,annodict,af_histogram1,af_histogram2,af_histogram3,tpmdict,allcounter,idtogenesymbol)
    writer.close()

def run():

    # vcf = "/mnt/share/ueda/RNA004/U87/U87_IVT/20231227_1535_MN32625_FAX73794_2cf3868f/result_filter.vcf"
    # tpm = "/mnt/share/ueda/RNA004/U87/U87_WT/20231215_1505_MN32625_FAX70236_31c89747/tpm.txt"
    # genetsv = "/mnt/share/ueda/Docker/source/gencode_v44.tsv"
    # out = "/mnt/share/ueda/RNA004/U87/U87_IVT/20231227_1535_MN32625_FAX73794_2cf3868f/U87ivt_summary_all.xlsx"
    # summary(vcf,tpm,genetsv,out)
    #
    # vcf = "/mnt/share/ueda/RNA004/U87/U87_WT/20231215_1505_MN32625_FAX70236_31c89747/result_filter.vcf"
    # tpm = "/mnt/share/ueda/RNA004/U87/U87_WT/20231215_1505_MN32625_FAX70236_31c89747/tpm.txt"
    # genetsv = "/mnt/share/ueda/Docker/source/gencode_v44.tsv"
    # out = "/mnt/share/ueda/RNA004/U87/U87_WT/20231215_1505_MN32625_FAX70236_31c89747/U87_summary_all.xlsx"
    # summary(vcf,tpm,genetsv,out)

    vcf = "/mnt/share/ueda/RNA004/U87/U87_STM2457/U87_STM2457/20240216_1602_MN32625_FAX70218_f0b7f679/result_filter.vcf"
    tpm = "/mnt/share/ueda/RNA004/U87/U87_STM2457/U87_STM2457/20240216_1602_MN32625_FAX70218_f0b7f679/tpm.txt"
    genetsv = "/mnt/share/ueda/Docker/source/gencode_v44.tsv"
    out = "/mnt/share/ueda/RNA004/U87/U87_STM2457/U87_STM2457/20240216_1602_MN32625_FAX70218_f0b7f679/U87_STM2457_summary_all.xlsx"
    summary(vcf,tpm,genetsv,out)

# run()