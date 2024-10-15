import mappy as mp
import pysam
import pandas as pd


def convertToGenomepos(x,refposs):

    if x < 0:
        return 0
    if x >= len(refposs):
        return 0
    conv =  refposs[x]
    if conv is not None:
        conv = conv+1
    return conv

def _find_index(lst, value):
    try:
        return lst.index(value)
    except ValueError:
        return -1

def find_index(lst, value):

    id1 = _find_index(lst, value)
    id2 = _find_index(lst, value-1)
    if id2 == -1:
        return id1
    if (id1-id2) == 1:
        return id1
    if (id1-id2) > 1:#insersion
        return id2+1
    return id1



def getDepth(pos,readlist):

    depth = 0
    for read, refposs in readlist:

        idx = find_index(refposs, pos)
        if idx >= 0:
            depth +=1
    return depth

import numpy as np
def calculate_normalization_factor(data, top_percentage=3):

    # Convert the data to a NumPy array
    data = np.array(data)

    # Get the value at the specified top percentage
    threshold_index = int(len(data) * (1 - top_percentage / 100))
    top_value = np.partition(data, threshold_index)[threshold_index]
    # print("top_value",top_value)

    # If the top value is less than 128, return 1
    if top_value < 127:
        return 1

    # Calculate the scaling factor
    scale_factor = 127 / top_value
    return scale_factor

import csv
def write_string_and_data_to_csv(file_path,data):

    # CSV
    with open(file_path, mode='w', newline='') as file:

        writer = csv.writer(file)
        for row_data in data:
            writer.writerow(row_data)

def takeStats(ref,bam,out,modnuc,refnuc):

    datadict = {}
    a = mp.Aligner(ref)
    for seqname in a.seq_names:


        genome = a.seq(seqname, end=0x7fffffff)
        seqlen = len(genome)
        print(seqlen,seqname)
        cnt=0
        readlist = []
        # try:
        bamfile = pysam.AlignmentFile(bam, "rb")
        for read in bamfile.fetch(seqname, 50, seqlen-50):

            refposs = read.get_reference_positions(full_length=True)
            readlist.append((read,refposs))
            if cnt > 1000:
                break
            cnt+=1
        bamfile.close()

        posmap = {}
        for read,refposs in readlist:

            if read.has_tag("MM"):

                modbase = read.modified_bases
                if modbase is not None:

                    modkeys = list(modbase.keys())
                    for m6Akey in modkeys:

                        if m6Akey[2] == modnuc:
                            modlist = modbase[m6Akey]
                            processed_tuples = [(convertToGenomepos(x, refposs), y) for x, y in modlist]
                            for tp in processed_tuples:

                                p, q = tp
                                if p is not None and q is not None:
                                    if p > 50 and p < seqlen-50:

                                        if p is None:
                                            continue
                                        modcnt = posmap.get(p,[])
                                        modcnt.append(q)
                                        posmap[p] = modcnt

        poskeys = sorted(posmap.keys())
        for pos in poskeys:

            depth = getDepth(pos,readlist)
            fiveMer = genome[pos-4:pos+2]
            if len(fiveMer)==6:
                middle = fiveMer[3]
                if middle == refnuc:
                    quallist = posmap[pos]

                    len1 = len(quallist)
                    len2 = depth
                    diff = len2-len1


                    if fiveMer in datadict:
                        counter = datadict[fiveMer]
                    else:
                        counter = np.zeros(256)

                    for q in quallist:
                        counter[q] = counter[q]+1

                    counter[0] = counter[0]+diff

                    datadict[fiveMer] = counter
        # except:
        #     print("err",cnt)
        #     pass
        # break

    alldata = []
    for key in datadict:

        counter = datadict[key]
        data = [key] + list(counter)
        print(data)
        alldata.append(data)
        # factor = calculate_normalization_factor(lis)
        # # print(key,factor,lis)

    write_string_and_data_to_csv(out, alldata)



ref = "/share/reference/IVTmix.fa"
bam = "/mnt/share/ueda/RNA004/Dorado0.8/bamout/BC1.bam"
out = "/mnt/share/ueda/RNA004/Dorado0.8/bamout/stats/m6A_FP.csv"
takeStats(ref,bam,out,"a","A")

out = "/mnt/share/ueda/RNA004/Dorado0.8/bamout/stats/Inosine_FP.csv"
takeStats(ref,bam,out,17596,"A")

out = "/mnt/share/ueda/RNA004/Dorado0.8/bamout/stats/m5C_FP.csv"
takeStats(ref,bam,out,"m","C")


out = "/mnt/share/ueda/RNA004/Dorado0.8/bamout/stats/Y_FP.csv"
takeStats(ref,bam,out,17802,"T")