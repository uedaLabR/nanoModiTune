import glob
import csv
import numpy as np
import math
import pysam
from numba import jit

def convert_to_cumulative_sum(data):
    """
    Convert the list of values into a cumulative sum from right to left.

    Args:
        data (list of str): List of string numbers to be converted.

    Returns:
        np.array: Array of cumulative sums.
    """
    data = [float(value) for value in data]

    cumulative_sum = []
    current_sum = 0
    for value in reversed(data):
        current_sum += value
        cumulative_sum.append(current_sum)

    return np.array(list(reversed(cumulative_sum)))


class Recalib:

    def loadStats(self, stats_dir):
        """
        Load modification statistics from CSV files and build a lookup table.

        Args:
            stats_dir (str): Directory path where the CSV files are located.
        """
        self.lookuptable = {}
        files = glob.glob(f"{stats_dir}/*.csv")

        for file in files:
            print("loading",file)
            modkey = ""
            if "m6A" in file:
                modkey = "a"
            if "m5C" in file:
                modkey = "m"
            if "Y" in file:
                modkey = "17802"
            if "Inosine" in file:
                modkey = "17596"

            with open(file, 'r') as f:
                reader = csv.reader(f)
                for row in reader:
                    key = row[0]
                    vals = row[1:]
                    vals = convert_to_cumulative_sum(vals)
                    pvalIdx = vals / vals[0]
                    kkey = modkey + "_" + key
                    self.lookuptable[kkey] = pvalIdx

    def getRecalibScore(self, pvaltable, originalscore):
        """
        Calculate the recalibration score based on the original score.

        Args:
            pvaltable (np.array): Array of p-values for recalibration.
            originalscore (int): Original score to be recalibrated.

        Returns:
            int: The recalibrated score.
        """
        if originalscore == 0:
            return 0

        pval = pvaltable[originalscore]
        # print(pvaltable)
        if pval ==0:
            return originalscore
        logP = -10 * math.log10(pval)
        recalib = int(((logP - 3) / 27) * 255)
        if recalib < 0:
            recalib = 0
        # Scale the score to a 0-255 range
        # print(originalscore, pval, logP,recalib)
        if originalscore < recalib:
            return originalscore
        # print("returnning recalib score",originalscore,recalib)
        return recalib

    def getModifiedScore(self, kmer,refnuc, mod, originalscore):

        if kmer[3]!=refnuc:
            return -1
        kkey = str(mod) + "_" + kmer
        if kkey in self.lookuptable:
            return self.getRecalibScore(self.lookuptable[kkey], originalscore)
        return originalscore


def convertToGenomepos(x,refposs):

    if x < 0:
        return 0
    if x >= len(refposs):
        return 0
    conv =  refposs[x]
    if conv is not None:
        conv = conv+1
    return conv

def sorfby(modkeys):

    lst = []
    lst.extend(modkeys)

    desired_order_first = ['A', 'T', 'C']
    desired_order_third = [17596, 'a', 17802, 'm']

    order_dict_first = {value: index for index, value in enumerate(desired_order_first)}
    order_dict_third = {value: index for index, value in enumerate(desired_order_third)}
    sorted_data = sorted(lst, key=lambda x: (
    order_dict_first.get(x[0], float('inf')), order_dict_third.get(x[2], float('inf'))))

    return sorted_data

import csv
def write_string_and_data_to_csv(file_path,data):

    # CSV
    with open(file_path, mode='w', newline='') as file:

        writer = csv.writer(file)
        for row_data in data:
            writer.writerow(row_data)


from numba import jit

# @jit(nopython=True)
def revcon(seq):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G','N':'N',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}

    revcom_seq = []

    # Iterate through the sequence in reverse
    for i in range(len(seq) - 1, -1, -1):
        base = seq[i]
        revcom_seq.append(complement[base])

    return ''.join(revcom_seq)


splicethres = 50
def  checkSplice(refposs,localpos,strand):

    if strand:
        s = localpos - 4
        e = localpos + 2
    else:
        s = localpos - 3
        e = localpos + 3

    if s < 0 or e >= len(refposs):
        return False
    ss = refposs[s]
    ee = refposs[e]
    if ss == 0 or ee == 0:
        return False
    if ee-ss > splicethres:
        return True
    return False

def getSpliceSeq(genome,localpos,refposs,strand):


    seqlist = []
    if strand:
        s = localpos - 4
        e = localpos + 2
    else:
        s = localpos - 3
        e = localpos + 3

    seq_range = range(s, e)
    for x in seq_range:
        gpos = convertToGenomepos(x, refposs)
        if gpos > 0:
            seqlist.append(genome[gpos])
    seq = "".join(seqlist)
    if not strand:
        seq = revcon(seq)

    return seq


def getSixMer(genome,read,pos,localpos,refposs):

    if read.is_unmapped:
        return None

    reverse = read.is_reverse

    if reverse:

        smer = genome[pos-3:pos+3].upper()
        smer = revcon(smer)
        #check splice
        if checkSplice(refposs,localpos,False):
            smer = getSpliceSeq(genome,localpos,refposs,False)

    else:

        smer = genome[pos-4:pos+2].upper()
        #check splice
        if checkSplice(refposs, localpos, True):
            smer = getSpliceSeq(genome,localpos,refposs,True)

    return smer



import copy
def run_recalib(inbam, outbam, refs, recalib_db, out_stats):
    """
    Perform recalibration on the BAM file using reference and statistics.

    """
    print("start recalib",inbam, outbam, refs, recalib_db, out_stats)
    datadict = {}
    recalibrator = Recalib()
    recalibrator.loadStats(recalib_db)
    fasta = pysam.FastaFile(refs)
    lastref = ""
    with pysam.AlignmentFile(inbam, "rb") as bam_in, pysam.AlignmentFile(outbam, "wb", template=bam_in) as bam_out:


        readcnt = 0
        recalibbase = 0
        unref = 0
        unchange = 0
        for read in bam_in:

            if not read.is_mapped:
                bam_out.write(read)
                continue

            referencename = bam_in.get_reference_name(read.reference_id)
            # print(readcnt,referencename,read)
            if referencename != lastref:
                genome = fasta.fetch(referencename,0,None)
                print("load",referencename)
            lastref = referencename

            if read.has_tag("MM"):

                refposs = read.get_reference_positions()
                modbase = read.modified_bases
                ML = read.get_tag("ML")
                orginal_array = copy.deepcopy(ML)

                index = 0
                if modbase is not None:
                    modkeys = list(modbase.keys())
                    modkeys = sorfby(modkeys)
                    for modkey in modkeys:

                        modlist = modbase[modkey]
                        processed_tuples = [(convertToGenomepos(x, refposs),x, y) for x, y in modlist]
                        refnuc = modkey[0]

                        for tp in processed_tuples:

                            pos,localpos,originalscore = tp

                            sixMer = getSixMer(genome,read,pos,localpos,refposs)
                            # print(referencename,not read.is_reverse,pos,sixMer)

                            if len(sixMer)==6:
                                recalibscore = recalibrator.getModifiedScore(sixMer,refnuc, modkey[2], originalscore)
                                unrefB = False
                                if originalscore == ML[index]:
                                    #if not something wrong
                                    if originalscore != recalibscore:

                                        if recalibscore<0:
                                            unref+=1
                                            unrefB = True
                                            recalibscore=0
                                        else:
                                            recalibbase+=1

                                        ML[index] = recalibscore

                                    else:
                                        unchange+=1

                                    if unrefB == False:
                                        kkey = str(modkey[2])
                                        if kkey in datadict:
                                            counter = datadict[kkey]
                                        else:
                                            counter = np.zeros(256)

                                        #
                                        if recalibscore > 0:
                                            counter[recalibscore] = counter[recalibscore]+1
                                            datadict[kkey] = counter

                            index+=1

                read.set_tag("ML",ML)
                read.set_tag("XM", orginal_array)

            readcnt+=1
            if readcnt%1000==0:
                print(readcnt,recalibbase,unref,unchange)

            # Write the modified read to the output BAM file
            bam_out.write(read)

        print("readcnt, recalibbase, nonref, unchange")
        print(readcnt, recalibbase, unref, unchange)
        alldata = []
        for key in datadict:
            counter = datadict[key]
            data = [key] + list(counter)
            print(data)
            alldata.append(data)
            # factor = calculate_normalization_factor(lis)
            # # print(key,factor,lis)

        alldata.sort(key=lambda x: x[0])
        write_string_and_data_to_csv(out_stats, alldata)

# Example execution
refs = "/mnt/ssdnas07/pipeline/rna_v08/source/mm10.fa"
inbam = "/mnt/ssdnas07/nanozero/rna/nanomoditune_v01/Adipocyte_1/Adipocyte_1/Adipocyte_1_recalib.bam"
outbam = "/mnt/share/ueda/RNA004/Dorado0.8/bamout/Adipocyte_2out_recalib.bam"
recalib_db = "/mnt/share/ueda/RNA004/Dorado0.8/bamout/stats"
out_stats = "/mnt/share/ueda/RNA004/Dorado0.8/bamout/stats/Adipocyte_2out_recalibstat.txt"

from functools import partial
import cProfile

# wrapped_function = partial(run_recalib, inbam, outbam, refs, recalib_db, out_stats)
# wrapped_function()
# cProfile.run('wrapped_function()')
