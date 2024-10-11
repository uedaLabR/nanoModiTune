import glob
import csv
import numpy as np
import math
import pysam

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

import copy
def recalib(refs, stats, bamfile, bamout):
    """
    Perform recalibration on the BAM file using reference and statistics.

    Args:
        refs (str): Path to the reference file.
        stats (str): Directory path containing the recalibration stats.
        bamfile (str): Input BAM file path.
        bamout (str): Output recalibrated BAM file path.
    """
    recalibrator = Recalib()
    recalibrator.loadStats(stats)
    fasta = pysam.FastaFile(refs)
    lastref = ""
    with pysam.AlignmentFile(bamfile, "rb") as bam_in, pysam.AlignmentFile(bamout, "wb", template=bam_in) as bam_out:

        readcnt = 0
        recalibbase = 0
        unref = 0
        unchange = 0
        for read in bam_in:

            if not read.is_mapped:
                bam_out.write(read)
                continue

            referencename = bam_in.get_reference_name(read.reference_id)
            if referencename != lastref:
                genome = fasta.fetch(referencename,0,None)

            if read.has_tag("MM"):

                refposs = read.get_reference_positions()
                modbase = read.modified_bases
                ML = read.get_tag("ML")
                orginal_array = copy.deepcopy(ML)

                index = -1
                if modbase is not None:
                    modkeys = list(modbase.keys())
                    modkeys = sorfby(modkeys)
                    for modkey in modkeys:

                        modlist = modbase[modkey]
                        processed_tuples = [(convertToGenomepos(x, refposs), y) for x, y in modlist]
                        refnuc = modkey[0]

                        for tp in processed_tuples:

                            recalibscore = 0
                            pos, originalscore = tp
                            sixMer = genome[pos-4:pos+2]
                            if len(sixMer)==6:
                                recalibscore = recalibrator.getModifiedScore(sixMer,refnuc, modkey[2], originalscore)

                            if originalscore == ML[index]:
                                #if not something wrong
                                if originalscore != recalibscore:

                                    if recalibscore<0:
                                        unref+=1
                                        recalibscore=0
                                    else:
                                        recalibbase+=1

                                    ML[index] = recalibscore

                                else:
                                    unchange+=1

                            index+=1

                read.set_tag("ML",ML)
                read.set_tag("XM", orginal_array)

            readcnt+=1
            if readcnt%10000==0:
                print(readcnt,recalibbase,unref,unchange)
            # Write the modified read to the output BAM file
            bam_out.write(read)

# Example execution
refs = "/share/reference/IVTmix.fa"
bamfile = "/mnt/share/ueda/RNA004/Dorado0.8/bamout/BC3.bam"
bamout = "/mnt/share/ueda/RNA004/Dorado0.8/bamout/BC3recalib_mod.bam"
stats = "/mnt/share/ueda/RNA004/Dorado0.8/bamout/stats"

recalib(refs, stats, bamfile, bamout)