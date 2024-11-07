import pysam

def getReferenceLength(read):

    refconsume = 0
    for cigaroprator, cigarlen in read.cigar:
        if cigaroprator == 0 or cigaroprator == 2: # match or del
            refconsume = refconsume + cigarlen
    return refconsume

# bamfile_name=bamfile_name, annotator=annotator,
#                           params=params, p_dict=p_dict, record=record,lowthres=lowthres

def takeSeq(neighborseqOrg, neighbor_seq):

    neighborseqOrg =  toSeqSeq(neighborseqOrg)
    neighbor_seq = toSeqSeq(neighbor_seq)
    if len(neighborseqOrg) == len(neighbor_seq):
        return neighbor_seq
    return neighborseqOrg

def toSeqSeq(seq):

    if seq is None:
        return ""
    if isinstance(seq, tuple):
        return seq[0]

    return str(seq)

def getNeighborSeq(record,chrom,gpos,strand,marjin = 20):

    # smer = fastafile.fetch(chrom, gpos - 1-marjin, gpos+marjin).upper()
    smer = str(record.seq[gpos - 1-marjin:gpos+marjin]).upper()
    if not strand:
        smer = reverse_complement(smer)
    post5 = smer[marjin:marjin+5]

    return smer,post5

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
def reverse_complement(seq):

    reverse_complement_seq = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement_seq

def convertToGenomepos(x,refposs):

    if x < 0:
        return 0
    if x >= len(refposs):
        return 0
    conv =  refposs[x]
    if conv is not None:
        conv = conv+1
    return conv

def initFilterFail(depth, ratio, modCount,params):

    mindepth = int(params.get('depth_min',10))
    ratio_min = float(params.get('ratio_min', 0.01))
    min_support_read = int(params.get('support_read_min', 4))

    if (depth < mindepth) or (ratio < ratio_min) or (modCount < min_support_read):
        return True

    return False

def count_intervals_simple(pthres,data, interval_ranges = [0, 63, 127, 191, 255]):

    interval_counts = [0] * (len(interval_ranges) - 1)
    cnt = 0
    for number in data:

        if number/255 >= pthres:
           cnt+=1

        for i in range(len(interval_ranges) - 1):
            if interval_ranges[i] <= number <= interval_ranges[i + 1]:
                interval_counts[i] += 1
                break
    return cnt,interval_counts

from collections import Counter
def most_frequent(list):
    counts = Counter(list)
    most_common_element = counts.most_common(1)[0][0]
    return most_common_element

def getSixmer(record,chrom,gpos,strand):

    if strand:
        # smer = fastafile.fetch(chrom,gpos-4,gpos+2).upper()
        smer = str(record.seq[gpos-4:gpos+2]).upper()
        mid5strand = smer[:5]
    else:
        # smer = fastafile.fetch(chrom, gpos - 3, gpos + 3).upper()
        smer = str(record.seq[gpos-4:gpos+2]).upper()
        mid5strand = smer[:5]
        smer = str(record.seq[gpos - 3:gpos + 3]).upper()
        smer = reverse_complement(smer)
    return smer,mid5strand


def calcP_null(q_list,pthres):

    qcnt = 0
    qsum = 0
    for q in q_list:
        qq = q/256
        if qq > pthres:
            qsum+=qq
            qcnt+=1
    if qcnt == 0:
        return 0.01
    qualmean = qsum/qcnt
    return 1- qualmean

def countHighMod(q_list,modqualthres=128):

    cnt=0
    for q in q_list:
        if q>=modqualthres:
            cnt+=1
    return cnt

from scipy.stats import binomtest
import math
def judgePval(params, depth,  q_list, pthres):

    p_null = float(params.get('p_null', 0.03))
    highmodCnt = countHighMod(q_list)
    if pthres > 0.4:
        highmodCnt = countHighMod(q_list,192)

    p_value = 0
    try:
        result = binomtest(highmodCnt, depth, p_null, alternative='greater')
        p_value = result.pvalue

    except ValueError as e:
        pass

    if p_value > 0:
        score = - math.log(p_value)
    elif p_value == 0:
        score = 1000

    print()
    return p_value,score


def refInConsistanse(modref,REF):
    REF = REF.upper()
    print(modref,REF)
    return modref != REF

def pileupMod(readlist,chrom,strand,start, end,record,annotator,params,p_dict):

    result = []
    posmap = {}

    for read, refposs in readlist:


        if read.has_tag("MM"):

            modbase = read.modified_bases
            if modbase is not None:

                modkeys = modbase.keys()
                mdict = {}
                for modkey in modkeys:
                    mdict[str(modkey[2])] = modkey[0]
                    modlist = modbase[modkey]
                    processed_tuples = [(convertToGenomepos(x, refposs), y) for x, y in modlist]
                    for tp in processed_tuples:
                        p,q = tp
                        if p is None:
                            continue
                        if q > 0:
                            key = str(modkey[2])
                            pdict = posmap.get(key,{})
                            posmap[key] = pdict
                            plist = pdict.get(p,[])
                            plist.append(q)
                            pdict[p] = plist
    depthCounter = {}
    for read, refposs in readlist:
        for genomepos in refposs:
            if genomepos is not None and genomepos >0:
                if genomepos in depthCounter:
                    depthCounter[genomepos] =  depthCounter[genomepos]+1
                else:
                    depthCounter[genomepos] = 1


    modkeys = sorted(posmap.keys())
    for mkey in modkeys:

        pthres = p_dict[mkey]
        for gpos in range(start, end):

            pdict = posmap[mkey]
            if gpos  not in pdict:
                continue

            q_list = pdict[gpos]

            modCnt,counts = count_intervals_simple(pthres,q_list)
            if gpos in depthCounter:
               depth = depthCounter[gpos]
            else:
               depth = len(q_list)

            ratio = 0
            if depth > 0:
                ratio = modCnt / depth
            initfilterFail = initFilterFail(depth, ratio,modCnt,params)
            if initfilterFail:
                continue
            fmer, mid5strand = getSixmer(record, chrom, gpos, strand)
            REF = fmer[3:4]
            if (mkey in mdict) and refInConsistanse(mdict[mkey],REF):
                continue

            # binomial test
            p_value,score = judgePval(params, depth,  q_list,pthres)
            statsTestOK = (p_value<0.01) or (p_value==0)
            neighborseq = None
            annoret = None
            if statsTestOK:

                annoret = annotator.annotate_genomic_position(chrom, strand, gpos)
                neighborseq = getNeighborSeq(record, chrom, gpos, strand)


            # depth, ratio, mostMismatchBase, mostMismatchCnt, diffcnt, base_counts, maints = ret
            formatted_af = f"{ratio:.4f}"
            info = "STRAND="+str(strand) +",AF=" + formatted_af + ",DP=" + str(depth) \
                   + ",MOD_Count=" + str(modCnt) +",MOD_Count_highbyQV="+str(counts)\
                   + ",PVAL=" + f"{p_value:.4f}"

            if neighborseq is not None:
                info = info + ",SEQ=" + neighborseq[0]
            if annoret is not None:
                info = info +"," + annoret


            ALT = mkey
            tp = (chrom, gpos, ".", REF, ALT, score, statsTestOK, info)
            print(tp)
            result.append(tp)

    return result

def pileup(interval,strand,bamfile_name,annotator,params,p_dict,record):

    chrom,start,end = interval.chrom,interval.start,interval.end
    # print(strand)
    strand = (strand == "p")

    #start pileup reads
    bamfile = pysam.AlignmentFile(bamfile_name, "rb")
    readlist = []
    for read in bamfile.fetch(chrom, start, end):

        if strand != (not read.is_reverse):
            continue

        refposs = read.get_reference_positions(full_length=True)
        readlist.append((read,refposs))

    result = pileupMod(readlist,chrom,strand,start, end,record,annotator,params,p_dict)
    bamfile.close()
    return result

from multiprocessing import Pool
from functools import partial
from annotation.GeneAnnotator import GenomeAnnotator
import pileup.PUtils as putils
import os
from Bio import SeqIO
def pileup_all(yamlf,recalib_stats,bamfile_name,outdir,ref, gtf_file,  stringtie_gtf,ncore=8):

    params, p_dict = loadPropFiles(yamlf,recalib_stats)
    print("threshold for mod call",p_dict)

    if not os.path.exists(bamfile_name):
        print("Could not find bam file",bamfile_name)
        sys.exit(1)
    if not os.path.exists(ref):
        print("Could not find ref file",ref)
        sys.exit(1)
    if not os.path.exists(recalib_stats):
        print("Could not find recalib stats file",recalib_stats)
        sys.exit(1)
    if not os.path.exists(recalib_stats):
        print("Could not find gtf_file file",gtf_file)
        sys.exit(1)


    if not os.path.exists(outdir):
        os.makedirs(outdir)

    print("load annotator")
    annotator = GenomeAnnotator(ref, gtf_file, stringtie_gtf,neibormargin=20)
    print("get interval")
    seq_index = SeqIO.index(ref, "fasta")

    autosomesOnly = params.get('autosomesOnly',True)
    vcfout = outdir + "/unfilter_result.vcf"

    if (os.path.isfile(vcfout)):
        os.remove(vcfout)

    p = Pool(ncore)
    intervalsByKey = annotator.getIntervalsByKey()

    for chrkey in intervalsByKey:

        chrom,strand = chrkey.split("-")
        print((chrom,strand))
        if autosomesOnly:
            if not putils.is_autosome_or_sex_chromosome(chrom):
                continue
        #
        if chrom in seq_index:
            record = seq_index[chrom]

        ivlist = intervalsByKey[chrkey]

        _pileup = partial(pileup,strand=strand,bamfile_name=bamfile_name, annotator=annotator,
                          params=params, p_dict=p_dict, record=record)

        retlist = p.map(_pileup,ivlist)

        outf = open(vcfout, 'a')
        for tp in retlist:

            result = tp
            for ret in result:
                line = '\t'.join(map(str, ret))
                print(line)
                outf.write(line + '\n')


        outf.close()




import csv
def csv_to_dict(filename, key_col, value_col,factor=1):

    dictionary = {}
    with open(filename, 'r') as file:
        reader = csv.reader(file)
        next(reader)  #skip title
        for row in reader:
            if row:
                key = row[key_col]
                value = float(row[value_col])*factor
                dictionary[key] = value
    return dictionary

def csv_to_dict_CplusDel(filename, factor=3):

    dictionary = {}
    with open(filename, 'r') as file:
        reader = csv.reader(file)
        next(reader)  #skip title
        for row in reader:
            if row:
                key = row[0]
                value = (float(row[4])+float(row[6]))*factor
                dictionary[key] = value
    return dictionary


import numpy as np
def calcparcentile(counter,parcentile):

    sum_p = np.sum(counter) / parcentile
    esum = 0
    for n in range(len(counter)):

        esum+=counter[n]
        if esum > sum_p:
            return (n / 256)
    return 0

import yaml
def loadPropFiles(yamlf,statfile):

    with open(yamlf, 'r') as file:
        params = yaml.safe_load(file)
    #
    pt = params.get("lowqual_parcentile",10)
    thres = params.get("min_mod_qual",0.4)

    p_dict = {}
    with open(statfile, 'r') as file:
        for line in file:
            line = line.split(",")
            key = line[0]
            counter = line[1:]
            counter = list(map(float, counter))
            pat = calcparcentile(counter, pt)
            p = thres
            if pat > thres:
                p = pat
            p_dict[key] = p

    return params,p_dict


# print("start")
ref = "/mnt/ssdnas07/pipeline/rna_v08/source/mm10.fa"
gtf_file = '/mnt/ssdnas07/pipeline/rna_v08/source/gencode.vM25.annotation.gff3'

# pileup_all(bamfile_name,bamfile_out,convertMatrix,out,gtf_file, ref,tpm,1)

def run():

    yamlf="/share/trna/project/nanoModiTune/nanoModiTune.yaml"
    recalib_stats = "/mnt/share/ueda/RNA004/Dorado0.8/bamout/stats/Adipocyte_2out_recalibstat.txt"
    bamfile_name = "/mnt/ssdnas07/nanozero/rna/nanomoditune_v01/Adipocyte_1/Adipocyte_1/Adipocyte_1_recalib.bam"
    stringtie_gtf = "/mnt/ssdnas07/nanozero/rna/nanomoditune_v01/Adipocyte_1/Adipocyte_1/Adipocyte_1.gtf"
    out = "/mnt/share/ueda/RNA004/Dorado0.8/bamout/test"
    ncore =12
    pileup_all(yamlf,recalib_stats,bamfile_name,out,ref,gtf_file,stringtie_gtf,ncore=ncore)


# run()