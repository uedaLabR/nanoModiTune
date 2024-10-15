import pysam

def getReferenceLength(read):

    refconsume = 0
    for cigaroprator, cigarlen in read.cigar:
        if cigaroprator == 0 or cigaroprator == 2: # match or del
            refconsume = refconsume + cigarlen
    return refconsume
def countGene(readlist,startEnd):

    index, chrom, strand, start, end = startEnd
    gekey = str(chrom)+":"+str(strand)+":"+str(start)+":"+str(end)
    geneCounter = {}
    for read, refposs in readlist:

        xstag = gekey
        if read.has_tag("XS"):
            xstag = read.get_tag("XS")
        # for TPM calclation
        refMaplen = getReferenceLength(read)
        if xstag not in geneCounter:
            geneCounter[xstag] = refMaplen
        else:
            geneCounter[xstag] = geneCounter[xstag] + refMaplen

    return geneCounter
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



from scipy.stats import binomtest
import math
def judgePval(params, depth, modCnt,  q_list,pthres):

    p_thres = float(params.get('p_val_thres', 0.01))
    p_null = calcP_null(q_list,pthres)
    p_value = 0
    try:
        result = binomtest(modCnt, depth, p_null, alternative='greater')
        p_value = result.pvalue

    except ValueError as e:
        pass

    if p_value > 0:
        score = - math.log(p_value)
    elif p_value == 0:
        score = 1000

    return p_null,p_value,score,(p_value <= p_thres)

def pileupMod(readlist,chrom,strand,start, end,record,annotator,params,p_dict,lowthres):

    result = []
    posmap = {}

    for read, refposs in readlist:

        xstag = "ge"
        if read.has_tag("XS"):
            xstag = read.get_tag("XS")

        if read.has_tag("MM"):

            modbase = read.modified_bases
            if modbase is not None:

                modkeys = modbase.keys()
                for modkey in modkeys:
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
                            plist.append((q,xstag))
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
        # print(p_dict,mkey)
        pthres = p_dict[mkey]
        for gpos in range(start, end):

            fmer, mid5strand = getSixmer(record, chrom, gpos, strand)
            refbase = fmer[2]

            pdict = posmap[mkey]
            if gpos  not in pdict:
                continue
            datalist = pdict[gpos]
            q_list = [x[0] for x in datalist]
            ts_list = [x[1] for x in datalist]


            modCnt,counts = count_intervals_simple(pthres,q_list)
            maints = most_frequent(ts_list)
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

            # binomial test
            p_null, p_value, score, judgeOK = judgePval(params, depth, modCnt,  q_list,pthres)
            # pass filter
            alternativeSplicingregion = False
            annotation = None
            neighborseqOrg = None
            if judgeOK:
                annoret = annotator.annotate_genomic_position(chrom, strand, gpos, maints,record)
                neighborseqOrg = getNeighborSeq(record, chrom, gpos, strand)
                # check motief and
                if annoret is not None:
                    neighborseq, alternativeSplicingregion, annotation = annoret
                    fmernew = neighborseq[17:23]
                    # reculculate if exon boundary
                    if fmer != fmernew:
                        fmer = fmernew

                    if annotation is not None:
                        (utrlabel1, gene_id, gene_symbol, dist, neighbor_seq, spliceJunction) = annotation
                        annotation = (utrlabel1, gene_id, gene_symbol, dist, spliceJunction)
                        neighborseqOrg = takeSeq(neighborseqOrg,neighbor_seq)


            # depth, ratio, mostMismatchBase, mostMismatchCnt, diffcnt, base_counts, maints = ret
            formatted_af = f"{ratio:.4f}"
            info = "STRAND="+str(strand) +",AF=" + formatted_af + ",DP=" + str(depth) \
                   + ",MOD_Count=" + str(modCnt) +",MOD_Count_byQV="+str(counts)\
                   + ",PVAL=" + f"{p_value:.4f}" + ",P_null=" + f"{p_null:.5f}" + ",MainTs=" + maints \
                   + ",ASR=" + str(alternativeSplicingregion) \
                   + ",SEQ=" + str(neighborseqOrg) + ", Anno=" + str(annotation)

            REF = fmer[3:4]
            ALT = mkey
            tp = (chrom, gpos, ".", REF, ALT, score, judgeOK, info)
            result.append(tp)

    return result

def pileup(startEnd,bamfile_name,annotator,params,p_dict,record,lowthres):

    index,chrom,strand,start,end = startEnd
    strand = (strand == "1")

    #start pileup reads
    bamfile = pysam.AlignmentFile(bamfile_name, "rb")
    readlist = []
    for read in bamfile.fetch(chrom, start, end):

        map_strand = True
        if read.is_reverse:
            map_strand = False
        if strand != map_strand:
            continue

        refposs = read.get_reference_positions(full_length=True)
        readlist.append((read,refposs))

    result = pileupMod(readlist,chrom,strand,start, end,record,annotator,params,p_dict,lowthres)
    geneCounter = countGene(readlist,startEnd)

    bamfile.close()
    return index,result,geneCounter

from multiprocessing import Pool
from functools import partial
from annotation.GeneAnnotation import GenomeAnnotator
import pileup.PUtils as putils
import os
from Bio import SeqIO
def pileup_all(yamlf,recalib_stats,bamfile_name,convertMatrix,outdir,gtf_file, ref,ncore=8,lowthres=True):

    params, p_dict = loadPropFiles(yamlf,recalib_stats)
    print(p_dict)

    if not os.path.exists(bamfile_name):
        print("Could not find file",bamfile_name)
        return

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    print("load annotator")
    annotator = GenomeAnnotator(gtf_file,convertMatrix, neibormargin=20)
    print("get interval")
    seq_index = SeqIO.index(ref, "fasta")

    autosomesOnly = params.get('autosomesOnly',True)
    vcfout = outdir + "/unfilter_result.vcf"

    if (os.path.isfile(vcfout)):
        os.remove(vcfout)

    tpm = outdir + "/tpm.txt"

    geneCounterAll = {}

    p = Pool(ncore)
    intervals = putils.getIntervals(convertMatrix)

    allout = []
    sampleout = []
    for intervalKey in intervals:

        data = intervalKey.split(":")
        chr = data[0]
        if autosomesOnly:
            if not putils.is_autosome_or_sex_chromosome(chr):
                continue

        strand = data[1]
        ivlist = intervals[intervalKey]

        if chr in seq_index:
            record = seq_index[chr]
        #
        ivlist = [((index,chr, strand) + x) for index, x in enumerate(ivlist)]
        #
        _pileup = partial(pileup, bamfile_name=bamfile_name, annotator=annotator,
                          params=params, p_dict=p_dict, record=record,lowthres=lowthres)

        retlist = p.map(_pileup,ivlist)
        retlist = sorted(retlist, key=lambda x: x[0])
        # result, bamTagInfo, geneCounter
        outf = open(vcfout, 'a')


        for tp in retlist:

            index,result,geneCounter = tp
            for ret in result:
                line = '\t'.join(map(str, ret))
                outf.write(line + '\n')


            for key, value in geneCounter.items():
                if key in geneCounterAll:
                    geneCounterAll[key] += value
                else:
                    geneCounterAll[key] = value


        outf.close()
        # break


    #calculate tpm
    putils.calcTpm(geneCounterAll, tpm, convertMatrix)
    #End
    # bintervallist = db.from_sequence(intervalsList)  # convert to dask bag

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


print("start")
ref = "/share/reference/hg38.fa"
gtf_file = '/mnt/share/ueda/RNA004/pipeline/gencode_v44.tsv'

# pileup_all(bamfile_name,bamfile_out,convertMatrix,out,gtf_file, ref,tpm,1)

def run():

    yamlf="/share/trna/project/nanoModiTune/nanoModiTune.yaml"
    recalib_stats = "/mnt/share/ueda/RNA004/Dorado0.8/bamout/stats/recalibstat.txt"
    bamfile_name = "/mnt/share/ueda/RNA004/nanoEvo/U87_inhibitor_bam/BC1.bam"
    convertMatrix = "/mnt/ssdnas/nanozero/rna/U87_inhibitors_DR14/U87_inhibitors_DR14/U87_inhibitors_DR14.matrix"
    out = "/mnt/share/ueda/RNA004/nanoEvo/test"
    ncore =12
    pileup_all(yamlf,recalib_stats,bamfile_name,convertMatrix,out,gtf_file, ref,ncore=ncore,lowthres=True)


# run()