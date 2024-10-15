def addOrExtend(il,start,end):

    def overlaps(interval):
        return not (interval[1] < start or interval[0] > end)

    for i, interval in enumerate(il):
        if overlaps(interval):
            new_start = min(interval[0], start)
            new_end = max(interval[1], end)
            il[i] = (new_start, new_end)
            return il

    il.append((start, end))
    return il

import csv
def getIntervals(file):

    intervals = {}
    print("start reading file")
    with open(file, encoding='utf-8', newline='') as f:

        for cols in csv.reader(f, delimiter='\t'):

            if cols[0].startswith("#"):
                continue

            chr = cols[1]
            strand = cols[2]
            poss = cols[3].split(",")
            start = int(poss[0].split(":")[0])
            end = int(poss[-1].split(":")[1])

            key = chr + ":" + str(strand)

            if key in intervals:
                il = intervals[key]
            else:
                il = []

            il = addOrExtend(il,start,end)
            intervals[key] = il

    return intervals

def getIntervalsList(file):

    ret = []
    intervals = getIntervals(file)
    keys = intervals.keys()
    #
    for intervalKey in keys:

        data = intervalKey.split(":")
        chr = data[0]
        strand = data[1]
        ivlist = intervals[intervalKey]
        ivlist = [((chr, strand) + x) for x in ivlist]
        ret.extend(ivlist)

    return ret



def getKey(qval,interval):

    if qval == 0:
        return "0"
    if qval == 255:
        return "255"

    div = qval//interval
    return str(div)

def addORupdateTag(read,mm,ml):


    try:
        current_value = read.get_tag("MM")
        new_value = current_value+";"+mm
    except KeyError:
        new_value = mm

    read.set_tag("MM", new_value, replace=True)

    addlist = (int(value) for value in current_value.split(','))
    try:
        current_value = read.get_tag("ML")
        current_value.extend(addlist)

    except KeyError:
        current_value = addlist

    read.set_tag("ML", current_value, replace=True)

    return read


import pysam
def updatebamtag(infile,output_bam,bamtagout):

    tag_dict = {}
    tsv_reader = csv.reader(bamtagout, delimiter='\t')
    for row in tsv_reader:
        if len(row) >= 3:
            readid = row[0]
            value = (row[1], row[2])
            tag_dict[readid] = value

    with pysam.AlignmentFile(infile, "rb") as input_bam:
        with pysam.AlignmentFile(output_bam, "wb", header=input_bam.header) as outfile:

            for read in input_bam:

                if read.query_name in tag_dict:
                    #add or update
                    mm,ml = tag_dict[read.query_name]
                    read = addORupdateTag(read,mm,ml)

                outfile.write(read)


def calcTpm(geneCounter,tpm,convertMatrix):

    genelist = []
    geneCounter2 = {}
    totalsam = 0

    with open(convertMatrix, encoding='utf-8', newline='') as f:

        for cols in csv.reader(f, delimiter='\t'):

            if cols[0].startswith("#"):
                continue

            id = cols[0]
            chrom = cols[1]
            strand = cols[2]
            exons = cols[3].split(",")
            txStart = exons[0].split(":")[0]
            txEnd = exons[-1].split(":")[1]
            genelist.append((id, chrom, strand, txStart, txEnd))
            poss = cols[3].split(",")
            geneLen = 0
            for posiv in poss:
                iv = posiv.split(":")
                geneLen += (int(iv[1]) - int(iv[0]))

            # print(id,id in geneCounter)
            if id in geneCounter and geneLen > 0:
                avedepth = geneCounter[id] / geneLen
                geneCounter2[id] = (avedepth, geneLen)
                totalsam += avedepth
            else:
                geneCounter2[id] = (0, geneLen)

    tplist = []
    tpmsum = 0
    for tp in genelist:

        id, chrom, strand, txStart, txEnd = tp
        tpmv, genelen = 0, 0
        if id in geneCounter2:
            avedepth, genelen = geneCounter2[id]
            tpmv = (avedepth / totalsam) * 1000000
            tp2 = id, tpmv, chrom, strand, txStart, txEnd,genelen
            tplist.append(tp2)
            tpmsum +=tpmv

    unit = (1000000 / tpmsum)
    outtpm = open(tpm, "w")
    outtpm.write("#ID\tTPM\tChrom\tSTRAND\ttxStart\ttxEnd\tGENELEN" + "\n")
    for tp in tplist:
        id, tpmv, chrom, strand, txStart, txEnd,genelen = tp
        tpmv = tpmv*unit
        outtpm.write(id + "\t" + str(tpmv) + "\t" + str(chrom) + "\t" + str(strand) + "\t" + str(txStart) + "\t" + str(
            txEnd) + "\t" + str(genelen) + "\n")
    outtpm.close()


def is_autosome_or_sex_chromosome(chromosome):

    if chromosome in ['chrX', 'chrY']:
        return True
    if chromosome.startswith('chr'):
        number_part = chromosome[3:]
        if number_part.isdigit() :
            return True

    return False


