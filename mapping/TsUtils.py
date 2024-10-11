import csv

class TsInfo:

    def __init__(self,cols,mode):

        if mode == 0:
            self.id = cols[0]
            self.chrom = cols[1]
            self.strand = cols[2]
            numExon = int(cols[7])
            exonStarts = cols[8].split(",")
            exonEnds = cols[9].split(",")
            self.exons = []
            for n in range(numExon):
                s = int(exonStarts[n])
                e = int(exonEnds[n])
                self.exons.append((s,e))
        else:

            self.id = cols[0]
            self.chrom = cols[1]
            self.strand = cols[2]
            exons = cols[3].split(",")
            # print(exons)
            self.exons = []
            for exon in exons:
                es = exon.split(":")
                # print("es",es)
                s = int(es[0])
                e = int(es[1])+1
                self.exons.append((s, e))


def getTsDict(file,mode):

    data = {}
    print("start reading file",mode)
    with open(file, encoding='utf-8', newline='') as f:

        for cols in csv.reader(f, delimiter='\t'):

            if cols[0].startswith("#"):
                continue

            id = cols[0]
            tsInfo = TsInfo(cols,mode)
            # print(tsInfo)
            data[id] = tsInfo

    return data

import mappy as mp
def getSeq(aligner,tsInfo):

    allseq = ""
    for s,e in tsInfo.exons:

        seq = aligner.seq(tsInfo.chrom,s,e)
        if seq is not None:
            allseq = allseq + seq

    if tsInfo.strand == "-":
        allseq = mp.revcomp(allseq)

    return allseq

def outputToFasta(ddict,ref,out):

    aligner = mp.Aligner(ref)
    f = open(out, 'w')
    for k in ddict:
        tsInfo = ddict[k]
        seq = getSeq(aligner,tsInfo)
        f.write(">"+k+"\n")
        f.write(seq+"\n")


    f.close()
