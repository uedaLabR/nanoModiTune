import glob

Flg_Error = 0
Flg_other = 1
Flg_m6A = 2
Flg_I = 3
Flg_m5C = 4
Flg_Y = 5

def loadKnownPos(knownPos,knowndir, genome):

    # loadEditingFile(knownPos, editingfile)
    counts = {}
    files = glob.glob(knowndir)
    print(files)
    for file in files:

        if genome in file:

            flg = getFlg(file)
            cnt = _loadKnownPos(knownPos, flg,file)
            counts[flg] = cnt

    return knownPos,counts

def getFlg(path):

    if "m5C" in path:
        return Flg_m5C
    if "m6A" in path:
        return Flg_m6A
    if "Pseudo" in path:
        return Flg_Y
    if "editing" in path:
        return Flg_I
    return -1

import pandas as pd
def _loadKnownPos(knownPos,flg,path):


    bed_df = pd.read_csv(path, header=None, sep='\t')
    cnt = 0
    for index, row in bed_df.iterrows():
        chr = row[0]
        pos = row[2]
        key = str(chr) + ":" + str(pos)
        knownPos[key] = flg
        cnt+=1
    return cnt

def stats(bed1,knownPos):

    cnt = 0
    scorethres = 20
    ratiothres = 10
    counter = {}
    counterknown = {}

    with open(bed1, "r") as bed_file:
        for line in bed_file:
            columns = line.strip().split("\t")
            alt = columns[3]
            score = int(columns[4])
            ratio = float(columns[10])
            if score >= scorethres and ratio > ratiothres:

                key0 = str(columns[0]) + ":" + str(columns[2])
                # print(key0,alt,score,ratio)
                inKnownDB = key0 in knownPos
                if inKnownDB:
                    if alt == "a":
                        flg = knownPos[key0]
                        if flg != Flg_m6A:
                            inKnownDB = False
                    if alt == "17596":
                        flg = knownPos[key0]
                        if flg != Flg_I:
                            inKnownDB = False

                if inKnownDB:
                    if alt in counterknown:
                        counterknown[alt] = counterknown[alt] + 1
                    else:
                        counterknown[alt] = 1

                if alt in counter:
                    counter[alt] = counter[alt] + 1
                else:
                    counter[alt] = 1
                cnt += 1

            if cnt % 1000 == 0:
                print(columns[0],columns[1],cnt, counter, counterknown)

        print(counter)
        print(counterknown)

import gzip
def loadEditingFile(knownPos,path):

    flg = Flg_I
    with gzip.open(path,'rt', encoding='utf-8') as inst:
        cnt = 0
        for line in inst:
            if cnt > 0:
                data = (line.split("\t"))
                key = str(data[1]) + ":" + str(data[2])
                knownPos[key] = flg
            cnt+=1
            if cnt % 10000==0:
                print(cnt,key)

def run():

    genome="hg38"
    knowndir = "/mnt/ssdnas07/pipeline/rna_v08/source/knownsites/human*.bed"
    # editingfile = "/mnt/share/ueda/RNA004/nanoModiTune/TABLE1_hg38_v3.txt.gz"

    knownPos = {}
    knownPos,counts = loadKnownPos(knownPos,knowndir, genome)
    # print(counts)

    #vcf1 = "/mnt/share/ueda/RNA004/hek293/result_filter.vcf"
    # vcf1 = "/mnt/ssdnas07/nanozero/rna/nanomoditune_v01/HEK293T_DR13/HEK293T_DR13/unfilter_result.vcf"
    vcf1 = "/share/ueda/nanoModiTune/Hek293pu.bed"
    stats(vcf1, knownPos)

    # vcf2 = "/mnt/share/ueda/RNA004/hek293/result_filter.vcf"
    # stats(vcf2, knownPos)


    # vcf3 = "/mnt/ssdnas07/nanozero/rna/nanomoditune_v01/U87_IVT/U87_IVT/unfilter_result.vcf"
    # stats(vcf3, knownPos)

run()