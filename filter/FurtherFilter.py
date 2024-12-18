import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import filter.NNModel as NNModel

Flg_Error = 0
Flg_other = 1
Flg_m6A = 2
Flg_I = 3
Flg_m5C = 4
Flg_Y = 5

import glob
def loadKnownPos(knownPosDir,genome):


    knownPos = {}
    files = glob.glob(knownPosDir)
    # print(files)
    for file in files:

        if genome in file:

            _loadKnownPos(knownPos, file)

    # loadEditingFile(knownPos,editingfile)
    return knownPos

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

def _loadKnownPos(knownPos,path):

    flg = getFlg(path)
    bed_df = pd.read_csv(path, header=None, sep='\t')
    for index, row in bed_df.iterrows():
        chr = row[0]
        pos = row[2]
        key = str(chr) + ":" + str(pos)
        knownPos[key] = flg

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
                print(cnt)


def printout_unchange(bed_df,output):

    edited_rows = []
    for index, row in bed_df.iterrows():

        edited_rows.append(row)

    title = ['#CHR','POS','ID','REF','ALT','QUAL','FILTER','INFO']
    outall = os.path.splitext(output)[0]+"_all.vcf"
    new_df = pd.DataFrame(edited_rows)
    new_df.columns = title
    new_df.to_csv(outall, sep='\t', index=False)

    new_df2 = new_df[new_df.iloc[:, 6] == True]
    new_df2.columns = title
    new_df2.to_csv(output, sep='\t', index=False)


def getCalledFlg(ALT):

    if ALT=="m":
        return Flg_m5C
    if ALT=="a":
        return Flg_m6A
    if ALT=="17802":
        return Flg_Y
    if ALT=="17596":
        return Flg_I

import numpy as np
from sklearn.model_selection import train_test_split
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.regularizers import L1L2
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras.utils import to_categorical
import pandas as pd
import os
import re

os.environ["CUDA_VISIBLE_DEVICES"] = "-1"


def isDrach(qseq):

    # Define the regular expression pattern for DRACH motif
    drach_pattern = re.compile(r"[AGUT][AG]AC[AUCUT]", re.IGNORECASE)

    # Search for the pattern in the input sequence
    return bool(drach_pattern.search(qseq))


def isRepresentitiveY(qseq):

    seqs = ["TTTTT", "GTTCG", "GTTCC", "GTTCA", "GTTCT",
            "TGTAG", "TCTAG", "GGTGC", "GCTGA", "CGTGA",
            "AGTGG", "TGTGG", "GGTGG"]


    return any(motif in qseq for motif in seqs)

def cgRich(qseq,threshold=0.7):

    qseq = qseq.upper()
    # Count the occurrences of C and G in the sequence
    cg_count = qseq.count('C') + qseq.count('G')
    # Calculate the CG content as a proportion of the total length
    cg_ratio = cg_count / len(qseq) if len(qseq) > 0 else 0

    # Check if the CG content meets or exceeds the threshold
    return cg_ratio >= threshold

def isMotief(called_flg,info):

    qseq = None
    ifs = info.split(',')
    for i in ifs:
        if "SEQ=" in i:
            print(i)
            qseq0 = i.replace("SEQ=", "")
            if len(qseq0) > 10:
                qseq = qseq0[17:23]
    if qseq == None:
        return False
    if called_flg == Flg_m6A:
        return isDrach(qseq)
    elif called_flg == Flg_Y:
        return isRepresentitiveY(qseq)
    elif called_flg == Flg_m5C:
        return cgRich(qseq)

    return False

def group_by_position(dataHolder, distance=2):
    """
    Groups items in dataHolder by chr, strand, and within a specified distance for pos.

    Parameters:
    dataHolder (list of tuples): List of tuples where each tuple represents
                                 (chr, strand, pos, called_flg, filteredout, predict_flg, knownAndFlgOk, known_motif, row).
    distance (int): The maximum distance (in bases) for grouping items by pos.

    Returns:
    list of lists: A list of grouped items, where each group is a list of tuples.
    """
    # Sort dataHolder by chr, strand, pos
    dataHolder = sorted(dataHolder, key=lambda x: (x[0], x[1], x[2]))

    grouped_data = []
    current_group = []
    currentset = set()

    for item in dataHolder:
        # Unpack item
        chr, strand, pos, called_flg, filteredout, predict_flg, knownAndFlgOk, known_motif, row = item

        if not current_group:
            # Start a new group
            current_group.append(item)
            currentset.add(called_flg)
        else:
            # Compare positions
            last_pos = current_group[-1][2]
            current_pos = pos

            # If conditions match, add to current group
            if (item[0] == current_group[-1][0] and item[1] == current_group[-1][1]
                    and (current_pos - last_pos) <= distance and called_flg not in currentset):
                current_group.append(item)
                currentset.add(called_flg)
            else:
                # Start a new group
                grouped_data.append(current_group)
                current_group = [item]
                currentset = {called_flg}  # Start fresh set with current flag

    # Append the last group if not empty
    if current_group:
        grouped_data.append(current_group)

    return grouped_data

def setfilter(group):
    """
    Sets the filter status for each item in the group based on known flags and motifs.

    Parameters:
    group (list of tuples): List of tuples representing each data item in the group,
                            including information about chr, strand, pos, called_flg,
                            filteredout, predict_flg, knownAndFlgOk, known_motif, and row.
    """
    knownidx = []
    knownmotiefL = []
    predictOK = []
    notfilteredout = []
    drachExist = False

    if len(group) ==1:
        return

    # Pass 1: Identify known flags, motifs, and check for drachExist
    for n, (chr, strand, pos, called_flg, filteredout, predict_flg, knownAndFlgOk, known_motif, row) in enumerate(group):
        if knownAndFlgOk:
            knownidx.append(n)
        if known_motif:
            knownmotiefL.append(n)
            if called_flg == Flg_m6A:
                drachExist = True
        if called_flg == predict_flg:
            predictOK.append(n)

    takeidx = []
    if len(knownidx) > 0:
        takeidx = knownidx
    elif len(predictOK) > 0:
        takeidx = predictOK
    elif len(knownmotiefL) > 0:
        takeidx = knownmotiefL


    if len(takeidx) == 0:
        #take min P value
        maxidx = 0
        scoremin = 0
        for n, (chr, strand, pos, called_flg, filteredout, predict_flg, knownAndFlgOk, known_motif, row) in enumerate(
                group):

            if row[6] == True:
                score = float(row[5])
                if score >  scoremin:
                    scoremin = score
                    maxidx = n

        takeidx.append(maxidx)


    # Pass 2: Set filter status based on conditions
    for n, (chr, strand, pos, called_flg, filteredout, predict_flg, knownAndFlgOk, known_motif, row) in enumerate(group):
        setTrue = False
        filter_true = row[6] == True
        # If drach exists, set m6A and neighboring Y to True
        if drachExist and called_flg == Flg_Y:
            if filter_true:
                setTrue = True

        # Set True if current index is in known indexes or known motifs
        if filter_true and n in takeidx:
            row[6] = True
            setTrue = True

        if filter_true and called_flg == Flg_I:
            setTrue = True

        # If none of the above conditions are met, set to False
        if not setTrue:
            row[6] = False




def checkFlgAndPoliNuc(called_flg,predict_flg,strand,info,filter_afthres):

    if called_flg != predict_flg:
        return True

    qseq = None
    ifs = info.split(',')
    af = 0
    for i in ifs:
        if "SEQ=" in i:
            # print(i)
            qseq0 = i.replace("SEQ=", "")
            if len(qseq0) > 10:
                qseq = qseq0[18:22]

        if "AF=" in i:
            afs = i.replace("AF=", "")
            if len(afs) > 0:
                af = float(af)


    if qseq == None:
        return False

    if af < filter_afthres:
        return True

    if  called_flg == Flg_m5C:
        cc = qseq.count('C')
        if cc >=2:
            return True
        gg = qseq.count('G')
        if gg >=2: #GCG
            return True

    if  called_flg == Flg_Y:
        tt = qseq.count('T')
        if tt >=2:
            return True
        aa = qseq.count('A')
        if aa >=2:
            return True


    return False



def classification(input,output,checkpoint_path,knowndir,genome="hg38"):

    model = NNModel.getModel()
    model.summary()
    model.compile()
    model.load_weights(checkpoint_path)
    # print(input,output,checkpoint_path,m5Cpath,psudepath)
    print("loading DB")
    knownPos = loadKnownPos(knowndir,genome)
    print(knownPos.keys())
    bed_df = pd.read_csv(input, header=None, sep='\t')

    tuple_list = []
    poslist = []
    counter = {}
    cnt = 0
    for index, row in bed_df.iterrows():

        chr = row[0]
        pos = row[1]
        alt = row[4]
        filter = row[6]
        info = row[7]
        if filter:
            if alt in counter:
                counter[alt] = counter[alt]+1
            else:
                counter[alt] = 1
            cnt+=1
        if cnt%1000==0:
            print(cnt,counter,len(tuple_list))

        seqexsist = False
        ifs = info.split(',')
        strand = True
        for i in ifs:
            if "SEQ=" in i:
                # print(i)
                qseq0 = i.replace("SEQ=", "")
                if len(qseq0) > 10:
                    qseq = qseq0[1:]
                    # print(len(qseq),qseq)
                    if len(qseq) == 40:

                        tuple_list.append((qseq,0))
                        key = str(chr) + ":" + str(pos)
                        poslist.append(key)



    numlist = NNModel.toNumberList(tuple_list)
    X_test = []
    for item in numlist:
        if item is not None:
            s = item[0]
            if len(s) ==38:
                X_test.append(s)


    X_test = np.array(X_test)
    print("X_len",len(X_test))

    # model = tf.keras.models.load_model(checkpoint_path)
    predict = model.predict(X_test)
    posdict = {}
    counter = {}
    for n in range(len(predict)):

        pre =np.argmax(predict[n])
        cc = counter.get(pre,0)
        cc = cc +1
        counter[pre] = cc
        key = poslist[n]
        posdict[key] = pre

    dataHolder = []
    for index, row in bed_df.iterrows():

        predict_flg = -1
        known_flg = -1
        filterflg = 0

        chr = row[0]
        pos = row[1]
        REF = row[3]
        ALT = row[4]
        info = row[7]

        strand = "STRAND=True" in info

        key0 = str(chr) + ":" + str(pos)
        if key0 in posdict:
            predict_flg = posdict[key0]

        inKnownDB = key0 in knownPos
        if inKnownDB:
            known_flg = knownPos[key0]
            inKnownDB = True

        called_flg = getCalledFlg(ALT)
        filteredout = False
        if known_flg < 0:
            if predict_flg == Flg_Error:
                row[6] = False # Filter out
                filteredout = True

            if predict_flg == Flg_other and called_flg != Flg_I:
                row[6] = False # Filter out
                filteredout = True


        filter_afthres = 0.2
        #For m5C, Y require flg match, and non polynuc
        if called_flg == Flg_m5C or called_flg == Flg_Y:
            checkNG = checkFlgAndPoliNuc(called_flg,predict_flg,strand,info,filter_afthres)
            if checkNG:
                row[6] = False  # Filter out
                filteredout = True

        #

        knownAndFlgOk = inKnownDB and (known_flg == called_flg)
        if knownAndFlgOk:
            row[6] = True
            row[7] = info + ",knownSites=True"
            if filteredout == True:
                row[7] = row[7] +",rescued=True"
            filteredout = False

        known_motif = isMotief(called_flg, info)
        if known_motif:
            row[6] = True
            row[7] = info + ",knownMotif=True"
            if filteredout == True:
                row[7] = row[7] +",rescued=True"
            filteredout = False



        if row[6] == True:
            dataHolder.append((chr,strand,pos,called_flg,filteredout,predict_flg,knownAndFlgOk,known_motif,row))

    #sort
    grouped_data = group_by_position(dataHolder)
    #filter neighbor
    edited_rows = []
    for group in grouped_data:

        setfilter(group)
        for data in group:
            row = data[8]
            edited_rows.append(row)


    title = ['#CHR','POS','ID','REF','ALT','QUAL','FILTER','INFO']
    # outall = os.path.splitext(output)[0]+"_all.vcf"
    new_df = pd.DataFrame(edited_rows)
    new_df.columns = title
    new_df.to_csv(output, sep='\t', index=False)
    # new_df2 = new_df[new_df.iloc[:, 6] == True]
    # new_df2.columns = title
    # new_df2.to_csv(output, sep='\t', index=False)


def run():

    checkpoint_path =  "/mnt/share/ueda/RNA004/resource/ntmodel.weights.h5"
    # input = "/mnt/share/ueda/RNA004/Dorado0.8/bamout/test_stats/unfilter_result.vcf"
    input = "/mnt/ssdnas07/nanozero/rna/nanomoditune_v02/HEK293T_DR13/HEK293T_DR13/unfilter_result.vcf"
    output = "/mnt/share/ueda/RNA004/hek293/result_filter.vcf"
    knowndir = "/mnt/ssdnas07/pipeline/rna_v08/source/knownsites/human*.bed"

    classification(input,output,checkpoint_path,knowndir,genome="hg38")

run()