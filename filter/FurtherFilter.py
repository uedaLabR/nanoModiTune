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
        pos = row[1]
        key = str(chr) + ":" + str(pos)
        knownPos[key] = flg



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

os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

def classification(input,output,checkpoint_path,knowndir,genome="hg38"):

    model = NNModel.getModel()
    model.summary()
    model.compile()
    model.load_weights('/mnt/share/ueda/RNA004/resource/test.weights.h5')
    # print(input,output,checkpoint_path,m5Cpath,psudepath)

    knownPos = loadKnownPos(knowndir,genome)
    print(knownPos.keys())

    bed_df = pd.read_csv(input, header=None, sep='\t')
    ccnt = 0
    tcnt = 0

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
                print(i)
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

    edited_rows = []
    for index, row in bed_df.iterrows():

        predict_flg = -1
        known_flg = -1

        chr = row[0]
        pos = row[1]
        REF = row[3]
        ALT = row[4]

        info = row[7]

        key0 = str(chr) + ":" + str(pos)
        if key0 in posdict:
            predict_flg = posdict[key0]

        inKnownDB = key0 in knownPos
        if inKnownDB:
            known_flg = knownPos[key0]
            inKnownDB = True

        filteredout = False
        if known_flg < 0:
            if predict_flg == Flg_Error:
                row[6] = False # Filter out
                filteredout = True
            if predict_flg == Flg_other:
                row[6] = False # Filter out
                filteredout = True
        called_flg = getCalledFlg(ALT)

        if inKnownDB and (known_flg ==called_flg) and filteredout == True :
            row[6] = True
            row[7] = info+",knownSites=True"

        if filter:
            edited_rows.append(row)

    title = ['#CHR','POS','ID','REF','ALT','QUAL','FILTER','INFO']
    # outall = os.path.splitext(output)[0]+"_all.vcf"
    new_df = pd.DataFrame(edited_rows)
    new_df.columns = title
    # new_df.to_csv(outall, sep='\t', index=False)
    new_df2 = new_df[new_df.iloc[:, 6] == True]
    new_df2.columns = title
    new_df2.to_csv(output, sep='\t', index=False)


def run():

    checkpoint_path =  "/mnt/share/ueda/RNA004/resource/ntmodel.weights.h5"
    input = "/mnt/ssdnas07/nanozero/rna/nanomoditune_v01/HEK293T_DR13/HEK293T_DR13/unfilter_result.vcf"
    output = "/mnt/share/ueda/RNA004/hek293/result_filter.vcf"
    knowndir = "/mnt/ssdnas07/pipeline/rna_v08/source/knownsites/human*.bed"
    classification(input,output,checkpoint_path,knowndir,genome="hg38")



run()