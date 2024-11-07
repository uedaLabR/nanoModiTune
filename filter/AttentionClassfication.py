
import numpy as np
from sklearn.model_selection import train_test_split
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.regularizers import L1L2
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras.utils import to_categorical
import pandas as pd
import filter.NNModel as NNModel
import tensorflow as tf

def train(data,weightpath,epoch,outhistory):

    # structlist = calc2ndaryStruct(data,calcstruct=False)
    # X = [item[0] for item in structlist]
    # y = [item[1] for item in structlist]

    numlist = NNModel.toNumberList(data)
    X = []
    y =[]
    for item in numlist:
        if item is not None:
            s = item[0]
            f = item[1]
            if len(s) ==38:
                X.append(s)
                y.append(f)

    X = np.array(X)
    y = np.array(y)
    print(X.shape)
    print(y.shape)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.15, random_state=42)

    np.save('/mnt/share/ueda/RNA004/resource/X_test.npy', X_test)
    np.save('/mnt/share/ueda/RNA004/resource/y_test.npy', y_test)

    y_train = to_categorical(y_train, num_classes=6)
    y_test = to_categorical(y_test, num_classes=6)

    X_train = np.array(X_train)
    y_train = np.array(y_train)


    print("X_train shape:", X_train.shape)
    print("y_train shape:", y_train.shape)
    print("X_test shape:", X_test.shape)
    print("y_test shape:", y_test.shape)

    print("X_train shape:", X_train.shape)
    print("y_train shape:", y_train.shape)
    print("X_train sample:", X_train[:5])
    print("y_train sample:", y_train[:5])

    model = NNModel.getModel()
    model.summary()

    checkpoint = ModelCheckpoint(filepath=weightpath,
                                 save_weights_only=True,
                                 save_best_only=True,
                                 monitor='val_acc',
                                 verbose=1)

    # Preparing callbacks.
    lr_schedule = tf.keras.optimizers.schedules.ExponentialDecay(
        initial_learning_rate=0.001,  # Starting learning rate
        decay_steps=10000,  # After how many steps to apply decay
        decay_rate=0.9,  # The decay rate to apply
        staircase=True)  # If True, decay the learning rate at discrete intervals

    # Initialize the Adam optimizer with the learning rate schedule
    adam = tf.keras.optimizers.Adam(
        learning_rate=lr_schedule,  # Use the learning rate schedule
        beta_1=0.9,
        beta_2=0.999,
        epsilon=1e-07,  # For default value, you can simply omit this or set it to 1e-7
        amsgrad=False)

    model.compile(optimizer=adam,
                  # loss='sparse_categorical_crossentropy',
                  loss='categorical_crossentropy',
                  metrics=['acc'])


    callbacks = [
        checkpoint
    ]

    # seed_everything(42)

    print(X_train)
    # Train the model.
    history = model.fit(x=X_train,
                        y= y_train,
                        batch_size=64,
                        epochs=epoch,
                        callbacks=callbacks,
                        shuffle=True,
                        validation_data=(X_test,y_test))

    history_df = pd.DataFrame(history.history)
    history_df.to_csv(outhistory, index=False)


# data = []
# data.append(("GCCTGCCCCCGCTAACCGGCTTTTTGCCCAAATGGGCCATT",1))
# data.append(("TTCATCGGCGTAAATCTAACTTTCTTCCCACAACACTTTCT",1))


import random


def fetch_random_sequences(fasta, chromosome, sequence_length=40, number_of_sequences=5):

    sequences = []
    chromosome_length = fasta.get_reference_length(chromosome)

    for _ in range(number_of_sequences*10):

        start = random.randint(0, chromosome_length - sequence_length)
        end = start + sequence_length
        sequence = fasta.fetch(chromosome, start, end)
        if "N" in sequence:
            continue

        mid = sequence[19]
        if mid == "T" or mid =="C":
            if mid == "C":
                if random.random() < 0.5:
                    continue
            sequences.append(sequence)

    return sequences


import pysam
import random
import re


def fetch_random_sequences(fasta_path, sequence_length=40, number_of_sequences=5):

    fasta = pysam.FastaFile(fasta_path)
    autosomes = [ref for ref in fasta.references if re.match(r'chr\d+$', ref)]
    non_autosomes = [ref for ref in fasta.references if not re.match(r'chr\d+$', ref)]
    if random.random() < 0.5:
        selected_chromosome = random.choice(autosomes)
    else:
        selected_chromosome = random.choice(non_autosomes)
    chromosome_length = fasta.get_reference_length(selected_chromosome)

    sequences = []
    for _ in range(number_of_sequences*10):
        start = random.randint(0, chromosome_length - sequence_length)
        end = start + sequence_length
        sequence = fasta.fetch(selected_chromosome, start, end).upper()
        if "N" in sequence:
            continue
        if len(sequences) ==  number_of_sequences:
            break
        if sequence[19] == "C" or sequence[19] == "T":
            sequences.append((sequence,Flg_other))

    fasta.close()

    return sequences

def reverse_complement(dna_seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_comp = ''.join([complement[base] for base in reversed(dna_seq)])
    return reverse_comp


def is_drach(sequence):

    # Check if the sequence length is exactly 5
    if len(sequence) != 5:
        return False

    # Define the allowed bases for each position in the DRACH motif
    d_bases = ['A', 'G', 'T']
    r_bases = ['A', 'G']
    a_base = ['A']
    c_base = ['C']
    h_bases = ['A', 'C', 'T']

    # Check each position against the allowed bases
    if sequence[0] in d_bases and sequence[1] in r_bases and sequence[2] in a_base and sequence[3] in c_base and \
            sequence[4] in h_bases:
        return True
    else:
        return False


Flg_Error = 0
Flg_other = 1
Flg_m6A = 2
Flg_I = 3
Flg_m5C = 4
Flg_Y = 5

def addData(data,file,flg,nuc,maxcnt):

    bed_df = pd.read_csv(file, header=None, sep='\t')
    column_19 = bed_df.iloc[:, 18]
    cnt = 0
    for x in column_19:
        pseq = x
        if x[19] == nuc:
            pseq = x[:-1]
        if len(pseq) == 40:
            data.append((pseq, flg))
            cnt += 1
            if cnt > maxcnt:
                break
    print("append",cnt)

import pandas as pd
import pysam
def getData(m6Apath,m5Cpath, psudepath,editingpath, fp_ivtpath,fasta_path,eachsize):

    maxcnt = eachsize
    data = []


    bed_df = pd.read_csv(fp_ivtpath, header=None, sep='\t')

    ccnt = 0
    tcnt = 0

    tuple_list2 = []
    for index, row in bed_df.iterrows():
        chr = row[0]
        pos = row[1]
        alt = row[4]
        info = row[7]
        filter = row[6]


        seqexsist = False
        ifs = info.split(',')
        if filter == True:
            strand = True
            for i in ifs:
                if "SEQ=" in i:
                    qseq0 = i.replace("SEQ=", "")
                    if len(qseq0) > 10:
                        qseq = qseq0[1:]
                        # print(qseq)
                        if len(qseq) == 40 and filter:
                            tuple_list2.append((qseq, Flg_Error))
                        seqexsist = True
                if "STRAND=False" in i:
                    strand = False



    data.extend(tuple_list2)
    print("fp data size=",len(tuple_list2))
    tuple_list3 = fetch_random_sequences(fasta_path, sequence_length=40, number_of_sequences=maxcnt*3)
    data.extend(tuple_list3)
    print("randomseq=", len(tuple_list3))
    addData(data, editingpath, Flg_I, 'T',maxcnt)
    addData(data, m6Apath, Flg_m6A, 'T', maxcnt)
    addData(data, m5Cpath, Flg_m5C, 'G', maxcnt)
    addData(data, psudepath, Flg_Y, 'A', maxcnt)
    return data

def trainNN(m6Apath,m5Cpath,psudepath,editingpath,fp_ivtpath,ref,weightpath,outhistory,eachsize =1000,epoch=200):

    data = getData(m6Apath,m5Cpath,psudepath,editingpath,fp_ivtpath,ref,eachsize)
    print("finish get data")
    train(data,weightpath,epoch,outhistory)


checkpoint_path =  "/mnt/share/ueda/RNA004/resource/ntmodel.weights.h5"
m6Apath = "/mnt/ssdnas07/pipeline/rna_v08/source/knownsites/human.hg38.m6A.result.col29.bed"
m5Cpath = "/mnt/ssdnas07/pipeline/rna_v08/source/knownsites/human.hg38.m5C.result.col29.bed"
psudepath = "/mnt/ssdnas07/pipeline/rna_v08/source/knownsites/human.hg38.Pseudo.result.col29.bed"
editingpath = "/mnt/ssdnas07/pipeline/rna_v08/source/knownsites/human.hg38.RNA-editing.result.col29.bed"

fp_ivtpath = "/mnt/ssdnas07/nanozero/rna/nanomoditune_v01/U87_IVT/U87_IVT/unfilter_result.vcf"
ref = "/mnt/share/ueda/RNA004/U87/U87_IVT/20231227_1535_MN32625_FAX73794_2cf3868f/bam_pass/U87ivt.fa"
outhistory = "/mnt/share/ueda/RNA004/resource/outhistory.csv"

# trainNN(m6Apath,m5Cpath,psudepath,editingpath,fp_ivtpath,ref,checkpoint_path,outhistory,eachsize =20000,epoch=100)