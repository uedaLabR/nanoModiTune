import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers


class TransformerBlock(layers.Layer):

    def __init__(self, embed_dim, num_heads, ff_dim, rate=0.1):
        super(TransformerBlock, self).__init__()
        self.att = layers.MultiHeadAttention(num_heads=num_heads, key_dim=embed_dim)
        self.ffn = keras.Sequential(
            [layers.Dense(ff_dim, activation="relu"), layers.Dense(embed_dim),]
        )
        self.layernorm1 = layers.LayerNormalization(epsilon=1e-6)
        self.layernorm2 = layers.LayerNormalization(epsilon=1e-6)
        self.dropout1 = layers.Dropout(rate)
        self.dropout2 = layers.Dropout(rate)

    def call(self, inputs, training):
        attn_output = self.att(inputs, inputs)
        attn_output = self.dropout1(attn_output, training=training)
        out1 = self.layernorm1(inputs + attn_output)
        ffn_output = self.ffn(out1)
        ffn_output = self.dropout2(ffn_output, training=training)
        return self.layernorm2(out1 + ffn_output)

class TokenAndPositionEmbedding(layers.Layer):

    def __init__(self, maxlen, vocab_size, embed_dim):
        super(TokenAndPositionEmbedding, self).__init__()
        self.token_emb = layers.Embedding(input_dim=vocab_size, output_dim=embed_dim)
        self.pos_emb = layers.Embedding(input_dim=maxlen, output_dim=embed_dim)

    def call(self, x):
        maxlen = tf.shape(x)[-1]
        positions = tf.range(start=0, limit=maxlen, delta=1)
        positions = self.pos_emb(positions)
        x = self.token_emb(x)
        return x + positions


def encode_dna_structure(dna, structure):

    # Mapping of DNA bases and secondary structure characters to numerical values
    base_mapping = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
    structure_mapping = {'.': 0, '(': 1, ')': 2}
    # Convert DNA sequence and secondary structure to numerical values
    dna_value = sum(base_mapping[base] * (4 ** i) for i, base in enumerate(reversed(dna)))
    structure_value = sum(structure_mapping[char] * (3 ** i) for i, char in enumerate(reversed(structure)))
    # Combine numerical values of DNA sequence and secondary structure to generate a unique integer
    # Here we use a coefficient larger than the maximum value of DNA sequence (4^3 - 1) to create an offset
    unique_value = dna_value + structure_value * (4 ** 3)

    return unique_value

def encode_dna(dna):

    base_mapping = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
    # Convert DNA sequence and secondary structure to numerical values
    dna_value = sum(base_mapping[base] * (4 ** i) for i, base in enumerate(reversed(dna)))
    return  dna_value

def tokonize(seq,dotb):

    ret = []
    for n in range(len(seq)-2):
        ps = seq[n:n+3]
        pdotb = dotb[n:n+3]
        token = encode_dna_structure(ps,pdotb)
        ret.append(token)

    return ret


def tokonizeN(seq):

    ret = []
    for n in range(len(seq)-2):
        ps = seq[n:n+3]
        token = encode_dna(ps)
        ret.append(token)

    return ret


def toNumberList(data):

    ret = []
    cnt = 0
    for seq,label in data:
        tokens = tokonizeN(seq)
        ret.append((tokens, label))
        cnt += 1

    return ret

import os
def seed_everything(seed):
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    tf.random.set_seed(seed)
    session_conf = tf.compat.v1.ConfigProto(
        intra_op_parallelism_threads=1,
        inter_op_parallelism_threads=1
    )
    sess = tf.compat.v1.Session(graph=tf.compat.v1.get_default_graph(), config=session_conf)
    tf.compat.v1.keras.backend.set_session(sess)


import numpy as np
from sklearn.model_selection import train_test_split
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.regularizers import L1L2
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras.utils import to_categorical
import pandas as pd

def train(data,weightpath,epoch,outhistory):

    # structlist = calc2ndaryStruct(data,calcstruct=False)
    # X = [item[0] for item in structlist]
    # y = [item[1] for item in structlist]

    structlist = toNumberList(data)
    X = []
    y =[]
    for item in structlist:
        if item is not None:
            s = item[0]
            f = item[1]
            if len(s) ==38:
                X.append(s)
                y.append(f)
    # X = [item[0] for item in structlist]
    # y = [item[1] for item in structlist]


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

    # X_train = X_train.reshape((-1,38,1))
    # X_test = X_test.reshape((-1, 38, 1))



    embed_dim = 64  # Embedding size for each token
    num_heads = 8  # Number of attention heads
    ff_dim = 128  # Hidden layer size in feed forward network inside transformer
    maxlen = 38
    # vocab_size = 1728 #12^3
    vocab_size = 4**3  # 4^3

    inputs = layers.Input(shape=(maxlen,))
    embedding_layer = TokenAndPositionEmbedding(maxlen, vocab_size, embed_dim)
    x = embedding_layer(inputs)
    transformer_block = TransformerBlock(embed_dim, num_heads, ff_dim)
    x = transformer_block(x,training=True)
    x = layers.GlobalAveragePooling1D()(x)
    x = layers.Dropout(0.1)(x)
    x = layers.Dense(64, activation='relu', kernel_regularizer=L1L2(l1=0.01, l2=0.01))(x)
    x = layers.Dropout(0.1)(x)
    outputs = layers.Dense(6, activation="softmax")(x)
    model = keras.Model(inputs=inputs, outputs=outputs)
    model.summary()

    checkpoint_path = weightpath
    checkpoint = ModelCheckpoint(filepath=checkpoint_path,
                                 save_weights_only=True,
                                 save_best_only=True,
                                 monitor='val_acc',
                                 verbose=1)

    # Preparing callbacks.
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


import pandas as pd
import pysam
def getData(m6Apath,m5Cpath, psudepath,editingpath, fp_ivtpath,fasta_path,eachsize):

    maxcnt = eachsize
    data = []
    addData(data, editingpath, Flg_I, 'T',maxcnt)
    addData(data, m6Apath, Flg_m6A, 'T', maxcnt)
    addData(data, m5Cpath, Flg_m5C, 'G', maxcnt)
    addData(data, psudepath, Flg_Y, 'A', maxcnt)

    bed_df = pd.read_csv(fp_ivtpath, header=None, sep='\t')

    ccnt = 0
    tcnt = 0

    tuple_list2 = []
    for index, row in bed_df.iterrows():
        chr = row[0]
        pos = row[1]
        alt = row[4]
        info = row[7]

        seqexsist = False
        ifs = info.split(',')
        strand = True
        for i in ifs:
            if "SEQ=" in i:
                qseq0 = i.replace("SEQ=", "")
                if len(qseq0) > 10:
                    qseq = qseq0[1:]
                    # print(qseq)
                    if len(qseq) == 40:
                        # print(qseq, qseq[19], 2, alt, info)
                        if "C" == qseq[19]:
                            ccnt+=1
                        if "T" == qseq[19]:
                            tcnt+=1

                        tuple_list2.append((qseq, Flg_Error))
                    seqexsist = True
            if "STRAND=False" in i:
                strand = False

        if len(tuple_list2) > maxcnt*5:
            break

    data.extend(tuple_list2)
    tuple_list3 = fetch_random_sequences(fasta_path, sequence_length=40, number_of_sequences=maxcnt)
    data.extend(tuple_list3)

    return data

def trainNN(m6Apath,m5Cpath,psudepath,editingpath,fp_ivtpath,ref,weightpath,outhistory,eachsize =15000,epoch=200):

    data = getData(m6Apath,m5Cpath,psudepath,editingpath,fp_ivtpath,ref,eachsize)
    print("finish get data")
    train(data,weightpath,epoch,outhistory)


weightpath = "/mnt/share/ueda/RNA004/resource/model.weights.h5"
m6Apath = "/mnt/ssdnas07/pipeline/rna_v08/source/knownsites/human.hg38.m6A.result.col29.bed"
m5Cpath = "/mnt/ssdnas07/pipeline/rna_v08/source/knownsites/human.hg38.m5C.result.col29.bed"
psudepath = "/mnt/ssdnas07/pipeline/rna_v08/source/knownsites/human.hg38.Pseudo.result.col29.bed"
editingpath = "/mnt/ssdnas07/pipeline/rna_v08/source/knownsites/human.hg38.RNA-editing.result.col29.bed"

fp_ivtpath = "/mnt/share/ueda/RNA004/resource/ivt_result.vcf"
ref = "/mnt/share/ueda/RNA004/U87/U87_IVT/20231227_1535_MN32625_FAX73794_2cf3868f/bam_pass/U87ivt.fa"
outhistory = "/mnt/share/ueda/RNA004/resource/outhistory.csv"

# trainNN(m6Apath,m5Cpath,psudepath,editingpath,fp_ivtpath,ref,weightpath,outhistory,eachsize =15000,epoch=200)