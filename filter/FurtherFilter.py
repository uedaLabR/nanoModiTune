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

    dna = dna.upper()
    # Mapping of DNA bases and secondary structure characters to numerical values
    base_mapping = {'A': 0, 'T': 1, 'C': 2, 'G': 3,'U':1,'N':0}
    structure_mapping = {'.': 0, '(': 1, ')': 2}
    # Convert DNA sequence and secondary structure to numerical values
    dna_value = sum(base_mapping[base] * (4 ** i) for i, base in enumerate(reversed(dna)))
    structure_value = sum(structure_mapping[char] * (3 ** i) for i, char in enumerate(reversed(structure)))
    # Combine numerical values of DNA sequence and secondary structure to generate a unique integer
    # Here we use a coefficient larger than the maximum value of DNA sequence (4^3 - 1) to create an offset
    unique_value = dna_value + structure_value * (4 ** 3)

    return unique_value

def encode_dna(dna):

    dna=dna.upper()
    base_mapping = {'A': 0, 'T': 1, 'C': 2, 'G': 3,'U':1,'N':0}
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
    for seq in data:
        tokens = tokonizeN(seq)
        ret.append(tokens)
        cnt += 1

    return ret


def loadKnownPos(knownPos,path,flg):

    bed_df = pd.read_csv(path, header=None, sep='\t')
    for index, row in bed_df.iterrows():
        chr = row[0]
        pos = row[1]
        key = str(chr) + ":" + str(pos)
        knownPos[key] = flg


def findShiftM5C(info):

    ifs = info.split(',')
    strandFalse = "STRAND=False" in info
    for i in ifs:
        if "SEQ=" in i:
            qseq = i.replace("SEQ=", "")
            if len(qseq) < 21:
                return 0
            if qseq[20] == "C":

                if strandFalse:
                    return -1
                else:
                    return 1

            if qseq[18] == "C":

                if strandFalse:
                    return 1
                else:
                    return -1
    if strandFalse:
        return -1
    else:
        return 1

def findShiftM6A(info):

    ifs = info.split(',')
    strandFalse = "STRAND=False" in info
    for i in ifs:
        if "SEQ=" in i:
            qseq = i.replace("SEQ=", "")
            # print(len(qseq))
            if len(qseq) < 21:
                return 0
            if qseq[20] == "A":

                if strandFalse:
                    return 1
                else:
                    return -1

            if qseq[18] == "A":

                if strandFalse:
                    return -1
                else:
                    return 1
    if strandFalse:
        return 1
    else:
        return -1

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



import numpy as np
from sklearn.model_selection import train_test_split
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.regularizers import L1L2
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras.utils import to_categorical
import pandas as pd
import os
def classification(input,output,checkpoint_path,knowndir):

    # print(input,output,checkpoint_path,m5Cpath,psudepath)
    Flg_m5C = 0
    Flg_Y = 1
    Flg_Error = 2
    Flg_other = 3
    Flg_nonDrachm6A = 4

    knownPos = loadKnownPos(knowndir)


    bed_df = pd.read_csv(input, header=None, sep='\t')
    ccnt = 0
    tcnt = 0

    tuple_list2 = []
    poslist = []
    for index, row in bed_df.iterrows():
        chr = row[0]
        pos = row[1]
        alt = row[4]
        info = row[7]
        if alt == "m6A":
            continue

        seqexsist = False
        ifs = info.split(',')
        strand = True
        for i in ifs:
            if "SEQ=" in i:
                print(i)
                qseq0 = i.replace("SEQ=", "")
                if len(qseq0) > 10:
                    qseq = qseq0[1:]
                    print(len(qseq),qseq)
                    if len(qseq) == 40:
                        # print(qseq, qseq[19], 2, alt, info)
                        if "C" == qseq[19]:
                            ccnt += 1
                        if "T" == qseq[19]:
                            tcnt += 1

                        tuple_list2.append(qseq)
                        key = str(chr) + ":" + str(pos)
                        poslist.append(key)

    if len(poslist) == 0:

        printout_unchange(bed_df,output)
        return

    # np.save('/mnt/share/ueda/RNA004/resource/X_test.npy', X_test)
    # np.save('/mnt/share/ueda/RNA004/resource/y_test.npy', y_test)
    X_test = toNumberList(tuple_list2)
    X_test = np.array(X_test)
    # print("X_test shape:", X_test.shape)


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
    x = transformer_block(x)
    x = layers.GlobalAveragePooling1D()(x)
    x = layers.Dropout(0.1)(x)
    x = layers.Dense(64, activation='relu', kernel_regularizer=L1L2(l1=0.01, l2=0.01))(x)
    x = layers.Dropout(0.1)(x)
    outputs = layers.Dense(5, activation="softmax")(x)
    model = keras.Model(inputs=inputs, outputs=outputs)
    model.summary()

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
        epsilon=None,  # For default value, you can simply omit this or set it to 1e-7
        amsgrad=False)

    model.compile(optimizer=adam,
                  # loss='sparse_categorical_crossentropy',
                  loss='categorical_crossentropy',
                  metrics=['acc'],run_eagerly=True)

    # callbacks = [
    #     checkpoint,EarlyStopping(patience=20)
    # ]
    # Train the model.
    model.load_weights(checkpoint_path)
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

        chr = row[0]
        pos = row[1]
        info = row[7]
        shift = 0

        key0 = str(chr) + ":" + str(pos)
        key1 = str(chr) + ":" + str(pos+1)
        key_1 = str(chr) + ":" + str(pos-1)
        shiftset = False
        inKnownDB = False
        if key0 in knownPos:
            flg = knownPos[key0]
        elif key1 in knownPos:
            flg = knownPos[key1]
            inKnownDB = True
            shift = 1
            shiftset = True
        elif key_1 in knownPos:
            flg = knownPos[key_1]
            inKnownDB = True
            shift = -1
            shiftset = True
        else:
            flg = posdict.get(key0,-1)

        if flg == Flg_m5C:

            if shiftset == False:
                shift = findShiftM5C(info)

            row[4] = "m5C"
            refbase = row[3]
            if refbase == "T":
                row[3] ="C"
                row[1] = int(row[1]) + shift
                row[7] = row[7]+",shift="+str(shift)

        if flg == Flg_nonDrachm6A:

            row[3] = "A"
            row[4] = "m6A"
            if shiftset == False:
                shift = findShiftM6A(info)

            row[1] = int(row[1]) + shift
            row[7] = row[7] + ",shift=" + str(shift)+",nonDRACH=True"

        if flg == Flg_Y:
            row[4] = "Y"
        if flg == Flg_Error:
            row[6] = False # Filter out
        if flg == Flg_other:
            row[6] = False # Filter out

        if inKnownDB:
            row[6] = True

        edited_rows.append(row)

    title = ['#CHR','POS','ID','REF','ALT','QUAL','FILTER','INFO']
    outall = os.path.splitext(output)[0]+"_all.vcf"
    new_df = pd.DataFrame(edited_rows)
    new_df.columns = title
    new_df.to_csv(outall, sep='\t', index=False)

    new_df2 = new_df[new_df.iloc[:, 6] == True]
    new_df2.columns = title
    new_df2.to_csv(output, sep='\t', index=False)


def run():

    checkpoint_path = "/mnt/share/ueda/RNA004/resource/model_weights_5.h5"
    # input = "/mnt/share/ueda/RNA004/U87/U87_WT/20231215_1505_MN32625_FAX70236_31c89747/unfilter_result.vcf"
    # output = "/mnt/share/ueda/RNA004/U87/U87_WT/20231215_1505_MN32625_FAX70236_31c89747/result_filter.vcf"
    m5Cpath = "/mnt/share/ueda/RNA004/resource/human.hg38.m5C.result.col29.bed"
    psudepath = "/mnt/share/ueda/RNA004/resource/human.hg38.Pseudo.result.col29.bed"

    # classification(input,output,checkpoint_path,m5Cpath,psudepath)

    input = "/mnt/share/ueda/RNA004/U87/U87_IVT/20231227_1535_MN32625_FAX73794_2cf3868f/unfilter_result.vcf"
    output = "/mnt/share/ueda/RNA004/U87/U87_IVT/20231227_1535_MN32625_FAX73794_2cf3868f/result_filter.vcf"

    classification(input,output,checkpoint_path)

    # input = "/mnt/share/ueda/RNA004/U87/U87_STM2457/U87_STM2457/20240216_1602_MN32625_FAX70218_f0b7f679/unfilter_result.vcf"
    # output = "/mnt/share/ueda/RNA004/U87/U87_STM2457/U87_STM2457/20240216_1602_MN32625_FAX70218_f0b7f679/result_filter.vcf"
    #
    # classification(input,output,checkpoint_path,m5Cpath,psudepath)

#run()