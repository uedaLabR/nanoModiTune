import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from sklearn.model_selection import train_test_split
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.regularizers import L1L2
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras.utils import to_categorical

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

    base_mapping = {'A': 0, 'T': 1, 'C': 2, 'G': 3,'N':0}
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

def getModel():

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

    return model