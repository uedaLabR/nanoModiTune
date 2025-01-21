import pandas as pd
import pyarrow.parquet as pq
import numpy as np
class DBAnnotator:

    def __init__(self, dirpath):

        self.dirpath = dirpath

    # def toDict(self,parquet_file_path):
    #
    #     df = pd.read_parquet(parquet_file_path)
    #     dictionary = pd.Series(df[df.columns[1]].values, index=df[df.columns[0]]).to_dict()
    #     return dictionary


    def lodDict(self,chr):

        filepath_strand = self.dirpath + '/' + chr + '_tc.parquet'
        filepath_anti_strand = self.dirpath + '/' + chr + '_ag.parquet'

        table1 = pq.read_table(filepath_strand)
        table2 = pq.read_table(filepath_anti_strand)

        column_keys = table1.column(0).to_pandas().astype(int)
        column_values = table1.column(1).to_pandas()
        column_keys2 = table2.column(0).to_pandas().astype(int)
        column_values2 = table2.column(1).to_pandas()

        combined_keys = np.concatenate([column_keys, column_keys2])
        combined_values = np.concatenate([column_values, column_values2])

        dictionary = dict(zip(combined_keys, combined_values))

        print("end loading " +str(len(dictionary)))

        self.dictionary = dictionary

    def isSNP(self,pos):

        return self.dictionary.get(pos,None)


