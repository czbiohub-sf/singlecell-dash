from scipy.sparse import *
import pandas as pd
import pickle

class sparse_df(object):


    def __init__(self):
        self.matrix = None
        self.rows = None
        self.columns = None

    def load_from_csv(self, filename):
        with open(filename, 'rb') as fp:
            df = pd.read_csv(fp, index_col = 0)
            self.load_from_dataframe(df)

    def load_from_dataframe(self, df):
        self.matrix = csc_matrix(df.as_matrix(), dtype='int32')
        self.rows = [str(i) for i in df.index]
        self.columns = [str(i) for i in df.columns]

    def save(self, filename):
        save_npz(filename + '.npz', self.matrix)
        with open(filename + '.p', 'wb') as fp:
            pickle.dump({"rows": self.rows, "columns": self.columns}, fp)

    def load(self, filename):
        self.matrix = load_npz(filename + '.npz')
        with open(filename + '.p', 'rb') as fp:
            d = pickle.load(fp)
            self.rows = d["rows"]
            self.columns = d["columns"]

    def _get_row(self, row):
        if isinstance(row, str):
            return self.rows.index(row)
        else:
            return row

    def _get_column(self, col):
        if isinstance(col, str):
            return self.columns.index(col)
        else:
            return col

    def _convert_item(self, item):

        if isinstance(item, tuple):
            row, col = item
            if isinstance(row, list):
                row = [self._get_row(i) for i in row]
            else:
                row = self._get_row(row)
            if isinstance(col, list):
                col = [self._get_column(i) for i in col]
            else:
                col = self._get_column(col)
            return (row, col)

        if isinstance(item, str):
            col = self._get_column(item)
            return slice(None, None, None), col
        elif isinstance(item, list):
            cols = [self._get_column(i) for i in item]
            return slice(None, None, None), cols
        else:
            col = item
            return slice(None, None, None), col

    def __getitem__(self, item):
        item = self._convert_item(item)
        return self.matrix[item]

    def __setitem__(self,item,value):
        item = self._convert_item(item)
        self.matrix[item] = value