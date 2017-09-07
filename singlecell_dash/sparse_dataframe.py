import os

import scipy.sparse

import pandas as pd
import pickle


class SparseDataFrame(object):
    """Custom built sparse DataFrame-like object for singlecell data"""

    def __init__(self, data_or_filename=None):
        self.matrix = None
        self.rows = None
        self.columns = None

        if isinstance(data_or_filename, pd.DataFrame):
            self.matrix, self.rows, self.columns = self.load_from_dataframe(
                    data_or_filename
            )
        elif isinstance(data_or_filename, str):
            if os.path.exists(data_or_filename):
                if data_or_filename.endswith('.csv'):
                    self.matrix, self.rows, self.columns = self.load_from_csv(
                            data_or_filename
                    )
                elif data_or_filename.endswith('.npz'):
                    self.matrix, self.rows, self.columns = self.load_from_npz(
                            data_or_filename
                    )
            else:
                raise ValueError(f'File {data_or_filename} not found')
        elif data_or_filename is not None:
            raise ValueError('Invalid value given for data_or_filename')

    def __repr__(self):
        n_cols = len(self.columns) if self.columns else 0
        n_rows = len(self.rows) if self.rows else 0

        return f'SparseDataFrame with {n_rows} rows and {n_cols} columns'

    @staticmethod
    def _strip_npz(filename):
        """Strip .npz extension if it is present, otherwise do nothing"""
        return filename[:-4] if filename.lower().endswith('.npz') else filename

    @staticmethod
    def load_from_dataframe(df):
        """Convert a dense DataFrame to sparse form"""
        matrix = scipy.sparse.csr_matrix(df.as_matrix(), dtype='int32')
        rows = [str(i) for i in df.index]
        columns = [str(i) for i in df.columns]
        return matrix, rows, columns

    @staticmethod
    def load_from_csv(filename, **df_kwargs):
        """Load a (dense) csv and convert to sparse format"""
        with open(filename, 'rb') as fp:
            df = pd.read_csv(fp, **df_kwargs)

        return self.load_from_dataframe(df)

    @staticmethod
    def load_from_npz(filename):
        matrix = scipy.sparse.load_npz(filename).tocsr()
        with open(SparseDataFrame._strip_npz(filename) + '.p', 'rb') as fp:
            d = pickle.load(fp)
            rows = d["rows"]
            columns = d["columns"]

        return matrix, rows, columns

    def save(self, filename):
        scipy.sparse.save_npz(filename, self.matrix)
        with open(self._strip_npz(filename) + '.p', 'wb') as fp:
            pickle.dump({"rows": self.rows, "columns": self.columns}, fp)

    def load(self, filename):
        self.matrix, self.rows, self.columns = self.load_from_npz(filename)

    def drop(self, columns, inplace=False):
        if isinstance(columns, str):
            drop_set = {columns}
        else:
            drop_set = set(columns)

        matrix = self[:, [c for c in self.columns
                          if c not in drop_set]]
        columns = [c for c in self.columns if c not in drop_set]

        if inplace:
            self.matrix, self.columns = matrix, columns
        else:
            sdf = SparseDataFrame()

            sdf.matrix = matrix
            sdf.columns = columns
            sdf.rows = self.rows[:]

            return sdf

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
            return row, col

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

    def __setitem__(self, item, value):
        item = self._convert_item(item)
        self.matrix[item] = value

    @property
    def shape(self):
        return self.matrix.shape
