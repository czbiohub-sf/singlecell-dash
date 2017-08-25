import datetime

import click
import numpy as np
import pandas as pd


def clean_plate_metadata(plate_metadata):
    """Fix column naming and date formatting for plate metadata"""

    # Remove whitespace from "tissue" column
    plate_metadata.tissue = plate_metadata.tissue.map(
        lambda x: x.strip() if isinstance(x, str) else x)

    # Add a column with both tissue and subtissue
    cleaned_subtissue = plate_metadata['subtissue'].map(
        lambda x: ': ' + x.strip() if isinstance(x, str) else '')
    plate_metadata['tissue_subtissue'] = plate_metadata['tissue'] \
                                         + cleaned_subtissue

    # Hard-coded column name of 21_55_F is actually the sample id column
    plate_metadata = plate_metadata.rename(
        columns={'mouse.id': 'Sample ID'})

    plate_metadata['Age (months)'] = plate_metadata['Sample ID'].map(
        lambda x: x.split('_')[0] if isinstance(x, str) else '')

    def parse_date(x):
        if isinstance(x, str):
            x = x.strip()
            if not x:
                return np.nan
            if x.endswith('/2017'):
                return datetime.datetime.strptime(x, '%m/%d/%Y')
            elif x.endswith('/17'):
                return datetime.datetime.strptime(x, '%m/%d/%y')
            else:
                return datetime.datetime.strptime(x, '%y%m%d')
        elif isinstance(x, float):
            return datetime.datetime.strptime(str(int(x)), '%y%m%d')
        else:
            raise TypeError

    for col in plate_metadata.columns:
        if 'date' in col.lower():
            plate_metadata[col] = plate_metadata[col].map(
                parse_date,
                na_action='ignore'
            )

    # Use only the metadata for the plates that have been sequenced
    plate_metadata = plate_metadata.dropna(how='all', axis=1)
    return plate_metadata


def read_and_clean_plate_metadata(filename):
    original_metadata = pd.read_csv(filename, index_col=0)
    cleaned_metadata = clean_plate_metadata(original_metadata)
    return cleaned_metadata
