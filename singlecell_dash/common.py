import glob
import os
import datetime

from collections import OrderedDict

import pandas as pd
import numpy as np

from sklearn.manifold import TSNE

import scipy.stats as stats
import scipy.sparse as sparse

from sparse_dataframe import SparseDataFrame


def combine_sdf_files(run_folder, folders, verbose=False, **kwargs):
    """function for concatenating SparseDataFrames together"""
    combined = SparseDataFrame()
    combined.rows = []
    columns = set()

    for folder in folders:
        filename = os.path.join(run_folder, folder, f'{folder}.mus.cell-gene.npz')
        if verbose:
            print(f'Reading {filename} ...')
        sdf = SparseDataFrame(filename)
        columns.add(tuple(sdf.columns))
        combined.rows.extend(sdf.rows)
        if combined.matrix is None:
            combined.matrix = sdf.matrix
        else:
            combined.matrix = sparse.vstack((combined.matrix, sdf.matrix),
                                            format='csr')

    assert len(columns) == 1

    combined.columns = columns.pop()

    return combined


def combine_csv_files(folder, globber, verbose=False, **kwargs):
    """generic function for concatentating a bunch of csv files into a single
    pandas Dataframe"""
    dfs = []

    for filename in glob.iglob(os.path.join(folder, globber)):
        if verbose:
            print(f'Reading {filename} ...')
        df = pd.read_csv(filename, **kwargs)
        dfs.append(df)
    combined = pd.concat(dfs)
    return combined


def maybe_to_numeric(series):
    try:
        return pd.to_numeric(series)
    except ValueError:
        return series


def clean_mapping_stats(mapping_stats_original, convert_to_percentage=None):
    """Remove whitespace from all values and convert to numbers"""

    if convert_to_percentage is None:
        convert_to_percentage = set()

    mapping_stats_original = mapping_stats_original.applymap(
        lambda x: (x.replace(',', '').strip().strip('%')
                   if isinstance(x, str) else x))

    numeric = mapping_stats_original.apply(maybe_to_numeric)

    numeric.columns = numeric.columns.map(str.strip)

    # for 10X mapping stats
    numeric.columns = numeric.columns.map(
            lambda x: ('Percent {}'.format(x.replace('Fraction ', ''))
                       if x in convert_to_percentage else x)
    )

    return numeric



def diff_exp(matrix, group1, group2, index):
    """Computes differential expression between group 1 and group 2
    for each column in the dataframe counts.

    Returns a dataframe of Z-scores and p-values."""

    g1 = matrix[group1, :]
    g2 = matrix[group2, :]

    g1mu = g1.mean(0)
    g2mu = g2.mean(0)

    mean_diff = np.asarray(g1mu - g2mu).flatten()
    # E[X^2] - (E[X])^2
    pooled_sd = np.sqrt(
        ((g1.power(2)).mean(0) - np.power(g1mu, 2)) / len(group1)
        + ((g2.power(2)).mean(0) - np.power(g2mu, 2)) / len(group2))
    pooled_sd = np.asarray(pooled_sd).flatten()

    z_scores = np.zeros_like(pooled_sd)
    nz = pooled_sd > 0
    z_scores[nz] = np.nan_to_num(mean_diff[nz] / pooled_sd[nz])

    # t-test
    p_vals = (1 - stats.norm.cdf(np.abs(z_scores))) * 2

    df = pd.DataFrame(OrderedDict([('z', z_scores), ('p', p_vals)]),
                      index=index)

    return df


class Plates(object):

    # Names of commonly accessed columns
    MEAN_READS_PER_CELL = 'Mean reads per well'
    MEDIAN_GENES_PER_CELL = 'Median genes per well'
    PERCENT_ERCC = 'Percent ERCC'
    PERCENT_MAPPED_READS = 'Percent mapped to genome'

    # maybe we should change this to the right thing
    SAMPLE_MAPPING = 'WELL_MAPPING'

    def __init__(self, data_folder, metadata, genes_to_drop='Rn45s',
                 verbose=False, nrows=None):

        plates_folder = os.path.join(data_folder, 'plates')

        counts = combine_csv_files(
            plates_folder, '*.htseq-count-by-cell.csv',
            index_col=[0, 1, 2, 3], verbose=verbose, nrows=nrows)
        mapping_stats = combine_csv_files(
            plates_folder, '*.log-by-cell.csv',
            index_col=[0, 1, 2, 3], verbose=verbose)
        self.genes, self.cell_metadata, self.mapping_stats = \
            self.clean_and_reformat(counts, mapping_stats)

        self.plate_summaries = self.calculate_plate_summaries()

        original_metadata = pd.read_csv(metadata, index_col=0)
        self.plate_metadata = self.clean_plate_metadata(original_metadata)
        self.plate_metadata = self.plate_metadata.loc[
            self.plate_summaries.index]

        if not os.path.exists(os.path.join(data_folder, 'coords')):
            os.mkdir(os.path.join(data_folder, 'coords'))

        self.bulk_smushed_cache_file = os.path.join(data_folder, 'coords',
                                                    'bulk_smushed.csv')
        self.cell_smushed_cache_file = os.path.join(data_folder, 'coords',
                                                    'cell_smushed.pickle')

        self.bulk_smushed = self.compute_bulk_smushing()
        self.cell_smushed = self.compute_cell_smushing()

        self.gene_names = sorted(self.genes.columns)
        self.plate_metadata_features = sorted(self.plate_metadata.columns)

        # Remove pesky genes
        self.genes = self.genes.drop(genes_to_drop, axis=1)

        # Get a counts per million rescaling of the genes
        self.counts_per_million = self.genes.divide(self.genes.sum(axis=1),
                                                    axis=0) * 1e6
        self.top_genes = self.compute_top_genes_per_cell()

        self.data = {'genes': self.genes,
                     'mapping_stats': self.mapping_stats,
                     'cell_metadata': self.cell_metadata,
                     'plate_metadata': self.plate_metadata,
                     'plate_summaries': self.plate_summaries}

    def __repr__(self):
        n_plates = self.plate_summaries.shape[0]
        n_barcodes = self.genes.shape[0]
        s = f'This is an object holding data for {n_plates} plates and ' \
            f'{n_barcodes} barcodes.\nHere are the accessible dataframes:\n'

        for name, df in self.data.items():
            s += f'\t"{name}" table dimensions: ' + str(df.shape) + '\n'
        return s

    @staticmethod
    def clean_and_reformat(counts, mapping_stats):
        """Move metadata information into separate dataframe and simplify ids

        Parameters
        ----------
        counts : pandas.DataFrame
            A (samples, genes) dataframe of integer number of reads that mapped
            to a gene in a cell, but also has extra columns of ERCC mapping and
            htseq-count output that we want to remove
        mapping_stats : pandas.DataFrame
            A (samples, mapping_statistics) dataframe of the time the alignment
            began, number of input reads, number of mapped reads, and other
            information output by STAR, but everything is a string instead of
            numbers which makes us sad

        Returns
        -------
        genes : pandas.DataFrame
            A (samples, genes) dataframe of integer number of reads that mapped
            to a gene in a cell
        cell_metadata : pandas.DataFrame
            A (samples, sample_features) dataframe of number of detected genes,
            total reads, ercc counts, and "WELL_MAPPING" (really,
            plate mapping)
        mapping_stats : pandas.DataFrame
            A (samples, mapping_statistics) dataframe of the time the alignment
            began, number of input reads, number of mapped reads, and other
            information output by STAR, with numbers properly formatted
        """
        mapping_stats = clean_mapping_stats(mapping_stats)
        cell_metadata = counts.index.to_frame()

        sample_ids = cell_metadata.index.droplevel([1, 2, 3])

        cell_metadata.index = sample_ids
        mapping_stats.index = sample_ids
        counts.index = sample_ids

        # Extract htseq-count outputs and save as separate files
        cols = [x for x in counts if x.startswith('__')]
        count_stats = counts[cols]
        count_stats.columns = [x.strip('_') for x in count_stats]

        # Separate spike-ins (ERCCs) and genes
        ercc_names = [col for col in counts.columns[3:] if 'ERCC-' in col]
        gene_names = [col for col in counts.columns[3:] if
                      'ERCC-' not in col and col[0] != '_']
        cell_metadata['total_reads'] = counts.sum(axis=1)

        # Separate counts of everything from genes-only
        genes = counts[gene_names]

        # Add mapping and ERCC counts to cell metadata
        cell_metadata['n_genes'] = (genes > 0).sum(axis=1)
        cell_metadata['mapped_reads'] = genes.sum(axis=1)
        cell_metadata['ercc'] = counts[ercc_names].sum(axis=1)
        cell_metadata = pd.concat([cell_metadata, count_stats], axis=1)

        # Remove not useful columns
        cell_metadata.drop(['too_low_aQual', 'not_aligned'], inplace=True,
                           axis=1)

        return genes, cell_metadata, mapping_stats

    def calculate_plate_summaries(self):
        """Get mean reads, percent mapping, etc summaries for each plate"""
        well_map = self.cell_metadata.groupby(Plates.SAMPLE_MAPPING)

        # these stats are from STAR mapping
        star_cols = ['Number of input reads', 'Uniquely mapped reads number']
        star_stats = self.mapping_stats[star_cols].groupby(
                self.cell_metadata[Plates.SAMPLE_MAPPING]).sum()

        total_reads = star_stats['Number of input reads']
        unique_reads = star_stats['Uniquely mapped reads number']

        percent_ercc = well_map.sum()['ercc'].divide(total_reads, axis=0)
        percent_mapped_reads = unique_reads / total_reads - percent_ercc

        plate_summaries = pd.DataFrame(OrderedDict([
            (Plates.MEAN_READS_PER_CELL, total_reads / well_map.size()),
            (Plates.MEDIAN_GENES_PER_CELL, well_map.median()['n_genes']),
            ('Percent not uniquely aligned', 100 * well_map.sum()['alignment_not_unique'].divide(total_reads, axis=0)),
            (Plates.PERCENT_MAPPED_READS, 100 * percent_mapped_reads),
            ('Percent no feature', 100 * well_map.sum()['no_feature'].divide(total_reads, axis=0)),
            ('Percent Rn45s', 100 * self.genes['Rn45s'].groupby(
                    self.cell_metadata[Plates.SAMPLE_MAPPING]).sum() / total_reads),
            (Plates.PERCENT_ERCC, 100 * percent_ercc),
            ('n_wells', well_map.size())
        ]))

        return plate_summaries

    @staticmethod
    def clean_plate_metadata(plate_metadata):
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

    def compute_bulk_smushing(self):
        """Get average signal from each plate ('bulk') and find 2d embedding"""

        grouped = self.genes.groupby(self.cell_metadata[self.SAMPLE_MAPPING])

        if os.path.exists(self.bulk_smushed_cache_file):
            smushed = pd.read_csv(self.bulk_smushed_cache_file, names=[0, 1],
                                  header=0, index_col=0)
            # if the set of plates hasn't changed, return the cached version
            if set(grouped.groups) == set(smushed.index):
                return smushed

        # if the cache was missing or invalid, compute a new projection
        medians = grouped.median()
        smusher = TSNE(random_state=0, perplexity=10, metric='cosine')
        smushed = pd.DataFrame(smusher.fit_transform(medians),
                               index=medians.index)

        smushed.to_csv(self.bulk_smushed_cache_file)

        return smushed

    def compute_cell_smushing(self):
        """Within each plate, find a 2d embedding of all cells"""
        grouped = self.genes.groupby(self.cell_metadata[self.SAMPLE_MAPPING])

        if os.path.exists(self.cell_smushed_cache_file):
            smusheds = pd.read_pickle(self.cell_smushed_cache_file)
            # if nothing is missing, return the cached version
            if not set(grouped.groups) - set(smusheds):
                return smusheds
        else:
            smusheds = {}

        for plate_name, genes_subset in grouped:
            if plate_name not in smusheds:
                cell_smusher = TSNE(metric='cosine', random_state=0)
                cell_smushed = pd.DataFrame(
                    cell_smusher.fit_transform(genes_subset),
                    index=genes_subset.index)
                smusheds[plate_name] = cell_smushed

        pd.to_pickle(smusheds, self.cell_smushed_cache_file)

        return smusheds

    def compute_top_genes_per_cell(self):
        """Get the most highly expressed genes in every cell

        Returns
        -------
        top_genes : pandas.Series
            A mapping of the cell barcode to a ranked list of the top 10 genes,
            where the first item is the most highly expressed (e.g. Rn45s)
        """
        ranks = self.genes.rank(axis=1, ascending=False)
        in_top10 = ranks[ranks <= 10]
        top_genes = in_top10.apply(
            lambda x: x.sort_values().dropna().index.tolist(), axis=1)
        return top_genes


class TenX_Runs(Plates):

    # Names of commonly accessed columns
    MEAN_READS_PER_CELL = 'Mean Reads per Cell'
    MEDIAN_GENES_PER_CELL = 'Median Genes per Cell'
    PERCENT_MAPPED_READS = 'Percent Reads Mapped Confidently to Transcriptome'

    SAMPLE_MAPPING = 'CHANNEL_MAPPING'

    COLUMNS_TO_CONVERT = {'Valid Barcodes',
                          'Reads Mapped Confidently to Transcriptome',
                          'Reads Mapped Confidently to Exonic Regions',
                          'Reads Mapped Confidently to Intronic Regions',
                          'Reads Mapped Confidently to Intergenic Regions',
                          'Reads Mapped Antisense to Gene',
                          'Sequencing Saturation',
                          'Q30 Bases in Barcode', 'Q30 Bases in RNA Read',
                          'Q30 Bases in Sample Index', 'Q30 Bases in UMI',
                          'Fraction Reads in Cells'}

    def __init__(self, data_folder, genes_to_drop='Rn45s',
                 verbose=False, nrows=None, tissue=None,
                 channels_to_use=None, tissue_folder='tissues'):

        run_folder = os.path.join(data_folder, '10x_data')

        self.plate_metadata = combine_csv_files(run_folder,
                                                'MACA_10X_P*.csv',
                                                index_col=0, nrows=nrows)

        if tissue is not None:
            tissues = tissue.split(',')
            folders = self.plate_metadata.index[self.plate_metadata['Tissue'].isin(tissues)]

        else:
            folders = self.plate_metadata.index

        folders = [f for f in folders if os.path.exists(os.path.join(run_folder, f))]

        if channels_to_use is not None:
            folders = [f for f in folders if f in channels_to_use]

        counts = combine_sdf_files(run_folder, folders,
                                   verbose=verbose)

        mapping_stats = self.combine_metrics_files(
                run_folder, folders)

        self.genes, self.cell_metadata, self.mapping_stats = \
            self.clean_and_reformat(counts, mapping_stats)

        self.plate_summaries = self.calculate_plate_summaries()

        self.plate_metadata = self.plate_metadata.loc[
            self.plate_summaries.index]

        self.cell_metadata = self.cell_metadata.join(self.plate_metadata,
                                                     on=self.SAMPLE_MAPPING)

        smushed_folder = os.path.join(run_folder, tissue_folder)

        if not os.path.exists(smushed_folder):
            os.mkdir(smushed_folder)

        self.cell_smushed = self.read_tissue_smushed(smushed_folder, verbose)

        self.gene_names = sorted(self.genes.columns)
        self.plate_metadata_features = sorted(self.plate_metadata.columns)

        # Remove pesky genes
        self.genes = self.genes.drop(genes_to_drop)

        # Get a counts per million rescaling of the genes
        # self.counts_per_million = self.genes.divide(self.genes.sum(axis=1),
        #                                             axis=0) * 1e6
        # self.top_genes = self.compute_top_genes_per_cell()

        self.data = {'genes': self.genes,
                     'mapping_stats': self.mapping_stats,
                     'cell_metadata': self.cell_metadata,
                     'plate_metadata': self.plate_metadata,
                     'plate_summaries': self.plate_summaries}

    def __repr__(self):
        n_channels = self.plate_summaries.shape[0]
        n_barcodes = len(self.genes.rows)
        s = f'This is an object holding data for {n_channels} 10X channels and ' \
            f'{n_barcodes} barcodes.\nHere are the accessible dataframes:\n'

        for name, df in self.data.items():
            s += f'\t"{name}" table dimensions: ' + str(df.shape) + '\n'
        return s

    @staticmethod
    def combine_cell_files(folder, globber, verbose=False):
        dfs = []

        for filename in glob.iglob(os.path.join(folder, globber)):
            if verbose:
                print(f'Reading {filename} ...')

            channel = os.path.basename(os.path.dirname(filename))

            df = pd.read_csv(filename, index_col=0)
            df.index = pd.MultiIndex.from_product(([channel], df.index),
                                                  names=['channel', 'cell_id'])
            dfs.append(df)

        combined = pd.concat(dfs)
        return combined

    @staticmethod
    def combine_metrics_files(run_folder, folders):
        dfs = []

        for folder in folders:
            filename = os.path.join(run_folder, folder, 'metrics_summary.csv')
            p_name = os.path.basename(os.path.dirname(filename))
            df = pd.read_csv(filename)
            df[TenX_Runs.SAMPLE_MAPPING] = p_name
            dfs.append(df)
        combined = pd.concat(dfs)
        combined.set_index(TenX_Runs.SAMPLE_MAPPING, inplace=True)
        return combined

    @staticmethod
    def clean_and_reformat(counts, mapping_stats):
        """Move metadata information into separate dataframe and simplify ids

        Parameters
        ----------
        counts : pandas.DataFrame
            A (samples, genes) dataframe of integer number of reads that mapped
            to a gene in a cell, but also has extra columns of ERCC mapping and
            htseq-count output that we want to remove
        mapping_stats : pandas.DataFrame
            A (samples, mapping_statistics) dataframe of the time the alignment
            began, number of input reads, number of mapped reads, and other
            information output by STAR, but everything is a string instead of
            numbers which makes us sad

        Returns
        -------
        genes : SparseDataFrame
            A (samples, genes) dataframe of integer number of reads that mapped
            to a gene in a cell
        cell_metadata : pandas.DataFrame
            A (samples, sample_features) dataframe of number of detected genes,
            total reads, ercc counts, and "WELL_MAPPING" (really,
            plate mapping)
        mapping_stats : pandas.DataFrame
            A (samples, mapping_statistics) dataframe of the time the alignment
            began, number of input reads, number of mapped reads, and other
            information output by CellRanger, with numbers properly formatted
        """
        # counts.sort_index(inplace=True)
        # channel_ids = counts.index.get_level_values(0)
        channel_ids = [c.rsplit('_', 1)[0] for c in counts.rows]

        mapping_stats = clean_mapping_stats(
                mapping_stats,
                convert_to_percentage=TenX_Runs.COLUMNS_TO_CONVERT
        )

        sample_ids = pd.Series(counts.rows)
        #         '{}_{}'.format(channel, index) for channel, index in
        #         counts.index
        # )

        cell_metadata = pd.DataFrame(
                index=sample_ids,
                data={TenX_Runs.SAMPLE_MAPPING: channel_ids}
        )

        counts.index = sample_ids

        # Separate spike-ins (ERCCs) and genes
        ercc_names = [col for col in counts.columns if col.startswith('ERCC-')]
        gene_names = [col for col in counts.columns if
                      not (col.startswith('ERCC-')
                           or col.endswith('_transgene'))]

        # Separate counts of everything from genes-only
        genes = SparseDataFrame()
        genes.matrix = counts[gene_names]
        genes.columns = gene_names[:]
        genes.rows = counts.rows[:]

        # Add mapping and ERCC counts to cell metadata
        cell_metadata['total_reads'] = counts.matrix.sum(axis=1)
        cell_metadata['n_genes'] = (genes.matrix > 0).sum(axis=1)
        cell_metadata['mapped_reads'] = genes.matrix.sum(axis=1)
        cell_metadata['ercc'] = counts[ercc_names].sum(axis=1)

        return genes, cell_metadata, mapping_stats

    def calculate_plate_summaries(self):
        """Get mean reads, percent mapping, etc summaries for each plate"""
        channel_map = self.cell_metadata.groupby(TenX_Runs.SAMPLE_MAPPING)

        total_reads = self.mapping_stats['Number of Reads']

        # percent_rn45s = pd.Series(self.genes['Rn45s'].todense()).groupby(
        #         self.cell_metadata[TenX_Runs.SAMPLE_MAPPING]
        # ).sum() / total_reads

        percent_ercc = channel_map['ercc'].sum().divide(total_reads, axis=0)

        plate_summaries = pd.concat(
                [self.mapping_stats,
                 pd.DataFrame(OrderedDict([
                     # ('Percent Rn45s', percent_rn45s),
                     (TenX_Runs.PERCENT_ERCC, percent_ercc),
                     ('n_barcodes', channel_map.size())
                 ]))], axis=1
        )

        return plate_summaries

    def read_tissue_smushed(self, folder, verbose=False):
        smusheds = {}
        glob.glob(os.path.join(folder, 'smushed-*'))
        for filename in glob.iglob(os.path.join(folder, 'smushed-*')):
            if verbose:
                print(f'Reading {filename} ...')
            tissue = filename.split('smushed-')[-1].split('.')[0]
            tissue = tissue.split('-')[0]
            df = pd.read_csv(filename, index_col=0)
            df.rename(columns={'0': 0, '1': 1}, inplace=True)
            smusheds[tissue] = df
            assert len(df.columns.difference([0, 1, 'cluster'])) == 0
        return smusheds
