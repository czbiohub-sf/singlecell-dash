import glob
import os
import datetime
import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
import scipy.stats as stats



def combine_cell_files(folder, globber, verbose=False, **kwargs):
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


def clean_mapping_stats(mapping_stats_original):
    """Remove whitespace from all values and convert to numbers"""
    mapping_stats_original = mapping_stats_original.applymap(
        lambda x: x.strip().strip('%') if isinstance(x, str) else x)
    numeric = mapping_stats_original.apply(maybe_to_numeric)
    return numeric

def diff_exp(counts, group1, group2):
    """Computes differential expression between group 1 and group 2
    for each column in the dataframe counts.

    Returns a dataframe of Z-scores and p-values."""

    mean_diff = counts.loc[group1].mean() - counts.loc[group2].mean()
    pooled_sd = np.sqrt(counts.loc[group1].var()/len(group1)
                        + counts.loc[group2].var()/len(group2))
    z_scores = mean_diff/pooled_sd
    z_scores = z_scores.fillna(0)

    # t-test
    p_vals = (1 - stats.norm.cdf(np.abs(z_scores)))*2

    df = pd.DataFrame({'z': z_scores})
    df['p'] = p_vals

    return df


class Plates(object):

    # Names of commonly accessed columns
    MEAN_READS_PER_CELL = 'mean_reads_per_cell'
    MEDIAN_GENES_PER_CELL = 'median_n_genes_per_cell'
    PERCENT_ERCC = 'percent_ercc'
    PERCENT_MAPPED_READS = 'percent_mapped_reads'

    def __init__(self, data_folder, metadata, verbose=False):

        plates_folder = os.path.join(data_folder, 'plates')

        counts = combine_cell_files(
            plates_folder, '*.htseq-count-by-cell.csv',
            index_col=[0, 1, 2, 3], verbose=verbose)
        mapping_stats = combine_cell_files(
            plates_folder, '*.log-by-cell.csv',
            index_col=[0, 1, 2, 3], verbose=verbose)
        self.genes, self.cell_metadata, self.mapping_stats = \
            self.clean_and_reformat(counts, mapping_stats)

        # Get a counts per million rescaling of the genes
        self.counts_per_million = self.genes.divide(self.genes.sum(axis=1),
                                                    axis=0) * 1e6

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
        means = self.cell_metadata.groupby('WELL_MAPPING').mean()
        means.columns = [f'mean_{x}' for x in means.columns]

        medians = self.cell_metadata.groupby('WELL_MAPPING').median()
        medians.columns = [f'median_{x}' for x in medians.columns]

        rn45s = self.genes['Rn45s'].groupby(
            self.cell_metadata['WELL_MAPPING']).mean()
        means_with_rn45s = pd.concat([means, rn45s], axis=1)

        percents = 100 * means_with_rn45s.divide(means['mean_total_reads'],
                                                 axis=0)
        percents.columns = [x.replace('mean', 'percent') for x in percents]
        percents = percents.rename(columns={'Rn45s': 'percent_Rn45s'})
        plate_summaries = pd.concat([means, medians, percents], axis=1)
        cols = ['mean_total_reads', 'median_n_genes', 'percent_mapped_reads',
                'percent_ercc',
                'percent_alignment_not_unique', 'percent_no_feature',
                'percent_Rn45s']
        plate_summaries = plate_summaries[cols]
        plate_summaries['n_cells'] = self.genes.groupby(
            self.cell_metadata['WELL_MAPPING']).size()

        # Rename columns to specify this is per cell fo each plate
        plate_summaries = plate_summaries.rename(
            columns={'mean_total_reads': 'mean_reads_per_cell',
                     'median_n_genes': 'median_n_genes_per_cell'})
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
            columns={'21_55_F': 'Sample ID'})

        plate_metadata['Age (months)'] = plate_metadata['Sample ID'].map(
            lambda x: x.split('_')[0] if isinstance(x, str) else '')

        def parse_date(x):
            if isinstance(x, str):
                x = x.strip()
                if not x:
                    return np.nan
                if x.endswith('/2017'):
                    return datetime.datetime.strptime(x,'%m/%d/%Y')
                elif x.endswith('/17'):
                    return datetime.datetime.strptime(x,'%m/%d/%y')
                else:
                    return datetime.datetime.strptime(x, '%y%m%d')
            elif isinstance(x,float):
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

        grouped = self.genes.groupby(self.cell_metadata['WELL_MAPPING'])

        if os.path.exists(self.bulk_smushed_cache_file):
            smushed = pd.read_csv(self.bulk_smushed_cache_file, names=[0, 1],
                                  header=0, index_col=0)
            # if the set of plates hasn't changed, return the cached version
            if set(grouped.groups) == set(smushed.index):
                return smushed

        # if the cache was missing or invalid, compute a new projection
        correlations = grouped.mean().T.rank().corr()
        smusher = TSNE(random_state=0)
        smushed = pd.DataFrame(smusher.fit_transform(correlations),
                               index=correlations.index)

        smushed.to_csv(self.bulk_smushed_cache_file)

        return smushed

    def compute_cell_smushing(self):
        """Within each plate, find a 2d embedding of all cells"""
        grouped = self.genes.groupby(self.cell_metadata['WELL_MAPPING'])

        if os.path.exists(self.cell_smushed_cache_file):
            smusheds = pd.read_pickle(self.cell_smushed_cache_file)
            # if nothing is missing, return the cached version
            if not set(grouped.groups) - set(smusheds):
                return smusheds
        else:
            smusheds = {}

        for plate_name, genes_subset in grouped:
            if plate_name not in smusheds:
                cell_correlations = genes_subset.T.rank().corr()
                cell_smusher = TSNE(metric='cosine', random_state=0)
                cell_smushed = pd.DataFrame(
                    cell_smusher.fit_transform(cell_correlations),
                    index=cell_correlations.index)
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
