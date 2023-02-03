# -*- coding: utf-8 -*-

""" Methods that concern processed files """
import logging
import os.path
import glob

import params


class ProcessingResults(object):
    """
    Container to hold the processing output folder and create individual objects based on the output file.
    """

    def __init__(self,
                 logger,
                 processed_files_dir,

                 ):
        """
        :param processed_files_folder: directory containing the processed output files
        :param logger: logger object
        """

        self.logger = logger

        self.dir = processed_files_dir

        list_of_file_names = [file.stem for file in self.dir.iterdir() if 'census-out' in file.stem]

        #self.rmats_mxe = self._read_rmats_mxe()
        self.list_of_file_names = list_of_file_names




    # def _read_rmats_mxe(self):
    #     """
    #     Read input data frame for Mutually Exclusive Exons (MXE)
    #     The slices should be anc-alt1-down and anc-alt2-down
    #     :return:
    #     """
    #
    #     df = pd.read_table(os.path.join(self.dir, 'MXE.MATS.JC.txt'), sep='\t', low_memory=False)
    #     df.columns = ['id', 'gene_id', 'gene_symbol', 'chr', 'strand',
    #                   'alt1_es', 'alt1_ee', 'alt2_es', 'alt2_ee', 'anc_es',
    #                   'anc_ee', 'down_es', 'down_ee', 'id0', 'ijc_s1', 'sjc_s1',
    #                   'ijc_s2', 'sjc_s2', 'inc_len', 'skp_len', 'p', 'fdr',
    #                   'inc_s1', 'inc_s2', 'inc_dif']
    #     df['jxn_type'] = 'MXE'
    #
    #     return df
    #
