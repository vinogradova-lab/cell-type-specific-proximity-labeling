# -*- coding: utf-8 -*-

"""jcast.main: Main function."""

import os
import argparse
from pathlib import Path
import params
from datacollection import ProcessingResults
from logger import get_logger

def main():
    """
    Main function for tiap
    :param args: parsed arguments
    :return: None
    """
    args = parse_arguments()

    # ---- Main logger setup ----
    logger = get_logger('tiap', args.out)
    logger.info(args)

    # ---- Read in input files ----
    print("hi")

    processing_results = ProcessingResults(logger=logger,
                                           processed_files_dir=Path(args.processed_files_folder),
                                           )

    print(processing_results.list_of_file_names)                                       
    # TODO: only read in the splice types that are needed

    # ---- Read the gtf file using  gtfpase then write as a pandas data frame. ----
    #gtf = ReadAnnotations(logger=logger,
    #                      path=args.gtf_file,
#                          )

    # ---- Read genome file into memory ----
#    genome = ReadGenome(logger=logger,
#                        path=args.genome,
#                        )
    return True



def parse_arguments():
    """ running main with parsed arguments from command line """

    import sys
    parser = argparse.ArgumentParser(description='TIAP - Turbo ID analysis pipeline')
    parser.add_argument('processed_files_folder',
                        help='path to folder storing input',
                        )
    parser.add_argument('TP/FP file',
                        help='path to TP/FP file',
                        type=argparse.FileType('r'),
                        )
    parser.add_argument('-o', '--out',
                        help='name of the output folder',
                        default='out')

    # print help message if no arguments are given
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    # parse all the arguments
    args = parser.parse_args()
    return(args)



if __name__ == '__main__':
    main()
