'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Andrew Lonsdale, 2017
License     : MIT
Maintainer  : Andrew Lonsdale
Portability : POSIX

The program reads one or more input FASTQ files. For each file it computes a
variety of statistics, and then prints a summary of the statistics as output.... in emoji.
'''

from __future__ import print_function
from argparse import ArgumentParser
from math import floor
import sys
from Bio import SeqIO
import logging
import pkg_resources
from pyemojify import emojify
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqIO import QualityIO
from . import fastqe_map as emaps # todo make maps illumin 1.9 specific etc


EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
EXIT_FASTA_FILE_ERROR = 3
DEFAULT_MIN_LEN = 0
DEFAULT_VERBOSE = False
#HEADER = 'FILENAME\tNUMSEQ\tTOTAL\tMIN\tAVG\tMAX'
HEADER = '# FASTQE sequence quality for:'

PROGRAM_NAME = "fastqe"


#:no_entry_sign:',
#2–9 6
#'#': ':skull:',
#10–19 15
#'+': ':poop:' ,
#20–24 22
#'5': ':rage:',
#25–29 27
#':': ':neutral_face:',
#30–34 33
#'?': ':smile:',
#35–39 37
#'D': ':sunglasses:',
#≥ 40 40
#'J': ':heart_eyes:',

def print_scale(full_quals,binned):
    count = 0
    if binned:
        print("#scale for fastqe (binned)")

        for i in full_quals:
            print("# ",count,i,emojify(emaps.fastq_emoji_map_binned.get(i,':heart_eyes:')))
            count = count +1
    else:
        print("#scale for fastqe")
        for i in full_quals:
            print("# ",count,i,emojify(emaps.fastq_emoji_map.get(i,':heart_eyes:')))
            count = count +1



try:
    PROGRAM_VERSION = pkg_resources.require(PROGRAM_NAME)[0].version
except pkg_resources.DistributionNotFound:
    PROGRAM_VERSION = "undefined_version"


def exit_with_error(message, exit_status):
    '''Print an error message to stderr, prefixed by the program name and 'ERROR'.
    Then exit program with supplied exit status.

    Arguments:
        message: an error message as a string.
        exit_status: a positive integer representing the exit status of the
            program.
    '''
    logging.error(message)
    print("{} ERROR: {}, exiting".format(PROGRAM_NAME, message), file=sys.stderr)
    sys.exit(exit_status)


def parse_args():
    '''Parse command line arguments.
    Returns Options object with command line argument values as attributes.
    Will exit the program on a command line error.
    '''
    parser = ArgumentParser(description='Read one or more FASTQ files, compute quality stats for each file, print as emoji... for some reason.')
    parser.add_argument(
        '--minlen',
        metavar='N',
        type=int,
        default=DEFAULT_MIN_LEN,
        help='Minimum length sequence to include in stats (default {})'.format(
            DEFAULT_MIN_LEN))
    parser.add_argument('--version',
        action='version',
        version='%(prog)s ' + PROGRAM_VERSION)
    parser.add_argument('--mean',
        default=True,
        action='store_true',
        help='show mean quality per position (DEFAULT)')
    parser.add_argument('--bin',
            action='store_true',
            help='use binned scores')
    parser.add_argument('--min',
        action='store_true',
        help='show minimum quality per position')
    parser.add_argument('--max',
        action='store_true',
        help='show maximum quality per position')
    parser.add_argument('--log',
        metavar='LOG_FILE',
        type=str,
        help='record program progress in LOG_FILE')
    parser.add_argument('--scale',
    action='store_true',
    help='show relevant scale in output')
    parser.add_argument('fasta_files',
        nargs='*',
        metavar='FASTA_FILE',
        type=str,
        help='Input FASTQ files')
    return parser.parse_args()


class FastaStats(object):
    '''Compute various statistics for a FASTA file:

    num_seqs: the number of sequences in the file satisfying the minimum
       length requirement (minlen_threshold).
    num_bases: the total length of all the counted sequences.
    min_len: the minimum length of the counted sequences.
    max_len: the maximum length of the counted sequences.
    average: the average length of the counted sequences rounded down
       to an integer.
    '''
    #pylint: disable=too-many-arguments
    def __init__(self,
                 num_seqs=None,
                 num_bases=None,
                 min_len=None,
                 max_len=None,
                 average=None,
                 quality_scores_mean=None,
                 quality_scores_min=None,
                 quality_scores_max=None):
        "Build an empty FastaStats object"
        self.num_seqs = num_seqs
        self.num_bases = num_bases
        self.min_len = min_len
        self.max_len = max_len
        self.average = average
        self.quality_score_mean = quality_scores_mean
        self.quality_score_min = quality_scores_min
        self.quality_score_max = quality_scores_max


    def __eq__(self, other):
        "Two FastaStats objects are equal iff their attributes are equal"
        if type(other) is type(self):
            return self.__dict__ == other.__dict__
        else:
            return False

    def __repr__(self):
        "Generate a printable representation of a FastaStats object"
        return "FastaStats(num_seqs={}, num_bases={}, min_len={}, max_len={}, " \
            "average={})".format(
                self.num_seqs, self.num_bases, self.min_len, self.max_len,
                self.average)

    def from_file(self, fasta_file, minlen_threshold=DEFAULT_MIN_LEN):
        '''Compute a FastaStats object from an input FASTA file.

        Arguments:
           fasta_file: an open file object for the FASTA file
           minlen_threshold: the minimum length sequence to consider in
              computing the statistics. Sequences in the input FASTA file
              which have a length less than this value are ignored and not
              considered in the resulting statistics.
        Result:
           A FastaStats object
        '''

        #FASTQ addition
        # works for up to 500bp reads
        means = np.zeros(500)
        mins = np.zeros(500)
        maxs = np.zeros(500)

        num_seqs = num_bases = 0
        min_len = max_len = None
        for seq in SeqIO.parse(fasta_file, "fastq"):

            # fastqe

            index=0
            for s in seq.letter_annotations["phred_quality"]:
                #maxs
                if s > maxs[index]:
                    maxs[index] = s

                #mins
                if mins[index] == 0:
                    mins[index] = s
                elif s < mins[index]:
                    mins[index] = s

                #means
                means[index] += s
                index = index+1


            # FASTA stat
            this_len = len(seq.seq)
            if this_len >= minlen_threshold:
                if num_seqs == 0:
                    min_len = max_len = this_len
                else:
                    min_len = min(this_len, min_len)
                    max_len = max(this_len, max_len)
                num_seqs += 1
                num_bases += this_len


        # after processing
        if num_seqs > 0:
            self.average = int(floor(float(num_bases) / num_seqs))
        else:
            self.average = None
        self.num_seqs = num_seqs
        self.num_bases = num_bases
        self.min_len = min_len
        self.max_len = max_len


        #fastq

        # create fake sequence for each type
        cleaned = np.trim_zeros(means)
        means_fp = cleaned/num_seqs
        fake_seq= ''.join(["a"]*len(means_fp.round()))

        record_mean = SeqRecord(Seq(fake_seq), id="test", name="mean scores",
                   description="example with mean fastq socres",
                   letter_annotations={'phred_quality':list(means_fp.round().astype(int))})

        self.quality_scores_mean = record_mean



        mins_trimmed = np.trim_zeros(mins)
        fake_seq_min= ''.join(["a"]*len(mins_trimmed))

        record_mins = SeqRecord(Seq(fake_seq_min), id="test", name="mean scores",
                   description="example with mean fastq socres",
                   letter_annotations={'phred_quality':list(mins_trimmed.astype(int))})

        self.quality_scores_mins= record_mins


        maxs_trimmed = np.trim_zeros(maxs)
        fake_seq_maxs= ''.join(["a"]*len(maxs_trimmed))

        record_maxs= SeqRecord(Seq(fake_seq_maxs), id="test", name="mean scores",
                   description="example with mean fastq socres",
                   letter_annotations={'phred_quality':list(maxs_trimmed.astype(int))})

        self.quality_scores_maxs= record_maxs



        return self

    def pretty(self, filename):
        '''Generate a pretty printable representation of a FastaStats object
        suitable for output of the program. The output is a tab-delimited
        string containing the filename of the input FASTA file followed by
        the attributes of the object. If 0 sequences were read from the FASTA
        file then num_seqs and num_bases are output as 0, and min_len, average
        and max_len are output as a dash "-".

        Arguments:
           filename: the name of the input FASTA file
        Result:
           A string suitable for pretty printed output
        '''



        return(filename)
        if self.num_seqs > 0:
            num_seqs = str(self.num_seqs)
            num_bases = str(self.num_bases)
            min_len = str(self.min_len)
            average = str(self.average)
            max_len = str(self.max_len)

        else:
            num_seqs = num_bases = "0"
            min_len = average = max_len = "-"
        return "\t".join([filename, num_seqs, num_bases, min_len, average,
                          max_len])


def process_files(options):
    '''Compute and print FastaStats for each input FASTA file specified on the
    command line. If no FASTA files are specified on the command line then
    read from the standard input (stdin).

    Arguments:
       options: the command line options of the program
    Result:
       None
    '''
    if options.fasta_files:
        for fasta_filename in options.fasta_files:
            logging.info("Processing FASTA file from {}".format(fasta_filename))
            try:
                fasta_file = open(fasta_filename)
            except IOError as exception:
                exit_with_error(str(exception), EXIT_FILE_IO_ERROR)
            else:
                with fasta_file:
                    stats = FastaStats().from_file(fasta_file, options.minlen)
                    #print(stats.pretty(fasta_filename))

                    if options.scale:
                        print_scale(emaps.all_qualities,options.bin)

                    #rewrite this
                    if options.bin:
                        logging.info("Binned calculations")
                        if options.max:
                            logging.info("Calculate max quality per position")
                            print(stats.pretty(fasta_filename),"max (binned)"," ".join([emojify(emaps.fastq_emoji_map_binned.get(s,':heart_eyes:')) for s in QualityIO._get_sanger_quality_str(stats.quality_scores_maxs)]),sep='\t')
                        logging.info("Calculate mean quality per position")
                        print(stats.pretty(fasta_filename),"mean (binned)"," ".join([emojify(emaps.fastq_emoji_map_binned.get(s,':heart_eyes:')) for s in QualityIO._get_sanger_quality_str(stats.quality_scores_mean)]),sep='\t')
                        if options.min:
                            logging.info("Calculate min quality per position")
                            print(stats.pretty(fasta_filename),"min (binned)"," ".join([emojify(emaps.fastq_emoji_map_binned.get(s,':heart_eyes:')) for s in QualityIO._get_sanger_quality_str(stats.quality_scores_mins)]),sep='\t')
                    else:
                        if options.max:
                            logging.info("Calculate max quality per position")
                            print(stats.pretty(fasta_filename),"max"," ".join([emojify(emaps.fastq_emoji_map[s]) for s in QualityIO._get_sanger_quality_str(stats.quality_scores_maxs)]),sep='\t')
                        logging.info("Calculate mean quality per position")
                        print(stats.pretty(fasta_filename),"mean"," ".join([emojify(emaps.fastq_emoji_map[s]) for s in QualityIO._get_sanger_quality_str(stats.quality_scores_mean)]),sep='\t')
                        if options.min:
                            logging.info("Calculate min quality per position")
                            print(stats.pretty(fasta_filename),"min"," ".join([emojify(emaps.fastq_emoji_map[s]) for s in QualityIO._get_sanger_quality_str(stats.quality_scores_mins)]),sep='\t')



                    #print("MAX:  "," ".join([s for s in QualityIO._get_sanger_quality_str(stats.quality_scores_maxs)]))
                    #print("MEAN: "," ".join([s for s in QualityIO._get_sanger_quality_str(stats.quality_scores_mean)]))
                    #print("MIN:  "," ".join([s for s in QualityIO._get_sanger_quality_str(stats.quality_scores_mins)]))


    else:
        logging.info("Processing FASTA file from stdin")
        stats = FastaStats().from_file(sys.stdin, options.minlen)
        print(stats.pretty("stdin"))


def init_logging(log_filename):
    '''If the log_filename is defined, then
    initialise the logging facility, and write log statement
    indicating the program has started, and also write out the
    command line from sys.argv

    Arguments:
        log_filename: either None, if logging is not required, or the
            string name of the log file to write to
    Result:
        None
    '''
    if log_filename is not None:
        logging.basicConfig(filename=log_filename,
            level=logging.DEBUG,
            filemode='w',
            format='%(asctime)s %(levelname)s - %(message)s',
            datefmt='%m-%d-%Y %H:%M:%S')
        logging.info('program started')
        logging.info('command line: {0}'.format(' '.join(sys.argv)))


def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    init_logging(options.log)
    #print(HEADER)
    process_files(options)


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
