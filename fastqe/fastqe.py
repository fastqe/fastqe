#!/usr/bin/env python

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
from argparse import ArgumentParser, FileType, Namespace
from math import floor
import sys
from Bio import SeqIO
import logging
from pyemojify import emojify
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import QualityIO
from . import fastqe_map as emaps # todo make maps illumin 1.9 specific etc
import os
import gzip
import ast
import binascii

from jinja2 import Environment, FileSystemLoader



EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
EXIT_FASTQ_FILE_ERROR = 3
DEFAULT_MIN_LEN = 0
DEFAULT_VERBOSE = False
HEADER = 'Filename\tStatistic\tQualities'
#HEADER = 'FILENAME\tNUMSEQ\tTOTAL\tMIN\tAVG\tMAX\tQUALITY'
#HEADER = '# FASTQE sequence quality for:'
DEFAULT_READ_LENGTH = 500
PROGRAM_NAME = "fastqe"

CASE_NEITHER = 0
CASE_MIN = 1
CASE_MAX = 2
CASE_BOTH = 3
CASE_HTML = 4
CASE_HTML_ESCAPE = 5

def print_scale(full_quals,mapping_dict,binned):
    count = 0
    print("#scale for fastqe")
    for i in full_quals:
        print("# ",count,i,emojify(mapping_dict.get(i,':heart_eyes:')))
        count = count +1



PROGRAM_VERSION = "0.5.2"


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
    parser = ArgumentParser(description='Read one or more FASTQ files, compute quality stats for each file, print as emoji... for some reason.'+emojify(":smile:"))
    parser.add_argument(
        '--minlen',
        metavar='N',
        type=int,
        default=DEFAULT_MIN_LEN,
        help='Minimum length sequence to include in stats (default {})'.format(
            DEFAULT_MIN_LEN))
    parser.add_argument('--scale',
                        action='store_true',
                        help='show relevant scale in output')
    parser.add_argument('--version',
        action='version',
        version='%(prog)s ' + PROGRAM_VERSION)
    parser.add_argument('--mean',
        default=True,
        action='store_true',
        help='show mean quality per position (DEFAULT)')
    parser.add_argument('--custom',
                        metavar='CUSTOM_DICT',
                        type=str,
                        help='use a mapping of custom emoji to quality in CUSTOM_DICT ('+emojify(":snake:")+emojify(":palm_tree:")+')')
    parser.add_argument('--bin',
                        action='store_true',
                        help='use binned scores ('+emojify(":no_entry_sign:")+emojify(":skull:")
                             +emojify(":poop:")+emojify(":warning:")+" "+emojify(":smile:")+emojify(":laughing:")+emojify(":sunglasses:")+emojify(":heart_eyes:")+")")
    parser.add_argument('--noemoji',
                        action='store_true',
                        help='use mapping without emoji (▁▂▃▄▅▆▇█)')
    parser.add_argument('--noheader',
        action='store_true',
        help='Hide the header before sample output')
    parser.add_argument('--html',
                        metavar='HTML_FILE',
                        type=FileType('w'),
                        help='output an additional HTML report in HTML_FILE')
    parser.add_argument(
        '--window',
        metavar='W',
        type=int,
        default=1,
        help='Window length to summarise reads in HTML report (default 1)')
    parser.add_argument('--html_escape',
                        action='store_true',
                        help='escape html within output, e.g. for Galaxy parsing')
    parser.add_argument('--min',
                        action='store_true',
                        help='show minimum quality per position')
    parser.add_argument('--max',
        action='store_true',
        help='show maximum quality per position')
    parser.add_argument('--output',
                        metavar='OUTPUT_FILE',
                        type=FileType('w'),
                        help = 'write output to OUTPUT_FILE instead of stdout'
                        )
    parser.add_argument('--long',
                        metavar='READ_LENGTH',
                        type=int,
                        help='set initial arrays to be READ_LENGTH bp for long reads')

    parser.add_argument('--log',
                        metavar='LOG_FILE',
                        type=str,
                        help='record program progress in LOG_FILE')
    parser.add_argument('fastq_files',
        nargs='*',
        metavar='FASTQ_FILE',
        type=str,
        help='Input FASTQ files')
    return parser.parse_args()


class FastqStats(object):
    '''Compute various statistics for a FASTQ file:

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
                 counts_per_position=None,
                 quality_scores_mean=None,
                 quality_scores_min=None,
                 quality_scores_max=None):
        "Build an empty FastqStats object"
        self.num_seqs = num_seqs
        self.num_bases = num_bases
        self.min_len = min_len
        self.max_len = max_len
        self.average = average
        self.counts_per_position =  counts_per_position
        self.quality_score_mean = quality_scores_mean
        self.quality_score_min = quality_scores_min
        self.quality_score_max = quality_scores_max


    def __eq__(self, other):
        "Two FastqStats objects are equal iff their attributes are equal"
        if type(other) is type(self):
            return self.__dict__ == other.__dict__
        else:
            return False

    def __repr__(self):
        "Generate a printable representation of a FastqStats object"
        return "FastqStats(num_seqs={}, num_bases={}, min_len={}, max_len={}, " \
            "average={})".format(
                self.num_seqs, self.num_bases, self.min_len, self.max_len,
                self.average)

    def from_file(self, fastq_file, read_size, minlen_threshold=DEFAULT_MIN_LEN):
        '''Compute a FastqStats object from an input FASTQ file.

        Arguments:
           fastq_file: an open file object for the FASTQ file
           minlen_threshold: the minimum length sequence to consider in
              computing the statistics. Sequences in the input FASTQ file
              which have a length less than this value are ignored and not
              considered in the resulting statistics.
        Result:
           A FastqStats object
        '''

        #FASTQ addition
        # works for up to 500bp reads
        means = np.zeros(read_size)
        mins = np.full(read_size,np.inf) # set minimum to inf so any value observed is smaller
        maxs = np.zeros(read_size)
        counts = np.zeros(read_size) # how often this position has a value - for long reads

        num_seqs = num_bases = 0
        min_len = max_len = None
        for seq in SeqIO.parse(fastq_file, "fastq"):

            # fastqe

            index=0
            for s in seq.letter_annotations["phred_quality"]:

                #assert(s < read_size) - this would never work
                # check if read is longer than the allocated statistics array(s), and as effiently as possible,
                # extend the arrays by the read size (default is 500bp)
                if (index >= read_size):
                # np.concatenate((a, b), axis=None)
                    means = np.concatenate((means,np.zeros(1)),axis=None)
                    mins = np.concatenate((mins,np.zeros(1)),axis=None)
                    maxs = np.concatenate((maxs,np.zeros(1)),axis=None)
                    # how often this position has a value - for long reads
                    counts = np.concatenate((counts,np.zeros(1)),axis=None)

                #maxs
                if s > maxs[index]:
                    maxs[index] = s

                #mins
                if mins[index] == 0:
                    mins[index] = s
                elif s < mins[index]:
                    mins[index] = s

                #counts
                counts[index] += 1

                #means
                means[index] += s
                index = index+1



            # FASTQ stat
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

        counts_cleaned = np.trim_zeros(counts, 'b' )
        self.counts_per_position = counts_cleaned


        # create fake sequence for each type
        cleaned = np.trim_zeros(means, 'b' )
        #means_fp = cleaned/num_seqs
        means_fp = cleaned/counts_cleaned # use counts to fix long read where not all seq same length
        fake_seq= ''.join(["a"]*len(means_fp.round()))

        record_mean = SeqRecord(Seq(fake_seq), id="test", name="mean scores",
                   description="example with mean fastq socres",
                   letter_annotations={'phred_quality':list(means_fp.round().astype(int))})

        self.quality_scores_mean = record_mean


        # mins can be 0 when mean and max are not, and so we need to check and make sure we only trim
        # as far as neccesary (perhaps TODO do this for all values )
        np.array([[0, 1], [2, 3]], order='C')
        mins.resize(len(means_fp))
        # mins_trimmed = np.trim_zeros(mins,'b')
        mins_trimmed = mins
        fake_seq_min= ''.join(["a"]*len(mins_trimmed))

        record_mins = SeqRecord(Seq(fake_seq_min), id="test", name="mean scores",
                   description="example with mean fastq socres",
                   letter_annotations={'phred_quality':list(mins_trimmed.astype(int))})

        self.quality_scores_mins= record_mins


        maxs_trimmed = np.trim_zeros(maxs,'b')
        fake_seq_maxs= ''.join(["a"]*len(maxs_trimmed))

        record_maxs= SeqRecord(Seq(fake_seq_maxs), id="test", name="mean scores",
                   description="example with mean fastq socres",
                   letter_annotations={'phred_quality':list(maxs_trimmed.astype(int))})

        self.quality_scores_maxs= record_maxs


        return self

    def pretty(self, filename):
        '''Generate a pretty printable representation of a FastqStats object
        suitable for output of the program. The output is a tab-delimited
        string containing the filename of the input FASTQ file followed by
        the attributes of the object. If 0 sequences were read from the FASTQ
        file then num_seqs and num_bases are output as 0, and min_len, average
        and max_len are output as a dash "-".

        Arguments:
           filename: the name of the input FASTQ file
        Result:
           A string suitable for pretty printed output
        '''



        return(filename)
        # original stats info TODO add option
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
    '''Compute and print FastqStats for each input FASTQ file specified on the
    command line. If no FASTQ files are specified on the command line then
    read from the standard input (stdin).

    Arguments:
       options: the command line options of the program
    Result:
       None
    '''

    #process options once for shared code between stdin and FASTQ files
    # set mapping to default
    mapping_dict = emaps.fastq_emoji_map
    mapping_text = ""
    mapping_default = ":heart_eyes:"
    mapping_spacer = " "

    if options.custom:
        with open(options.custom) as f:
            mapping_dict = ast.literal_eval(f.read())
            mapping_text = " (custom)"
            logging.info("Custom emoji map:", options.custom)

    elif options.noemoji:
            # list of tuples to process ("min", stats.quality_scores_min) then loop and print results of map_scores with each
            mapping_dict = emaps.fastq_noemoji_map
            mapping_text = " (no-emoji)"
            mapping_default = '█'
            mapping_spacer = ""
            logging.info("Use no-emoji map")
    elif options.bin:
        # list of tuples to process ("min", stats.quality_scores_min) then loop and print results of map_scores with each
        mapping_dict = emaps.fastq_emoji_map_binned
        mapping_text = " (binned)"
        logging.info("Binned emoji map")

    else:
        logging.info("Default emoji map")

    OUTPUT_OPTIONS = CASE_NEITHER
    if options.html:
        if options.html_escape:
            OUTPUT_OPTIONS = CASE_HTML_ESCAPE
            logging.info("Output mean, min, max quality per position in html (escape)")
        else:
            OUTPUT_OPTIONS = CASE_HTML
            logging.info("Output mean, min, max quality per position in html")

    elif options.max and options.min:
        OUTPUT_OPTIONS = CASE_BOTH
        logging.info("Calculate max quality per position")
        logging.info("Calculate min quality per position")
    elif options.max:
        OUTPUT_OPTIONS = CASE_MAX
        logging.info("Calculate max quality per position")
    elif options.min:
        OUTPUT_OPTIONS = CASE_MIN
        logging.info("Calculate min quality per position")
    else:
        OUTPUT_OPTIONS = CASE_NEITHER

    if options.output:
        output_file = options.output
        logging.info("Redirect output to:", options.output)
    else:
        output_file = sys.stdout

    if options.long:
        read_size = options.long
        logging.info("Long read support up for bp up to:",options.long)
    else:
        read_size = DEFAULT_READ_LENGTH


    # before file processing
    # scale - print scale first before output, with lines starting with #
    if options.scale:
        print_scale(emaps.all_qualities, mapping_dict, options.bin)


    if options.html:
        logging.info("Collating outputs for HTML report")
        jinja_data = []
        template_dir = os.path.join(os.path.dirname(__file__), "templates")

        # Setup Jinja2 environment
        env = Environment(loader=FileSystemLoader(template_dir))
        template = env.get_template("index.html")

    if options.fastq_files:
        for fastq_filename in options.fastq_files:
            logging.info("Processing FASTQ file from {}".format(fastq_filename))
            try:
                if fastq_filename.endswith(".gz"):
                    fastq_file = gzip.open(fastq_filename, 'rt')
                else:
                    fastq_file = open(fastq_filename)

            except IOError as exception:
                exit_with_error(str(exception), EXIT_FILE_IO_ERROR)
            else:
                with fastq_file:

                    stats = FastqStats().from_file(fastq_file, read_size, options.minlen)
                    print_output(stats,fastq_filename, mapping_dict, mapping_text, mapping_default, output_file, OUTPUT_OPTIONS, spacer = mapping_spacer, window=options.window)
                    if options.html:
                        jinja_data.append(
                                html_output(stats,fastq_filename, mapping_dict, mapping_text, mapping_default, output_file, OUTPUT_OPTIONS, spacer = mapping_spacer, window=options.window))


    else:
        logging.info("Processing FASTQ file from stdin")
        print("Processing FASTQ file from stdin, press Ctrl-C to cancel. fastqe --help for more info", file = sys.stderr)

        # peek to see if gzipped
        #print(binascii.hexlify(sys.stdin.buffer.peek(1)[:2]) == b'1f8b')
        if (binascii.hexlify(sys.stdin.buffer.peek(1)[:2]) == b'1f8b'):
            #print("zipped")
            stdin_file = gzip.open(sys.stdin.buffer, 'rt')
        else:
            stdin_file = sys.stdin

        stats = FastqStats().from_file(stdin_file, read_size, options.minlen)
        print_output(stats, "-", mapping_dict, mapping_text, mapping_default, output_file, OUTPUT_OPTIONS, spacer = mapping_spacer,window = options.window)
        if options.html:
             jinja_data.append(
                                html_output(stats,"-", mapping_dict, mapping_text, mapping_default, output_file, OUTPUT_OPTIONS, spacer = mapping_spacer, window=options.window))
    

    if options.html:
        logging.info("Merging data into template for HTML report")
        html_content = template.render(items=jinja_data)

        output_path =  options.html
        #output_path = os.path.join(os.getcwd(), "output.html")
        #with open(output_path, "w", encoding="utf-8") as f:
        #    f.write(html_content)
        options.html.write(html_content)
        logging.info("Rendered HTML to file")


def html_output(stats_object,fastq_filename, mapping_dict, mapping_text,mapping_default, output_file ,output_type, sep = "\t",spacer = " ",window = 1 ):
    '''
    :param stats_object:
    :param filename:
    :param mapping_dict:
    :param mapping_text:
    :param mapping_default:
    :param seperator:
    :param output_file:
    :param output_type:
    :param sep:
    :param window:
    :return:
    '''

    # Returns a seq object for jinja template that summarises a file

    #print(stats_object)

    if output_type == CASE_HTML or output_type == CASE_HTML_ESCAPE:

        seq_object = {
                "name": stats_object.pretty(fastq_filename),
                "numseq": stats_object.num_seqs ,
                "max": stats_object.max_len ,
                "min":  stats_object.min_len ,
                "mean":  stats_object.average ,
                "bases":  stats_object.num_bases ,
                "reads": []
                }

        # sliding window - summarise in groups of N nucluotides, useful for long reads
        seq_len = len(stats_object.quality_scores_mean)+1
        window_range = range(1,seq_len,window )
        #print(list(window_range))

        mean_windowed = np.nanmean(np.r_[stats_object.quality_scores_mean.letter_annotations["phred_quality"], np.nan + np.zeros((-(seq_len-1) % window,))].reshape(-1, window), axis=-1)
        min_windowed = np.nanmean(np.r_[stats_object.quality_scores_mins.letter_annotations["phred_quality"], np.nan + np.zeros((-(seq_len-1) % window,))].reshape(-1, window), axis=-1)
        max_windowed = np.nanmean(np.r_[stats_object.quality_scores_maxs.letter_annotations["phred_quality"], np.nan + np.zeros((-(seq_len-1) % window,))].reshape(-1, window), axis=-1)


        # mean added twice so that it is shown on HTML as well as tool tip
        interleaved_scores = np.dstack((max_windowed,mean_windowed,min_windowed,mean_windowed
                                        )).flatten()

        # print(interleaved_scores)
        fake_seq_interleaved = ''.join(["i"] * len(interleaved_scores))
        #print(fake_seq_interleaved,interleaved_scores)

        interleaved_scores_seq = SeqRecord(Seq(fake_seq_interleaved), id="test", name="mean scores",
                                description="example with mean fastq socres",
                                letter_annotations={'phred_quality': list(interleaved_scores.astype(int))})


        #html currently supports a HTML output
        just_emojis = interleaved_scores

        # print([str(i) + "," + str(i + 3) for i in range(0, len(just_emojis), 3)])
        just_emojis_3 = [just_emojis[i:i + 3] for i in range(0, len(just_emojis), 3)]
        just_emojis_4 = [just_emojis[i:i + 4] for i in range(0, len(just_emojis), 4)]
        # print(just_emojis_3)
        

        html_seqs = []

        for ind, (max_e,mean_e,min_e,mean_e_2) in enumerate(just_emojis_4):
            # print(max_e)
            # print(np.array(max_e))
            interleaved_scores_subset = np.array([max_e,mean_e,min_e,mean_e_2])
            print_seq = SeqRecord(Seq("mmmm"), id="test", name="max mean min mean scores",
                                  description="example with max mean min mean fastq socres",
                                  letter_annotations={'phred_quality': list(interleaved_scores_subset.astype(int))})

            # print_values = str(ind+1) +" " + map_scores(print_seq, mapping_dict=mapping_dict,default_value=mapping_default,
            #                                             spacer=' ')
            (seq_max, seq_mean, seq_min,_) = map_scores_html(print_seq, mapping_dict=mapping_dict,default_value=mapping_default,
                                                        spacer='')
            # print(print_values)
            # print(' '.join(print_values))
            if window == 1:
                paging_text = ""
            else:
                paging_text ="-"+str(min(window_range[ind]+window-1,seq_len-1))
            
            seq_object["reads"].append( {
                    "pos": str(window_range[ind])+paging_text,
                    "max": seq_max,
                    "mean": seq_mean,
                    "min": seq_min
                    } )


        # print(stats_object.pretty(fastq_filename), "(html)" + mapping_text,
        #  html_header + ''.join(html_text) + html_footer
        #       ,
        #       sep=sep,file=output_file)

        # print(stats_object)


        return seq_object


def print_output(stats_object,fastq_filename, mapping_dict, mapping_text,mapping_default, output_file ,output_type, sep = "\t",spacer = " ",window = 1 ):
    '''
    :param stats_object:
    :param filename:
    :param mapping_dict:
    :param mapping_text:
    :param mapping_default:
    :param seperator:
    :param output_file:
    :param output_type:
    :param sep:
    :param window:
    :return:
    '''
    if output_type == CASE_BOTH or output_type == CASE_MAX:
            print(stats_object.pretty(fastq_filename), "max" + mapping_text,
              map_scores(stats_object.quality_scores_maxs, mapping_dict=mapping_dict, default_value = mapping_default,spacer=spacer), sep=sep, file=output_file)

    print(stats_object.pretty(fastq_filename), "mean" + mapping_text,map_scores(stats_object.quality_scores_mean, mapping_dict=mapping_dict,default_value = mapping_default,spacer=spacer), sep=sep,file=output_file)

    if output_type == CASE_BOTH or output_type == CASE_MIN:
            print(stats_object.pretty(fastq_filename), "min" + mapping_text,
              map_scores(stats_object.quality_scores_mins, mapping_dict=mapping_dict, default_value = mapping_default,spacer=spacer), sep=sep, file=output_file)



def nomap(input):

    return(input)


def map_scores_html(sequence,
               mapping_dict = emaps.fastq_emoji_map,
               default_value = ":heart_eyes:",
               mapping_function = emojify,
              spacer = " "):
    '''
    :param sequence:
    :param mapping_dict:
    :param default_value:
    :param mapping_function:
    :param spacer:
    :return:
    '''
    # print(sequence)

    mapped_values = spacer.join([mapping_function(mapping_dict.get(s, default_value)) for s in QualityIO._get_sanger_quality_str(sequence)])
    return(mapped_values)

def map_scores(sequence,
               mapping_dict = emaps.fastq_emoji_map,
               default_value = ":heart_eyes:",
               mapping_function = emojify,
              spacer = " "):
    '''
    :param sequence:
    :param mapping_dict:
    :param default_value:
    :param mapping_function:
    :param spacer:
    :return:
    '''
    # print(sequence)
    mapped_values = spacer.join([mapping_function(mapping_dict.get(s, default_value)) for s in QualityIO._get_sanger_quality_str(sequence)])
    return(mapped_values)






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

def run_fastqe(fastq_files,minlen=0,scale=False,version=False,
               mean = True,custom=None,noemoji=False,min=False,max=False,
               output=None,long=None,log=None,bin=False,html=False,html_escape=False,noheader=False, window=1):
    options = Namespace(fastq_files=fastq_files, bin=bin, custom=custom,log=log, long=long, max=max, mean=mean, min=min, minlen=minlen, noemoji=noemoji, output=output, scale=scale, version=version,html=html,html_escape=html_escape,noheader=noheader, window=window)
    if options.version:
        print(PROGRAM_NAME,PROGRAM_VERSION)
        return
    process_files(options)

def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    init_logging(options.log)
    if (not options.noheader):
        if(options.output):
            print(HEADER,file=options.output)
        else:
            print(HEADER)
    process_files(options)


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
