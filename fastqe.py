from pyemojify import emojify
import numpy as np
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqIO import QualityIO


seq_emoji_map = {
    'A': ':apple:',
    'C': ':corn:',
    'T': ':tomato:',
    'G': ':grapes:',
    'N': ':question:'
}


#From https://en.wikipedia.org/wiki/FASTQ_format
fastq_emoji_map = {
    '!': ':skull:',
    '\"': ':poop:',
    '#': ':no_entry_sign:',
    '$': ':alien:',
    '%': ':ghost:',
    '&': ':japanese_ogre:',
    '\'': ':space_invader:',
    '(': ':japanese_goblin:',
    ')': ':imp:',
    '*': ':smiling_imp:',
    '+': ':zap:',
    ',': ':cyclone:',
    '-': ':umbrella:',
    '.': ':foggy:',
    '/': ':cloud:',
    '0': ':rat:',
    '1': ':snowman:',
    '2': ':snowflake:',
    '3': ':volcano:',
    '4': ':anger:',
    '5': ':fire:',
    '6': ':bomb:',
    '7': ':pray:',
    '8': ':broken_heart:',
    '9': ':thumbsdown:',
    ':': ':x:',
    ';': ':zzz:',
    'A': ':point_down:',
    'B': ':scream_cat:',
    'D': ':crying_cat_face:',
    'J': ':pouting_cat:',
    '@': ':no_good:',
    '<': ':rage:',
    '=': ':angry:',
    'C': ':warning:',
    '>': ':scream:',
    'E': ':grimacing:',
    'F': ':see_no_evil:',
    'G': ':hear_no_evil:',
    'H': ':speak_no_evil:',
    'I': ':monkey_face:',
    '?': ':fearful:',
    'K': ':cold_sweat:',
    'L': ':sweat_smile:',
    'M': ':sweat:',
    'N': ':weary:',
    'O': ':tired_face:',
    'P': ':anguished:',
    'Q': ':neutral_face:',
    'R': ':sob:',
    'S': ':cry:',
    'T': ':hushed:',
    'U': ':innocent:',
    'V': ':worried:',
    'W': ':expressionless:',
    'X': ':no_mouth:',
    'Y': ':unamused:',
    'Z': ':confused:',
    '[': ':frowning:',
    '\\': ':dizzy_face:',
    ']': ':disappointed:',
    '^': ':persevere:',
    '_': ':pensive:',
    '`': ':disappointed_relieved:',
    'a': ':confounded:',
    'b': ':sleeping:',
    'c': ':sleepy:',
    'd': ':mask:',
    'e': ':open_mouth:',
    'f': ':smirk:',
    'g': ':satisfied:',
    'h': ':yum:',
    'i': ':astonished:',
    'j': ':smile:',
    'k': ':laughing:',
    'l': ':triumph:',
    'm': ':relieved:',
    'n': ':joy:',
    'o': ':grin:',
    'p': ':flushed:',
    'q': ':stuck_out_tongue:',
    'r': ':stuck_out_tongue_closed_eyes:',
    's': ':stuck_out_tongue_winking_eye:',
    't': ':kissing_smiling_eyes:',
    'u': ':kissing:',
    'v': ':kissing_closed_eyes:',
    'w': ':kissing_heart:',
    'x': ':heart_eyes:',
    'y': ':wink:',
    'z': ':relaxed:',
    '{': ':blush:',
    '|': ':grinning:',
    '}': ':smiley:',
    '~': ':sunglasses:'
}


def mean_emoji(filename):
    # works for up to 500bp reads
    means = np.zeros(500)
    seq_count = 0
    
    for r in SeqIO.parse(filename, "fastq"):
        index=0
        for s in r.letter_annotations["phred_quality"]:
            means[index] += s
            index = index+1
        seq_count = seq_count + 1
        
    cleaned = np.trim_zeros(means)
    means_fp = cleaned/seq_count


    fake_seq= ''.join(["a"]*len(means_fp.round()))
                   
    record = SeqRecord(Seq(fake_seq), id="test", name="mean scores", 
                   description="example with mean fastq socres",
                   letter_annotations={'phred_quality':list(means_fp.round().astype(int))}

    )
    
    print("".join([emojify(fastq_emoji_map[s]) for s in QualityIO._get_sanger_quality_str(record)]))



def main():
    mean_emoji(sys.argv[1])

if __name__ == "__main__":
    main()

