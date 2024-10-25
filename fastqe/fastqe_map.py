#!/usr/bin/env python

seq_emoji_map = {
    'A': ':apple:',  # avocado? differnet colours?
    'C': ':corn:',
    'T': ':tomato:',
    'G': ':grapes:',
    'N': ':question:'
}

all_qualities = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"

#From https://en.wikipedia.org/wiki/FASTQ_format
# note order not exact here
fastq_emoji_map = {
    '!': ':no_entry_sign:',
    '\"': ':x:',
    '#': ':japanese_goblin:',
    '$': ':broken_heart:',
    '%': ':no_good:',
    '&': ':space_invader:',
    '\'': ':imp:',
    '(': ':skull:',
    ')': ':ghost:',
    '': ':pouting_cat:',
    '*': ':see_no_evil:',
    '+': ':hear_no_evil:',
    ',': ':speak_no_evil:',
    '/': ':pouting_cat:',
    '-': ':monkey_face:',
    '.': ':crying_cat_face:',
    '0': ':scream_cat:',
    '1': ':bomb:',
    '2': ':fire:',
    '3': ':rage:',
    '4': ':poop:',
    '5': "🚨", # #':warning:',
    '6': ':grinning:',
    '7': ':sweat_smile:',
    '8': ':smirk:',
    '9': ':blush:',
    ':': ':kissing_smiling_eyes:',
    ';': ':kissing:',
    '<': ':kissing_closed_eyes:',
    '>': ':kissing_heart:',
    '@': ':smile:',
    '=': ':smiley:',
    '?': ':laughing:',
    'A': ':yum:',
    'B': ':relieved:', # was :relaxed:
    'D': ':stuck_out_tongue:',
    'C': ':stuck_out_tongue_closed_eyes:',
    'E': ':stuck_out_tongue_winking_eye:',
    'G': ':grin:',
    'H': ':smile:',
    'I': ':sunglasses:',
    'J': ':heart_eyes:',
    'F': ':wink:'



}

# only use 0 to 50  use same emoji


# binning - i.e https://www.illumina.com/documents/products/technotes/technote_understanding_quality_scores.pdf
fastq_emoji_map_binned= {
#N (no call) N (no call)
'!': ':no_entry_sign:',
'"': ':no_entry_sign:',

#2–9 6
'#': ':skull:',
'$': ':skull:',
'%': ':skull:',
'&': ':skull:',
'\'': ':skull:',
'(': ':skull:',
')': ':skull:',
'*': ':skull:',

#10–19 15
'+': ':poop:' ,
',': ':poop:' ,
'-': ':poop:' ,
'.': ':poop:' ,
'/': ':poop:' ,
'0': ':poop:' ,
'1': ':poop:' ,
'2': ':poop:' ,
'3': ':poop:' ,
'4': ':poop:' ,

#20–24 22
'5': '🚨',#:warning:',
'6': '🚨',#:warning:',
'7': '🚨',#:warning:',
'8': '🚨',#:warning:',
'9': '🚨',#:warning:',

#25–29 27
':': ':smile:',
';': ':smile:',
'<': ':smile:',
'=': ':smile:',
'>': ':smile:',


#30–34 33
'?': ':laughing:',
'@': ':laughing:',
'A': ':laughing:',
'B': ':laughing:',
'C': ':laughing:',

#35–39 37
'D': ':sunglasses:',
'E': ':sunglasses:',
'F': ':sunglasses:',
'G': ':sunglasses:',
'H': ':sunglasses:',

#≥ 40 40
'I': ':heart_eyes:',
'J': ':heart_eyes:',

}


# binning - i.e https://www.illumina.com/documents/products/technotes/technote_understanding_quality_scores.pdf
fastq_noemoji_map = {
#N (no call) N (no call)
'!': '▁',
'"': '▁',

#2–9 6
'#': '▂',
'$': '▂',
'%': '▂',
'&': '▂',
'\'': '▂',
'(': '▂',
')': '▂',
'*': '▂',

#10–19 15
'+': '▃' ,
',': '▃' ,
'-': '▃' ,
'.': '▃' ,
'/': '▃' ,
'0': '▃' ,
'1': '▃' ,
'2': '▃' ,
'3': '▃' ,
'4': '▃' ,

#20–24 22
'5': '▄',
'6': '▄',
'7': '▄',
'8': '▄',
'9': '▄',


#25–29 27
':': '▅',
';': '▅',
'<': '▅',
'=': '▅',
'>': '▅',


#30–34 33
'?': '▆',
'@': '▆',
'A': '▆',
'B': '▆',
'C': '▆',

#35–39 37
'D': '▇',
'E': '▇',
'F': '▇',
'G': '▇',
'H': '▇',

#≥ 40 40
'I': '█',
'J': '█',

}
