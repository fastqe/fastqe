import fastqe.fastqe
import fastqe.fastqe_map as emaps
from pyemojify import emojify


def test_fastq_scores_mean():
  '''Test calculation of mean phred quality scores from FASTQ file'''
  filename = 'tests/data/test.fastq'
  expected_scores = [29, 30, 33, 33, 32, 36, 36, 36, 32, 36, 36, 
    37, 38, 38, 38, 38, 38, 37, 37, 38, 38, 38, 39, 38, 38, 
    39, 39, 38, 39, 38, 38, 38, 39, 39, 39, 39, 38, 39, 37, 
    38, 39, 38, 38, 39, 39, 38, 38, 39, 37, 38, 38, 34, 34, 
    38, 38, 33, 38, 37, 38, 38, 38, 38, 38, 34, 38, 38, 38, 
    38, 34, 38, 38, 37, 38, 38, 37, 38, 37, 28, 37, 38, 38, 
    38, 39, 38, 39, 39, 38, 39, 38, 34, 38, 39, 39, 38, 37, 
    38, 37, 38, 38, 38, 38, 39, 39, 38, 38, 33, 37, 36, 38, 
    38, 38, 38, 38, 37, 38, 36, 38, 34, 38, 38, 38, 38, 38, 
    38, 38, 38, 38, 37, 39, 36, 38, 39, 38, 39, 38, 39, 38, 
    38, 36, 36, 36, 36, 37, 36, 37, 38, 37, 32, 36, 35, 37, 
    36, 36, 36, 37, 37, 37, 38, 38, 38, 38, 38, 32, 35, 36, 
    36, 38, 35, 31, 32, 35, 34, 37, 37, 31, 33, 36, 36, 35, 
    30, 35, 37, 36, 36, 36, 36, 32, 34, 32, 36, 37, 37, 31, 
    37, 34, 34, 35, 32, 36, 36, 36, 31, 30, 36, 34, 37, 36, 
    32, 37, 36, 31, 33, 33, 36, 29, 34, 30, 31, 36, 36, 36, 
    37, 35, 31, 33, 25, 24, 31, 34, 30, 36, 37, 36, 31, 31, 
    33, 33, 36, 37, 36, 37, 30, 30, 35, 25, 28, 25, 29, 34, 25, 18]
  read_length = 500
  min_length = 0
  stats = fastqe.fastqe.FastqStats().from_file(
     filename, read_length, min_length)
  calculated_scores = stats.quality_scores_mean.letter_annotations['phred_quality']
  assert calculated_scores==expected_scores

def test_fastq_mapping():
  '''Test calculation of emoji-mapped mean scores from FASTQ file'''
  filename = 'tests/data/test.fastq'
  expected_emojis = "😘 😆 😌 😌 😋 😜 😜 😜 😋 😜 😜 😉 😁 😁 😁 😁 😁 "+\
    "😉 😉 😁 😁 😁 😄 😁 😁 😄 😄 😁 😄 😁 😁 😁 😄 😄 😄 😄 "+\
    "😁 😄 😉 😁 😄 😁 😁 😄 😄 😁 😁 😄 😉 😁 😁 😝 😝 😁 😁 "+\
    "😌 😁 😉 😁 😁 😁 😁 😁 😝 😁 😁 😁 😁 😝 😁 😁 😉 😁 😁 "+\
    "😉 😁 😉 😃 😉 😁 😁 😁 😄 😁 😄 😄 😁 😄 😁 😝 😁 😄 😄 "+\
    "😁 😉 😁 😉 😁 😁 😁 😁 😄 😄 😁 😁 😌 😉 😜 😁 😁 😁 😁 "+\
    "😁 😉 😁 😜 😁 😝 😁 😁 😁 😁 😁 😁 😁 😁 😁 😉 😄 😜 😁 "+\
    "😄 😁 😄 😁 😄 😁 😁 😜 😜 😜 😜 😉 😜 😉 😁 😉 😋 😜 😛 "+\
    "😉 😜 😜 😜 😉 😉 😉 😁 😁 😁 😁 😁 😋 😛 😜 😜 😁 😛 😄 "+\
    "😋 😛 😝 😉 😉 😄 😌 😜 😜 😛 😆 😛 😉 😜 😜 😜 😜 😋 😝 😋 "+\
    "😜 😉 😉 😄 😉 😝 😝 😛 😋 😜 😜 😜 😄 😆 😜 😝 😉 😜 😋 😉 "+\
    "😜 😄 😌 😌 😜 😘 😝 😆 😄 😜 😜 😜 😉 😛 😄 😌 😙 😊 😄 😝 😆 "+\
    "😜 😉 😜 😄 😄 😌 😌 😜 😉 😜 😉 😆 😆 😛 😙 😃 😙 😘 😝 😙 😡"
  read_length = 500
  min_length = 0
  stats = fastqe.fastqe.FastqStats().from_file(
     filename, read_length, min_length)
  mapped = fastqe.fastqe.map_scores(stats.quality_scores_mean)
  assert mapped==expected_emojis
