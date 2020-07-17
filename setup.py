from setuptools import setup
LONG_DESCRIPTION = \
'''The program reads one or more input FASTQ files.
For each file it computes the minimum, maximum and mean FASTQ quality score at each position across all reads in a file.

For some reason, it then represents these as emoji.
'''


setup(
  name = 'fastqe',
  packages = ['fastqe'],
  version = '0.2.6',
  license = 'MIT',
  description = 'A emoji based bioinformatics command line tool',
  long_description=(LONG_DESCRIPTION),
  author = 'Andrew Lonsdale',
  author_email = 'andrew.lonsdale@lonsbio.com.au',
  url = 'https://github.com/lonsbio/fastqe',
  download_url = 'https://github.com/lonsbio/fastqe/archive/fastqe-0.2.5.tar.gz',
  keywords = ['emoji', 'bioinformatics', 'next-generation sequencing'],
  classifiers = [],
  install_requires=["biopython>=1.66",'pyemojify'],
  entry_points={
        'console_scripts': ['fastqe = fastqe.fastqe:main']
    },

)
