---
title: 'FASTQE: sequence quality visualisation with Emoji'
tags:
  - emoji
  - bioinformatics
  - visualisation
  - next-generation sequencing
  - education
authors:
  - name: Andrew Lonsdale
    orcid: 0000-0002-0292-2880
    affiliation: "1, 2, 3" # (Multiple affiliations must be quoted)
affiliations:
 - name: Peter MacCallum Cancer Centre, Melbourne, VIC, Australia
   index: 1
 - name: Sir Peter MacCallum Department of Oncology, University of Melbourne, Parkville, VIC, Australia
   index: 2
 - name: Murdoch Children’s Research Institute, Parkville, VIC, Australia;
   index: 3
date: 25 October 2024
bibliography: paper.bib
---

# Summary

Bioinformatics is the science 🧑🏻‍🔬 of understanding and analysing biological information 💡, such as the genetic information contained in DNA 🧬. Modern scientific equipment can measure biological sequences with high accuracy 🎯, yet no technology is flawless and called nucleotide bases can be wrong ❌. Data quality issues need to be assessed and addressed to ensure confidence in downstream interpretation 📝 in medicine and science.

FASTQE 🤔 is a utility for viewing the quality of biological sequence data by using emoji 😂. It summarises the average quality score for each position in a set of biological sequence measurements, and transcribes that average quality into a corresponding emoji to see the good 😍, the bad 💩,and the ugly 💀 of sequencing data.  When invoked from the command line it can also display the minimum 📉 and maximum 📈 quality scores per position, and bin quality 🗑️ scores into a reduced set of emoji. Custom emoji can also be used 🐍 🌵 👍.

FASTQE can be used to rapidly 🏃 assess the quality of sequence data. It also helps transform complex 🤯 bioinformatics data into engaging, emoji-based visualisations 📊, making bioinformatics concepts more accessible 😌 and adding an element of fun 🤪 to scientific education 📚 and communication 🗣️.

# Statement of need

Assessing the quality of high-throughput sequencing data is critical for ensuring medical and scientific conclusions drawn from that data are correct. Accounting for technical limitations, contamination, and incorrect data can prevent costly mistakes or incorrect results.  Existing tools for quality assessment often use a GUI or require in-depth knowledge to assess. FASTQC [@fastqc] is the classic tool for sequence quality analysis. Experienced bioinformaticians are familiar with its graphical output format, and education in bioinformatics will include it as elementary content. Recent alternatives that are used where speed and efficiency are critical include fastp [@fastp], falco [@falco] and others, but FASTQC remains the de facto standard.  It is currently the  most effective way to communicate quality scores for sequencing data to other scientists.

FASTQE has the same goal of expressing sequence quality, but fills a need in accessibility and simplicity. As a command-line-first tool, it is designed to be run alongside other bioinformatics tools in environments such as high performance computing clusters. It gives insight into sequence quality in an informal manner. Its feature set is similar to but not as detailed as FASTQC, and focuses on core concepts of sequence quality analysis with an engaging visualization format.

The interpretability of emoji in education and outreach is also important, especially when considering non-English speakers, and FASTQE improves the accessibility and diversity of the audience that biomedical researchers can communicate to. FASTQE is a Python package that can be used in programming exercises, as well as a command line tool. It can be installed both via PyPI and Bioconda. FASTQE offers a fun gateway into the use of the command line, which may also benefit users used to graphical interfaces.


# Usage

The utility of FASTQE can easily be seen by comparing before and after quality filtering on sequencing data. For some (compressed) data in the FASTQ format [@cock_sanger_2010], FASTQE will produce by default an emoji for the mean score at each base position. This data clearly has quality issues that need investigating:

```
$ fastqe sample.50bp.fastq.gz
sample.50bp.fastq	mean	😀 😀 🚨 🚨 🚨 🚨 🚨 🚨 🚨 🚨 🚨 💩 🚨 🚨 🚨 🚨 💩 💩 💩 🚨 💩 💩 💩 💩 🚨 💩 😡 💩 💩 💩 💩 💩 😡 💩 💩 💩 💩 💩 😡 😡 💩 😡 😡 😡 😡 😡 😡 😡 😡 😿
```

After the removal of low-quality sequences, for example with `fastp` or `Trimmomatic` [@trimmomatic], the remaining files can be read in with FASTQE to see the effect:

```
$ sample.50bp.filtered.fastq.gz
sample.50bp.filtered.fastq.gz        mean    😝 😌 😌 😌 😌 😌 😝 😌 😌 😋 😌 😄 😋 😌 😝 😌 😄 😄 😄 😋 😋 😄 😄 😄 😋 😋 😆 😆 😄 😄 😆 😄 😄 😄 😆 😄 😄 😆 😆 😆 😆 😄 😄 😆 😆 😆 😆 😘 😘 😡
```

# Design

FASTQE employs a simple algorithm to parse FASTQ format using BioPython [@biopython] utilities. For each read, we extract a list of numerical quality scores at each position, and maintain a record of the minimum, maximum and average quality for each position. The numerical scores are originally encoded in the FASTQ format as ASCII-encoded Phred summary scores, where the score relates to probability of an error. BioPython interprets these numerically, however FASTQE then re-creates these ASCII encodings for the summarised minimum, average and maximum scores, and then applies a mapping from quality score ASCII to emoji.

The default mappings with Phred score, ASCII character and emoji are as follows:

```
  0 ! 🚫	  1 " ❌	  2  👺	      3 $ 💔	  4 % 🙅	  5 & 👾
  6 ' 👿	  7 ( 💀	  8 ) 👻	  9 * 🙈	  10 + 🙉	  11 , 🙊
  12 - 🐵	  13 . 😿	  14 / 😾	  15 0 🙀	  16 1 💣	  17 2 🔥
  18 3 😡	  19 4 💩	  20 5 🚨	  21 6 😀	  22 7 😅	  23 8 😏
  24 9 😊	  25 : 😙	  26 ; 😗	  27 < 😚	  28 = 😃	  29 > 😘
  30 ? 😆	  31 @ 😄	  32 A 😋	  33 B 😌	  34 C 😝	  35 D 😛
  36 E 😜	  37 F 😉	  38 G 😁	  39 H 😄	  40 I 😎     41 J 😍
```

Binning into simplified emoji is also available with the `--bin` option to improve impact:

```
  0 ! 🚫	  1 " 🚫	  2  💀	      3 $ 💀	  4 % 💀	  5 & 💀
  6 ' 💀	  7 ( 💀	  8 ) 💀	  9 * 💀	  10 + 💩	  11 , 💩
  12 - 💩	  13 . 💩	  14 / 💩	  15 0 💩	  16 1 💩	  17 2 💩
  18 3 💩	  19 4 💩	  20 5 🚨	  21 6 🚨	  22 7 🚨	  23 8 🚨
  24 9 🚨	  25 : 😄	  26 ; 😄	  27 < 😄	  28 = 😄	  29 > 😄
  30 ? 😆	  31 @ 😆	  32 A 😆	  33 B 😆	  34 C 😆	  35 D 😎
  36 E 😎	  37 F 😎	  38 G 😎	  39 H 😎	  40 I 😎     41 J 😍
```

For those who can't (or won't, cowards) use emoji, ASCII boxes can be used to proportionally indicate sequence quality with `--noemoji`:

```
$ fastqe --noemoji sample.50bp.fastq.gz	
sample.50bp.fastq.gz	mean (no-emoji)	▄▄▄▄▄▄▄▄▄▄▄▃▄▄▄▄▃▃▃▄▃▃▃▃▄▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
```

```
$ fastqe --noemoji sample.50bp.filtered.fastq.gz
sample.50bp.filtered.fastq.gz	mean (no-emoji)	▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▆▅▅▃
```

Users can provide a custom mapping of custom emoji to quality in a text file. FASTQE is designed with `pyemojify` to use emoji aliases, e.g. `:crying_cat_face:`, however direct use of emoji in the dictionary is also supported. Revisiting the sample data before and after quality filtering demonstrates the visual narratives possible with custom emoji, such as in this case, turning silver into gold.

```
$ fastqe --custom custom.txt  sample.50bp.fastq.gz
sample.50bp.fastq	mean (custom)	🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥈 🥇 🥇 🥇 🥇 🥈 🥈 🥈 🥇 🥈 🥈 🥈 🥈 🥇 🥈 🥈 🥈 🥈 🥈 🥈 🥈 🥈 🥈 🥈 🥈 🥈 🥈 🥈 🥈 🥈 🥈 🥈 🥈 🥈 🥈 🥈 🥈 🥈 🥈
```

```
$ fastqe --custom custom.txt  sample.50bp.filtered.fastq.gz
sample.50bp.filtered.fastq.gz	mean (custom)	🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥇 🥈
```

# Conclusion

FASTQE is a software tool that serves both a practical and educational purpose. It can be adapted for many purposes. A resource on sequence quality suitable for high school and undergraduate students has also been developed [@jacques], and it has also been used in short courses and bioinformatics training [@batut]. Everyone knows that using emoji to visualise biological sequencing data is a silly idea. We have shown here that maybe it isn’t as silly as it sounds.

# Acknowledgements

Thanks to the users of FASTQE.  We would also like to thank Ray Enke and the NIBLES incubator team for developing teaching materials using FASTQE.  We would like to acknowledge contributions to FASTQE, including during the BCC2020 CoFest: Björn Grüning, Catherine Bromhead, Clare Sloggett, Clarissa Womack, Helena Rasche, Maria Doyle, Michael Franklin, Nicola Soranzo, Phil Ewels. Additional thanks to Clare Slogget for manuscript feedback, and the PyCon Au 2024 academic team Maia Sauren and Alan Rubin. 

# References
