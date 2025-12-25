#!/usr/bin/env python

from setuptools import setup, find_packages

# Remember to change when making a new release
version = '0.5.0'
dl_version = 'master' if 'dev' in version else 'v{}'.format(version)

with open('README.md') as f:
    readme = f.read()

setup(
    name = 'fastqe',
    version = version,
    description = 'A emoji based bioinformatics command line tool',
    long_description = readme,
    long_description_content_type = 'text/markdown',
    keywords = ['emoji', 'bioinformatics', 'next-generation sequencing'],
    author = 'Andrew Lonsdale',
    author_email = 'andrew.lonsdale@lonsbio.com.au',
    url = 'https://github.com/fastqe/fastqe',
    download_url = 'https://github.com/fastqe/fastqe/tarball/{}'.format(dl_version),
    license = 'BSD-3-Clause',
    entry_points = {'console_scripts': ['fastqe = fastqe.fastqe:main']},
    install_requires = ["biopython>=1.66",'pyemojify','setuptools>=38.6'],
    setup_requires = ['twine>=1.11.0', 'setuptools>=38.6'],
    packages = find_packages(exclude=('test', 'docs')),
    classifiers = [
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)
