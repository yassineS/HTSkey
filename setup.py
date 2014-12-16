from distutils.core import setup
from setuptools import find_packages

__version__ = 'a 0.1'

README = open('README.md').read()

setup(name='GenomeKey2',
            version=__version__,
            description = "Next-generation Sequencing Analysis Pipeline for COSMOS 2.0",
            license='Research only',
            long_description=README,
            packages=find_packages(),
            scripts=['bin/genomekey'],
            install_requires=['pysam','ipdb']
      )
