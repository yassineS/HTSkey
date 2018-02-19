GenomeKey: for COSMOS 2.0
==========

GenomeKey is a Whole Genome Analysis pipeline, that can call variants from FASTQ or BAM files, as well as massively
annotate VCF files.  It is implemented and made possible by the Cosmos workflow management system.

Components include:

* *BWA + GATK Best Practices v4* Cosmos workflow
* *AnnovarExtensions annotation* Cosmos workflow

Install
=======

1) Install Cosmos using virtualenvwrapper

2) Clone git@github.com:LPM-HMS/htsKey.git

3) Activate Cosmos virtualenv

    $ workon cosmos

4) Add GenomeKey to your PYTHONPATH when you're in the cosmos virtualenv

    add2virtualenv /path/to/htsKey


Configuration
=============

If you're running things on Orchestra or AWS, htsKey does not need any configuration, and the rest of this
section is only for educational purposes.

GenomeKey is configured in ``htskey/htskey/wga_settings.py`` where it points to the correct paths to the
GATK bundle, reference genome, and binaries.  It chooses these paths based on the ``cosmos.ini`` ``server_name``
setting.  If ``server_name`` is set to ``orchestra``, it will point to ``/scratch/esg21/WGA`` where all the files such as
annotation databases and binaries for GATK, BWA, AnnovarExtensions, etc. are located.

AnnovarExtensions is configured in WGA/annovarext_data/config.ini which may need to be edited if you are using an
installation of
of the WGA folder that is not ``/scratch/esg21/WGA`` (for ex, you copied it to AWS)

Usage
======

Inside the htsKey directory, execute:

$ bin/htskey -h

From BAM
+++++++++

    htskey bam -n "My Workflow from BAM" -i /path/to/bam1

    htskey bam -n "My Multi-BAM Workflow" -il /path/to/bam.list

From FASTQ
++++++++++

    genomekey json -n "My workflow from a JSON file" '/path/to/json'

    json file should be of the format:

.. code-block:: json

    [
        {
            'chunk': 001,
            'library': 'LIB-1216301779A',
            'sample_name': '1216301779A',
            'platform': 'ILLUMINA',
            'platform_unit': 'C0MR3ACXX.001'
            'pair': 0, #0 or 1
            'path': '/path/to/fastq'
        },
        {..}
    ]

.. note::
    I have GenomeKey set to launch you into an ipdb post mortem debugging session on any exceptions.  That behavior is
    set in bin/genomekey.  To quit enter **q** then enter.

Testing
========

**-test** will inform htsKey you are running a test dataset.  It will only analyse chr20, and
drmaa_native_specification() will be adjusted accordingly automatically for Orchestra, so that requests are sent to
the mini queue with a cpu_requirement of 1.  htsKey comes with some test data, so you can just
run this from the htsKey directory:

.. code-block:: bash

    $ htskey -t bam -n 'Test GK' -il htskey/test/bams.list

Issues
======

* If there are unpaired reads when converting a BAM to FASTQ, they're not used in the re-alignment
