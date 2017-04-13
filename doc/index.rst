#######################
CpG Dinucleotides Stats
#######################

**cpgStats** is a C application to parse CpG Dinucleotides file and get a summary of statistics and bed file annotation.
It can also compare two dinucleotide files to get intersection stats.

*******************
cpgStats QuickStart
*******************

.. note::

    Before starting to run **cpgStats** please make sure input files are coordinate sorted.

0. Download
===========

**cpgStats** can be downloaded from: `github cpgStats`_

.. _github cpgStats: https://github.com/MarcosFernandez/cpgStats/

1. Compile
==========

1.1. Release
------------

Download and extract the code and compile it from ``cpgStats/Release`` path.

::

    cd cpgStats/Release
    make clean
    make 

1.2. Debug
----------

For cases of code testing it is possible to compile a debug version code.

::

    cd cpgStats/Debug
    make clean
    make 

2. Run cpgStats
===============

2.1. Get Summary Stats
----------------------

2.1.1 From a zipped file
^^^^^^^^^^^^^^^^^^^^^^^^

::

    cpgStats -i file_cpg.txt.gz -z

2.1.2 From a unzipped file
^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    cpgStats -i file_cpg.txt

2.1.3 Save stats in json file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    cpgStats -i file_cpg.txt.gz -z -o stats.json


2.2. Annotate Bed File
----------------------

::

    cpgStats -i file_cpg.txt.gz -z -b promoters.sorted.bed -a promoters.annotated.bed

2.3. Intersection Stats
-----------------------

::

    cpgStats -x first_file_cpg.txt.gz -y second_file_cpg.txt.gz -g

    -g In case dinucleotide files are compressed.

3. Format files
===============

3.1. Input CpG Dinucleotides Format. **(Should be coordinate sorted.)**
-----------------------------------------------------------------------

+----------+-----------+------------+---------------+-------------+--------------------+------------------+---------------+---------------+------------+-------------+-----------+------------+
| CONTIG   |POSITION   |CONTEXT     |CONTEXT        | Phred       |Methylation         |Methylation       |Info Cs & Ts.  |Info Ts & Cs.  |Info Ts,Cs, |All Info and | MC8       | MC8 next   |
|          |           |(Reference) |(Genotype Call)| Score       |(Mode Beta Distrib.)|Deviation         |For Meth.      |For Meth.      |As non Info |non Info. For|           | by base    |
|          |           |            |               |             |                    |(SD Beta Distrib.)|               |               |for Meth.   |methylation. |           |            |  
+==========+===========+============+===============+=============+====================+==================+===============+===============+============+=============+===========+============+
|Contig    |Position   |Dinucleotide|Dinucleotide   |Phred scaled |Methylation value   |Standard Deviation|Informative Cs |Informative Ts |Informative |All base     |Base counts|Base counts |
|Chromosome|in the     |base call   |base according |probability  |according to Mode of|of the methylation|for methylation|for methylation|Ts,Cs and   |counts.      |non        |non         |
|name.     |reference  |according to|to             |that REF/ALT |Beta Distribution   |value.            |or Cs and Ts in|or Ts and As in|non         |             |informative|informative |
|          |sequence.  |reference.  |genotype call. |polymorphism |estimated from      |St. Dev of Beta   |case of        |case of        |informative |             |for        |for         |
|          |           |            |               |exists.      |informative         |Distribution.     |CG/CC/GG.      |CG/CC/GG.      |Cs for      |             |methylation|methylation |
|          |           |            |               |             |methylation counts. |                  |               |               |methylation |             |(ACGT)     |(ACGT)      |
|          |           |            |               |             |                    |                  |               |               |or          |             |followed   |followed    |
|          |           |            |               |             |                    |                  |               |               |informative |             |by         |by          |
|          |           |            |               |             |                    |                  |               |               |Cs,Ts,As    |             |informative|informative |
|          |           |            |               |             |                    |                  |               |               |and non     |             |for        |for         |
|          |           |            |               |             |                    |                  |               |               |informative |             |methylation|methylation |
|          |           |            |               |             |                    |                  |               |               |Cs,Gs in    |             |(ACGT).    |(ACGT)      |
|          |           |            |               |             |                    |                  |               |               |case of     |             |           |next by     |
|          |           |            |               |             |                    |                  |               |               |CG/CC/GG.   |             |           |base.       |
+----------+-----------+------------+---------------+-------------+--------------------+------------------+---------------+---------------+------------+-------------+-----------+------------+


3.2. Input BED File Format. **(Should be coordinate sorted.)**
--------------------------------------------------------------

+----------+----------+---------+-------------+-----+--------------+
|CONTIG    |START     |END      |EXTRA FIELD 1|...  |EXTRA FIELD N |
+==========+==========+=========+=============+=====+==============+
|Contig    |Start     |Last     |Extra        |     |Extra         |
|Chromosome|Position  |Position |optional     |...  |optional      |
|name.     |          |         |field one    |     |field N       |
+----------+----------+---------+-------------+-----+--------------+

3.3. Output **annotated** BED
----------------------------- 

+----------+----------+---------+-------------+-----+--------------+----------------+------------------+--------------------+-------------+-------------+
|CONTIG    |START     |END      |EXTRA FIELD 1|...  |EXTRA FIELD N |MEAN            |MEDIAN            |STANDARD DEVIATION  |CpG          |HETEROZYGOUS |
|          |          |         |             |     |              |METHYLATION     |METHYLATION       |METHYLATION         |DINUCLEOTIDES|DINUCLEOTIDES|
+==========+==========+=========+=============+=====+==============+================+==================+====================+=============+=============+
|Contig    |Start     |Last     |Extra        |     |Extra         |Mean Methylation|Median Methylation|Standard Deviation  |Number of    |Number of CpG|
|Chromosome|Position  |Position |optional     |...  |optional      |value of CpG    |value of CpG      |methylation value   |CpG          |Dinucleotides|
|name.     |          |         |field one    |     |field N       |Dinucleotides   |Dinucleotides     |of CpG Dinucleotides|Dinucleotides|heterozygous.|
|          |          |         |             |     |              |present in the  |present in the    |present in the      |in the       |             |
|          |          |         |             |     |              |window.         |window.           |window.             |window.      |             |
+----------+----------+---------+-------------+-----+--------------+----------------+------------------+--------------------+-------------+-------------+



            


