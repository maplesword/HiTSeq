HiTSeq
======

A Java package to process data from high-throughput sequencing assays

======

Detailed description

Based on Picard, the Java package of SamTools, HiTSeq processes high-throughput sequencing data, especially RNA-Seq.

Several commands are included in HiTSeq, each of which is designed for different usage:

info         Extract mapping information of high-throughput sequencing data, including the number of reads,
             mapped reads, uniquely mapped reads, and reads with different distance to reference.
         
uniq         Extract uniquely mapped reads of high-throughput sequencing data.

count/rpkm   Given genome features (e.g. gene annotation, peak intervals), count reads which are overlapping with
             each feature, followed by normalization using RPKM if necessary.
             

======

Package structure

HiTSeq.java                   The main class.
|
|-Annotation.java             The class of annotation, i.e. genome features.
|
|-ReadCounter.java            The class of read counter, to count reads overlapping with each genome feature.
| |
| |-SAMRecordProcessor.java   The class of record processor, to identify genome features that overlap with the
|                             current read.
|
|-MappingProcessor.java       The class of SAM/BAM file processor, for processing which sees the whole SAM/BAM file
                              as a whole, e.g. mapping information extraction, unique reads extraction.


======

Version

0.1
