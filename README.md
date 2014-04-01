HiTSeq
======

A Java package to process data from high-throughput sequencing assays



### Detailed description ###

Based on Picard, the Java package of SamTools, HiTSeq processes high-throughput sequencing data, especially RNA-Seq.

Several commands are included in HiTSeq, each of which is designed for different usage:

<table>
	<tr>
		<td><strong>Command</strong></td>
		<td><strong>Description</strong></td>
	<tr>
		<td><strong>info</strong></td>
		<td>Extract mapping information of high-throughput sequencing data, including the number of reads, mapped reads, uniquely mapped reads, and reads with different distance to reference.</td>
	</tr>
         
	<tr>
		<td><strong>uniq</strong></td>
		<td>Extract uniquely mapped reads of high-throughput sequencing data.</td>
	</tr>

	<tr>
		<td><strong>count/rpkm</strong></td>
		<td>Given genome features (e.g. gene annotation, peak intervals), count reads which are overlapping with each feature, followed by normalization using RPKM if necessary.</td>
	</tr>
</table>             


### Package structure ###

* <strong>HiTSeq.java</strong> | The main class.<br />
* <strong>Annotation.java</strong> | The class of annotation, i.e. genome features.<br/>
* <strong>ReadCounter.java</strong> | The class of read counter, to count reads overlapping with each genome feature.<br />
* <strong>SAMRecordProcessor.java</strong> | The class of record processor, to identify genome features that overlap with the<br />
* <strong>MappingProcessor.java</strong> | The class of SAM/BAM file processor, for processing which sees the whole SAM/BAM file as a whole, e.g. mapping information extraction, unique reads extraction.<br />



### Version ###

0.1



### Usage ###
<code>java -jar HiTSeq.jar <command> <...></code>
<br />	Or:<br />
<code>HiTSeq.sh <command> <...></code>
