HiTSeq
======

A Java package to process data from high-throughput sequencing assays



### Detailed description ###

Based on htsjdk, the Java library to manipulate SAM/BAM file, HiTSeq processes high-throughput sequencing data, especially RNA-Seq.

Several commands are included in HiTSeq, each of which is designed for different usage:
<table>
	<tr>
		<td><strong>Command</strong></td>
		<td><strong>Description</strong></td>
	</tr>
	<tr>
		<td><strong>count/rpkm</strong></td>
		<td>Given genome features (e.g. gene annotation, peak intervals), count reads which are overlapping with each feature, followed by normalization using RPKM if necessary.</td>
	</tr>
	<tr>
		<td><strong>bias</strong></td>
		<td>Given genome features (e.g. gene annotation, peak intervals), count reads covering each percentile of genes from 5'-end to 3'-end, to estimate degree of 3'-bias in the data as a proxy of RNA quality.</td>
	</tr>
	<tr>
		<td><strong>demultiplex</strong></td>
		<td>Given a list of sites with varied bases in different species/groups, scan for unique UMIs covering the sites and assign them to the species. Useful when cells/nuclei from close species (human and chimp) are pooled for single cell RNA-seq.</td>
	</tr>
	<tr>
		<td><strong>tostruc</strong></td>
		<td>Convert GTF formatted transcriptome annotation into the simplified struc format by discarding the isoform information. Much smaller and faster to read.</td>
	</tr>
	<tr>
		<td><strong>tojuncs</strong></td>
		<td>Combine and convert junction lists into a single junction list in 'juncs' format as provided by Tophat2.</td>
	</tr>
	<tr>
		<td><strong>toevents</strong></td>
		<td>Construct alternative splicing events based on the provided junction list.</td>
	</tr>
	<tr>
		<td><strong>countjunc</strong></td>
		<td>Given a junction or event list, count the number of reads representing each junction in the provided SAM/BAM file.</td>
	</tr>
	<tr>
		<td><strong>gui</strong> (under development)</td>
		<td>Open the GUI of HiTSeq. Only provide read counting function.</td>
	</tr>
	<tr>
		<td><strong>correct</strong> (inefficient)</td>
		<td>Correct the 'properpaired' flag for paired-ended RNA-seq data, by resetting two mates from the same transcript but far away from each other due to presence of huge introns as properly paired.</td>
	</tr>
	<tr>
		<td><strong>info</strong> (deprecated)</td>
		<td>Extract mapping information of high-throughput sequencing data, including the number of reads, mapped reads, uniquely mapped reads, and reads with different distance to reference.</td>
	</tr>
	<tr>
		<td><strong>uniq (deprecated)</strong></td>
		<td>Extract uniquely mapped reads of high-throughput sequencing data.</td>
	</tr>
</table>             

The list of commands can be seen by revoking the help page of HiTSeq using command: `java -jar HiTSeq.jar`. The help page of each command is also available by using `java -jar HiTSeq.jar <command> -h`


### Version ###

0.2


### Usage ###
<code>java -jar HiTSeq.jar <command> <...></code>
<br />	Or:<br />
<code>HiTSeq.sh <command> <...></code>
