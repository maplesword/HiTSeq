/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hitseq;

import hitseq.annotation.Annotation;
import hitseq.annotation.Junction;
import hitseq.annotation.JunctionSet;
import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

/**
 *
 * @author hezhisong
 */
public class ReadCounter {

    private File inputFile;
    private Annotation annotation;
    private int strandSpecific;
    private int modeForMultiGenesOverlap;
    private boolean calculatedRPKM;
    private double totalNumReads;
    private HashMap<String, Double> counts;
    private HashMap<String, Double> rpkm;
    private JunctionSet junctions;
    private HashMap<Junction, Double> junctionCounts;

    /**
     * Generate a new ReadCounter object with given annotation
     *
     * @param file File object of the input BAM/SAM file
     * @param annotation Annotation object
     * @param strandSpecific The strandness of the input BAM/SAM file. 0 for no
     * strand information, 1 for with the same strand, -1 for with the opposite
     * strand.
     * @param modeForMultiGenesOverlap The mode to process ambiguous reads.
     */
    public ReadCounter(File file, Annotation annotation, int strandSpecific, int modeForMultiGenesOverlap) {
        if (!file.exists()) {
            System.err.println("Cannot find input file: " + file.getAbsolutePath());
            System.exit(1);
        }
        this.inputFile = file;
        this.annotation = annotation;
        this.strandSpecific = strandSpecific;
        this.modeForMultiGenesOverlap = modeForMultiGenesOverlap;

        totalNumReads = 0;
        this.calculatedRPKM = false;
        counts = new HashMap<>();
        rpkm = new HashMap<>();
        for (String gene : annotation.getGeneSet()) {
            counts.put(gene, 0.0);
            rpkm.put(gene, 0.0);
        }

        junctions = new JunctionSet(annotation);
        junctionCounts = new HashMap<>();
        for (String chrom : junctions.getJunctions().keySet()) {
            for (Junction junction : junctions.getJunctions().get(chrom)) {
                junctionCounts.put(junction, 0.0);
            }
        }
    }

    /**
     * Generate a ReadCounter object with the given JunctionSet annotation.
     * Notice: ReadCounter objects generated using this method can only be used
     * to quantify read count of junctions.
     *
     * @param file File object of the input BAM/SAM file
     * @param junctions JunctionSet annotation.
     * @param strandSpecific The strandness of the input BAM/SAM file. 0 for no
     * strand information, 1 for with the same strand, -1 for with the opposite
     * strand.
     */
    public ReadCounter(File file, JunctionSet junctions, int strandSpecific) {
        if (!file.exists()) {
            System.err.println("Cannot find input file: " + file.getAbsolutePath());
            System.exit(1);
        }
        this.inputFile = file;
        this.strandSpecific = strandSpecific;

        this.junctions = junctions;
        this.junctionCounts = new HashMap<>();
        for (String chrom : junctions.getJunctions().keySet()) {
            for (Junction junction : junctions.getJunctions().get(chrom)) {
                junctionCounts.put(junction, 0.0);
            }
        }
    }

    public double getTotalNumReads() {
        return (totalNumReads);
    }

    public HashMap<String, Double> getCounts() {
        return (counts);
    }

    HashMap<String, Double> getRPKM() {
        if (calculatedRPKM) {
            return (rpkm);
        } else {
            return (null);
        }
    }

    void initAnnotationPointers() {
        annotation.resetPointer();
    }

    /**
     * FunName: replaceInputFile. Description: Replace the input file with the
     * given new one, and initiate all the calculation results.
     *
     * @param newFile The new input BAM/SAM file
     * @param strandSpecific The strandness of the new BAM/SAM file
     * @return true if the new file exists.
     */
    boolean replaceInputFile(File newFile, int strandSpecific) {
        if (!newFile.exists()) {
            return (false);
        } else {
            this.inputFile = newFile;
            for (String gene : counts.keySet()) {
                counts.put(gene, 0.0);
                rpkm.put(gene, 0.0);
            }
            this.calculatedRPKM = false;
            this.totalNumReads = 0;
            this.strandSpecific = strandSpecific;
            return (true);
        }
    }

    /**
     * FunName: replaceAnnotation. Description: Replace the annotation with the
     * given new one, and initiate all the calculation results.
     *
     * @param newAnnotation The object of the new annotation.
     */
    void replaceAnnotation(Annotation newAnnotation) {
        this.annotation = newAnnotation;
        counts = new HashMap<>();
        rpkm = new HashMap<>();
        for (String gene : annotation.getGeneSet()) {
            counts.put(gene, 0.0);
            rpkm.put(gene, 0.0);
        }
        this.calculatedRPKM = false;
        this.totalNumReads = 0;
    }

    void setStrandSpecific(int newStrandSpecific) {
        this.strandSpecific = newStrandSpecific;
    }

    /**
     * FunName: estimateCounts. Description: Count reads which are overlapping
     * with each feature (gene, peak, etc) in the annotation.
     *
     * @param considerNHAttrib If true, use NH attribute in the SAM record to
     * separate uniquely and multiply mapped reads.
     * @param readCollapse If true, collapse the plausible PCR duplicates, i.e.
     * reads with the same mapping chromosome, strand, coordinate and cigar,
     * into one read.
     * @param convergLimit The upper limit of iteration time if the counter is
     * set to use iteration to better assign ambiguous reads to genes.
     */
    public void estimateCounts(boolean considerNHAttrib, boolean readCollapse, int convergLimit) {
        if (modeForMultiGenesOverlap < 3) {
            estimateCountsSimply(considerNHAttrib, readCollapse);
        } else if (modeForMultiGenesOverlap == 3) {
            estimateCountsIteratively(considerNHAttrib, readCollapse, convergLimit);
        }
    }

    /**
     * FunName: estimateCountsSimply. Description: count reads which are
     * overlapping with each feature (gene, peak, etc) in the annotation, no
     * iteration.
     *
     * @param considerNHAttrib If true, use NH attribute in the SAM record to
     * separate uniquely and multiply mapped reads.
     * @param readCollapse If true, use NH attribute in the SAM record to
     * separate uniquely and multiply mapped reads.
     * @param modeForMultiGenesOverlap The mode to process ambiguous reads. 0
     * for simply abandon, 1 for equally assign, 2 for assign with equal
     * probability.
     * @param markAmbiguous If true, output the ambiguous reads to a temp BAM
     * file.
     */
    private void estimateCountsSimply(boolean considerNHAttrib, boolean readCollapse, int modeForMultiGenesOverlap, boolean markAmbiguous) {
        totalNumReads = 0;
        int numNoFeature = 0, numAmbiguous = 0;
        boolean notifyPairedCollapse = false;

        try (SamReader inputSam = SamReaderFactory.makeDefault().open(inputFile)) {
            inputSam.getFileHeader().setSortOrder(SAMFileHeader.SortOrder.coordinate); // for the possible SAM/BAM file for the ambiguous reads when necessary
            SAMRecordIterator iterator=inputSam.iterator();
            iterator.assertSorted(SAMFileHeader.SortOrder.coordinate); // the reads should be sorted by coordinate, for a exception will be thrown out
            
            File outAmbiguous;
            SAMFileWriter outputSam = null;
            if (markAmbiguous) {
                outAmbiguous = new File(inputFile.getAbsolutePath() + "-ambiguous.bam");
                outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(inputSam.getFileHeader(), false, outAmbiguous);
            }

            String chromLast = "";
            String strandLast = "";
            String cigarLast = "";
            int alignmentStartLast = -1;
            HashMap<String, HashSet<String>> overlapGenesForPair = new HashMap<>();
            HashMap<String, SAMRecord> firstRecordOfPair = new HashMap<>();

            while(iterator.hasNext()){
                SAMRecord record=iterator.next();
                if (record.getReadUnmappedFlag()) // skip if this read is unmapped
                    continue;
                if (record.getReadPairedFlag() && (!record.getProperPairFlag())) // skip if the read if paired but not in the proper paired mapping
                    continue;

                List<AlignmentBlock> hitsList = record.getAlignmentBlocks();
                ArrayList<AlignmentBlock> hits = new ArrayList<>();
                for (AlignmentBlock block : hitsList) {
                    hits.add(block);
                }
                int alignmentStart = hits.get(0).getReferenceStart();

                String chrom = record.getReferenceName();
                String strand = record.getReadNegativeStrandFlag() ? "-" : "+";
                String cigar = record.getCigarString();

                // remove PCR artifact if necessary (only work for single-ended data)
                if (readCollapse && !record.getReadPairedFlag()) {
                    if (chromLast.equals(chrom) && strandLast.equals(strand) && cigarLast.equals(cigar) && alignmentStartLast == alignmentStart) {
                        continue;
                    } else {
                        chromLast = chrom;
                        strandLast = strand;
                        cigarLast = cigar;
                        alignmentStartLast = alignmentStart;
                    }
                } else if (readCollapse && record.getReadPairedFlag() && !notifyPairedCollapse) {
                    System.err.println("WARNING: Redundant reads removal is not supported for paired-ended RNA-seq data yet.");
                    notifyPairedCollapse = true;
                }

                // count total reads
                if ((record.getReadPairedFlag() && record.getSecondOfPairFlag()) || !record.getReadPairedFlag()) {
                    if (considerNHAttrib && record.getIntegerAttribute("NH") != null) {
                        totalNumReads += 1.0 / record.getIntegerAttribute("NH");
                    } else {
                        totalNumReads++;
                    }
                    if (java.lang.Math.ceil(totalNumReads) % 1000000 == 0) {
                        System.err.println("reading reads " + Double.valueOf(java.lang.Math.ceil(totalNumReads)).longValue() + "...");
                    }
                }

                // skip if no gene exists in this chromosome in annotation
                if (!annotation.chromIsExisted(chrom)) {
                    numNoFeature++;
                    continue;
                }

                // get overlapping genes for the record
                SAMRecordProcessor recordProcessor = new SAMRecordProcessor(record, annotation);
                ArrayList<String> overlappedGenes = recordProcessor.getOverlapGenes(strandSpecific);

                if (record.getReadPairedFlag() && record.getFirstOfPairFlag()) { // if this is the first segment of a pair, denote its overlapping genes, and the record as well if needed
                    String id = record.getReadName();
                    id = id.replaceAll("[12]$", "");
                    id = id + " " + record.getMateReferenceName() + ":" + Integer.toString(record.getMateAlignmentStart());

                    if (!overlappedGenes.isEmpty()) {
                        HashSet<String> genes = new HashSet<>();
                        genes.addAll(overlappedGenes);
                        overlapGenesForPair.put(id, genes);
                    }
                    if (markAmbiguous) {
                        firstRecordOfPair.put(id, record);
                    }
                } else if (record.getReadPairedFlag() && record.getSecondOfPairFlag()) { // if this is the second segment of a pair, combine the overlapping genes of the two segments
                    String id = record.getReadName();
                    id = id.replaceAll("[12]$", "");
                    id = id + " " + record.getReferenceName() + ":" + Integer.toString(record.getAlignmentStart());

                    if (overlapGenesForPair.containsKey(id)) {
                        HashSet<String> genes = overlapGenesForPair.get(id);
                        genes.addAll(overlappedGenes);
                        overlappedGenes.clear();
                        overlappedGenes.addAll(genes);

                        overlapGenesForPair.remove(id);
                    }
                }

                // add count to the genes
                double add = 1;
                if (considerNHAttrib && record.getIntegerAttribute("NH") != null) {
                    add = add / record.getIntegerAttribute("NH");
                }
                if (overlappedGenes.isEmpty()) {
                    numNoFeature++;
                    if (record.getReadPairedFlag() && record.getFirstOfPairFlag()) {
                        numNoFeature--;
                    }
                } else if (overlappedGenes.size() == 1 || modeForMultiGenesOverlap == 1) { // Mode 1: For the multi-genes hits, equally assign 1/n to each gene
                    add /= overlappedGenes.size();
                    for (String gene : overlappedGenes) {
                        counts.put(gene, counts.get(gene) + add);
                    }
                } else if (overlappedGenes.size() > 1 && modeForMultiGenesOverlap == 2) { // Mode 2: For the multi-genes hits, randomly assign to one of the gene
                    java.util.Random random = new java.util.Random();
                    int selected = java.lang.Math.abs(random.nextInt()) % overlappedGenes.size();
                    counts.put(overlappedGenes.get(selected), counts.get(overlappedGenes.get(selected)) + add);
                } else if (modeForMultiGenesOverlap == 0) { // Mode 0: See the multi-genes hits as ambiguous hits
                    numAmbiguous++;
                    if (record.getReadPairedFlag() && record.getFirstOfPairFlag()) {
                        numAmbiguous--;
                    }
                }

                if (markAmbiguous) { // Mode 3: For the multi-genes hits, save them for the iterative counting
                    if (record.getReadPairedFlag() && record.getSecondOfPairFlag()) {
                        String id = record.getReadName();
                        id = id.replaceAll("[12]$", "");
                        id = id + " " + record.getReferenceName() + ":" + Integer.toString(record.getAlignmentStart());
                        if (overlappedGenes.size() > 1) {
                            outputSam.addAlignment(firstRecordOfPair.get(id));
                            outputSam.addAlignment(record);
                        }
                        firstRecordOfPair.remove(id);
                    } else if (!record.getReadPairedFlag() && overlappedGenes.size() > 1) {
                        outputSam.addAlignment(record);
                    }
                }
            }
            iterator.close();
            CloserUtil.close(inputSam);
            if (markAmbiguous) {
                outputSam.close();
            }
        } catch (Exception e) {
            System.err.println("ERROR! " + e + " in estimateCountsSimply at ReadCounter!\n");
            System.exit(1);
        }

        System.err.println(inputFile.getAbsolutePath() + ":");
        System.err.printf("%45s|          %d\n", "Number of mapped reads", (int) totalNumReads);
        System.err.printf("%45s|          %d\n", "Number of reads with assigned feature", (int) totalNumReads - numNoFeature - numAmbiguous);
        System.err.printf("%45s|          %d\n", "Number of reads with no feature", numNoFeature);
        System.err.printf("%45s|          %d\n", "Number of ambiguous reads", numAmbiguous);
    }

    /**
     * FunName: estimateRPKM. Description: According to read counting result for
     * features in the annotation, calculate RPKM.
     *
     * @param modeForMultiGenesOverlap The mode to process ambiguous reads.
     */
    private void estimateRPKM(int modeForMultiGenesOverlap) {
        for (String gene : counts.keySet()) {
            int length;
            if (modeForMultiGenesOverlap == 0) {
                if (annotation.getExclusiveGeneLengthNoStrand(gene) == -1) {
                    annotation.estimateAmbiguousGeneRegions();
                }
                length = strandSpecific == 0 ? annotation.getExclusiveGeneLengthNoStrand(gene) : annotation.getExclusiveGeneLength(gene);
            } else {
                length = annotation.getGeneLength(gene);
            }

            if (length != -1 && length != 0) {
                rpkm.put(gene, counts.get(gene) * 1000 * 1000000 / totalNumReads / length);
            } else if (length == 0) {
                rpkm.put(gene, 0.0);
            }
        }
        calculatedRPKM = true;
    }

    /**
     * FunName: estimateCountsSimply. Description: Count reads which are
     * overlapping with each feature (gene, peak, etc) in the annotation, no
     * iteration. The mode to process ambiguous reads is the one set to the
     * ReadCounter object
     *
     * @param considerNHAttrib If true, use NH attribute in the SAM record to
     * separate uniquely and multiply mapped reads.
     * @param readCollapse If true, collapse the plausible PCR duplicates, i.e.
     * reads with the same mapping chromosome, strand, coordinate and cigar,
     * into one read.
     */
    void estimateCountsSimply(boolean considerNHAttrib, boolean readCollapse) {
        estimateCountsSimply(considerNHAttrib, readCollapse, this.modeForMultiGenesOverlap, false);
    }

    /**
     * FunName: estimateRPKM. Description: According to read counting result for
     * features in the annotation, calculate RPKM. The set mode in the
     * ReadCounter object will be used to process ambiguous reads.
     */
    void estimateRPKM() {
        estimateRPKM(this.modeForMultiGenesOverlap);
    }

    /**
     * FunName: estimateCountsIteratively. Description: Count reads which are
     * overlapping with each feature (gene, peak, etc) in the annotation, and
     * use iteration to better assign the ambiguous reads to one of its
     * overlapping feature.
     *
     * @param considerNHAttrib If true, use NH attribute in the SAM record to
     * separate uniquely and multiply mapped reads.
     * @param readCollapse If true, collapse the plausible PCR duplicates, i.e.
     * reads with the same mapping chromosome, strand, coordinate and cigar,
     * into one read.
     * @param convergLimit The upper limit of iteration time if the counter is
     * set to use iteration to better assign ambiguous reads to genes.
     */
    private void estimateCountsIteratively(boolean considerNHAttrib, boolean readCollapse, int convergLimit) {
        estimateCountsSimply(considerNHAttrib, readCollapse, 0, true);
        estimateRPKM(0);

        boolean continueIterate = true;
        int iterationTime = 0;
        HashMap<String, Double> countsBackup = new HashMap<>();
        while (continueIterate) {
            System.err.println("start iteration: round " + String.valueOf(iterationTime + 1));

            HashMap<String, Double> influencedGenesLastRPKM = new HashMap<>();
            annotation.resetPointer();

            try (SamReader inputSam = SamReaderFactory.makeDefault().open(new File(inputFile.getAbsolutePath() + "-ambiguous.bam"))) {
                inputSam.getFileHeader().setSortOrder(SAMFileHeader.SortOrder.coordinate);

                HashMap<String, HashSet<String>> overlapGenesForPair = new HashMap<>(); // overlapping genes of the first segment of a pair
                
                for (SAMRecord record : inputSam) {
                    SAMRecordProcessor recordProcessor = new SAMRecordProcessor(record, annotation);
                    ArrayList<String> overlappedGenes = recordProcessor.getOverlapGenes(strandSpecific);

                    // if paired-ended data, save if it's the first segment, union the ones of the first segment if it's the second one
                    if(record.getReadPairedFlag() && record.getFirstOfPairFlag() && !overlappedGenes.isEmpty()){
                        String id=record.getReadName();
                        id = id.replaceAll("[12]$", "");
                        id = id + " " + record.getMateReferenceName() + ":" + Integer.toString(record.getMateAlignmentStart());
                        HashSet<String> genes=new HashSet<>();
                        genes.addAll(overlappedGenes);
                        overlapGenesForPair.put(id, genes);
                    } else if(record.getReadPairedFlag() && record.getSecondOfPairFlag()){
                        String id=record.getReadName();
                        id = id.replaceAll("[12]$", "");
                        id = id + " " + record.getReferenceName() + ":" + Integer.toString(record.getAlignmentStart());
                        if(overlapGenesForPair.containsKey(id)){
                            HashSet<String> genes=overlapGenesForPair.get(id);
                            genes.addAll(overlappedGenes);
                            overlappedGenes.clear();
                            overlappedGenes.addAll(genes);
                            overlapGenesForPair.remove(id);
                        }
                    }
                    
                    if(!record.getReadPairedFlag() || record.getSecondOfPairFlag()){
                        ArrayList<Double> cdfOverlappedGenes = new ArrayList<>();
                        for (String gene : overlappedGenes) {
                            if (counts.get(gene).intValue() > 0) {
                                influencedGenesLastRPKM.put(gene, rpkm.get(gene));
                            }

                            if (cdfOverlappedGenes.isEmpty()) {
                                cdfOverlappedGenes.add(rpkm.get(gene));
                            } else {
                                cdfOverlappedGenes.add(cdfOverlappedGenes.get(cdfOverlappedGenes.size() - 1) + rpkm.get(gene));
                            }

                            if (!countsBackup.containsKey(gene)) {
                                countsBackup.put(gene, new Double(counts.get(gene).doubleValue()));
                            }
                        }

                        if (cdfOverlappedGenes.get(cdfOverlappedGenes.size() - 1) > 0) {
                            double add = 1;
                            if (considerNHAttrib && record.getIntegerAttribute("NH") != null) {
                                add = add / record.getIntegerAttribute("NH");
                            }

                            java.util.Random random = new java.util.Random();
                            double randomNum = random.nextDouble() * cdfOverlappedGenes.get(cdfOverlappedGenes.size() - 1);
                            int selected = 0;
                            while (randomNum > cdfOverlappedGenes.get(selected) && selected < cdfOverlappedGenes.size()) {
                                selected++;
                            }

                            counts.put(overlappedGenes.get(selected), counts.get(overlappedGenes.get(selected)) + add);
                        }
                    }
                }

                inputSam.close();
            } catch (Exception e) {
                System.err.println("ERROR! " + e + " in estimateCountsIteratively at ReadCounter!\n");
                System.exit(1);
            }

            estimateRPKM(2);

            int numConverg = 0;
            for (String gene : influencedGenesLastRPKM.keySet()) {
                Double rpkmThis = rpkm.get(gene);
                if (java.lang.Math.abs(rpkmThis - influencedGenesLastRPKM.get(gene)) < 0.1 * influencedGenesLastRPKM.get(gene)) {
                    numConverg++;
                }
            }
            System.err.println("Converged genes: " + numConverg + " out of " + influencedGenesLastRPKM.size());

            iterationTime++;
            if (numConverg > influencedGenesLastRPKM.keySet().size() * 0.9 || iterationTime >= convergLimit) {
                continueIterate = false;
            } else {
                for (String gene : countsBackup.keySet()) {
                    counts.put(gene, new Double(countsBackup.get(gene).doubleValue()));
                }
                //System.err.println("ENSG00000000457.9\t"+counts.get("ENSG00000000457.9").intValue());
            }
        }

        File file = new File(inputFile.getAbsolutePath() + "-ambiguous.bam");
        file.delete();
    }

    /**
     * FunName: estimateJunctionCounts. Description: Count the number of
     * junction reads supporting each junction in the junction list.
     *
     * @param considerNHAttrib If true, use NH attribute in the SAM record to
     * separate uniquely and multiply mapped reads.
     * @param readCollapse If true, collapse the plausible PCR duplicates, i.e.
     * reads with the same mapping chromosome, strand, coordinate and cigar,
     * into one read.
     */
    void estimateJunctionCounts(boolean considerNHAttrib, boolean readCollapse) {
        int numJuncReads = 0;
        boolean notifyPairedCollapse=false;
        HashMap<String, HashSet<Junction>> juncListChrom = junctions.getJunctions();
        try (SamReader inputSam = SamReaderFactory.makeDefault().open(inputFile)) {
            String chromLast = "";
            String strandLast = "";
            String cigarLast = "";
            int alignmentStartLast = -1;

            for (SAMRecord record : inputSam) {
                if (record.getReadUnmappedFlag()) // skip if this read is unmapped
                {
                    continue;
                }
                if (record.getReadPairedFlag() && (!record.getProperPairFlag())) // skip if the read if paired but not in the proper paired mapping
                {
                    continue;
                }

                List<AlignmentBlock> hitsList = record.getAlignmentBlocks();
                ArrayList<AlignmentBlock> hits = new ArrayList<>();
                for (AlignmentBlock block : hitsList) {
                    hits.add(block);
                }

                // remove PCR artifact if necessary
                String chrom = record.getReferenceName();
                String strand = record.getReadNegativeStrandFlag() ? "-" : "+";
                String cigarString = record.getCigarString();
                int alignmentStart = hits.get(0).getReferenceStart();
                if (readCollapse && ! record.getReadPairedFlag()) {
                    if (chromLast.equals(chrom) && strandLast.equals(strand) && cigarLast.equals(cigarString) && alignmentStartLast == alignmentStart) {
                        continue;
                    } else {
                        chromLast = chrom;
                        strandLast = strand;
                        cigarLast = cigarString;
                        alignmentStartLast = alignmentStart;
                    }
                } else if(readCollapse && record.getReadPairedFlag() && !notifyPairedCollapse){
                    System.err.println("WARNING: Redundant reads removal is not supported for paired-ended RNA-seq data yet.");
                    notifyPairedCollapse = true;
                }

                // count total reads
                double add = 1.0;
                if (considerNHAttrib && record.getIntegerAttribute("NH") != null) {
                    add /= record.getIntegerAttribute("NH");
                }
                totalNumReads += add;
                if (java.lang.Math.ceil(totalNumReads) % 1000000 == 0) {
                    System.err.println("reading reads " + (int) java.lang.Math.ceil(totalNumReads) + "...");
                }

                // determine junction strand
                strand = "*";
                if (strandSpecific == 1) {
                    strand = record.getReadNegativeStrandFlag() ? "-" : "+";
                } else if (strandSpecific == -1) {
                    strand = record.getReadNegativeStrandFlag() ? "+" : "1";
                }

                // get cigar
                Cigar cigar = record.getCigar();
                if (cigar.numCigarElements() == 1) {
                    continue;
                }

                // determine whether it is junction read, and count the junction
                int lastBlockEnd = hits.get(0).getReferenceStart() - 1;
                boolean containJunction = false;

                for (int i = 0; i < cigar.numCigarElements(); i++) {
                    CigarElement element = cigar.getCigarElement(i);
                    if (element.getOperator().consumesReferenceBases()) {
                        int thisBlockEnd = lastBlockEnd + element.getLength();
                        if (element.getOperator().equals(CigarOperator.N)) {
                            containJunction = true;

                            // skip reads to unannotated chromosomes
                            if (!juncListChrom.containsKey(chrom)) {
                                break;
                            }

                            Junction junc = new Junction(chrom, strand, lastBlockEnd, thisBlockEnd + 1);
                            // count it only the junction is in the junction list with the correct strand.
                            if (juncListChrom.get(chrom).contains(junc)) {
                                if (junctionCounts.containsKey(junc)) {
                                    junctionCounts.put(junc, junctionCounts.get(junc) + add);
                                } else {
                                    junctionCounts.put(junc, add);
                                }
                            } else if (strand.equals("*")) {
                                /*
                                 * If the library has no strand information,
                                 * i.e. the strand of the junction cannot be
                                 * determined, and there is no existed junction
                                 * with no strand information has the same
                                 * coordinates, then count it if either an
                                 * existed junction in the list has the same
                                 * coordinates and at the forward strand, or at
                                 * the reverse strand, but not both.
                                 */
                                Junction juncPos = new Junction(chrom, "+", lastBlockEnd, thisBlockEnd + 1);
                                Junction juncNeg = new Junction(chrom, "-", lastBlockEnd, thisBlockEnd + 1);
                                if (juncListChrom.get(chrom).contains(juncPos) && !juncListChrom.get(chrom).contains(juncNeg)) {
                                    if (junctionCounts.containsKey(juncPos)) {
                                        junctionCounts.put(juncPos, junctionCounts.get(juncPos) + add);
                                    } else {
                                        junctionCounts.put(juncPos, add);
                                    }
                                } else if (juncListChrom.get(chrom).contains(juncNeg) && !juncListChrom.get(chrom).contains(juncPos)) {
                                    if (junctionCounts.containsKey(juncNeg)) {
                                        junctionCounts.put(juncNeg, junctionCounts.get(juncNeg) + add);
                                    } else {
                                        junctionCounts.put(juncNeg, add);
                                    }
                                }
                            }
                        }
                        lastBlockEnd = thisBlockEnd;
                    }
                }

                if (containJunction) {
                    numJuncReads += add;
                    //System.err.println(numJuncReads);
                } else {
                    //System.err.println(record.getReadName()+"\t"+record.getCigarString());
                }
            }

            CloserUtil.close(inputSam);
        } catch (Exception e) {
            System.err.println("Error in ReadCounter when estimate junction count: " + e);
            System.exit(1);
        }

        System.err.println(inputFile.getAbsolutePath() + ":");
        System.err.printf("%45s|          %d\n", "Number of mapped reads", (int) totalNumReads);
        System.err.printf("%45s|          %d\n", "Number of junction reads", (int) numJuncReads);
    }

    HashMap<Junction, Double> getJunctionCounts() {
        return (junctionCounts);
    }
}
