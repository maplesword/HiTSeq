/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hitseq;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import net.sf.samtools.*;


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
    
    ReadCounter(File file, Annotation annotation, int strandSpecific, int modeForMultiGenesOverlap){
        if(! file.exists()){
            System.err.println("Cannot find input file: "+file.getAbsolutePath());
            System.exit(1);
        }
        this.inputFile=file;
        this.annotation=annotation;
        this.strandSpecific=strandSpecific;
        this.modeForMultiGenesOverlap=modeForMultiGenesOverlap;
        
        totalNumReads=0;
        this.calculatedRPKM=false;
        counts=new HashMap<>();
        rpkm=new HashMap<>();
        for(String gene : annotation.getGeneSet()){
            counts.put(gene, 0.0);
            rpkm.put(gene, 0.0);
        }
    }
    
    double getTotalNumReads(){
        return(totalNumReads);
    }
    
    HashMap<String, Double> getCounts(){
        return(counts);
    }
    
    HashMap<String, Double> getRPKM(){
        if(calculatedRPKM)
            return(rpkm);
        else
            return(null);
    }
    
    void initAnnotationPointers(){
        annotation.resetPointer();
    }
    
    boolean replaceInputFile(File newFile, int strandSpecific){
        if(! newFile.exists())
            return(false);
        else{
            this.inputFile=newFile;
            for(String gene : counts.keySet()){
                counts.put(gene, 0.0);
                rpkm.put(gene, 0.0);
            }
            this.calculatedRPKM=false;
            this.totalNumReads=0;
            this.strandSpecific=strandSpecific;
            return(true);
        }
    }
    
    void replaceAnnotation(Annotation newAnnotation){
        this.annotation=newAnnotation;
        counts=new HashMap<>();
        rpkm=new HashMap<>();
        for(String gene : annotation.getGeneSet()){
            counts.put(gene, 0.0);
            rpkm.put(gene, 0.0);
        }
        this.calculatedRPKM=false;
    }
    
    void setStrandSpecific(int newStrandSpecific){
        this.strandSpecific=newStrandSpecific;
    }
    
    void estimateCounts(boolean considerNHAttrib, boolean readCollapse, int convergLimit){
        if(modeForMultiGenesOverlap<3)
            estimateCountsSimply(considerNHAttrib, readCollapse);
        else if(modeForMultiGenesOverlap==3)
            estimateCountsIteratively(considerNHAttrib, readCollapse, convergLimit);
    }
    
    private void estimateCountsSimply(boolean considerNHAttrib, boolean readCollapse, int modeForMultiGenesOverlap, boolean markAmbiguous){
        totalNumReads=0;
        int numNoFeature=0, numAmbiguous=0;
        try (SAMFileReader inputSam = new SAMFileReader(inputFile)) {
            File outAmbiguous;
            SAMFileWriter outputSam = null;
            if(markAmbiguous){
                outAmbiguous=new File(inputFile.getAbsolutePath()+"-ambiguous.bam");
                outputSam=new SAMFileWriterFactory().makeSAMOrBAMWriter(inputSam.getFileHeader(), true, outAmbiguous);
            }
            
            String chromLast="";
            String strandLast="";
            String cigarLast="";
            int alignmentStartLast=-1;
            
            for(SAMRecord record : inputSam){
                if(record.getReadUnmappedFlag()) // skip if this read is unmapped
                    continue;
                if(record.getReadPairedFlag() && (! record.getProperPairFlag())) // skip if the read if paired but not in the proper paired mapping
                    continue;
                
                List<AlignmentBlock> hitsList=record.getAlignmentBlocks();
                ArrayList<AlignmentBlock> hits=new ArrayList<>();
                for(AlignmentBlock block : hitsList)
                    hits.add(block);
                int alignmentStart=hits.get(0).getReferenceStart();
                
                String chrom=record.getReferenceName();
                String strand=record.getReadNegativeStrandFlag() ? "-" : "+";
                String cigar=record.getCigarString();
                
                // remove PCR artifact if necessary
                if(readCollapse){
                    if(chromLast.equals(chrom) && strandLast.equals(strand) && cigarLast.equals(cigar) && alignmentStartLast==alignmentStart){
                        continue;
                    }
                    else{
                        chromLast=chrom;
                        strandLast=strand;
                        cigarLast=cigar;
                        alignmentStartLast=alignmentStart;
                    }
                }
                
                // count total reads
                if(considerNHAttrib && record.getIntegerAttribute("NH")!=null)
                    totalNumReads+=1.0/record.getIntegerAttribute("NH");
                else
                    totalNumReads++;
                if(java.lang.Math.ceil(totalNumReads)%1000000==0)
                    System.err.println("reading reads "+Double.valueOf(java.lang.Math.ceil(totalNumReads)).longValue()+"...");
                
                // skip if no gene exists in this chromosome in annotation
                if(! annotation.chromIsExisted(chrom)){
                    numNoFeature++;
                    continue;
                }
                
                // get overlapping genes for the record, add count to the genes
                SAMRecordProcessor recordProcessor=new SAMRecordProcessor(record, annotation);
                ArrayList<String> overlappedGenes=recordProcessor.getOverlapGenes(strandSpecific);
                
                double add=1;
                if(considerNHAttrib && record.getIntegerAttribute("NH")!=null)
                    add=add/record.getIntegerAttribute("NH");
                if(overlappedGenes.isEmpty()){
                    numNoFeature++;
                } else if(overlappedGenes.size()==1 || modeForMultiGenesOverlap==1){ // Mode 1: For the multi-genes hits, equally assign 1/n to each gene
                    add/=overlappedGenes.size();
                    for(String gene : overlappedGenes)
                        counts.put(gene, counts.get(gene)+add);
                } else if(overlappedGenes.size()>1 && modeForMultiGenesOverlap==2){ // Mode 2: For the multi-genes hits, randomly assign to one of the gene
                    java.util.Random random=new java.util.Random();
                    int selected=java.lang.Math.abs(random.nextInt())%overlappedGenes.size();
                    counts.put(overlappedGenes.get(selected), counts.get(overlappedGenes.get(selected))+add);
                } else if(modeForMultiGenesOverlap==0){
                    numAmbiguous++;
                }
                
                if(markAmbiguous && overlappedGenes.size()>1)
                    outputSam.addAlignment(record);
            }
            inputSam.close();
            if(markAmbiguous)
                outputSam.close();
        }
        catch(Exception e){
            System.err.println("ERROR! "+e+" in estimateCountsSimply at ReadCounter!\n");
            System.exit(1);
        }
        
        System.err.println(inputFile.getAbsolutePath()+":");
        System.err.printf("%45s|          %d\n","Number of mapped reads",(int)totalNumReads);
        System.err.printf("%45s|          %d\n","Number of reads with assigned feature",(int)totalNumReads-numNoFeature-numAmbiguous);
        System.err.printf("%45s|          %d\n","Number of reads with no feature",numNoFeature);
        System.err.printf("%45s|          %d\n","Number of ambiguous reads",numAmbiguous);
    }
    
    private void estimateRPKM(int modeForMultiGenesOverlap){
        for(String gene : counts.keySet()){
            int length;
            if(modeForMultiGenesOverlap==0){
                if(annotation.getExclusiveGeneLengthNoStrand(gene)==-1)
                    annotation.estimateExclusiveGeneLength();
                length=strandSpecific==0 ? annotation.getExclusiveGeneLengthNoStrand(gene) : annotation.getExclusiveGeneLength(gene);
            } else
                length=annotation.getGeneLength(gene);
            
            if(length!=-1 && length!=0)
                rpkm.put(gene, counts.get(gene)*1000*1000000/totalNumReads/length);
            else if(length==0)
                rpkm.put(gene, 0.0);
        }
        calculatedRPKM=true;
    }
    
    void estimateCountsSimply(boolean considerNHAttrib, boolean readCollapse){
        estimateCountsSimply(considerNHAttrib, readCollapse, this.modeForMultiGenesOverlap, false);
    }
    
    void estimateRPKM(){
        estimateRPKM(this.modeForMultiGenesOverlap);
    }
    
    private void estimateCountsIteratively(boolean considerNHAttrib, boolean readCollapse, int convergLimit){
        estimateCountsSimply(considerNHAttrib, readCollapse, 0, true);
        estimateRPKM(0);
        
        boolean continueIterate=true;
        int iterationTime=0;
        HashMap<String, Double> countsBackup=new HashMap<>();
        while(continueIterate){
            System.err.println("start iteration: round "+String.valueOf(iterationTime+1));

            HashMap<String, Double> influencedGenesLastRPKM=new HashMap<>();
            annotation.resetPointer();
            
            try (SAMFileReader inputSam = new SAMFileReader(new File(inputFile.getAbsolutePath()+"-ambiguous.bam"))) {
                for(SAMRecord record : inputSam){
                    SAMRecordProcessor recordProcessor=new SAMRecordProcessor(record, annotation);
                    ArrayList<String> overlappedGenes=recordProcessor.getOverlapGenes(strandSpecific);

                    ArrayList<Double> cdfOverlappedGenes=new ArrayList<>();
                    for(String gene : overlappedGenes){
                        if(counts.get(gene).intValue()>0)
                            influencedGenesLastRPKM.put(gene, rpkm.get(gene));

                        if(cdfOverlappedGenes.isEmpty())
                            cdfOverlappedGenes.add(rpkm.get(gene));
                        else
                            cdfOverlappedGenes.add(cdfOverlappedGenes.get(cdfOverlappedGenes.size()-1)+rpkm.get(gene));
                        
                        if(!countsBackup.containsKey(gene))
                            countsBackup.put(gene, new Double(counts.get(gene).doubleValue()));
                    }

                    if(cdfOverlappedGenes.get(cdfOverlappedGenes.size()-1)>0){
                        double add=1;
                        if(considerNHAttrib && record.getIntegerAttribute("NH")!=null)
                            add=add/record.getIntegerAttribute("NH");

                        java.util.Random random=new java.util.Random();
                        double randomNum=random.nextDouble()*cdfOverlappedGenes.get(cdfOverlappedGenes.size()-1);
                        int selected=0;
                        while(randomNum>cdfOverlappedGenes.get(selected) && selected<cdfOverlappedGenes.size())
                            selected++;

                        counts.put(overlappedGenes.get(selected), counts.get(overlappedGenes.get(selected))+add);
                    }
                }

                inputSam.close();
            }
            catch(Exception e){
                System.err.println("ERROR! "+e+" in estimateCountsIteratively at ReadCounter!\n");
                System.exit(1);
            }
            
            estimateRPKM(2);
            
            int numConverg=0;
            for(String gene : influencedGenesLastRPKM.keySet()){
                Double rpkmThis=rpkm.get(gene);
                if(java.lang.Math.abs(rpkmThis-influencedGenesLastRPKM.get(gene))<0.1*influencedGenesLastRPKM.get(gene))
                    numConverg++;
            }
            System.err.println("Converged genes: "+numConverg+" out of "+influencedGenesLastRPKM.size());
            
            iterationTime++;
            if(numConverg>influencedGenesLastRPKM.keySet().size()*0.9 || iterationTime>=convergLimit)
                continueIterate=false;
            else{
                for(String gene : countsBackup.keySet())
                    counts.put(gene, new Double(countsBackup.get(gene).doubleValue()));
                //System.err.println("ENSG00000000457.9\t"+counts.get("ENSG00000000457.9").intValue());
            }
        }
        
        File file=new File(inputFile.getAbsolutePath()+"-ambiguous.bam");
        file.delete();
    }
}
