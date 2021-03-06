/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hitseq;

import hitseq.annotation.Annotation;
import hitseq.annotation.Gene;
import hitseq.annotation.Transcript;
import hitseq.annotation.Exon;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import java.util.ArrayList;
import java.util.List;
/**
 *
 * @author hezhisong
 */
public class SAMRecordProcessor {
    SAMRecord record;
    Annotation annotation;
    
    SAMRecordProcessor(SAMRecord record, Annotation annotation){
        this.record=record;
        this.annotation=annotation;
    }
    
    void replaceNewSAMRecord(SAMRecord newRecord){
        this.record=newRecord;
    }
    
    void replaceNewAnnotation(Annotation newAnnotation){
        this.annotation=newAnnotation;
    }
    
    SAMRecord getSAMRecord(){
        return(this.record);
    }
    
    Annotation getAnnotation(){
        return(this.annotation);
    }
    
    String recordToString(){
        String readName = record.getReadName();
        readName = readName.replaceAll("[12]$", "");
        String readRef = record.getReferenceName();
        String readStart = Integer.toString(record.getAlignmentStart());
        String readStrand = Boolean.toString(! record.getReadNegativeStrandFlag());
        return(readName+"|"+readRef+"|"+readStart+"|"+readStrand);
    }
    
    String recordMateToString(){
        String readName = record.getReadName();
        readName = readName.replaceAll("[12]$", "");
        String mateRef = record.getMateReferenceName();
        String mateStart = Integer.toString(record.getMateAlignmentStart());
        String mateStrand = Boolean.toString(! record.getMateNegativeStrandFlag());
        return(readName+"|"+mateRef+"|"+mateRef+"|"+mateStrand);
    }
    
    ArrayList<String> getOverlapGenes(int strandSpecific){
        ArrayList<String> overlappedGenes=new ArrayList<>();
        List<AlignmentBlock> hitsList=record.getAlignmentBlocks();
        ArrayList<AlignmentBlock> hits=new ArrayList<>();
        for(AlignmentBlock block : hitsList)
            hits.add(block);
        int alignmentStart=hits.get(0).getReferenceStart();
        int alignmentEnd=hits.get(hits.size()-1).getReferenceStart()+hits.get(hits.size()-1).getLength()-1;

        String chrom=record.getReferenceName();
        String strand=record.getReadNegativeStrandFlag() ? "-" : "+";
        
        int pointer=annotation.getCurrentGenePointer(chrom);
        if(pointer!=-1){
            String currentGene=annotation.getGene(chrom, pointer);
            int currentGeneEnd=annotation.getGeneEnd(currentGene);
            while(currentGeneEnd < alignmentStart){
                pointer=annotation.movePointerGene(chrom);
                if(pointer==-1)
                    break;
                currentGene=annotation.getGene(chrom, annotation.getCurrentGenePointer(chrom));
                currentGeneEnd=annotation.getGeneEnd(currentGene);
            }
            int currentGeneStart=annotation.getGeneStart(currentGene);
            int additionalIndex=0;
            while(currentGeneStart<=alignmentEnd){ 
                // this gene regions is overlapping with the read region
                if((alignmentStart <= currentGeneStart && currentGeneStart <= alignmentEnd) || (currentGeneStart <= alignmentStart && alignmentStart <= currentGeneEnd))
                    if(record.getReadPairedFlag()){
                        if(strandSpecific==0
                                || (strandSpecific==1 && record.getFirstOfPairFlag() && strand.equals(annotation.getGeneStrand(currentGene)))
                                || (strandSpecific==1 && record.getSecondOfPairFlag() && (! strand.equals(annotation.getGeneStrand(currentGene))))
                                || (strandSpecific==-1 && record.getFirstOfPairFlag() && (! strand.equals(annotation.getGeneStrand(currentGene))))
                                || (strandSpecific==-1 && record.getSecondOfPairFlag() && strand.equals(annotation.getGeneStrand(currentGene))))
                                overlappedGenes.add(currentGene);
                    } else if((strandSpecific==0) 
                            || (strandSpecific==1 && strand.equals(annotation.getGeneStrand(currentGene))) 
                            || (strandSpecific==-1 && (!strand.equals(annotation.getGeneStrand(currentGene)))))
                        overlappedGenes.add(currentGene);

                additionalIndex++;
                currentGene=annotation.getGene(chrom, annotation.getCurrentGenePointer(chrom)+additionalIndex);
                if(currentGene==null)
                    break;
                currentGeneStart=annotation.getGeneStart(currentGene);
                currentGeneEnd=annotation.getGeneEnd(currentGene);
            }
        }

        // for each gene in the candidate set, check whether the alignments are overlapping with the exon regions
        ArrayList<String> notExonicOverlap=new ArrayList<>();
        for(String gene : overlappedGenes){
            //System.out.println(record.getReadName()+"\t"+gene);
            int currentExonIndex=annotation.getCurrentExonIndex(gene);
            int additionalIndexExon=0;
            int indexAlignmentBlock=0;

            boolean overlap=false;
            while(! overlap && currentExonIndex!=-1 && indexAlignmentBlock < hits.size()){
                int currentExonStart=annotation.getNonRedundantExonStart(gene, currentExonIndex+additionalIndexExon);
                int currentExonEnd=annotation.getNonRedundantExonEnd(gene, currentExonIndex+additionalIndexExon);
                int currentAlignmentBlockStart=hits.get(indexAlignmentBlock).getReferenceStart();
                int currentAlignmentBlockEnd=hits.get(indexAlignmentBlock).getReferenceStart()+hits.get(indexAlignmentBlock).getLength()-1;

                if(currentExonStart==-1 || currentExonEnd==-1)
                    break;
                else if(currentExonEnd < alignmentStart) // the current exon is in front of all the alignment blocks
                    currentExonIndex=annotation.movePointerExon(gene);
                else if(currentExonEnd < currentAlignmentBlockStart) // the current exon is in front of the current alignment block, while the current alignment block is not the first one
                    additionalIndexExon++;                        
                else if(currentAlignmentBlockEnd < currentExonStart) // the current exon is after the current alignment block
                    indexAlignmentBlock++;
                else if(currentExonStart <= currentAlignmentBlockStart && currentAlignmentBlockStart <= currentExonEnd) // overlapping
                    overlap=true;
                else if(currentAlignmentBlockStart <= currentExonStart && currentExonStart <= currentAlignmentBlockEnd) // overlapping
                    overlap=true;
            }
            
            if(!overlap)
                notExonicOverlap.add(gene);
        }
        overlappedGenes.removeAll(notExonicOverlap);
        return(overlappedGenes);
    }
    
    double[] getGeneBodyQuantile(Gene gene, int strandSpecific){
        double[] quantiles;
        String chrom = record.getReferenceName();
        String strand = record.getReadNegativeStrandFlag() ? "-" : "+";
        if(strandSpecific == -1)
            strand = strand.equals("+") ? "-" : "+";
        
        if(! gene.getChrom().equals(chrom) ||
                (strandSpecific != 0 && ! strand.equals(gene.getStrand()))){
            quantiles = new double[]{-1, -1};
        } else{
            int recordStart = record.getAlignmentStart();
            int recordEnd = record.getAlignmentEnd();
            Transcript recordTranscript = new Transcript(record.getReadName(), chrom, gene.getStrand(), recordStart, recordEnd);
            for(int i = 0; i < record.getAlignmentBlocks().size(); i++){
                AlignmentBlock block = record.getAlignmentBlocks().get(i);
                Exon exonBlock = new Exon(chrom, gene.getStrand(), block.getReferenceStart(), block.getReferenceStart() + block.getLength() - 1);
                recordTranscript.addExon(exonBlock);
            }
            
            ArrayList<Exon> overlap = gene.getNonredundantTranscript().getOverlappingRegions(recordTranscript, strandSpecific != 0);
            if(! overlap.isEmpty()){
                int first = overlap.get(0).getStart();
                int last = overlap.get(overlap.size()-1).getEnd();
                quantiles = new double[]{gene.getNonredundantTranscript().getQuantileAtPosition(chrom, gene.getStrand(), first), gene.getNonredundantTranscript().getQuantileAtPosition(chrom, gene.getStrand(), last)};
                if(gene.getStrand() == "-"){
                    double temp = quantiles[1];
                    quantiles[1] = quantiles[0];
                    quantiles[0] = temp;
                }
            } else{
                quantiles = new double[]{-1, -1};
            }
        }
        
        return(quantiles);
    }
}
