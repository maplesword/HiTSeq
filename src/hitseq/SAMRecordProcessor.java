/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hitseq;

import java.util.ArrayList;
import java.util.List;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMRecord;
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
                    if((strandSpecific==0) || (strandSpecific==1 && strand.equals(annotation.getGeneStrand(currentGene))) || (strandSpecific==-1 && (!strand.equals(annotation.getGeneStrand(currentGene)))))
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
}
