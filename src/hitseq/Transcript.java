/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hitseq;

import java.util.ArrayList;

/**
 *
 * @author hezhisong
 */
public class Transcript {
    private String id;
    private String chrom;
    private String strand;
    private int start;
    private int end;
    private ArrayList<Exon> exons;
    private int length;
    
    Transcript(String id, String chrom, String strand, int start, int end){
        this.id=id;
        this.chrom=chrom;
        this.strand=strand;
        this.start=start;
        this.end=end;
        this.exons=new ArrayList<>();
        this.length=0;
    }
    
    Transcript(String id, String chrom, String strand){
        this(id, chrom, strand, -1, -1);
    }
    
    void addExon(Exon newExon){
        exons.add(newExon);
        if(start==-1 || newExon.getStart()<start)
            start=newExon.getStart();
        if(end==-1 || newExon.getEnd()>end)
            end=newExon.getEnd();
        
        this.length+=newExon.getLength();
    }
    
    void addExon(int exonStart, int exonEnd){
        Exon newExon=new Exon(chrom, strand, exonStart, exonEnd);
        addExon(newExon);
    }
    
    void addExons(ArrayList<Exon> newExons){
        for(Exon newExon : newExons)
            addExon(newExon);
    }
    
    void sortExons(){
        if(exons.size()>1){
            for(int i=1; i<exons.size(); i++){
                int j=0;
                while(j<i && exons.get(i).compareTo(exons.get(j))>0)
                    j++;
                exons.add(j, exons.get(i));
                exons.remove(i+1);
            }
        }
        
        this.start=exons.get(0).getStart();
        this.end=exons.get(exons.size()-1).getEnd();
    }
    
    void mergeExons(){
        if(exons.size()>1){
            sortExons();
            
            for(int i=0; i<exons.size()-1; i++){
                int exonEndThis=exons.get(i).getEnd();
                int exonStartNext=exons.get(i+1).getStart();
                int exonEndNext=exons.get(i+1).getEnd();

                if(exonStartNext <= exonEndThis){ // these two exons are overlapping
                    int newExonStart=exons.get(i).getStart();
                    int newExonEnd=exonEndThis > exonEndNext ? exonEndThis : exonEndNext;
                    Exon newExon=new Exon(chrom, strand, newExonStart, newExonEnd);
                    exons.add(i, newExon);
                    exons.remove(i+1);
                    exons.remove(i+1);
                    i--;
                }
                else if(exonStartNext == exonEndThis+1){
                    int newExonStart=exons.get(i).getStart();
                    int newExonEnd=exons.get(i+1).getEnd();
                    Exon newExon=new Exon(chrom, strand, newExonStart, newExonEnd);
                    exons.add(i, newExon);
                    exons.remove(i+1);
                    exons.remove(i+1);
                    i--;
                }
            }
            
            length=0;
            for(Exon exon : exons)
                length+=exon.getLength();
        }
    }
    
    String getID(){
        return(id);
    }
    
    String getChrom(){
        return(chrom);
    }
    
    String getStrand(){
        return(strand);
    }
    
    int getStart(){
        return(start);
    }
    
    int getEnd(){
        return(end);
    }
    
    ArrayList<Exon> getExons(){
        return(exons);
    }
    
    int getTotalExonLength(){
        return(length);
    }
    
    ArrayList<Junction> getJunctions(){
        sortExons();
        ArrayList<Junction> junctions=new ArrayList<>();
        if(exons.size()>1){
            for(int i=0; i<exons.size()-1; i++){
                if(exons.get(i).getEnd()<exons.get(i+1).getStart()){
                    Junction newJunction=new Junction(chrom, strand, exons.get(i).getEnd(), exons.get(i+1).getStart());
                    junctions.add(newJunction);
                }
            }
        }
        
        if(junctions.size()!=exons.size()-1)
            System.err.println("Warning: Junction number incompatible for transcript: "+id+" (should be "+(exons.size()-1)+" but "+junctions.size()+")");
        return(junctions);
    }
    
    ArrayList<Exon> getOverlappingRegions(Transcript transcript, boolean considerStrand){
        ArrayList<Exon> answer=new ArrayList<>();
        
        boolean sameStrand=strand.equals(transcript.getStrand());
        int idx1=0, idx2=0;
        ArrayList<Exon> thisExons=exons, thatExons=transcript.getExons();
        
        while(idx1<thisExons.size() && idx2<thatExons.size()){
            if(thisExons.get(idx1).getStart()<=thatExons.get(idx2).getStart() && thatExons.get(idx2).getStart()<=thisExons.get(idx1).getEnd()){
                if(thisExons.get(idx1).getEnd()<thatExons.get(idx2).getEnd()){
                    Exon overlap=new Exon(chrom, strand, thatExons.get(idx2).getStart(), thisExons.get(idx1).getEnd());
                    if((!considerStrand) || sameStrand)
                        answer.add(overlap);
                    idx1++;
                } else if(thisExons.get(idx1).getEnd()==thatExons.get(idx2).getEnd()){
                    Exon overlap=new Exon(chrom, strand, thatExons.get(idx2).getStart(), thisExons.get(idx1).getEnd());
                    if((!considerStrand) || sameStrand)
                        answer.add(overlap);
                    idx1++; idx2++;
                } else{
                    Exon overlap=new Exon(chrom, strand, thatExons.get(idx2).getStart(), thatExons.get(idx2).getEnd());
                    if((!considerStrand) || sameStrand)
                        answer.add(overlap);
                    idx2++;
                }
            } else if(thatExons.get(idx2).getStart()<=thisExons.get(idx1).getStart() && thisExons.get(idx1).getStart()<=thatExons.get(idx2).getEnd()){
                if(thatExons.get(idx2).getEnd()<thisExons.get(idx1).getEnd()){
                    Exon overlap=new Exon(chrom, strand, thisExons.get(idx1).getStart(), thatExons.get(idx2).getEnd());
                    if((!considerStrand) || sameStrand)
                        answer.add(overlap);
                    idx2++;
                } else if(thatExons.get(idx2).getEnd()==thisExons.get(idx1).getEnd()){
                    Exon overlap=new Exon(chrom, strand, thisExons.get(idx1).getStart(), thatExons.get(idx2).getEnd());
                    if((!considerStrand) || sameStrand)
                        answer.add(overlap);
                    idx2++; idx1++;
                } else{
                    Exon overlap=new Exon(chrom, strand, thisExons.get(idx1).getStart(), thisExons.get(idx1).getEnd());
                    if((!considerStrand) || sameStrand)
                        answer.add(overlap);
                    idx1++;
                }
            } else{
                if(thisExons.get(idx1).getStart()<thatExons.get(idx2).getStart())
                    idx1++;
                else
                    idx2++;
            }
        }
        
        return(answer);
    }
}
