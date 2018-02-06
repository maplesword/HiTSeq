/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hitseq.annotation;

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
    private String type;
    
    public Transcript(String id, String chrom, String strand, int start, int end, String type){
        this.id=id;
        this.chrom=chrom;
        this.strand=strand;
        this.start=start;
        this.end=end;
        this.exons=new ArrayList<>();
        this.length=0;
        this.type = type;
    }
    
    public Transcript(String id, String chrom, String strand, int start, int end){
        this(id, chrom, strand, start, end, "others");
    }
    
    public Transcript(String id, String chrom, String strand, String type){
        this(id, chrom, strand, -1, -1, type);
    }
    
    public Transcript(String id, String chrom, String strand){
        this(id, chrom, strand, -1, -1);
    }
    
    public void addExon(Exon newExon){
        exons.add(newExon);
        if(start==-1 || newExon.getStart()<start)
            start=newExon.getStart();
        if(end==-1 || newExon.getEnd()>end)
            end=newExon.getEnd();
        
        this.length+=newExon.getLength();
    }
    
    public void addExon(int exonStart, int exonEnd){
        Exon newExon=new Exon(chrom, strand, exonStart, exonEnd);
        addExon(newExon);
    }
    
    public void addExons(ArrayList<Exon> newExons){
        for(Exon newExon : newExons)
            addExon(newExon);
    }
    
    public void sortExons(){
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
    
    public void mergeExons(){
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
        
        int numOfExons=exons.size();
        if(numOfExons>0){
            this.start=exons.get(0).getStart();
            this.end=exons.get(numOfExons-1).getEnd();
        }
    }
    
    public String getID(){
        return(id);
    }
    
    public String getChrom(){
        return(chrom);
    }
    
    public String getStrand(){
        return(strand);
    }
    
    public int getStart(){
        return(start);
    }
    
    public int getEnd(){
        return(end);
    }
    
    public String getType(){
        return(type);
    }
    
    public Exon getExon(int idx){
        if(idx<0 || idx>exons.size()-1)
            return(null);
        else
            return(exons.get(idx));
    }
    
    public int getExonNumber(){
        return(exons.size());
    }
    
    public ArrayList<Exon> getExons(){
        return((ArrayList<Exon>) exons.clone());
    }
    
    public int getTotalExonLength(){
        return(length);
    }
    
    public ArrayList<Junction> getJunctions(){
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
    
    public ArrayList<Exon> getOverlappingRegions(Transcript transcript, boolean considerStrand){
        ArrayList<Exon> answer=new ArrayList<>();
        
        boolean sameStrand=strand.equals(transcript.getStrand());
        int idx1=0, idx2=0;
        ArrayList<Exon> thisExons=exons, thatExons=transcript.exons;
        
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
    
    public ArrayList<Exon> exclusiveRegions(Transcript transcript, boolean considerStrand){
        ArrayList<Exon> answer = new ArrayList<>();
        if(considerStrand && ! strand.equals(transcript.getStrand())){
            for(Exon exon : exons){
                Exon thisExon = new Exon(exon.getID(), exon.getChrom(), exon.getStrand(), exon.getStart(), exon.getEnd());
                answer.add(thisExon);
            }
            return answer;
        }
        
        ArrayList<Exon> thisExons = (ArrayList<Exon>)exons.clone();
        ArrayList<Exon> thatExons = transcript.getExons();
        int idx1=0, idx2=0;
        while(idx1<thisExons.size() && idx2<thatExons.size()){
            if(thisExons.get(idx1).getEnd() < thatExons.get(idx2).getStart()){
                answer.add(new Exon(thisExons.get(idx1).getChrom(), thisExons.get(idx1).getStrand(), thisExons.get(idx1).getStart(), thisExons.get(idx1).getEnd()));
                idx1++;
            } else if(thatExons.get(idx2).getEnd() < thisExons.get(idx1).getStart())
                idx2++;
            else if(thisExons.get(idx1).getStart() < thatExons.get(idx2).getStart()){
                answer.add(new Exon(thisExons.get(idx1).getChrom(), thisExons.get(idx1).getStrand(), thisExons.get(idx1).getStart(), thatExons.get(idx2).getStart()-1));
                if(thisExons.get(idx1).getEnd() < thatExons.get(idx2).getEnd()){
                    idx1++;
                } else if(thisExons.get(idx1).getEnd() == thatExons.get(idx2).getEnd()){
                    idx1++;
                    idx2++;
                } else{
                    thisExons.set(idx1, new Exon(thisExons.get(idx1).getChrom(), thisExons.get(idx1).getStrand(), thatExons.get(idx2).getEnd()+1, thisExons.get(idx1).getEnd()));
                    idx2++;
                }
            }
            else if(thatExons.get(idx2).getStart() <= thisExons.get(idx1).getStart()){
                if(thisExons.get(idx1).getEnd() <= thatExons.get(idx2).getEnd())
                    idx1++;
                else{
                    thisExons.set(idx1, new Exon(thisExons.get(idx1).getChrom(), thisExons.get(idx1).getStrand(), thatExons.get(idx2).getEnd()+1, thisExons.get(idx1).getEnd()));
                    idx2++;
                }
            }
        }
        if(idx1 < thisExons.size())
            for(int i = idx1; i < thisExons.size(); i++)
                answer.add(thisExons.get(i));
        
        return answer;
    }
    
    public double getQuantileAtPosition(String chrom, String strand, int pos){
        if(this.chrom.equals(chrom) && this.strand.equals(strand) && this.start <= pos && this.end >= pos){
            double quantile = 0;
            for (Exon exon : exons) {
                if (pos < exon.getStart()) {
                    return -1;
                } else if (pos <= exon.getEnd()) {
                    quantile += pos - exon.getStart() + 1;
                    break;
                } else {
                    quantile += exon.getLength();
                }
            }
            quantile /= this.getTotalExonLength();
            if(strand == "-")
                quantile = 1-quantile;
            return quantile;
        }
        return -1;
    }
}
