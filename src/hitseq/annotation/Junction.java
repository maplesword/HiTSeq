/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hitseq.annotation;

import java.util.HashSet;
import java.util.Objects;

/**
 *
 * @author hezhisong
 */
public class Junction {
    private String id;
    private String chrom;
    private String strand;
    private int startSite;
    private int endSite;
    private HashSet<Gene> annotatedGenes;
    private HashSet<Transcript> annotatedTranscripts;
    
    public Junction(String id, String chrom, String strand, int startSite, int endSite, HashSet<Gene> annotatedGenes, HashSet<Transcript> annotatedTranscripts){
        this.id=id;
        this.chrom=chrom;
        this.strand=strand;
        this.startSite=startSite;
        this.endSite=endSite;
        this.annotatedGenes=annotatedGenes;
        this.annotatedTranscripts=annotatedTranscripts;
    }
    
    public Junction(String id, String chrom, String strand, int startSite, int endSite){
        this(id,chrom,strand,startSite, endSite, new HashSet<Gene>(), new HashSet<Transcript>());
    }
    
    public Junction(String chrom, String strand, int startSite, int endSite){
        this("undefined",chrom,strand,startSite,endSite,new HashSet<Gene>(), new HashSet<Transcript>());
    }
    
    public void addAnnotatedGene(Gene gene){
        annotatedGenes.add(gene);
    }
    
    public void addAnnotatedGeneSet(java.util.Collection<Gene> newAnnotatedGenes){
        annotatedGenes.addAll(newAnnotatedGenes);
    }
    
    public void addAnnotatedTranscript(Transcript transcript){
        annotatedTranscripts.add(transcript);
    }
    
    public void addAnnotatedTranscriptSet(java.util.Collection<Transcript> newAnnotatedTranscripts){
        annotatedTranscripts.addAll(newAnnotatedTranscripts);
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
    
    public int getStartSite(){
        return(startSite);
    }
    
    public int getEndSite(){
        return(endSite);
    }
    
    public int getSize(){
        return(endSite-startSite+1);
    }
    
    public HashSet<Gene> getAnnotatedGenes(){
        return((HashSet<Gene>) annotatedGenes.clone());
    }
    
    public HashSet<Transcript> getAnnotatedTranscripts(){
        return((HashSet<Transcript>) annotatedTranscripts.clone());
    }

    @Override
    public int hashCode() {
        int hash = 3;
        hash = 73 * hash + Objects.hashCode(this.chrom);
        hash = 73 * hash + Objects.hashCode(this.strand);
        hash = 73 * hash + this.startSite;
        hash = 73 * hash + this.endSite;
        return hash;
    }
    
    @Override
    public boolean equals(Object junc2){
        if(!(junc2 instanceof Junction))
            return false;
        Junction junc=(Junction) junc2;
        return(chrom.equals(junc.getChrom()) && strand.equals(junc.getStrand()) && startSite==junc.getStartSite() && endSite==junc.getEndSite());
    }
    
    public int compareTo(Junction junc2){
        if(!chrom.equals(junc2.chrom))
            return(chrom.compareTo(junc2.chrom));
        else if(startSite!=junc2.startSite)
            return((new Integer(startSite)).compareTo(new Integer(junc2.startSite)));
        else if(endSite!=junc2.endSite)
            return((new Integer(endSite)).compareTo(new Integer(junc2.endSite)));
        else
            return(strand.compareTo(junc2.strand));
    }
    
    public boolean touch(Junction junc2){
        if(!chrom.equals(junc2.chrom))
            return(false);
        else if(!strand.equals(junc2.strand) && !strand.equals("*") & !junc2.strand.equals("*"))
            return(false);
        else if(startSite!=junc2.startSite && endSite!=junc2.endSite)
            return(false);
        else
            return(true);
    }
    
    public boolean cross(Junction junc2){
        if(!chrom.equals(junc2.chrom))
            return(false);
        else if(!strand.equals(junc2.strand) && !strand.equals("*") & !junc2.strand.equals("*"))
            return(false);
        else if(startSite < junc2.startSite && junc2.startSite < endSite && endSite < junc2.endSite)
            return(true);
        else if(junc2.startSite < startSite && startSite < junc2.endSite && junc2.endSite < endSite)
            return(true);
        else
            return(false);
    }
    
    public boolean withStartSite(int donor){
        return(donor==startSite);
    }
    
    public boolean withEndSite(int acceptor){
        return(acceptor==endSite);
    }
    
    @Override
    public String toString(){
        String answer;
        if(id.equals("undefined"))
            answer=chrom+":"+strand+":"+startSite+"-"+endSite;
        else
            answer=id;
        return(answer);
    }
}
