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
    private int donorSite;
    private int acceptorSite;
    private HashSet<Gene> annotatedGenes;
    private HashSet<Transcript> annotatedTranscripts;
    
    public Junction(String id, String chrom, String strand, int donorSite, int acceptorSite, HashSet<Gene> annotatedGenes, HashSet<Transcript> annotatedTranscripts){
        this.id=id;
        this.chrom=chrom;
        this.strand=strand;
        this.donorSite=donorSite;
        this.acceptorSite=acceptorSite;
        this.annotatedGenes=annotatedGenes;
        this.annotatedTranscripts=annotatedTranscripts;
    }
    
    public Junction(String id, String chrom, String strand, int donorSite, int acceptorSite){
        this(id,chrom,strand,donorSite, acceptorSite, new HashSet<Gene>(), new HashSet<Transcript>());
    }
    
    public Junction(String chrom, String strand, int donorSite, int acceptorSite){
        this("undefined",chrom,strand,donorSite,acceptorSite,new HashSet<Gene>(), new HashSet<Transcript>());
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
    
    public int getDonorSite(){
        return(donorSite);
    }
    
    public int getAcceptorSite(){
        return(acceptorSite);
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
        hash = 73 * hash + this.donorSite;
        hash = 73 * hash + this.acceptorSite;
        return hash;
    }
    
    @Override
    public boolean equals(Object junc2){
        if(!(junc2 instanceof Junction))
            return false;
        Junction junc=(Junction) junc2;
        return(chrom.equals(junc.getChrom()) && strand.equals(junc.getStrand()) && donorSite==junc.getDonorSite() && acceptorSite==junc.getAcceptorSite());
    }
    
    public int compareTo(Junction junc2){
        if(!chrom.equals(junc2.chrom))
            return(chrom.compareTo(junc2.chrom));
        else if(donorSite!=junc2.donorSite)
            return((new Integer(donorSite)).compareTo(new Integer(junc2.donorSite)));
        else if(acceptorSite!=junc2.acceptorSite)
            return((new Integer(acceptorSite)).compareTo(new Integer(junc2.acceptorSite)));
        else
            return(strand.compareTo(junc2.strand));
    }
    
    public boolean touch(Junction junc2){
        if(!chrom.equals(junc2.chrom))
            return(false);
        else if(!strand.equals(junc2.strand) && !strand.equals("*") & !junc2.strand.equals("*"))
            return(false);
        else if(donorSite!=junc2.donorSite && acceptorSite!=junc2.acceptorSite)
            return(false);
        else
            return(true);
    }
    
    public boolean cross(Junction junc2){
        if(!chrom.equals(junc2.chrom))
            return(false);
        else if(!strand.equals(junc2.strand) && !strand.equals("*") & !junc2.strand.equals("*"))
            return(false);
        else if(donorSite < junc2.donorSite && junc2.donorSite < acceptorSite && acceptorSite < junc2.acceptorSite)
            return(true);
        else if(junc2.donorSite < donorSite && donorSite < junc2.acceptorSite && junc2.acceptorSite < acceptorSite)
            return(true);
        else
            return(false);
    }
    
    public boolean withDonorSite(int donor){
        return(donor==donorSite);
    }
    
    public boolean withAcceptorSite(int acceptor){
        return(acceptor==acceptorSite);
    }
}
