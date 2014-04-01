/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hitseq;

/**
 *
 * @author hezhisong
 */
public class Junction {
    String id;
    String chrom;
    String strand;
    int donorSite;
    int acceptorSite;
    
    Junction(String id, String chrom, String strand, int donorSite, int acceptorSite){
        this.id=id;
        this.chrom=chrom;
        this.strand=strand;
        this.donorSite=donorSite;
        this.acceptorSite=acceptorSite;
    }
    
    Junction(String chrom, String strand, int donorSite, int acceptorSite){
        this.id="undefined";
        this.chrom=chrom;
        this.strand=strand;
        this.donorSite=donorSite;
        this.acceptorSite=acceptorSite;
    }
    
    String getChrom(){
        return(chrom);
    }
    
    String getStrand(){
        return(strand);
    }
    
    int getDonorSite(){
        return(donorSite);
    }
    
    int getAcceptorSite(){
        return(acceptorSite);
    }
}
