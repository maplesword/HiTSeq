/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hitseq.annotation;

/**
 *
 * @author hezhisong
 */
public class Exon{
    private String id;
    private String chrom;
    private String strand;
    private int start;
    private int end;
    
    public Exon(String id, String chrom, String strand, int start, int end){
        if(start>end){
            System.err.println("Error when creating Exon: the start coordinate should not be larger than the end coordinate.");
            System.exit(1);
        }
        
        this.id=id;
        this.chrom=chrom;
        this.strand=strand;
        this.start=start;
        this.end=end;
    }
    
    public Exon(String chrom, String strand, int start, int end){
        this("undefined", chrom, strand, start, end);
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
    
    public static boolean atSameChrom(Exon exon1, Exon exon2){
        return(exon1.chrom.equals(exon2.chrom));
    }
    
    public static boolean atSameStrand(Exon exon1, Exon exon2){
        return(exon1.strand.equals(exon2.strand));
    }
    
    public int startCompareTo(Exon exon2){
        int compare;
        if(this.start>exon2.start)
            compare=1;
        else if(this.start<exon2.start)
            compare=-1;
        else
            compare=0;
        
        return(compare);
    }
    
    public int endCompareTo(Exon exon2){
        int compare;
        if(this.end>exon2.end)
            compare=1;
        else if(this.end<exon2.end)
            compare=-1;
        else
            compare=0;
        
        return(compare);
    }
    
    public int compareTo(Exon exon2){
        int compare=startCompareTo(exon2);
        if(compare==0)
            compare=endCompareTo(exon2);
        return(compare);
    }
    
    public int getLength(){
        return(end-start+1);
    }
}
