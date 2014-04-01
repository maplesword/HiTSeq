/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hitseq;

import java.util.ArrayList;
import java.util.HashMap;
        
/**
 *
 * @author hezhisong
 */
public class Gene {
    private String id;
    private String chrom;
    private String strand;
    private int start;
    private int end;
    private HashMap<String, Transcript> transcripts;
    private Transcript nonredundantTranscript;
    private ArrayList<Exon> segments;
    private ArrayList<Junction> junctions;
    
    
    Gene(String id, String chrom, String strand, int start, int end){
        this.id=id;
        this.chrom=chrom;
        this.strand=strand;
        this.start=start;
        this.end=end;
        this.transcripts=new HashMap<>();
        this.nonredundantTranscript=new Transcript("Structure", chrom, strand);
        this.segments=new ArrayList<>();
        this.junctions=new ArrayList<>();
    }
    
    Gene(String id, String chrom, String strand){
        this(id, chrom, strand, -1, -1);
    }
    
    void addTranscript(Transcript newTranscript){
        String transcriptID=newTranscript.getID();
        transcripts.put(transcriptID, newTranscript);
        if(start==-1 || newTranscript.getStart()<start)
            start=newTranscript.getStart();
    }
    
    void generateNonredundantTranscript(){
        nonredundantTranscript=new Transcript("Structure", chrom, strand);
        for(String transcriptID : transcripts.keySet())
            nonredundantTranscript.addExons(transcripts.get(transcriptID).getExons());
        nonredundantTranscript.mergeExons();
    }
    
    void generateAllJunctions(){
        for(Transcript transcript : transcripts.values())
            junctions.addAll(transcript.getJunctions());
    }
    
    void generateSegments(){
        
    }
    
    boolean containTranscript(String transcriptID){
        return(transcripts.containsKey(transcriptID));
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
    
    Transcript getTranscript(String transcriptID){
        if(containTranscript(transcriptID))
            return(transcripts.get(transcriptID));
        else
            return(null);
    }
    
    Transcript getNonredundantTranscript(){
        return(nonredundantTranscript);
    }
    
    ArrayList<Junction> getAllJunctions(){
        return(junctions);
    }
    
    ArrayList<Exon> getSegments(){
        return(segments);
    }
    
    int getTotalExonLength(){
        if(nonredundantTranscript.getStart()==-1)
            generateNonredundantTranscript();
        return(nonredundantTranscript.getTotalExonLength());
    }
    
    ArrayList<Exon> getOverlappingRegion(Gene gene, boolean considerStrand){
        ArrayList<Exon> answer=nonredundantTranscript.getOverlappingRegions(gene.getNonredundantTranscript(), considerStrand);
        return(answer);
    }
}
