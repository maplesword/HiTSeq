/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hitseq.annotation;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
        
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
    private HashSet<Junction> junctions;
    
    private Transcript ambiguousRegions;
    private Transcript ambiguousRegionsIgnoreStrand;
    
    /**
     * Generate a new Gene object with given gene ID, chromosome, strand, gene start and gene end coordinates.
     * @param id Gene ID.
     * @param chrom Chromosome.
     * @param strand Strand.
     * @param start Gene start coordinates. Count from 1.
     * @param end Gene end coordinates. Count from 1.
     */
    public Gene(String id, String chrom, String strand, int start, int end){
        this.id=id;
        this.chrom=chrom;
        this.strand=strand;
        this.start=start;
        this.end=end;
        this.transcripts=new HashMap<>();
        this.nonredundantTranscript=new Transcript("Structure", chrom, strand);
        this.segments=new ArrayList<>();
        this.junctions=new HashSet<>();
        this.ambiguousRegions=new Transcript("Ambiguous", chrom, strand);
        this.ambiguousRegionsIgnoreStrand=new Transcript("Ambiguous", chrom, strand);
    }
    
    /**
     * Generate a new Gene object without gene start and end coordinates.
     * @param id Gene ID.
     * @param chrom Chromosome.
     * @param strand Strand.
     */
    public Gene(String id, String chrom, String strand){
        this(id, chrom, strand, -1, -1);
    }
    
    /**
     * FunName: addTranscript.
     * Description: Add a new transcript to the gene.
     * @param newTranscript The new transcript.
     */
    public void addTranscript(Transcript newTranscript){
        String transcriptID=newTranscript.getID();
        transcripts.put(transcriptID, newTranscript);
        if(start==-1 || newTranscript.getStart()<start)
            start=newTranscript.getStart();
    }
    
    /**
     * FunName: generateNonredundantTranscript.
     * Description: Generate the non-redundant transcript which covers all the possible exonic regions of this gene
     * according its transcript. The generated transcript will be represented as nonredundantTranscript.
     */
    public void generateNonredundantTranscript(){
        nonredundantTranscript=new Transcript("Structure", chrom, strand);
        for(String transcriptID : transcripts.keySet())
            for(int i=0; i<transcripts.get(transcriptID).getExonNumber(); i++)
                nonredundantTranscript.addExon(transcripts.get(transcriptID).getExon(i));
        nonredundantTranscript.mergeExons();
    }
    
    /**
     * FunName: generateAllJunctions.
     * Description: According to the transcripts of the gene, identify all the annotated exon-exon junctions.
     */
    public void generateAllJunctions(){
        for(Transcript transcript : transcripts.values())
            junctions.addAll(transcript.getJunctions());
        for(Junction junc : junctions)
            junc.addAnnotatedGene(this);
    }
    
    public void generateSegments(){
        
    }
    
    /**
     * FunName: containTranscript.
     * Description: Judge whether this gene contains a transcript with the given transcript ID.
     * @param transcriptID The transcript ID of the transcript to be determined.
     * @return true if the transcript exists.
     */
    public boolean containTranscript(String transcriptID){
        return(transcripts.containsKey(transcriptID));
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
    
    /**
     * FunName: getTranscript.
     * Description: Given the transcript ID, return the corresponding Transcript object if the transcript exists.
     * @param transcriptID The transcript ID.
     * @return The corresponding transcript object.
     */
    public Transcript getTranscript(String transcriptID){
        if(containTranscript(transcriptID))
            return(transcripts.get(transcriptID));
        else
            return(null);
    }
    
    public Transcript getNonredundantTranscript(){
        return(nonredundantTranscript);
    }
    
    public Transcript getAmbiguousRegions(){
        return(ambiguousRegions);
    }
    
    public Transcript getAmbiguousRegionsIgnoringStrand(){
        return(ambiguousRegionsIgnoreStrand);
    }
    
    public HashSet<Junction> getAllJunctions(){
        return((HashSet<Junction>)junctions.clone());
    }
    
    public boolean refleshJunction(Junction newJunc){
        if(junctions.contains(newJunc)){
            junctions.remove(newJunc);
            junctions.add(newJunc);
            return(true);
        }
        else
            return(false);
    }
    
    public ArrayList<Exon> getSegments(){
        return(segments);
    }
    
    public int getTotalExonLength(){
        if(nonredundantTranscript.getStart()==-1)
            generateNonredundantTranscript();
        return(nonredundantTranscript.getTotalExonLength());
    }
    
    public ArrayList<Exon> getOverlappingRegion(Gene gene, boolean considerStrand){
        ArrayList<Exon> answer=nonredundantTranscript.getOverlappingRegions(gene.getNonredundantTranscript(), considerStrand);
        return(answer);
    }
    
    @Override
    public String toString(){
        return(id);
    }
}
