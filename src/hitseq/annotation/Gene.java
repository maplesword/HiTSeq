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
    private String type;
    
    private HashMap<String, Transcript> transcripts;
    private Transcript nonredundantTranscript;
    private ArrayList<Exon> segments;
    private HashSet<Junction> junctions;
    private HashSet<Integer> transcriptTSS;
    
    private Transcript ambiguousRegions;
    private Transcript ambiguousRegionsIgnoreStrand;
    
    private Transcript exclusiveRegions;
    private Transcript exclusiveRegionsIgnoreStrand;
    
    /**
     * Generate a new Gene object with given gene ID, chromosome, strand, gene start, gene end coordinates and gene type.
     * @param id Gene ID.
     * @param chrom Chromosome.
     * @param strand Strand.
     * @param start Gene start coordinates. Count from 1.
     * @param end Gene end coordinates. Count from 1.
     * @param type Gene type.
     */
    public Gene(String id, String chrom, String strand, int start, int end, String type){
        this.id=id;
        this.chrom=chrom;
        this.strand=strand;
        this.start=start;
        this.end=end;
        this.transcripts=new HashMap<>();
        this.nonredundantTranscript=new Transcript("Structure", chrom, strand);
        this.segments=new ArrayList<>();
        this.junctions=new HashSet<>();
        this.transcriptTSS=new HashSet<>();
        this.type=type;
        this.ambiguousRegions=new Transcript("Ambiguous", chrom, strand);
        this.ambiguousRegionsIgnoreStrand=new Transcript("Ambiguous", chrom, strand);
        this.exclusiveRegions=new Transcript("Exclusive", chrom, strand);
        this.exclusiveRegionsIgnoreStrand=new Transcript("Exclusive", chrom, strand);
    }
    
    /**
     * Generate a new Gene object with only given gene ID, chromosome, strand, gene start and gene end coordinates.
     * @param id Gene ID.
     * @param chrom Chromosome.
     * @param strand Strand.
     * @param start Gene start coordinates. Count from 1.
     * @param end Gene end coordinates. Count from 1.
     */
    public Gene(String id, String chrom, String strand, int start, int end){
        this(id, chrom, strand, start, end, "others");
    }
    
    /**
     * Generate a new Gene object without gene start, end coordinates as well as gene type.
     * @param id Gene ID.
     * @param chrom Chromosome.
     * @param strand Strand.
     */
    public Gene(String id, String chrom, String strand){
        this(id, chrom, strand, -1, -1);
    }
    
    /**
     * Generate a new Gene object without gene start and end coordinates
     * @param id Gene ID.
     * @param chrom Chromosome.
     * @param strand Strand.
     * @param type Gene type.
     */
    public Gene(String id, String chrom, String strand, String type){
        this(id, chrom, strand, -1, -1, type);
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
        if(end==-1 || newTranscript.getEnd()>end)
            end = newTranscript.getEnd();
        
        if(strand.equals("+"))
            transcriptTSS.add(newTranscript.getStart());
        else
            transcriptTSS.add(newTranscript.getEnd());
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
        refreshStartAndEnd();
        if(strand.equals("+"))
            transcriptTSS.add(nonredundantTranscript.getStart());
        else
            transcriptTSS.add(nonredundantTranscript.getEnd());
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
    
    public String getType(){
        return(type);
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
    
    public Transcript getExclusiveRegions(){
        return(exclusiveRegions);
    }
    
    public Transcript getExclusiveRegionsIgnoringStrand(){
        return(exclusiveRegionsIgnoreStrand);
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
    
    public HashSet<Integer> getTSS(){
        return(transcriptTSS);
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
    
    private void refreshStartAndEnd(){
        start=nonredundantTranscript.getStart();
        end=nonredundantTranscript.getEnd();
    }
}
