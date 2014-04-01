/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hitseq;

import java.io.File;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author hezhisong
 */
public class Annotation {
    private HashMap<String, ArrayList<String>> genesInChrom;
    private HashMap<String, String> strandOfGene;
    private HashMap<String, Integer> startOfGene;
    private HashMap<String, Integer> endOfGene;
    private HashMap<String, Integer> lengthOfGene;
    private HashMap<String, Integer> exclusiveLengthOfGene;
    private HashMap<String, Integer> exclusiveLengthOfGeneNoStrand;
    private HashMap<String, ArrayList<Integer>> exonStartOfGene;
    private HashMap<String, ArrayList<Integer>> exonEndOfGene;
    
    private HashMap<String, Integer> pointerOfChrom;
    private HashMap<String, Integer> pointerOfGene;
    
    public Annotation(){
        reset();
    }
    
    public Annotation(File file, String fileType){
        reset();
        addAdditionalAnnotations(file, fileType);
    }
    
    public Annotation(File file){
        this(file, "struc");
    }
    
    final void reset(){
        resetAnnotation();
        resetPointer();
    }
    
    final void resetAnnotation(){
        genesInChrom=new HashMap<>();
        strandOfGene=new HashMap<>();
        lengthOfGene=new HashMap<>();
        exclusiveLengthOfGene=new HashMap<>();
        exclusiveLengthOfGeneNoStrand=new HashMap<>();
        startOfGene=new HashMap<>();
        endOfGene=new HashMap<>();
        exonStartOfGene=new HashMap<>();
        exonEndOfGene=new HashMap<>();
    }
    
    final void resetPointer(){
        pointerOfChrom=new HashMap<>();
        pointerOfGene=new HashMap<>();
        if(! genesInChrom.isEmpty())
            for(String chrom : genesInChrom.keySet())
                pointerOfChrom.put(chrom, 0);
        if(! startOfGene.isEmpty())
            for(String gene : startOfGene.keySet())
                pointerOfGene.put(gene, 0);
    }
    
    final void addAdditionalAnnotations(File file, String fileType){
        boolean emptyAtFirst=genesInChrom.isEmpty() ? true : false; // whether the annotation set is empty at the beginning
        
        try{
            RandomAccessFile fileIn=new RandomAccessFile(file,"r");
            
            // Read the annotation file
            int numLines=0;
            String line;
            while((line=fileIn.readLine()) != null){ // deal with each line separately
                if(line.startsWith("#"))
                    continue;
                String[] elements=line.split("\t");
                switch (fileType) {
                    case "struc":
                        // File format: gene structure (my customized file type)
                        // Fields: geneID, chrom, strand, exonLength, exonStructure
                        
                        if(! genesInChrom.containsKey(elements[1])){
                            genesInChrom.put(elements[1], new ArrayList<String>());
                            pointerOfChrom.put(elements[1], 0);
                        }
                        genesInChrom.get(elements[1]).add(elements[0]);
                        
                        Integer geneStart=-1;
                        Integer geneEnd=-1;
                        Pattern geneStartPattern=Pattern.compile("^\\d+");
                        Pattern geneEndPattern=Pattern.compile("\\d+$");
                        Matcher geneStartMatcher=geneStartPattern.matcher(elements[4]);
                        Matcher geneEndMatcher=geneEndPattern.matcher(elements[4]);
                        if(geneStartMatcher.find() && geneEndMatcher.find()){
                            geneStart=Integer.valueOf(geneStartMatcher.group());
                            geneEnd=Integer.valueOf(geneEndMatcher.group());
                        }
                        startOfGene.put(elements[0], geneStart);
                        endOfGene.put(elements[0], geneEnd);
                        pointerOfGene.put(elements[0], 0);
                        switch (elements[2]) {
                            case "1":
                                elements[2]="+";
                                break;
                            case "-1":
                                elements[2]="-";
                                break;
                        }
                        strandOfGene.put(elements[0], elements[2]);
                        lengthOfGene.put(elements[0], Integer.valueOf(elements[3]));
                        
                        String[] exons=elements[4].split(",");
                        exonStartOfGene.put(elements[0], new ArrayList<Integer>());
                        exonEndOfGene.put(elements[0], new ArrayList<Integer>());
                        for(String exon : exons){
                            String[] coord=exon.split("\\.\\.");
                            exonStartOfGene.get(elements[0]).add(Integer.valueOf(coord[0]));
                            exonEndOfGene.get(elements[0]).add(Integer.valueOf(coord[1]));
                        }
                        break;
                    case "gtf":
                        // File format: GTF file from Ensembl
                        if(elements[2].equals("exon")){
                            if(! genesInChrom.containsKey(elements[0])){
                                genesInChrom.put(elements[0], new ArrayList<String>());
                                pointerOfChrom.put(elements[0], 0);
                            }
                            Integer exonStart=Integer.valueOf(elements[3]);
                            Integer exonEnd=Integer.valueOf(elements[4]);
                            
                            Pattern geneIdPattern=Pattern.compile("gene_id \"[\\w\\.]+\";");
                            Matcher geneIdMatcher=geneIdPattern.matcher(elements[8]);
                            if(geneIdMatcher.find()){
                                String geneId=Pattern.compile("gene_id \"").matcher(geneIdMatcher.group()).replaceAll("");
                                geneId=Pattern.compile("\";").matcher(geneId).replaceAll("");
                                
                                if(!genesInChrom.get(elements[0]).contains(geneId))
                                    genesInChrom.get(elements[0]).add(geneId);
                                strandOfGene.put(geneId, elements[6]);
                                
                                if(! exonStartOfGene.containsKey(geneId)){
                                    exonStartOfGene.put(geneId, new ArrayList<Integer>());
                                    exonEndOfGene.put(geneId, new ArrayList<Integer>());
                                    startOfGene.put(geneId, exonStart);
                                    endOfGene.put(geneId, exonEnd);
                                    pointerOfGene.put(geneId, 0);
                                }
                                exonStartOfGene.get(geneId).add(exonStart);
                                exonEndOfGene.get(geneId).add(exonEnd);
                                if(exonStart < startOfGene.get(geneId))
                                    startOfGene.put(geneId, exonStart);
                                if(exonEnd > endOfGene.get(geneId))
                                    endOfGene.put(geneId, exonEnd);
                            }
                        }
                        break;
                    case "bed":
                        // File format: BED4/BED8
                        // Fields: chrom, chromStart, chromEnd, name, score, strand; the other fields will not be considered.
                        if(elements.length<3)
                            break;
                        
                        if(! genesInChrom.containsKey(elements[0])){
                            genesInChrom.put(elements[0], new ArrayList<String>());
                            pointerOfChrom.put(elements[0], 0);
                        }
                        Integer exonStart=Integer.valueOf(elements[1])+1;
                        Integer exonEnd=Integer.valueOf(elements[2]);
                        String name=elements.length>3 ? elements[3] : "Interval."+numLines;
                        String strand=elements.length>5 ? elements[5] : "*";
                        if(strandOfGene.containsKey(name)){
                            int idx=1;
                            while(strandOfGene.containsKey(name+String.valueOf(idx)))
                                idx++;
                            name=name+String.valueOf(idx);
                        }
                        
                        genesInChrom.get(elements[0]).add(name);
                        strandOfGene.put(name, strand);
                        lengthOfGene.put(name, Integer.valueOf(exonEnd-exonStart+1));
                        startOfGene.put(name, exonStart);
                        endOfGene.put(name, exonEnd);
                        exonStartOfGene.put(name, new ArrayList<Integer>());
                        exonEndOfGene.put(name, new ArrayList<Integer>());
                        exonStartOfGene.get(name).add(exonStart);
                        exonEndOfGene.get(name).add(exonEnd);
                        break;
                }
                numLines++;
                if(numLines%100000==0)
                    System.err.println("read annotation (in "+fileType+") "+String.valueOf(numLines)+" lines");
            }
            
            System.err.println("finish reading annotation file (in "+fileType+" format)");
            
            // sort the gene list for each chromosome
            for(String chrom : genesInChrom.keySet()){
                if(genesInChrom.get(chrom).size()==1) // if there is only one gene in the set, sorting is not necessary
                    continue;
                
                for(int i=0; i<genesInChrom.get(chrom).size()-1; i++){ // i: last index of ordered area
                    int j=i+1; // j: first index of unordered area
                    String geneFirstUnordered=genesInChrom.get(chrom).get(j);
                    
                    while(i>=0 && (startOfGene.get(genesInChrom.get(chrom).get(i)) > startOfGene.get(geneFirstUnordered) || (startOfGene.get(genesInChrom.get(chrom).get(i)).intValue() == startOfGene.get(geneFirstUnordered).intValue() && endOfGene.get(genesInChrom.get(chrom).get(i)) > endOfGene.get(geneFirstUnordered)))){
                        i--;
                    }
                    genesInChrom.get(chrom).add(i+1, geneFirstUnordered);
                    genesInChrom.get(chrom).remove(j+1);
                    i=j-1;
                }
            }
            
            if(fileType.equals("gtf") || ! emptyAtFirst){
                // If the annotation is in GTF format, or the operation is to add additional annotation, the overlapping exons of the same gene should be merged
                
                System.err.println("start sorting gene structure");
                System.err.println("in total "+exonStartOfGene.keySet().size()+" genes");
                for(String gene : exonStartOfGene.keySet()){
                    if(exonStartOfGene.get(gene).size()==1){ // if there is only one exon of this gene, only the exonic length needs to be calculated
                        int length=exonEndOfGene.get(gene).get(0) - exonStartOfGene.get(gene).get(0) + 1;
                        lengthOfGene.put(gene, length);
                        continue;
                    }
                    // first, sort the exon start and end coordinate hash table
                    for(int i=0; i<exonStartOfGene.get(gene).size()-1; i++){
                        int j=i+1;
                        int exonStartFirstUnordered=exonStartOfGene.get(gene).get(j);
                        int exonEndFirstUnordered=exonEndOfGene.get(gene).get(j);
                        while(i>=0 && (exonStartOfGene.get(gene).get(i) > exonStartFirstUnordered || (exonStartOfGene.get(gene).get(i) == exonStartFirstUnordered && exonEndOfGene.get(gene).get(i) > exonEndFirstUnordered))){
                            i--;
                        }
                        exonStartOfGene.get(gene).add(i+1, exonStartFirstUnordered);
                        exonEndOfGene.get(gene).add(i+1, exonEndFirstUnordered);
                        exonStartOfGene.get(gene).remove(j+1);
                        exonEndOfGene.get(gene).remove(j+1);
                        i=j-1;
                    }
                    
                    // second, merge the neighboring exons if they are overlapping with each other
                    for(int i=0; i<exonStartOfGene.get(gene).size()-1; i++){
                        int exonEndThis=exonEndOfGene.get(gene).get(i);
                        int exonStartNext=exonStartOfGene.get(gene).get(i+1);
                        int exonEndNext=exonEndOfGene.get(gene).get(i+1);
                        
                        if(exonStartNext <= exonEndThis){ // these two exons are overlapping
                            exonStartOfGene.get(gene).remove(i+1);
                            if(exonEndThis > exonEndNext)
                                exonEndOfGene.get(gene).remove(i+1);
                            else
                                exonEndOfGene.get(gene).remove(i);
                            i--;
                        }
                        else if(exonStartNext == exonEndThis+1){
                            exonStartOfGene.get(gene).remove(i+1);
                            exonEndOfGene.get(gene).remove(i);
                            i--;
                        }
                    }
                    
                    // third, calculate exon length for each gene
                    int length=0;
                    for(int i=0; i<exonStartOfGene.get(gene).size(); i++)
                        length+=exonEndOfGene.get(gene).get(i) - exonStartOfGene.get(gene).get(i) +1;
                    lengthOfGene.put(gene, length);
                }
            }
        }
        catch(java.io.IOException e){
            System.err.println("IO error: cannot open file"+file.getPath());
        }
    }
    
    void outputInStruc(){
        // Fields: geneID, chrom, strand, exonLength, exonStructure
        for(String chr : genesInChrom.keySet()){
            for(String gene : genesInChrom.get(chr)){
                String strand=strandOfGene.get(gene);
                strand=strand.equals("+") ? "1" : "-1";
                Integer length=lengthOfGene.get(gene);
                //System.err.println(gene+"\t"+length);
                
                String struc="";
                for(int i=0; i<exonStartOfGene.get(gene).size(); i++){
                    if(! struc.isEmpty())
                        struc=struc+",";
                    String start=exonStartOfGene.get(gene).get(i).toString();
                    String end=exonEndOfGene.get(gene).get(i).toString();
                    struc=struc+start+".."+end;
                }
                
                System.out.println(gene+"\t"+chr+"\t"+strand+"\t"+length.toString()+"\t"+struc);
            }
        }
    }
    
    void outputInStruc(File file){
        try{
            RandomAccessFile fileOut=new RandomAccessFile(file,"w");
            for(String chr : genesInChrom.keySet()){
                for(String gene : genesInChrom.get(chr)){
                    String strand=strandOfGene.get(gene);
                    strand=strand.equals("+") ? "1" : "-1";
                    Integer length=lengthOfGene.get(gene);
                    
                    String struc="";
                    for(int i=0; i<exonStartOfGene.get(gene).size(); i++){
                        if(! struc.isEmpty())
                            struc=struc+",";
                        String start=exonStartOfGene.get(gene).get(i).toString();
                        String end=exonEndOfGene.get(gene).get(i).toString();
                        struc=struc+start+".."+end;
                    }
                    
                    fileOut.writeUTF(gene+"\t"+chr+"\t"+strand+"\t"+length.toString()+"\t"+struc+System.getProperty("line.separator"));
                }
            }
        } catch(java.io.IOException e){
            System.err.println("IO error: cannot open file"+file.getPath());
        }
    }
    
    boolean chromIsExisted(String chrom){
        return genesInChrom.containsKey(chrom);
    }
    
    String[] getAvailableChromosomes(){
        return (String[])genesInChrom.keySet().toArray();
    }
    
    String[] getGeneSet(){
        Object[] results=startOfGene.keySet().toArray();
        String[] answer=new String[results.length];
        for(int i=0; i<results.length; i++)
            answer[i]=(String) results[i];
        return answer;
    }
    
    String[] getGeneSet(String chrom){
        if(genesInChrom.containsKey(chrom))
            return (String[]) genesInChrom.get(chrom).toArray();
        else
            return null;
    }
    
    int getCurrentGenePointer(String chrom){
        if(genesInChrom.containsKey(chrom))
            return pointerOfChrom.get(chrom);
        else
            return -1;
    }
    
    String getGene(String chrom, int pointer){
        if(genesInChrom.containsKey(chrom) && genesInChrom.get(chrom).size()>pointer)
            return genesInChrom.get(chrom).get(pointer);
        else
            return null;
    }
    
    String getGeneStrand(String gene){
        if(strandOfGene.containsKey(gene))
            return strandOfGene.get(gene);
        else
            return null;
    }
    
    int getGeneLength(String gene){
        if(lengthOfGene.containsKey(gene))
            return lengthOfGene.get(gene);
        else
            return -1;
    }
    
    int getExclusiveGeneLength(String gene){
        if(exclusiveLengthOfGene.containsKey(gene))
            return exclusiveLengthOfGene.get(gene);
        else
            return -1;
    }
    
    int getExclusiveGeneLengthNoStrand(String gene){
        if(exclusiveLengthOfGeneNoStrand.containsKey(gene))
            return exclusiveLengthOfGeneNoStrand.get(gene);
        else
            return -1;
    }
    
    int getGeneStart(String gene){
        if(startOfGene.containsKey(gene))
            return startOfGene.get(gene);
        else
            return -1;
    }
    
    int getGeneEnd(String gene){
        if(endOfGene.containsKey(gene))
            return endOfGene.get(gene);
        else
            return -1;
    }
    
    int getCurrentExonIndex(String gene){
        if(pointerOfGene.containsKey(gene))
            return pointerOfGene.get(gene);
        else
            return -1;
    }
    
    int getExonStart(String gene, int index){
        if(index < exonStartOfGene.get(gene).size() && index >=0)
            return exonStartOfGene.get(gene).get(index);
        else
            return -1;
    }
    
    int getExonEnd(String gene, int index){
        if(index < exonEndOfGene.get(gene).size() && index >=0)
            return exonEndOfGene.get(gene).get(index);
        else
            return -1;
    }
    
    int movePointerGene(String chrom){
        if(pointerOfChrom.get(chrom)>=genesInChrom.get(chrom).size())
            pointerOfChrom.put(chrom, -1);
        else
            pointerOfChrom.put(chrom, pointerOfChrom.get(chrom)+1);
        return pointerOfChrom.get(chrom);
    }
    
    int movePointerExon(String gene){
        pointerOfGene.put(gene, pointerOfGene.get(gene)+1);
        return pointerOfGene.get(gene);
    }
    
    void estimateExclusiveGeneLength(boolean outputPairs){
        for(String gene : lengthOfGene.keySet()){
            exclusiveLengthOfGene.put(gene, lengthOfGene.get(gene));
            exclusiveLengthOfGeneNoStrand.put(gene, lengthOfGene.get(gene));
        }
        
        for(String chrom : genesInChrom.keySet()){
            ArrayList<String> genesInThisChrom=genesInChrom.get(chrom);
            for(int i=0; i<genesInThisChrom.size()-1; i++){
                String thisGene=genesInThisChrom.get(i);
                int geneStartThis=startOfGene.get(thisGene);
                int geneEndThis=endOfGene.get(thisGene);

                ArrayList<String> overlapGenes=new ArrayList<>();
                for(int j=i+1; j<genesInThisChrom.size(); j++){
                    int geneStartThat=startOfGene.get(genesInThisChrom.get(j));
                    if(geneStartThis <= geneStartThat && geneStartThat <= geneEndThis)
                        overlapGenes.add(genesInThisChrom.get(j));
                    else if(geneEndThis < geneStartThat)
                        break;
                }
                
                //if(overlapGenes.size()>0) System.out.println("overlap candidates for "+thisGene+": "+overlapGenes);
                for(String gene : overlapGenes){
                    boolean sameStrand=strandOfGene.get(thisGene).equals(strandOfGene.get(gene));
                    int overlapSize=0;
                    int overlapSizeNoStrand=0;
                    
                    ArrayList<Integer> exonStartThis=exonStartOfGene.get(thisGene);
                    ArrayList<Integer> exonEndThis=exonEndOfGene.get(thisGene);
                    ArrayList<Integer> exonStartThat=exonStartOfGene.get(gene);
                    ArrayList<Integer> exonEndThat=exonEndOfGene.get(gene);
                    
                    int idx1=0, idx2=0; // pointers
                    while(idx1<exonStartThis.size() && idx2<exonStartThat.size()){
                        if(exonStartThis.get(idx1).intValue()<=exonStartThat.get(idx2).intValue() && exonStartThat.get(idx2).intValue()<=exonEndThis.get(idx1).intValue()){
                            if(exonEndThis.get(idx1).intValue()<exonEndThat.get(idx2).intValue()){
                                overlapSizeNoStrand+=exonEndThis.get(idx1).intValue()-exonStartThat.get(idx2).intValue()+1;
                                if(sameStrand) overlapSize+=exonEndThis.get(idx1).intValue()-exonStartThat.get(idx2).intValue()+1;
                                idx1++;
                            } else if(exonEndThis.get(idx1).intValue()==exonEndThat.get(idx2).intValue()){
                                overlapSizeNoStrand+=exonEndThis.get(idx1).intValue()-exonStartThat.get(idx2).intValue()+1;
                                if(sameStrand) overlapSize+=exonEndThis.get(idx1).intValue()-exonStartThat.get(idx2).intValue()+1;
                                idx1++; idx2++;
                            } else{
                                overlapSizeNoStrand+=exonEndThat.get(idx2).intValue()-exonStartThat.get(idx2).intValue()+1;
                                if(sameStrand) overlapSize+=exonEndThat.get(idx2).intValue()-exonStartThat.get(idx2).intValue()+1;
                                idx2++;
                            }
                        } else if(exonStartThat.get(idx2).intValue()<=exonStartThis.get(idx1) && exonStartThis.get(idx1).intValue()<=exonEndThat.get(idx2).intValue()){
                            if(exonEndThat.get(idx2).intValue()<exonEndThis.get(idx1).intValue()){
                                overlapSizeNoStrand+=exonEndThat.get(idx2).intValue()-exonStartThis.get(idx1).intValue()+1;
                                if(sameStrand) overlapSize+=exonEndThat.get(idx2).intValue()-exonStartThis.get(idx1).intValue()+1;
                                idx2++;
                            } else if(exonEndThat.get(idx2).intValue()==exonEndThis.get(idx1).intValue()){
                                overlapSizeNoStrand+=exonEndThat.get(idx2).intValue()-exonStartThis.get(idx1).intValue()+1;
                                if(sameStrand) overlapSize+=exonEndThat.get(idx2).intValue()-exonStartThis.get(idx1).intValue()+1;
                                idx2++; idx1++;
                            } else{
                                overlapSizeNoStrand+=exonEndThis.get(idx1).intValue()-exonStartThis.get(idx1).intValue()+1;
                                if(sameStrand) overlapSize+=exonEndThis.get(idx1).intValue()-exonStartThis.get(idx1).intValue()+1;
                                idx1++;
                            }
                        } else{
                            if(exonStartThis.get(idx1).intValue()<exonStartThat.get(idx2).intValue())
                                idx1++;
                            else
                                idx2++;
                        }
                        //System.out.println(thisGene+"\t"+gene+"\t"+idx1+"\t"+idx2);
                    }
                    
                    exclusiveLengthOfGene.put(thisGene, exclusiveLengthOfGene.get(thisGene).intValue()-overlapSize);
                    exclusiveLengthOfGeneNoStrand.put(thisGene, exclusiveLengthOfGeneNoStrand.get(thisGene).intValue()-overlapSizeNoStrand);
                    exclusiveLengthOfGene.put(gene, exclusiveLengthOfGene.get(gene).intValue()-overlapSize);
                    exclusiveLengthOfGeneNoStrand.put(gene, exclusiveLengthOfGeneNoStrand.get(gene).intValue()-overlapSizeNoStrand);

                    if(outputPairs && overlapSizeNoStrand>0)
                        System.out.println(thisGene+"\t"+gene+"\t"+sameStrand+"\t"+overlapSizeNoStrand+"\t"+overlapSize);
                }
                //System.out.println("done "+i+" genes in "+chrom);
            }
        }
    }
    
    void estimateExclusiveGeneLength(){
        estimateExclusiveGeneLength(false);
    }
}
