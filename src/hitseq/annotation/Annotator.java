/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hitseq.annotation;
import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 *
 * @author Chih-sung
 */
public class Annotator {
    Annotation annotation;
    File file;
    ArrayList<AnnotatedIntervals> annotatedIntervals;
    
    public Annotator(Annotation annotation, File file){
        this.annotation = annotation;
        annotatedIntervals = new ArrayList<>();
        try{
            if(file.exists() && file.isFile()){
                this.file = file;
            } else{
                throw(new java.io.FileNotFoundException());
            }
        } catch(java.io.FileNotFoundException e){
            System.err.println("The given file: " + file.getPath() + " cannot be processed.");
            System.exit(1);
        }
    }
    
    public Annotator(Annotation annotation, String filePath){
        this(annotation, new File(filePath));
    }
    
    public void annotateTSS(){
        HashMap<String, HashMap<Integer, HashSet<String>>> structuredTSS = annotation.getStructuredTSS();
        HashMap<String, ArrayList<Integer>> sortedTSS = new HashMap<>();
        HashMap<String, Integer> pointer = new HashMap<>();
        for(String chrom : structuredTSS.keySet()){
            ArrayList<Integer> allTSSChrom = new ArrayList<>();
            allTSSChrom.addAll(structuredTSS.get(chrom).keySet());
            java.util.Collections.sort(allTSSChrom);
            
            sortedTSS.put(chrom, allTSSChrom);
            pointer.put(chrom, 0);
        }
        
        try{
            RandomAccessFile fileIn=new RandomAccessFile(file,"r");
            
            // Read BED file
            int numLines=0;
            String line;
            while((line=fileIn.readLine()) != null){ // deal with each line separately
                if(line.startsWith("#"))
                    continue;
                String[] elements=line.split("\t");
                
                if(elements.length<3)
                    break;
                String chrom = elements[0];
                Integer exonStart = Integer.valueOf(elements[1])+1;
                Integer exonEnd = Integer.valueOf(elements[2]);
                String name=elements.length>3 ? elements[3] : "Interval."+numLines;
                String strand=elements.length>5 ? elements[5] : "*";
                Exon interval = new Exon(name, chrom, strand, exonStart, exonEnd);
                
                
                
                numLines++;
                if(numLines%100000==0)
                    System.err.println("read annotation "+String.valueOf(numLines)+" lines");
            }
            
        } catch(Exception e){
            System.err.println("IO error: cannot open file"+file.getPath());
        }
    }
    
    public class AnnotatedIntervals {
        private final Exon interval;
        private final Gene gene;
        private final String type;
        
        public AnnotatedIntervals(Exon interval, Gene gene, String type){
            this.interval = interval;
            this.gene = gene;
            this.type = type;
        }
        
        public Exon getInterval(){
            return interval;
        }
        
        public Gene getGene(){
            return gene;
        }
        
        public String getType(){
            return type;
        }
    }
}
