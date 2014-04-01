/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hitseq;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;

/**
 * The class for SAM/BAM processing, including information collecting and uniquely mapped reads extraction
 * @author hezhisong
 */
public class MappingProcessor {
    private File inputFile;
    private ArrayList<Integer> info;
    
    MappingProcessor(File file){
        if(! file.exists()){
            System.err.println("Cannot find input file: "+file.getAbsolutePath());
            System.exit(1);
        }
        this.inputFile=file;
        info=new ArrayList<>();
    }
    
    ArrayList<Integer> getMappingInfo(){
        if(info.isEmpty())
            return(null);
        else
            return(info);
    }
    
    void collectMappingInformation(){
        HashMap<String, Integer> numTimesOfReads=new HashMap<>();
        java.util.HashSet<String> unmappedReads=new java.util.HashSet<>();
        HashMap<Integer, ArrayList<String>> readsWithMismatches=new HashMap<>();
            
        try (SAMFileReader inputSam = new SAMFileReader(inputFile)) {
            for(SAMRecord record : inputSam){
                String readName=record.getReadName();
                
                if(numTimesOfReads.containsKey(readName))
                    numTimesOfReads.put(readName, numTimesOfReads.get(readName)+1);
                else
                    numTimesOfReads.put(readName, 1);
                
                if(record.getReadUnmappedFlag() || (record.getReadPairedFlag() && ! record.getProperPairFlag()))
                    unmappedReads.add(readName);
                else{
                    if(record.getIntegerAttribute("NM")!=null){
                        int numMismatches=record.getIntegerAttribute("NM");
                        if(!readsWithMismatches.containsKey(numMismatches))
                            readsWithMismatches.put(numMismatches, new ArrayList<String>());
                        readsWithMismatches.get(numMismatches).add(readName);
                    } else if(record.getIntegerAttribute("nM")!=null){
                        int numMismatches=record.getIntegerAttribute("nM");
                        if(!readsWithMismatches.containsKey(numMismatches))
                            readsWithMismatches.put(numMismatches, new ArrayList<String>());
                        readsWithMismatches.get(numMismatches).add(readName);
                    }
                }
                
                if(numTimesOfReads.size()%1000000==0)
                    System.err.println("finish reading "+numTimesOfReads.size()+" reads.");
            }
            inputSam.close();
        } catch(Exception e){
            System.err.println(e);
        }
        
        int numTotalReads=numTimesOfReads.keySet().size();
        int numUnmappedReads=unmappedReads.size();
        int numMappedReads=numTotalReads-numUnmappedReads;
        int numUniquelyMapped=0;
        for(String read : numTimesOfReads.keySet())
            if(numTimesOfReads.get(read).equals(1) && ! unmappedReads.contains(read))
                numUniquelyMapped++;
        info.add(numTotalReads);
        info.add(numMappedReads);
        info.add(numUniquelyMapped);
        
        java.util.TreeSet<Integer> sortedMismatchesNum=new java.util.TreeSet<>(new java.util.Comparator<Integer>() {
            @Override
            public int compare(Integer o1, Integer o2) {
                return o1.compareTo(o2);
            }
        });
        for(Integer numMistaches : readsWithMismatches.keySet())
            sortedMismatchesNum.add(numMistaches);
        for (Iterator<Integer> it = sortedMismatchesNum.iterator(); it.hasNext();) {
            Integer numMismatches=it.next();
            int numReadsWithThisNum=readsWithMismatches.get(numMismatches).size();
            int numReadsWithThisNumUniq=0;
            for(String read : readsWithMismatches.get(numMismatches))
                if(numTimesOfReads.get(read).equals(1))
                    numReadsWithThisNumUniq++;
            info.add(numReadsWithThisNum);
            info.add(numReadsWithThisNumUniq);
        }
    }
    
    int extractUniquelyMappedReads(File output){
        int numUnique=0;
        
        SAMFileReader inputSam=new SAMFileReader(inputFile);
        SAMFileWriter outputSam=new SAMFileWriterFactory().makeSAMOrBAMWriter(inputSam.getFileHeader(), true, output);
        boolean withNHTag=false;
        HashMap<String, Integer> numHits=new HashMap<>();
        
        for(SAMRecord record : inputSam){
            if(record.getIntegerAttribute("NH")!=null){ // there is NH tag, which shows the number of hits of this read
                withNHTag=true;
                if(record.getIntegerAttribute("NH").equals(1)){
                    outputSam.addAlignment(record);
                    numUnique++;
                }
            }
            else{
                if(numHits.containsKey(record.getReadName()))
                    numHits.put(record.getReadName(), numHits.get(record.getReadName())+1);
                else
                    numHits.put(record.getReadName(), 1);
            }
        }
        inputSam.close();
        
        if(! withNHTag){
            inputSam=new SAMFileReader(inputFile);
            for(SAMRecord record : inputSam)
                if(numHits.get(record.getReadName()) == 1){
                    record.setAttribute("NH", 1);
                    outputSam.addAlignment(record);
                    numUnique++;
                }
        }
        
        inputSam.close();
        outputSam.close();
        return(numUnique);
    }
}
