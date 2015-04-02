/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hitseq;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

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
            
        try (SamReader inputSam = SamReaderFactory.makeDefault().open(inputFile)) {
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
    
    int extractUniquelyMappedReads(File output, boolean attemptNH){
        // try the simple way to extract, but only work for single-ended data with NH tag
        if(attemptNH){
            System.err.println("Attempting to use the NH tag...");
            int numUnique=0;
            SamReader inputSam=SamReaderFactory.makeDefault().open(inputFile);
            try (SAMFileWriter outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(inputSam.getFileHeader(), true, output)) {
                for(SAMRecord record : inputSam){
                    if(record.getReadPairedFlag() || record.getIntegerAttribute("NH")==null){
                        System.err.println("NH attemption failed: The input SAM/BAM file is paired-ended or contains no NH tag.");
                        System.err.println("Continue the extraction using the typical way.");
                        CloserUtil.close(inputSam);
                        outputSam.close();
                        
                        numUnique=this.extractUniquelyMappedReads(output, false);
                        return(numUnique);
                    } else if(record.getIntegerAttribute("NH").equals(1)){
                        outputSam.addAlignment(record);
                        numUnique++;
                    }
                }
                
                CloserUtil.close(inputSam);
            }
            return(numUnique);
            
        // the time-consuming way, firstly sort it with read name, and then extract the needed reads, and re-sort them by coordinate again
        } else{
            int numUnique=0;
            
            // firstly, sort the reads by read IDs
            System.err.println("Start sorting by QueryName...");
            boolean tempGenerated=true;
            File tempFile=new File("tmp_sorting.bam");
            SamReader inputSam=SamReaderFactory.makeDefault().open(inputFile);
            if(inputSam.getFileHeader().getSortOrder().equals(SAMFileHeader.SortOrder.queryname)){
                tempFile=inputFile;
                tempGenerated=false;
            } else{ // only do the sorting when the input is not sorted by read name
                inputSam.getFileHeader().setSortOrder(SAMFileHeader.SortOrder.queryname);

                try (SAMFileWriter outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(inputSam.getFileHeader(), false, tempFile)) {
                    for (SAMRecord record : inputSam) {
                        outputSam.addAlignment(record);
                    }
                    CloserUtil.close(inputSam);
                    outputSam.close();
                } catch (Exception e) {
                    System.err.println("Error in the first sorting: extractUniquelyMappedReads of MappingProcessor.");
                    System.err.println(e);
                    return (-1);
                }
            }
            
            // second, extract the uniquely mapped reads or pairs of reads and then re-sort them into sorted-by-coordinate
            System.err.println("Start extraction...");
            inputSam=SamReaderFactory.makeDefault().open(tempFile);
            inputSam.getFileHeader().setSortOrder(SAMFileHeader.SortOrder.coordinate);
            try (SAMFileWriter outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(inputSam.getFileHeader(), false, output);
                    SAMRecordIterator iterator=inputSam.iterator()) {
                boolean paired=false;
                String readName=null; // the name of the last read
                ArrayList<SAMRecord> recordsTheRead=new ArrayList<>(); // the records from the reads with their names the same as the last read
                
                boolean theLastLoop=false;
                SAMRecord record=null;
                while(iterator.hasNext() || theLastLoop){
                    if(!theLastLoop)
                        record = iterator.next();

                    if (record.getReadPairedFlag()) {
                        paired = true;
                    } else if (paired) {
                        throw new Exception("Error: The input SAM/BAM file should only contain either single-ended or paired-ended reads.");
                    }
                    String id = record.getReadName();
                    if (paired) {
                        id = id.replaceAll("[12]$", "");
                    }
                    
                    // the current read has the same ID as the previous one, then save it
                    if((readName==null || id.equals(readName)) && !theLastLoop){
                        readName=id;
                        recordsTheRead.add(record);
                        if(!iterator.hasNext())
                            theLastLoop=true;
                    }
                    
                    // the current read has a different ID, indicating that all the reads with the previous read ID have been saved. start the processing
                    if(!iterator.hasNext() || (readName!=null  && !id.equals(readName))){
                        // single-ended
                        if(!paired && recordsTheRead.size()==1){
                            outputSam.addAlignment(recordsTheRead.get(0));
                            numUnique++;
                            
                        // paired-ended
                        } else if(paired){
                            HashMap<String, ArrayList<SAMRecord>> properPairs=new HashMap<>();
                            for(SAMRecord thisRecord : recordsTheRead){
                                if(thisRecord.getProperPairFlag()){ // only consider the proper pairs
                                    String pairID=null;
                                    if(thisRecord.getFirstOfPairFlag()){
                                        pairID=thisRecord.getReferenceName()+":"+Integer.toString(thisRecord.getAlignmentStart())+" "
                                            +thisRecord.getMateReferenceName()+":"+Integer.toString(thisRecord.getMateAlignmentStart());
                                    } else if(thisRecord.getSecondOfPairFlag()){
                                        pairID=thisRecord.getMateReferenceName()+":"+Integer.toString(thisRecord.getMateAlignmentStart())+" "
                                            +thisRecord.getReferenceName()+":"+Integer.toString(thisRecord.getAlignmentStart());
                                    }
                                    if(pairID!=null){
                                        if (!properPairs.containsKey(pairID))
                                            properPairs.put(pairID, new ArrayList<SAMRecord>());
                                        properPairs.get(pairID).add(thisRecord);
                                    }
                                }
                            }
                            
                            // check the number of proper pairs
                            if(properPairs.size()==1){
                                numUnique++;
                                //System.out.println(id);
                                for(SAMRecord outputRecord : properPairs.get(properPairs.keySet().iterator().next()))
                                    outputSam.addAlignment(outputRecord);
                            }
                        }
                        
                        recordsTheRead.clear();
                        recordsTheRead.add(record);
                        readName=id;
                        if(!theLastLoop && !iterator.hasNext()){
                            theLastLoop=true;
                        } else if(!iterator.hasNext()){
                            theLastLoop=false;
                        }
                    }
                }
                
                CloserUtil.close(inputSam);
                outputSam.close();
                if(tempGenerated)
                    tempFile.delete();
                return(numUnique);
            }
            catch (Exception e){
                System.err.println(e);
                System.err.println("Error in the extraction process: extractUniquelyMappedReads of MappingProcessor.");
                return(-1);
            }
        }
    }
    int extractUniquelyMappedReads(File output){
        int numUnique=0;
        boolean withNHTag = false;
        boolean paired = false;
        
        HashMap<String, Integer> numHits = new HashMap<>(); // for single-ended RNA-seq without NH tag
        HashMap<String, HashSet<String>> properPairs = new HashMap<>(); // for paired-ended RNA-seq
        
        SamReader inputSam=SamReaderFactory.makeDefault().open(inputFile);
        SAMFileWriter outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(inputSam.getFileHeader(), true, output);
        
        try(SAMRecordIterator iterator=inputSam.iterator()){
            while(iterator.hasNext()) {
                SAMRecord record=iterator.next();
                
                if(!record.getReadPairedFlag()){ // Single-ended RNA-seq
                    if(paired){
                        System.err.println("Error: Records in one SAM/BAM file should be all unpaired or paired.");
                        throw new Exception();
                    }
                    
                    if (record.getIntegerAttribute("NH") != null) { // there is NH tag, which shows the number of hits of this read
                        withNHTag = true;
                        if (record.getIntegerAttribute("NH").equals(1)) {
                            outputSam.addAlignment(record);
                            numUnique++;
                        }
                    } else { // there is no NH tag
                        if (numHits.containsKey(record.getReadName())) {
                            numHits.put(record.getReadName(), numHits.get(record.getReadName()) + 1);
                        } else {
                            numHits.put(record.getReadName(), 1);
                        }
                    }
                } else if(record.getProperPairFlag() && record.getFirstOfPairFlag()){ // Paired-ended RNA-seq, but only consider the properly paired ones
                    paired=true;
                    
                    String id=record.getReadName();
                    id=id.replaceAll("{1,2}$", "");
                    String idWithCoord=record.getReferenceName()+":"+Integer.toString(record.getAlignmentStart())
                            +" "+record.getMateReferenceName()+":"+Integer.toString(record.getMateAlignmentStart());
                    
                    properPairs.put(id, new HashSet<String>());
                    properPairs.get(id).add(idWithCoord);
                }
            }
            iterator.close();
        }
        catch(Exception e){
            System.err.println("ERROR! "+e+" in extractUniquelyMappedReads at MappingProcessor!\n");
            System.exit(1);
            return(-1);
        }
        
        // if the input SAM/BAM has no NH tag, or it is paired-ended, a re-scanning of the SAM/BAM file is needed
        if (!withNHTag || paired) {
            System.err.println("Start re-scanning the SAM/BAM file...");
            if(! paired){
                CloserUtil.close(inputSam);
                inputSam = SamReaderFactory.makeDefault().open(inputFile);
                for (SAMRecord record : inputSam) {
                    if (numHits.get(record.getReadName()) == 1) {
                        record.setAttribute("NH", 1);
                        outputSam.addAlignment(record);
                        numUnique++;
                    }
                }
            } else{
                CloserUtil.close(inputSam);
                inputSam = SamReaderFactory.makeDefault().open(inputFile);
                for (SAMRecord record : inputSam){
                    String id=record.getReadName();
                    id=id.replaceAll("{1,2}$", "");
                    if(properPairs.containsKey(id))
                        if(record.getFirstOfPairFlag()){
                            String idWithCoord=record.getReferenceName()+":"+Integer.toString(record.getAlignmentStart())
                                    +" "+record.getMateReferenceName()+":"+Integer.toString(record.getMateAlignmentStart());
                            if(properPairs.get(id).contains(idWithCoord))
                                outputSam.addAlignment(record);
                        } else if(record.getSecondOfPairFlag()){
                            String idWithCoord=record.getMateReferenceName()+":"+Integer.toString(record.getMateAlignmentStart())
                                    +" "+record.getReferenceName()+":"+Integer.toString(record.getAlignmentStart());
                            if(properPairs.get(id).contains(idWithCoord)){
                                outputSam.addAlignment(record);
                                numUnique++;
                            }
                        }
                }
            }
        }

        CloserUtil.close(inputSam);
        CloserUtil.close(outputSam);
        return(numUnique);
    }
}
