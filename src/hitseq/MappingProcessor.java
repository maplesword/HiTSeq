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
import java.util.regex.*;
import hitseq.annotation.*;
import java.io.OutputStream;
import java.io.PrintStream;
import java.nio.file.Files;

/**
 * The class for SAM/BAM processing, including information collecting and uniquely mapped reads extraction
 * @author hezhisong
 */
public class MappingProcessor {
    private File inputFile;
    private ArrayList<Integer> info;
    
    public MappingProcessor(File file){
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
    
    public int extractUniquelyMappedReads(File output, boolean attemptNH, boolean onlyChr){
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
                        
                        numUnique=this.extractUniquelyMappedReads(output, false, onlyChr);
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
            File tempFile=new File(System.getProperty("java.io.tmpdir")+"/tmp_sorting.bam");
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
                        Pattern pattern=Pattern.compile("\\W[12]$");
                        Matcher matcher=pattern.matcher(id);
                        if(matcher.find())
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
                                if(thisRecord.getProperPairFlag() || (onlyChr && thisRecord.getReferenceName().equals(thisRecord.getMateReferenceName()))){ 
                                // onlyChr is false: only consider the proper pairs
                                // onlyChr is true: also consider other pairs with both mates at the same chromosome
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
                                String pairID=properPairs.keySet().iterator().next();
                                if(properPairs.get(pairID).size()==2){
                                    numUnique++;
                                    for (SAMRecord outputRecord : properPairs.get(properPairs.keySet().iterator().next()))
                                        outputSam.addAlignment(outputRecord);
                                }
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
    
    public int reassignProperPairFlag(File output, int stranded, Annotation annotation){
        File tempFileProper=new File(System.getProperty("java.io.tmpdir")+"/tmp_proper.bam");
        File tempFileImproper=new File(System.getProperty("java.io.tmpdir")+"/tmp_improper.bam");
        File tempFileOther=new File(System.getProperty("java.io.tmpdir")+"/tmp_others.bam");
        File tempFileImproperCorrected=new File(System.getProperty("java.io.tmpdir")+"/tmp_improper_corrected.bam");
        boolean isPaired=true;
        int numCorrected = 0;
        
        SamReader inputSam=SamReaderFactory.makeDefault().open(inputFile);
        try{
            SAMFileHeader headerProper=inputSam.getFileHeader();
            headerProper.setSortOrder(SAMFileHeader.SortOrder.queryname);
            SAMFileWriter properSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(headerProper, false, tempFileProper);
            SAMFileWriter otherSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(headerProper, false, tempFileOther);
            SAMFileHeader headerImproper=inputSam.getFileHeader();
            headerImproper.setSortOrder(SAMFileHeader.SortOrder.coordinate);
            SAMFileWriter improperSam = new SAMFileWriterFactory().makeBAMWriter(headerImproper, false, tempFileImproper);
            SAMRecordIterator iterator=inputSam.iterator();
            int numImproper = 0;
            
            while(iterator.hasNext()){
                SAMRecord record = iterator.next();
                if(! record.getReadPairedFlag()){
                    isPaired = false;
                    break;
                }
                
                if(record.getProperPairFlag())
                    properSam.addAlignment(record);
                else if((!(record.getReadUnmappedFlag() || record.getMateUnmappedFlag())) && record.getReferenceName().equals(record.getMateReferenceName())){
                    improperSam.addAlignment(record);
                    numImproper++;
                } else{
                    otherSam.addAlignment(record);
                }
            }
            CloserUtil.close(inputSam);
            CloserUtil.close(properSam);
            CloserUtil.close(improperSam);
            CloserUtil.close(otherSam);
            
            if(isPaired && numImproper>0){
                SAMFileWriter outputCorrectedSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(headerProper, false, tempFileImproperCorrected);
                inputSam = SamReaderFactory.makeDefault().open(tempFileImproper);
                iterator = inputSam.iterator();
                int recordNum = 0;
                
                // 1st scan the improper reads, determine the additional proper pairs
                HashMap<String, ArrayList<ArrayList<String>>> readGenes=new HashMap<>();
                HashMap<String, ArrayList<Integer>> readIdx=new HashMap<>();
                HashSet<Integer> properReadIdx = new HashSet<>();
                while(iterator.hasNext()){
                    SAMRecord record = iterator.next();
                    SAMRecordProcessor proccessor = new SAMRecordProcessor(record, annotation);
                    String recordID = proccessor.recordToString()+"|"+Boolean.toString(record.getFirstOfPairFlag());
                    String mateID = proccessor.recordMateToString()+"|"+Boolean.toString(! record.getFirstOfPairFlag());
                    
                    if(readIdx.containsKey(mateID)){
                        ArrayList<String> genesThis = proccessor.getOverlapGenes(stranded);
                        ArrayList<String> genesMate = readGenes.get(mateID).get(0);
                        HashSet<String> allGenes=new HashSet<>();
                        for(String gene : genesThis)
                            allGenes.add(gene);
                        for(String gene : genesMate)
                            allGenes.add(gene);
                        if(allGenes.size() < genesThis.size()+genesMate.size()  // have overlapping genes
                                && (!(record.getReadNegativeStrandFlag() && record.getMateNegativeStrandFlag()))){ // one forward one reverse
                            properReadIdx.add(recordNum);
                            properReadIdx.add(readIdx.get(mateID).get(0));
                            numCorrected++;
                        }
                        
                        readGenes.get(mateID).remove(0);
                        if(readGenes.get(mateID).isEmpty())
                            readGenes.remove(mateID);
                        readIdx.get(mateID).remove(0);
                        if(readIdx.get(mateID).isEmpty())
                            readIdx.remove(mateID);
                        
                        
                    } else{
                        if(! readIdx.containsKey(recordID))
                            readIdx.put(recordID, new ArrayList<Integer>());
                        readIdx.get(recordID).add(recordNum);
                        
                        ArrayList<String> genes = proccessor.getOverlapGenes(stranded);
                        if(! readGenes.containsKey(recordID))
                            readGenes.put(recordID, new ArrayList<ArrayList<String>>());
                        readGenes.get(recordID).add(genes);
                    }
                    
                    recordNum++;
                }
                
                // 2nd scan of the improper reads, re-assign the proper pair flag to the additional ones
                // write the corrected reads to the file; if no change is needed, simply copy the file
                if(properReadIdx.size()>0){
                    iterator = inputSam.iterator();
                    recordNum = 0;
                    while(iterator.hasNext()){
                        SAMRecord record = iterator.next();
                        if(properReadIdx.contains(recordNum))
                            record.setProperPairFlag(true);
                        outputCorrectedSam.addAlignment(record);
                    }
                } else{
                    Files.copy(tempFileImproper.toPath(), tempFileImproperCorrected.toPath(), java.nio.file.StandardCopyOption.REPLACE_EXISTING);
                }
                Files.delete(tempFileImproper.toPath());
                CloserUtil.close(inputSam);
                CloserUtil.close(outputCorrectedSam);
            }
            
            // output the corrected SAM files and sort by queryname 
            SamReader inputProper = SamReaderFactory.makeDefault().open(tempFileProper);
            SamReader inputImproperCorrected = SamReaderFactory.makeDefault().open(tempFileImproperCorrected);
            SamReader inputOther = SamReaderFactory.makeDefault().open(tempFileOther);
            SAMFileWriter outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(headerProper, false, output);
            for(SAMRecord record : inputProper)
                outputSam.addAlignment(record);
            for(SAMRecord record : inputImproperCorrected)
                outputSam.addAlignment(record);
            for(SAMRecord record : inputOther)
                outputSam.addAlignment(record);
            CloserUtil.close(inputProper);
            CloserUtil.close(inputImproperCorrected);
            CloserUtil.close(inputOther);
            CloserUtil.close(outputSam);
        }
        catch(Exception e){
            System.err.println(e);
        }
        
        return(numCorrected);
    }
}
