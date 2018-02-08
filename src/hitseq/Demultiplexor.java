/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hitseq;

import java.io.File;
import java.util.HashSet;
import java.util.HashMap;
import java.util.ArrayList;
import htsjdk.samtools.*;
import hitseq.annotation.DivergentSiteSet;
import hitseq.annotation.DivergentSite;

/**
 *
 * @author Chih-sung
 */
public class Demultiplexor {
    final private File samFile;
    final private DivergentSiteSet divergentSites;
    final private HashSet<String> barcodes;
    final boolean onlyUnique;
    final boolean filter;
    
    public Demultiplexor(File file, DivergentSiteSet divergentSites, HashSet<String> barcodes, boolean onlyUnique, boolean filter){
        this.samFile = file;
        this.divergentSites = divergentSites;
        this.barcodes = barcodes;
        this.onlyUnique = onlyUnique;
        this.filter = filter;
    }
    
    public Demultiplexor(File file, DivergentSiteSet divergentSites, HashSet<String> barcodes, boolean onlyUnique){
        this(file, divergentSites, barcodes, onlyUnique, true);
    }
    
    public Demultiplexor(File file, DivergentSiteSet divergentSites, HashSet<String> barcodes){
        this(file, divergentSites, barcodes, true);
    }
    
    public Demultiplexor(File file, DivergentSiteSet divergentSites){
        this(file, divergentSites, null);
    }
    
    public HashMap<String, double[]> demultiplex(boolean onlyUniqueSites){
        HashMap<String, double[]> results = new HashMap<>(); // KEY: CB; VALUE: array of total match scores
        try(SamReader inputSam = SamReaderFactory.makeDefault().open(samFile)) {
            if(! inputSam.getFileHeader().getSortOrder().equals(SAMFileHeader.SortOrder.coordinate))
                throw new IllegalArgumentException("Critical error: the given SAM/BAM file is not sorted by coordinate");
            
            ArrayList<SAMRecord> candidates = new ArrayList<>();
            DivergentSite currentSite = null;
            for (SAMRecord record : inputSam) {
                if(record.getStringAttribute("CB") == null) // cellular barcode is a must
                    continue;
                if(barcodes != null && ! barcodes.contains(record.getStringAttribute("CB"))) // if provided, only consider CBs which are unlikely to be background
                    continue;
                if(onlyUnique && record.getIntegerAttribute("NH") != null && ! record.getIntegerAttribute("NH").equals(1)) // when requiring uniquely mapped reads, filter based on the NH tag
                    continue;
                if(filter && (record.getReadUnmappedFlag() || record.getNotPrimaryAlignmentFlag()|| record.getDuplicateReadFlag() || record.getReadFailsVendorQualityCheckFlag()))
                    continue;
                
                String chrom=record.getReferenceName();
                if(currentSite == null && (currentSite = divergentSites.next(chrom)) == null){
                    continue;
                }
                boolean chromChanged = !currentSite.getChromosome().equals(chrom);
                
                // this record goes into a different chromosome or passes the current site
                while(! candidates.isEmpty() && (chromChanged || currentSite.getCoordinate() < record.getAlignmentStart())){
                    // scan all the potential records, get the ones covering the current site, and calculate the matching scores
                    HashMap<String, HashMap<String, double[]>> transcriptScores = new HashMap<>();
                    for(SAMRecord candidate : candidates){
                        int readOffset = candidate.getReadPositionAtReferencePosition(currentSite.getCoordinate());
                        if(readOffset != 0){ // this read covers the divergent site
                            String id = candidate.getStringAttribute("UB") == null ? candidate.getReadName() : candidate.getReadString() + "-" + candidate.getStringAttribute("UB");
                            String cb = candidate.getStringAttribute("CB");
                            if(! transcriptScores.containsKey(cb))
                                transcriptScores.put(cb, new HashMap<String, double[]>());
                            
                            if(! transcriptScores.get(cb).containsKey(id)){
                                String base = (new String(new byte[]{candidate.getReadBases()[readOffset - 1]})).toUpperCase();
                                double[] scores = currentSite.matchScore(base, onlyUniqueSites);
                                transcriptScores.get(cb).put(id, scores);
                                //System.out.println(candidate.getReadName() + "\t" + currentSite.getCoordinate() + "\t" + readOffset);
                            }
                        }
                    }
                    
                    // sum up the matching scores for each cellular barcode (each cell/nucleus)
                    for(String cb : transcriptScores.keySet()){
                        if(! results.containsKey(cb)){
                            results.put(cb, new double[currentSite.getReferenceBases().size() + 1]);
                            for(int i = 0; i < currentSite.getReferenceBases().size() + 1; i++)
                                results.get(cb)[i] = 0;
                        }
                        for(String id : transcriptScores.get(cb).keySet())
                            for(int i = 0; i < transcriptScores.get(cb).get(id).length; i++)
                                results.get(cb)[i] += transcriptScores.get(cb).get(id)[i];
                    }
                    
                    // move to the next site, remove current record candidates which have been passed
                    currentSite = divergentSites.next(currentSite.getChromosome());
                    if(currentSite == null){
                        candidates = new ArrayList<>();
                        //System.out.println("Clear");
                    } else
                        while(candidates.size() > 0){
                            SAMRecord candidate = candidates.get(0);
                            if(currentSite.getCoordinate() > candidate.getAlignmentEnd()){
                                //System.out.println("remove (" + candidates.size() + "): " + candidate.getReadName() + "\t<-\t" + currentSite.getCoordinate() + "\t" + candidate.getAlignmentStart() + "\t" + candidate.getCigarString());
                                candidates.remove(0);
                            } else
                                break;
                        }
                }
                
                // put this record into the candidate list
                if(chromChanged){
                    candidates = new ArrayList<>();
                    currentSite = divergentSites.next(chrom);
                    //System.err.println("...start to process chromosome " + chrom);
                }
                while(currentSite != null && currentSite.getCoordinate() < record.getAlignmentStart()){
                    currentSite = divergentSites.next(chrom);
                }
                if(currentSite != null){
                    candidates.add(record);
                    //System.out.println("adding (" + candidates.size() + "): " + record.getReadName() + "\t->\t" + currentSite.getCoordinate() + "\t" + record.getAlignmentStart() + "\t" + record.getCigarString());
                }
            }
            
            // if there are candidates remained after reading the whole BAM file, process them to finalize the whole analysis
            while (currentSite != null && candidates.size() > 0) {
                // scan all the potential records, get the ones covering the current site, and calculate the matching scores
                HashMap<String, HashMap<String, double[]>> transcriptScores = new HashMap<>();
                for (SAMRecord candidate : candidates) {
                    int readOffset = candidate.getReadPositionAtReferencePosition(currentSite.getCoordinate());
                    if (readOffset != 0) { // this read covers the divergent site
                        String id = candidate.getStringAttribute("UB") == null ? candidate.getReadName() : candidate.getReadString() + "-" + candidate.getStringAttribute("UB");
                        String cb = candidate.getStringAttribute("CB");
                        if (!transcriptScores.containsKey(cb))
                            transcriptScores.put(cb, new HashMap<String, double[]>());

                        if (!transcriptScores.get(cb).containsKey(id)) {
                            String base = (new String(new byte[]{candidate.getReadBases()[readOffset - 1]})).toUpperCase();
                            double[] scores = currentSite.matchScore(base, onlyUniqueSites);
                            transcriptScores.get(cb).put(id, scores);
                            //System.out.println(candidate.getReadName() + "\t" + currentSite.getCoordinate() + "\t" + readOffset);
                        }
                    }
                }

                // sum up the matching scores for each cellular barcode (each cell/nucleus)
                for (String cb : transcriptScores.keySet()) {
                    if (!results.containsKey(cb)) {
                        results.put(cb, new double[currentSite.getReferenceBases().size()+1]);
                        for (int i = 0; i < currentSite.getReferenceBases().size()+1; i++)
                            results.get(cb)[i] = 0;
                    }
                    for (String id : transcriptScores.get(cb).keySet())
                        for (int i = 0; i < transcriptScores.get(cb).get(id).length; i++)
                            results.get(cb)[i] += transcriptScores.get(cb).get(id)[i];
                }
                
                // move to the next site, remove current record candidates which have been passed
                currentSite = divergentSites.next(currentSite.getChromosome());
                if(currentSite == null)
                    candidates = new ArrayList<>();
                else
                    while(candidates.size() > 0){
                        SAMRecord candidate = candidates.get(0);
                        if(currentSite.getCoordinate() > candidate.getAlignmentEnd())
                            candidates.remove(0);
                        else
                            break;
                    }
            }
        } catch (Exception e) {
            System.err.println("Critical error when demultiplexing cells: Method demultiplex in Class Demultiplexor");
            e.printStackTrace();
            System.exit(1);
        }
        
        return results;
    }
    
    public static String output(HashMap<String, double[]> matches, String[] refnames){       
        int numRefs = -1;
        String results = "";
        ArrayList<String> cellBarcodes = new ArrayList<>();
        cellBarcodes.addAll(matches.keySet());
        cellBarcodes.sort(null);
        for(String cb : cellBarcodes){
            if(numRefs == -1)
                numRefs = matches.get(cb).length;
            results += cb;
            for(double score : matches.get(cb))
                results += "\t" + score;
            results += System.getProperty("line.separator");
        }
        
        String header = "Barcode";
        if(refnames != null){
            for(String ref : refnames)
                header += "\t" + ref;
        } else
            for(int i=0; i<numRefs; i++)
                header += "\tC" + i;
        header += System.getProperty("line.separator");

        String output = header + results;
        return output;
    }
}
