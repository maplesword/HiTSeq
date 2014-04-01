/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hitseq;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * @author hezhisong
 */
public class JunctionSet {
    HashMap<String,ArrayList<Junction>> juncInChrom;
    
    JunctionSet(File file, String fileType){
        if(! file.exists()){
            System.err.println("Cannot find file "+file.getAbsolutePath()+".");
            System.exit(1);
        }
        juncInChrom=new HashMap();
        
        try{
            RandomAccessFile fileIn=new RandomAccessFile(file,"r");
            int numLines=0;
            String line;
            while((line=fileIn.readLine()) != null){ // deal with each line separately
                if(line.startsWith("#"))
                    continue;
                numLines++;
                
                String[] elements=line.split("\t");
                switch (fileType) {
                    case "bed": // Fields: chrom, donor(3'end of last exon), acceptor(5'end of next exon), id, score, strand
                        String chrom=elements[0];
                        int donorSite=Integer.parseInt(elements[1])+1;
                        int acceptorSite=Integer.parseInt(elements[2]);
                        String id=elements.length>3 ? elements[3] : "junc"+numLines;
                        String strand=elements.length>5 ? elements[5] : "*";
                        
                        if(!juncInChrom.containsKey(chrom))
                            juncInChrom.put(chrom, new ArrayList<Junction>());
                        juncInChrom.get(chrom).add(new Junction(id,chrom,strand,donorSite,acceptorSite));
                        
                        break;
                    case "gtf":
                        break;
                    case "junc":
                        break;
                    case "bam":
                        break;
                }
            }
        }
        catch(IOException | NumberFormatException e){
            System.err.println("Error in JunctionSet.java: "+e);
        }
    }
}
