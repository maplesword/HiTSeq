/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hitseq.annotation;

import java.util.HashMap;
import java.util.ArrayList;
import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;

/**
 *
 * @author Chih-sung
 */
public class DivergentSiteSet {
    private HashMap<String, ArrayList<DivergentSite>> sitesInChromosomes;
    private HashMap<String, Integer> pointers;
    private int processed;
    
    public DivergentSiteSet(){
        sitesInChromosomes = new HashMap<>();
        pointers = new HashMap<>();
        processed = -1;
    }
    
    public DivergentSiteSet(File file){
        this();
        addSitesInFile(file);
    }
    
    public final void addSitesInFile(File file){
        try{
            RandomAccessFile reader = new RandomAccessFile(file, "r");
            String line;
            int numRefs = -1;
            while((line = reader.readLine()) != null) {
                String[] elem = line.split("\t");
                if(elem.length <= 3)
                    throw new IllegalArgumentException("Each divergent site should contain at least two alternative bases");
                
                String chr = elem[0];
                int coord = Integer.parseInt(elem[1]);
                ArrayList<String> refBases = new ArrayList<>();
                for(int i = 2; i < elem.length; i++)
                    refBases.add(elem[i]);
                if(numRefs == -1)
                    numRefs = refBases.size();
                if(numRefs != refBases.size())
                    throw new IllegalArgumentException("Every divergent site should have the same number of alternative bases");
                
                DivergentSite newSite = new DivergentSite(chr, coord, refBases);
                
                if(! sitesInChromosomes.containsKey(chr))
                    sitesInChromosomes.put(chr, new ArrayList<DivergentSite>());
                sitesInChromosomes.get(chr).add(newSite);
            }
            reader.close();
        } catch(IOException e){
            System.err.println(e);
            System.err.println("Critical I/O error: cannot read the divergent site list file: " + file.getPath());
        } catch (IllegalArgumentException e){
            System.err.println(e);
        }
        
        for(String chr : sitesInChromosomes.keySet()){
            sitesInChromosomes.get(chr).sort(null);
            pointers.put(chr, -1);
        }
    }
    
    public DivergentSite next(String chr){
        if(! sitesInChromosomes.containsKey(chr))
            return null;
        
        pointers.put(chr, pointers.get(chr) + 1);
        if(pointers.get(chr) >= sitesInChromosomes.get(chr).size())
            return null;
        else{
            processed++;
            if(processed != 0 && processed % 100000 == 0)
                System.err.println("...processed " + processed + " sites");
            return sitesInChromosomes.get(chr).get(pointers.get(chr));
        }
    }
    
    public boolean containChromosome(String chr){
        return sitesInChromosomes.containsKey(chr);
    }
    
    public void reset(){
        for(String chr : pointers.keySet())
            pointers.put(chr, -1);
        processed = -1;
    }
}
