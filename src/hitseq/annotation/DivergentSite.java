/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hitseq.annotation;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Objects;

/**
 *
 * @author Chih-sung
 */
public class DivergentSite implements Comparable<DivergentSite>{
    private final String chr;
    private final int coord;
    private final ArrayList<String> refBases;
    
    public DivergentSite(String chr, int coord, ArrayList<String> refBases){
        this.chr = chr;
        this.coord = coord;
        this.refBases = refBases;
        
        HashSet<String> bases = new HashSet<>();
        for(String base : refBases)
            bases.add(base);
        if(bases.size() == 1)
            throw new IllegalArgumentException("There should be at least one species with a different base!");
    }
    
    public String getChromosome(){
        return chr;
    }
    
    public int getCoordinate(){
        return coord;
    }
    
    public ArrayList<String> getReferenceBases(){
        return (ArrayList<String>)refBases.clone();
    }
    
    public boolean[] match(String base){
        boolean[] match = new boolean[refBases.size() + 1];
        boolean matched = false;
        for(int i = 0; i < refBases.size(); i++){
            match[i] = refBases.get(i).equals(base);
            matched = matched || match[i];
        }
        match[match.length-1] = ! matched;
        
        return match;
    }
    
    public double[] matchScore(String base, boolean onlyUniqueSite){
        boolean[] matches = match(base);
        int numMatches = 0;
        for(boolean match : matches)
            if(match)
                numMatches++;
        
        double[] scores = new double[matches.length];
        for(int i = 0; i < matches.length; i++){
            if(onlyUniqueSite)
                scores[i] = matches[i] & numMatches == 1 ? 1 : 0;
            else
                scores[i] = matches[i] ? 1.0/numMatches : 0;
        }
        return scores;
    }
    
    public double[] matchScore(String base){
        return matchScore(base, true);
    }
    
    public boolean conflict(DivergentSite site2){
        return this.chr.equals(site2.chr) && this.coord == site2.coord && ! this.refBases.equals(site2.refBases);
    }
    
    @Override
    public int compareTo(DivergentSite site2){
        int compareChr = this.chr.compareTo(site2.getChromosome());
        int compareCoord = new Integer(this.coord).compareTo(new Integer(site2.getCoordinate()));
        int compare = compareChr == 0 ? compareCoord : compareChr;
        return compare;
    }
    
    @Override
    public boolean equals(Object site2){
        if(site2 == this)
            return true;
        if(!(site2 instanceof DivergentSite))
            return false;
        
        boolean equalChr = this.chr.equals(((DivergentSite) site2).getChromosome());
        boolean equalCoord = this.coord == ((DivergentSite) site2).getCoordinate();
        boolean equalRefBases = this.refBases.equals(((DivergentSite) site2).getReferenceBases());
        return equalChr && equalCoord && equalRefBases;
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 23 * hash + Objects.hashCode(this.chr);
        hash = 23 * (23 * hash + this.coord);
        for(String base : refBases)
            hash += Objects.hashCode(base);
        return hash;
    }
}
