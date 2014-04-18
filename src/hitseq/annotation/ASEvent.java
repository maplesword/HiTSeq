/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hitseq.annotation;

import java.util.Collection;
import java.util.HashSet;
import java.util.Objects;

/**
 *
 * @author hezhisong
 */
public class ASEvent {
    public static final int SKIPP_EXON=0;
    public static final int MUTUAL_EXCL=1;
    public static final int ALT_END=2;
    
    private int type;
    private HashSet<Junction> inclusiveJunctionsP1;
    private HashSet<Junction> inclusiveJunctionsP2;
    private HashSet<Junction> exclusiveJunctionsP1;
    private HashSet<Junction> exclusiveJunctionsP2;
    
    public ASEvent(int type, Collection<Junction> inclusiveJunctionsP1, Collection<Junction> inclusiveJunctionsP2, Collection<Junction> exclusiveJunctionsP1, Collection<Junction> exclusiveJunctionsP2){
        if(type!=1){
            System.err.println("AS Event type error. This constructor is only for mutual exclusive exons {1}.");
            System.exit(1);
        }
        this.inclusiveJunctionsP1=new HashSet<>();
        this.inclusiveJunctionsP2=new HashSet<>();
        this.exclusiveJunctionsP1=new HashSet<>();
        this.exclusiveJunctionsP2=new HashSet<>();
        this.inclusiveJunctionsP1.addAll(inclusiveJunctionsP1);
        this.inclusiveJunctionsP2.addAll(inclusiveJunctionsP2);
        this.exclusiveJunctionsP1.addAll(exclusiveJunctionsP1);
        this.exclusiveJunctionsP2.addAll(exclusiveJunctionsP2);
    }
    
    public ASEvent(int type, Collection<Junction> inclusiveJunctionsP1, Collection<Junction> inclusiveJunctionsP2, Collection<Junction> exclusiveJunctions){
        if(type!=0){
            System.err.println("AS Event type error. This constructor is only for skipped exons {0}.");
            System.exit(1);
        }
        this.inclusiveJunctionsP1=new HashSet<>();
        this.inclusiveJunctionsP2=new HashSet<>();
        this.exclusiveJunctionsP1=new HashSet<>();
        this.exclusiveJunctionsP2=new HashSet<>();
        this.inclusiveJunctionsP1.addAll(inclusiveJunctionsP1);
        this.inclusiveJunctionsP2.addAll(inclusiveJunctionsP2);
        this.exclusiveJunctionsP1.addAll(exclusiveJunctions);
    }
    
    public ASEvent(int type, Collection<Junction> inclusiveJunctions, Collection<Junction> exclusiveJunctions){
        if(type!=2){
            System.err.println("AS Event type error. This constructor is only for alternative 3'/5' end {2}.");
            System.exit(1);
        }
        this.inclusiveJunctionsP1=new HashSet<>();
        this.inclusiveJunctionsP2=new HashSet<>();
        this.exclusiveJunctionsP1=new HashSet<>();
        this.exclusiveJunctionsP2=new HashSet<>();
        this.inclusiveJunctionsP1.addAll(inclusiveJunctions);
        this.exclusiveJunctionsP1.addAll(exclusiveJunctions);
    }
    
    public HashSet<Junction> getInclusiveJunctions(){
        HashSet<Junction> incluJuncsReturn=new HashSet<>();
        incluJuncsReturn.addAll(inclusiveJunctionsP1);
        if(!inclusiveJunctionsP2.isEmpty()) incluJuncsReturn.addAll(inclusiveJunctionsP2);
        return(incluJuncsReturn);
    }
    
    public HashSet<Junction> getExclusiveJunctions(){
        HashSet<Junction> excluJuncsReturn=new HashSet<>();
        excluJuncsReturn.addAll(exclusiveJunctionsP1);
        if(!exclusiveJunctionsP2.isEmpty()) excluJuncsReturn.addAll(exclusiveJunctionsP2);
        return(excluJuncsReturn);
    }
    
    public int getType(){
        return(type);
    }
    
    public boolean containAll(Collection<Junction> juncs){
        boolean answer=true;
        for(Junction junc : juncs){
            if(!inclusiveJunctionsP1.contains(junc) && !exclusiveJunctionsP1.contains(junc) && !inclusiveJunctionsP2.contains(junc) && !exclusiveJunctionsP2.contains(junc)){
                answer=false;
                break;
            }
        }
        return(answer);
    }
    
    @Override
    public boolean equals(Object event2){
        if(!(event2 instanceof ASEvent))
            return(false);
        ASEvent asEvent=(ASEvent) event2;
        if(asEvent.type!=type)
            return(false);
        if(asEvent.inclusiveJunctionsP1.size()!=inclusiveJunctionsP1.size())
            return(false);
        if(asEvent.exclusiveJunctionsP1.size()!=exclusiveJunctionsP1.size())
            return(false);
        if(asEvent.inclusiveJunctionsP2.size()!=inclusiveJunctionsP2.size())
            return(false);
        if(asEvent.exclusiveJunctionsP2.size()!=exclusiveJunctionsP2.size())
            return(false);
        for(Junction junc : asEvent.inclusiveJunctionsP1)
            if(!inclusiveJunctionsP1.contains(junc))
                return(false);
        for(Junction junc : asEvent.exclusiveJunctionsP1)
            if(!exclusiveJunctionsP1.contains(junc))
                return(false);
        for(Junction junc : asEvent.inclusiveJunctionsP2)
            if(!inclusiveJunctionsP2.contains(junc))
                return(false);
        for(Junction junc : asEvent.exclusiveJunctionsP2)
            if(!exclusiveJunctionsP2.contains(junc))
                return(false);
        return(true);
    }

    @Override
    public int hashCode() {
        int hash = 3;
        hash = 47 * hash + this.type;
        hash = 47 * hash + Objects.hashCode(this.inclusiveJunctionsP1);
        hash = 47 * hash + Objects.hashCode(this.exclusiveJunctionsP1);
        return hash;
    }
}
