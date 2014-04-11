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
    public static final int ALT_FIVE=2;
    public static final int ALT_THREE=3;
    
    private int type;
    private HashSet<Junction> inclusiveJunctionsP1;
    private HashSet<Junction> inclusiveJunctionsP2;
    private HashSet<Junction> exclusiveJunctions;
    
    public ASEvent(int type, Collection<Junction> inclusiveJunctionsP1, Collection<Junction> inclusiveJunctionsP2, Collection<Junction> exclusiveJunctions){
        if(type<0 || type>3){
            System.err.println("AS Event type error. It should be an integer in [0,3].");
            System.exit(1);
        }
        this.inclusiveJunctionsP1=new HashSet<>();
        this.inclusiveJunctionsP2=new HashSet<>();
        this.exclusiveJunctions=new HashSet<>();
        this.inclusiveJunctionsP1.addAll(inclusiveJunctionsP1);
        this.inclusiveJunctionsP2.addAll(inclusiveJunctionsP2);
        this.exclusiveJunctions.addAll(exclusiveJunctions);
    }
    
    public ASEvent(int type, Collection<Junction> inclusiveJunctions, Collection<Junction> exclusiveJunctions){
        if(type!=2 && type!=3){
            System.err.println("AS Event type error. This constructor is only for alternative 3'/5' end {2,3}.");
            System.exit(1);
        }
        this.inclusiveJunctionsP1=new HashSet<>();
        this.inclusiveJunctionsP2=new HashSet<>();
        this.exclusiveJunctions=new HashSet<>();
        this.inclusiveJunctionsP1.addAll(inclusiveJunctionsP1);
        this.exclusiveJunctions.addAll(exclusiveJunctions);
    }
    
    public void addInclusiveJunction(Junction newIncluJunc){
        inclusiveJunctionsP1.add(newIncluJunc);
    }
    
    public void addExclusiveJunction(Junction newExcluJunc){
        exclusiveJunctions.add(newExcluJunc);
    }
    
    public HashSet<Junction> getInclusiveJunctions(){
        HashSet<Junction> incluJuncsReturn=new HashSet<>();
        incluJuncsReturn.addAll(inclusiveJunctionsP1);
        return(incluJuncsReturn);
    }
    
    public HashSet<Junction> getExclusiveJunctions(){
        HashSet<Junction> excluJuncsReturn=new HashSet<>();
        excluJuncsReturn.addAll(exclusiveJunctions);
        return(excluJuncsReturn);
    }
    
    public int getType(){
        return(type);
    }
    
    @Override
    public boolean equals(Object event2){
        if(!(event2 instanceof ASEvent))
            return(false);
        ASEvent asEvent=(ASEvent) event2;
        if(asEvent.type!=type)
            return(false);
        for(Junction junc : asEvent.inclusiveJunctionsP1)
            if(!inclusiveJunctionsP1.contains(junc))
                return(false);
        for(Junction junc : asEvent.exclusiveJunctions)
            if(!exclusiveJunctions.contains(junc))
                return(false);
        return(true);
    }

    @Override
    public int hashCode() {
        int hash = 3;
        hash = 47 * hash + this.type;
        hash = 47 * hash + Objects.hashCode(this.inclusiveJunctionsP1);
        hash = 47 * hash + Objects.hashCode(this.exclusiveJunctions);
        return hash;
    }
}
