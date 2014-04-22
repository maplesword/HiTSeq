/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hitseq.annotation;

import java.util.*;

/**
 *
 * @author hezhisong
 */
public class ASEvent {
    public static final int SKIPP_EXON=0;
    public static final int MUTUAL_EXCL=1;
    public static final int ALT_END=2;
    
    private String group;
    private int type;
    private HashSet<Junction> inclusiveJunctionsP1;
    private HashSet<Junction> inclusiveJunctionsP2;
    private HashSet<Junction> exclusiveJunctionsP1;
    private HashSet<Junction> exclusiveJunctionsP2;
    
    public ASEvent(String group, int type, Collection<Junction> inclusiveJunctionsP1, Collection<Junction> inclusiveJunctionsP2, Collection<Junction> exclusiveJunctionsP1, Collection<Junction> exclusiveJunctionsP2){
        if(type!=1){
            System.err.println("AS Event type error. This constructor is only for mutual exclusive exons {1}.");
            System.exit(1);
        }
        this.group=group;
        this.type=type;
        this.inclusiveJunctionsP1=new HashSet<>();
        this.inclusiveJunctionsP2=new HashSet<>();
        this.exclusiveJunctionsP1=new HashSet<>();
        this.exclusiveJunctionsP2=new HashSet<>();
        this.inclusiveJunctionsP1.addAll(inclusiveJunctionsP1);
        this.inclusiveJunctionsP2.addAll(inclusiveJunctionsP2);
        this.exclusiveJunctionsP1.addAll(exclusiveJunctionsP1);
        this.exclusiveJunctionsP2.addAll(exclusiveJunctionsP2);
    }
    
    public ASEvent(String group, int type, Collection<Junction> inclusiveJunctionsP1, Collection<Junction> inclusiveJunctionsP2, Collection<Junction> exclusiveJunctions){
        if(type!=0){
            System.err.println("AS Event type error. This constructor is only for skipped exons {0}.");
            System.exit(1);
        }
        this.group=group;
        this.type=type;
        this.inclusiveJunctionsP1=new HashSet<>();
        this.inclusiveJunctionsP2=new HashSet<>();
        this.exclusiveJunctionsP1=new HashSet<>();
        this.exclusiveJunctionsP2=new HashSet<>();
        this.inclusiveJunctionsP1.addAll(inclusiveJunctionsP1);
        this.inclusiveJunctionsP2.addAll(inclusiveJunctionsP2);
        this.exclusiveJunctionsP1.addAll(exclusiveJunctions);
    }
    
    public ASEvent(String group, int type, Collection<Junction> inclusiveJunctions, Collection<Junction> exclusiveJunctions){
        if(type!=2){
            System.err.println("AS Event type error. This constructor is only for alternative 3'/5' end {2}.");
            System.exit(1);
        }
        this.group=group;
        this.type=type;
        this.inclusiveJunctionsP1=new HashSet<>();
        this.inclusiveJunctionsP2=new HashSet<>();
        this.exclusiveJunctionsP1=new HashSet<>();
        this.exclusiveJunctionsP2=new HashSet<>();
        this.inclusiveJunctionsP1.addAll(inclusiveJunctions);
        this.exclusiveJunctionsP1.addAll(exclusiveJunctions);
    }
    
    public HashSet<Junction> getAllInclusiveJunctions(){
        HashSet<Junction> incluJuncsReturn=new HashSet<>();
        incluJuncsReturn.addAll(inclusiveJunctionsP1);
        if(!inclusiveJunctionsP2.isEmpty()) incluJuncsReturn.addAll(inclusiveJunctionsP2);
        return(incluJuncsReturn);
    }
    
    public HashSet<Junction> getAllExclusiveJunctions(){
        HashSet<Junction> excluJuncsReturn=new HashSet<>();
        excluJuncsReturn.addAll(exclusiveJunctionsP1);
        if(!exclusiveJunctionsP2.isEmpty()) excluJuncsReturn.addAll(exclusiveJunctionsP2);
        return(excluJuncsReturn);
    }
    
    public HashSet<Junction> getLeftInclusiveJunctions(){
        HashSet<Junction> incluJuncsReturn=new HashSet<>();
        incluJuncsReturn.addAll(inclusiveJunctionsP1);
        return(incluJuncsReturn);
    }
    
    public HashSet<Junction> getRightInclusiveJunctions(){
        HashSet<Junction> incluJuncsReturn=new HashSet<>();
        incluJuncsReturn.addAll(inclusiveJunctionsP2);
        return(incluJuncsReturn);
    }
    
    public HashSet<Junction> getLeftExclusiveJunctions(){
        HashSet<Junction> excluJuncsReturn=new HashSet<>();
        excluJuncsReturn.addAll(exclusiveJunctionsP1);
        return(excluJuncsReturn);
    }
    
    public HashSet<Junction> getRightExclusiveJunctions(){
        HashSet<Junction> excluJuncsReturn=new HashSet<>();
        excluJuncsReturn.addAll(exclusiveJunctionsP2);
        return(excluJuncsReturn);
    }
    
    public int getType(){
        return(type);
    }
    
    public String getGroup(){
        return(group);
    }
    
    public boolean containsAll(Collection<Junction> juncs){
        boolean answer=true;
        for(Junction junc : juncs){
            if(!inclusiveJunctionsP1.contains(junc) && !exclusiveJunctionsP1.contains(junc) && !inclusiveJunctionsP2.contains(junc) && !exclusiveJunctionsP2.contains(junc)){
                answer=false;
                break;
            }
        }
        return(answer);
    }
    
    public ArrayList<Double> countInclusion(HashMap<Junction, Double> counts){
        ArrayList<Double> answer=new ArrayList<>();
        double inclu1=0; double inclu2=0; double exclu1=0; double exclu2=0;
        for(Junction junc : inclusiveJunctionsP1)
            if(counts.containsKey(junc))
                inclu1+=counts.get(junc);
        for(Junction junc : inclusiveJunctionsP2)
            if(counts.containsKey(junc))
                inclu2+=counts.get(junc);
        for(Junction junc : exclusiveJunctionsP1)
            if(counts.containsKey(junc))
                exclu1+=counts.get(junc);
        for(Junction junc : exclusiveJunctionsP2)
            if(counts.containsKey(junc))
                exclu2+=counts.get(junc);
        
        inclu1=(inclusiveJunctionsP2.isEmpty() || inclu1<inclu2) ? inclu1 : inclu2;
        exclu1=(exclusiveJunctionsP2.isEmpty() || exclu1<exclu2) ? exclu1 : exclu2;
        answer.add(inclu1); answer.add(exclu1);
        
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
    
    public int compareTo(ASEvent event2){
        if(!group.equals(event2.group))
            return(group.compareTo(event2.group));
        else if(type!=event2.type)
            return(Integer.valueOf(type).compareTo(Integer.valueOf(event2.type)));
        else{
            TreeSet<Junction> sortedJuncs=new TreeSet<>(new Comparator<Junction>(){
                @Override
                public int compare(Junction arg0, Junction arg1) {
                    return(arg0.compareTo(arg1));
                }
            });
            
            sortedJuncs.addAll(inclusiveJunctionsP1);
            Junction firstIncluThis=sortedJuncs.first();
            sortedJuncs.clear();
            sortedJuncs.addAll(event2.inclusiveJunctionsP1);
            Junction firstIncluThat=sortedJuncs.first();
            
            if(!firstIncluThis.equals(firstIncluThat))
                return(firstIncluThis.compareTo(firstIncluThat));
            
            sortedJuncs.clear();
            sortedJuncs.addAll(exclusiveJunctionsP1);
            Junction firstExcluThis=sortedJuncs.first();
            sortedJuncs.clear();
            sortedJuncs.addAll(event2.exclusiveJunctionsP1);
            Junction firstExcluThat=sortedJuncs.first();
            
            return(firstExcluThis.compareTo(firstExcluThat));
        }
    }
    
    @Override
    public String toString(){
        String[] eventTypes={"SKIPPED_EXON","MUTUAL_EXCLU_EXON","ALT_END"};
        String info=group+"\t"+eventTypes[type];

        TreeSet<Junction> sortedIncluJuncs=new TreeSet<>(new Comparator<Junction>(){
            @Override
            public int compare(Junction arg0, Junction arg1) {
                return(arg0.compareTo(arg1));
            }
        });
        sortedIncluJuncs.addAll(inclusiveJunctionsP1);
        String inclus="";
        for(Iterator<Junction> it=sortedIncluJuncs.iterator(); it.hasNext();){
            Junction junc=it.next();
            inclus+=junc.toString()+";";
        }
        sortedIncluJuncs.clear();
        sortedIncluJuncs.addAll(inclusiveJunctionsP2);
        if(!sortedIncluJuncs.isEmpty()){
            inclus+="\t";
            for(Iterator<Junction> it=sortedIncluJuncs.iterator(); it.hasNext();){
                Junction junc=it.next();
                inclus+=junc.toString()+";";
            }
        } else
            inclus+="\tNA";

        TreeSet<Junction> sortedExcluJuncs=new TreeSet<>(new Comparator<Junction>(){
        @Override
        public int compare(Junction arg0, Junction arg1){
            return(arg0.compareTo(arg1));
        }
        });
        sortedExcluJuncs.addAll(exclusiveJunctionsP1);
        String exclus="";
        for(Iterator<Junction> it=sortedExcluJuncs.iterator(); it.hasNext();){
            Junction junc=it.next();
            exclus+=junc.toString()+";";
        }
        sortedExcluJuncs.clear();
        sortedExcluJuncs.addAll(exclusiveJunctionsP2);
        if(!sortedExcluJuncs.isEmpty()){
            exclus+="\t";
            for(Iterator<Junction> it=sortedExcluJuncs.iterator(); it.hasNext();){
                Junction junc=it.next();
                exclus+=junc.toString()+";";
            }
        } else
            exclus+="\tNA";

        info+="\t"+inclus+"\t"+exclus;
        return(info);
    }
}
