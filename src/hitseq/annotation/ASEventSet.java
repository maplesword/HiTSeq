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
public class ASEventSet {
    private HashSet<ASEvent> events;
    private ArrayList<ArrayList<ASEvent>> eventsInType;
    
    public ASEventSet(JunctionSet junctions){
        if(junctions.getJunctionGroups()==null)
            junctions.groupJuncSet();
        HashMap<String, HashSet<Junction>> juncGroups=junctions.getJunctionGroups();
        
        events=new HashSet<>();
        eventsInType=new ArrayList<>();
        
        for(String group : juncGroups.keySet()){
            if(juncGroups.get(group).size()<2)
                continue;
            
            // sort junctions in the group
            TreeSet<Junction> sortedJunctionsGroupTreeSet=new TreeSet<>(new Comparator<Junction>(){
                @Override
                public int compare(Junction arg0, Junction arg1) {
                    return(arg0.compareTo(arg1));
                }
            });
            sortedJunctionsGroupTreeSet.addAll(juncGroups.get(group));
            ArrayList<Junction> sortedJunctionsGroup=new ArrayList<>();
            for(Iterator<Junction> it=sortedJunctionsGroupTreeSet.iterator(); it.hasNext();){
                Junction junc=it.next();
                sortedJunctionsGroup.add(junc);
            }
            
            ArrayList<ArrayList<ASEvent>> eventsInTypeThisGroup=new ArrayList<>();
            // Look for skipped exon events
            eventsInType.add(new ArrayList<ASEvent>());
            eventsInTypeThisGroup.add(new ArrayList<ASEvent>());
            if(sortedJunctionsGroup.size()>=3){
                for(int i=1; i<sortedJunctionsGroup.size()-1; i++){
                    Junction exclu=sortedJunctionsGroup.get(i);
                    ArrayList<Junction> excluList=new ArrayList<>();
                    excluList.add(exclu);
                    ArrayList<Junction> incluLeft=new ArrayList<>();
                    ArrayList<Junction> incluRight=new ArrayList<>();

                    for(int j=i-1; j>=0; j--){
                        Junction left=sortedJunctionsGroup.get(j);
                        if(exclu.withDonorSite(left.getDonorSite()))
                            incluLeft.add(left);
                    }
                    if(incluLeft.isEmpty()) continue;
                    for(int j=i+1; j<sortedJunctionsGroup.size(); j++){
                        Junction right=sortedJunctionsGroup.get(j);
                        if(exclu.withAcceptorSite(right.getAcceptorSite()))
                            incluRight.add(right);
                    }
                    if(incluRight.isEmpty()) continue;
                    
                    ASEvent event=new ASEvent(ASEvent.SKIPP_EXON, incluLeft, incluRight, excluList);
                    events.add(event);
                    eventsInType.get(ASEvent.SKIPP_EXON).add(event);
                    eventsInTypeThisGroup.get(ASEvent.SKIPP_EXON).add(event);
                }
            }
            
            // Look for mutual exclusive exons event
            eventsInType.add(new ArrayList<ASEvent>());
            eventsInTypeThisGroup.add(new ArrayList<ASEvent>());
            if(sortedJunctionsGroup.size()>=4){
                
            }
            
            // Look for alternative 5' end
            
            // Look for alternative 3' end
            
        }
    }
}
