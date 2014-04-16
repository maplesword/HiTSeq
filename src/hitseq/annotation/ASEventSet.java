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
            // if there is only one junction in the group, skip it
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
            ArrayList<Junction> sortedJunctionsGroup=new ArrayList<>(); // the sorted list of junctions in this group
            ArrayList<Integer> sortedJunctionCoord=new ArrayList<>(); // the sorted sites of junction site (exon boundary); the overlapping sites will be collapsed
            for(Iterator<Junction> it=sortedJunctionsGroupTreeSet.iterator(); it.hasNext();){
                Junction junc=it.next();
                sortedJunctionsGroup.add(junc);
                
                if(sortedJunctionCoord.isEmpty()){
                    sortedJunctionCoord.add(junc.getStartSite());
                    sortedJunctionCoord.add(junc.getEndSite());
                } else{
                    int i=sortedJunctionCoord.size()-1;
                    while(i>0 && sortedJunctionCoord.get(i)>junc.getEndSite())
                        i--;
                    if(!sortedJunctionCoord.get(i).equals(junc.getEndSite()))
                        sortedJunctionCoord.add(i+1, junc.getEndSite());
                    while(i>0 && sortedJunctionCoord.get(i)>junc.getStartSite())
                        i--;
                    if(!sortedJunctionCoord.get(i).equals(junc.getStartSite()))
                        sortedJunctionCoord.add(i+1, junc.getStartSite());
                }
            }
            
            // if there is no touched junction in the group, skip it
            if(sortedJunctionCoord.size()==sortedJunctionsGroup.size()*2) // if there is touched junction, because of the collapse of junction sites, the number of coordinates should be smaller than 2 folds of the number of junctions
                continue;
            
            ArrayList<ArrayList<ASEvent>> eventsInTypeThisGroup=new ArrayList<>();
            // Look for skipped exon events
            eventsInType.add(new ArrayList<ASEvent>());
            eventsInTypeThisGroup.add(new ArrayList<ASEvent>());
            if(sortedJunctionsGroup.size()>=3){
                for(int i=1; i<sortedJunctionsGroup.size()-1; i++){
                    /*
                     * assume this junction is the exclusive junction of a skipped exon event,
                     * try to identify its corresponding inclusive junction pairs.
                     * if succeed, add a skipped exon event
                     */
                    Junction exclu=sortedJunctionsGroup.get(i); 
                    ArrayList<Junction> excluList=new ArrayList<>();
                    excluList.add(exclu);
                    ArrayList<Junction> incluLeft=new ArrayList<>();
                    ArrayList<Junction> incluRight=new ArrayList<>();

                    for(int j=i-1; j>=0; j--){
                        Junction left=sortedJunctionsGroup.get(j);
                        if(exclu.withStartSite(left.getStartSite()))
                            incluLeft.add(left);
                    }
                    if(incluLeft.isEmpty()) continue;
                    for(int j=i+1; j<sortedJunctionsGroup.size(); j++){
                        Junction right=sortedJunctionsGroup.get(j);
                        if(exclu.withEndSite(right.getEndSite()))
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
                ArrayList<ArrayList<Junction>> eventCandidates=new ArrayList<>();
                for(int i=0; i<sortedJunctionsGroup.size(); i++){
                    Junction thisJunc=sortedJunctionsGroup.get(i);
                    if(i>0){
                        if(!eventCandidates.isEmpty()){
                            ArrayList<ArrayList<Junction>> existedCandidates=new ArrayList<>();
                            existedCandidates.addAll(eventCandidates);
                            
                            for(ArrayList<Junction> candidate : existedCandidates){
                                if(candidate.size()==2 && candidate.get(1).cross(thisJunc)){
                                    /* 
                                     * only if the end site of the first junction is the closest boundary of the start site of this junction,
                                     * this junction will be added to possible candidate set
                                     */
                                    if(sortedJunctionCoord.indexOf(candidate.get(0).getEndSite()) == sortedJunctionCoord.indexOf(thisJunc.getStartSite()-1)){
                                        ArrayList<Junction> newCandidate=new ArrayList<>();
                                        newCandidate.addAll(candidate);
                                        newCandidate.add(thisJunc);
                                        eventCandidates.add(newCandidate);
                                    }
                                } else if(candidate.size()==3 && candidate.get(2).getEndSite()==thisJunc.getEndSite()){
                                    /* 
                                     * only if the end site of the second junction is the closest boundary of the start site of this junction,
                                     * this junction will be combined with the other three junctions to form a MXE event
                                     */
                                    if(sortedJunctionCoord.indexOf(candidate.get(1).getEndSite()) == sortedJunctionCoord.indexOf(thisJunc.getStartSite()-1)){
                                        ArrayList<Junction> newEventIncluP1=new ArrayList<>();
                                        ArrayList<Junction> newEventIncluP2=new ArrayList<>();
                                        ArrayList<Junction> newEventExcluP1=new ArrayList<>();
                                        ArrayList<Junction> newEventExcluP2=new ArrayList<>();
                                        newEventIncluP1.add(candidate.get(0));
                                        newEventIncluP2.add(candidate.get(2));
                                        newEventExcluP1.add(candidate.get(1));
                                        newEventExcluP2.add(thisJunc);
                                        ASEvent newEvent=new ASEvent(ASEvent.MUTUAL_EXCL,newEventIncluP1,newEventIncluP2,newEventExcluP1,newEventExcluP2);
                                        events.add(newEvent);
                                        eventsInType.get(ASEvent.MUTUAL_EXCL).add(newEvent);
                                        eventsInTypeThisGroup.get(ASEvent.MUTUAL_EXCL).add(newEvent);
                                    }
                                }
                            }
                        }
                        
                        // firstly, if there are two touched junctions, add them as the initial candidate
                        int j=i-1;
                        while(j>=0 && sortedJunctionsGroup.get(j).getStartSite()==thisJunc.getStartSite()){ 
                            ArrayList<Junction> candidate=new ArrayList<>();
                            candidate.add(sortedJunctionsGroup.get(j));
                            candidate.add(thisJunc);
                            eventCandidates.add(candidate);
                            j--;
                        }
                    }
                }
            }
            
            // Look for alternative 5'/3' end
            eventsInType.add(new ArrayList<ASEvent>());
            eventsInTypeThisGroup.add(new ArrayList<ASEvent>());
            if(juncGroups.get(group).size()>1 && sortedJunctionsGroup.size()*2!=sortedJunctionCoord.size()){
                // first check start-shared junctions
                int i=0;
                ArrayList<Junction> juncShareStart=new ArrayList<>();
                while(i<sortedJunctionsGroup.size()){
                    if(juncShareStart.isEmpty()){
                        juncShareStart.add(sortedJunctionsGroup.get(i));
                        i++;
                        continue;
                    }
                    else if(juncShareStart.get(0).withStartSite(sortedJunctionsGroup.get(i).getStartSite())){
                        juncShareStart.add(sortedJunctionsGroup.get(i));
                        i++;
                        if(i!=sortedJunctionsGroup.size())
                            continue;
                        else
                            i--;
                    }
                    
                    if(juncShareStart.size()>1)
                        for(int j1=0; j1<juncShareStart.size()-1; j1++)
                            for(int j2=j1+1; j2<juncShareStart.size(); j2++){
                                ArrayList<Junction> candidates=new ArrayList<>();
                                candidates.add(juncShareStart.get(j1));candidates.add(juncShareStart.get(j2));

                                // check whether these two junctions are already involved in other SE/MXE events
                                boolean involved=false;
                                for(ASEvent event : eventsInTypeThisGroup.get(ASEvent.SKIPP_EXON))
                                    if(event.containAll(candidates)){
                                        involved=true;
                                        break;
                                    }
                                if(!involved)
                                    for(ASEvent event : eventsInTypeThisGroup.get(ASEvent.MUTUAL_EXCL))
                                        if(event.containAll(candidates)){
                                            involved=true;
                                            break;
                                        }

                                if(!involved){ // if they are not involved in other events, add them as an alternative 5'/3' end events
                                    ArrayList<Junction> newEventInclu=new ArrayList<>();
                                    ArrayList<Junction> newEventExclu=new ArrayList<>();
                                    newEventInclu.add(candidates.get(0));
                                    newEventExclu.add(candidates.get(1));
                                    ASEvent newEvent=new ASEvent(ASEvent.ALT_END, newEventInclu, newEventExclu);
                                    events.add(newEvent);
                                    eventsInType.get(ASEvent.ALT_END).add(newEvent);
                                    eventsInTypeThisGroup.get(ASEvent.ALT_END).add(newEvent);
                                }
                            }

                    juncShareStart.clear();
                    i++;
                }
                
                // second check end-shared junctions
                HashMap<Integer, ArrayList<Junction>> endOfJuncs=new HashMap<>();
                for(i=0; i<sortedJunctionsGroup.size(); i++){
                    if(!endOfJuncs.containsKey(sortedJunctionsGroup.get(i).getEndSite()))
                        endOfJuncs.put(sortedJunctionsGroup.get(i).getEndSite(), new ArrayList<Junction>());
                    endOfJuncs.get(sortedJunctionsGroup.get(i).getEndSite()).add(sortedJunctionsGroup.get(i));
                }
                for(ArrayList<Junction> candidates : endOfJuncs.values()){
                    boolean involved=false;
                    for(ASEvent event : eventsInTypeThisGroup.get(ASEvent.SKIPP_EXON))
                        if(event.containAll(candidates)){
                            involved=true;
                            break;
                        }
                    if(!involved)
                        for(ASEvent event : eventsInTypeThisGroup.get(ASEvent.MUTUAL_EXCL))
                            if(event.containAll(candidates)){
                                involved=true;
                                break;
                            }

                    if(!involved){ // if they are not involved in other events, add them as an alternative 5'/3' end events
                        ArrayList<Junction> newEventInclu=new ArrayList<>();
                        ArrayList<Junction> newEventExclu=new ArrayList<>();
                        newEventInclu.add(candidates.get(0));
                        newEventExclu.add(candidates.get(1));
                        ASEvent newEvent=new ASEvent(ASEvent.ALT_END, newEventInclu, newEventExclu);
                        events.add(newEvent);
                        eventsInType.get(ASEvent.ALT_END).add(newEvent);
                        eventsInTypeThisGroup.get(ASEvent.ALT_END).add(newEvent);
                    }
                }
            }
        }
    }
    
    public void outputASEventSet(){
        String[] typeName={"SKIPPED_EXON","MUTUAL_EXCLU_EXON","ALT_END"};
        for(int i=0; i<eventsInType.size(); i++)
            for(ASEvent event : eventsInType.get(i)){
                HashSet<Junction> incluJuncs=event.getInclusiveJunctions();
                String inclus="";
                for(Junction junc : incluJuncs)
                    inclus+=junc.toString()+";";
                
                HashSet<Junction> excluJuncs=event.getExclusiveJunctions();
                String exclus="";
                for(Junction junc : excluJuncs)
                    exclus+=junc.toString()+";";
                
                System.out.println(typeName+"\t"+inclus+"\t"+exclus);
            }
    }
}
