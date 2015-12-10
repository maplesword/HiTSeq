/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hitseq.gui;

import hitseq.MappingProcessor;
import hitseq.ReadCounter;
import hitseq.annotation.Annotation;
import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.regex.*;

/**
 *
 * @author Chih-sung
 */
public class CommandRunner extends Thread{
    final private String cmd; // Command to run, including count, uniq
    final private ParameterSet parameters; // Parameters
    
    private HashMap<String, Double> totalNumMappedReads;
    private HashMap<String, HashMap<String, Double>> readCount;
    
    final private javax.swing.JProgressBar progressBar;
    final private javax.swing.JButton[] buttons;
    final private javax.swing.JButton buttonCancel;
    final private javax.swing.JTabbedPane pane;
    
    private Annotation annotation;
    
    CommandRunner(){
        cmd=null;
        parameters=null;
        progressBar=null;
        buttons=null;
        pane=null;
        totalNumMappedReads=null;
        readCount=null;
        buttonCancel=null;
        annotation=null;
    }
    CommandRunner(String cmd, ParameterSet parameters){
        this.cmd=cmd;
        this.parameters=parameters;
        progressBar=null;
        buttons=null;
        pane=null;
        totalNumMappedReads=null;
        readCount=null;
        buttonCancel=null;
        annotation=null;
    }
    CommandRunner(String cmd, ParameterSet parameters, javax.swing.JProgressBar progressBar, javax.swing.JButton[] buttons, javax.swing.JButton cancel, javax.swing.JTabbedPane pane){
        this.cmd=cmd;
        this.parameters=parameters;
        this.progressBar=progressBar;
        this.buttons=buttons;
        this.pane=pane;
        totalNumMappedReads=null;
        readCount=null;
        buttonCancel=cancel;
        annotation=null;
    }
    
    HashMap<String, Double> getTotalNumMappedReads(){
        return totalNumMappedReads;
    }
    
    HashMap<String, HashMap<String, Double>> getReadCount(){
        return readCount;
    }
    
    @Override
    public void run(){
        switch (cmd) {
            case "count":
                runCount();
                break;
            case "uniq":
                runUniq();
                break;
        }
        
        if(progressBar!=null)
            progressBar.setString("Finish");
        if(buttons!=null)
            for(javax.swing.JButton button : buttons)
                button.setEnabled(true);
        if(buttonCancel!=null)
            buttonCancel.setText("Close");
        if(pane!=null){
            pane.setEnabled(true);
            pane.setFocusable(true);
        }
    }
    
    void runCount(){
        int strandSpecific = parameters.getStrandedness();
        boolean considerNH = parameters.getConsiderNHTag();
        boolean onlyUnique = parameters.getOnlyUnique();
        boolean readCollapse = parameters.getReadCollapseTag();
        
        int numMappingFiles=parameters.getMappingFiles().size();
        int stepLength=100/(numMappingFiles+1);
        
        // read annotation and initial data vectors
        annotation = new Annotation(parameters.getAnnotFile(), parameters.getAnnotFormat());
        if (!this.isInterrupted()) {
            if (progressBar != null)
                progressBar.setValue(progressBar.getValue() + stepLength / 2);
 
            annotation.estimateAmbiguousGeneRegions();
            totalNumMappedReads = new HashMap<>();
            readCount = new HashMap<>();
            for (String gene : annotation.getGeneSet()) {
                readCount.put(gene, new HashMap<String, Double>());
            }
            
            if (progressBar != null) 
                progressBar.setValue(progressBar.getValue() + stepLength / 2);
        }
        
        // read counting
        for(File mappingFile : parameters.getMappingFiles()){
            if(this.isInterrupted())
                break;
            annotation.resetPointer();
            ReadCounter counter = new ReadCounter(mappingFile, annotation, strandSpecific, 0, false);
            counter.estimateCounts(considerNH, onlyUnique, readCollapse, stepLength);
            HashMap<String, Double> count = counter.getCounts();
            
            for (String gene : count.keySet()) {
                readCount.get(gene).put(mappingFile.getAbsolutePath(), count.get(gene));
            }
            totalNumMappedReads.put(mappingFile.getAbsolutePath(), counter.getTotalNumReads());
            
            if(progressBar!=null)
                progressBar.setValue(progressBar.getValue()+stepLength);
        }
    }
    
    void runUniq(){
        boolean considerNH = parameters.getConsiderNHTag();
        File inputSam = parameters.getMappingFiles().get(0);
        File outputSam = parameters.getMappingFiles().get(1);
        MappingProcessor processor = new MappingProcessor(inputSam);
        int numUniqueReads = processor.extractUniquelyMappedReads(outputSam, considerNH, false);
        System.err.println("Total Number of Uniquely Mapped Reads of " + inputSam.getAbsolutePath() + ": " + numUniqueReads);
    }
    
    String generateOutputString(){
        // output header
        String output = "GENE_ID\tLENGTH";
        for (File mappingFile : parameters.getMappingFiles()) {
            output = output + "\t" + mappingFile.getAbsolutePath();
        }
        System.out.println(output);
        
        // output the total number of mapped reads
        String totalReads = String.valueOf(totalNumMappedReads.get(parameters.getMappingFiles().get(0).getAbsolutePath()).intValue());
        if (totalNumMappedReads.size() > 1) {
            for (int i = 1; i < totalNumMappedReads.size(); i++) {
                totalReads = totalReads + "\t" + String.valueOf(totalNumMappedReads.get(parameters.getMappingFiles().get(i).getAbsolutePath()).intValue());
            }
        }
        output+="TOTAL_READS\tNA\t" + totalReads + "\n";
        System.out.println("TOTAL_READS\tNA\t" + totalReads);
        
        // output read count for genes
        java.util.TreeSet<String> sortedGeneNames = new java.util.TreeSet<>(new java.util.Comparator<String>() {
            @Override
            public int compare(String o1, String o2) {
                return o1.compareTo(o2);
            }
        });

        for (String gene : readCount.keySet()) {
            sortedGeneNames.add(gene);
        }
        for (String gene : sortedGeneNames) {
            int geneLength;
            if (parameters.getStrandedness() == 0) // no strand information
            {
                geneLength = annotation.getExclusiveGeneLengthNoStrand(gene);
            } else {
                geneLength = annotation.getExclusiveGeneLength(gene);
            }

            String readNum = "";
            for (int i = 0; i < parameters.getMappingFiles().size(); i++) {
                readNum = readNum + "\t" + String.valueOf(readCount.get(gene).get(parameters.getMappingFiles().get(i).getAbsolutePath()).intValue());
            }
            //output += gene + "\t" + geneLength + readNum + "\n";
            System.out.println(gene + "\t" + geneLength + readNum);
        }
        
        return output;
    }
}

class ParameterSet {
    private int strandSpecific;
    private boolean considerNH;
    private boolean onlyUnique;
    private boolean readCollapse;
    
    private String annotFormat;
    private File annotFile;
    final private ArrayList<File> mappingFiles;
    
    ParameterSet(){
        strandSpecific=0;
        considerNH=false;
        onlyUnique=false;
        readCollapse=false;
        
        annotFormat=null;
        annotFile=null;
        mappingFiles=new ArrayList<>();
    }
    
    void setAnnotFormat(String format){
        if(format.toLowerCase().equals("bed") || format.toLowerCase().equals("gtf") || format.toLowerCase().equals("struc") || format.toLowerCase().equals("juncs"))
            annotFormat=format;
    }
    
    void setAnnotFile(File file){
        String path=file.getAbsolutePath().toLowerCase();
        if(path.endsWith(".bed") || path.endsWith(".gtf") || path.endsWith(".struc") || path.endsWith(".juncs")){
            annotFile=file;
            Pattern pattern=Pattern.compile("(bed|gtf|struc|juncs)$");
            Matcher matcher=pattern.matcher(path);
            matcher.find();
            setAnnotFormat(matcher.group());
        }
    }
    
    void setConsiderNHTag(boolean consider){
        considerNH=consider;
    }
    
    void setOnlyUnique(boolean consider){
        onlyUnique=consider;
    }
    
    void setReadCollapse(boolean consider){
        readCollapse=consider;
    }
    
    boolean setStrandness(int s){
        if(s==0 || s==-1 || s==1){
            strandSpecific=s;
            return(true);
        } else
            return(false);
    }
    
    void addMappingFile(File file){
        mappingFiles.add(file);
    }
    
    void addMappingFiles(Collection<File> collection){
        mappingFiles.addAll(collection);
    }
    
    String getAnnotFormat(){
        return(annotFormat);
    }
    
    File getAnnotFile(){
        return(annotFile);
    }
    
    ArrayList<File> getMappingFiles(){
        return(mappingFiles);
    }
    
    int getStrandedness(){
        return (strandSpecific);
    }

    boolean getConsiderNHTag() {
        return (considerNH);
    }
    
    boolean getOnlyUnique(){
        return (onlyUnique);
    }

    boolean getReadCollapseTag() {
        return (readCollapse);
    }
}