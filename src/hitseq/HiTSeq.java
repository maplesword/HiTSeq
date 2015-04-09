/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hitseq;

import hitseq.annotation.*;
import htsjdk.samtools.SAMFileHeader;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

/**
 *
 * @author hezhisong
 */
public class HiTSeq {
    /**
     * Show the help page
     */
    private static void showHelp(){
        System.err.println("\nThis is the help page of HiTSeq.\n"
                    + "Usage: java -jar HiTSeq.jar <command> <...>\n"
                    + "   Or: HiTSeq.sh <command> <...>\n\n"
                    + "Commands:  help       This help page\n"
                    + "           info       Output information of mapping, including total number of reads, mapped reads, etc\n"
                    + "           count      Given annotation in struc/gtf/bed format, output table of read count of each gene for the input alignments\n"
                    + "           rpkm       Given annotation in struc/gtf/bed format, output table of RPKM of each gene for the input alignment\n"
                    + "           uniq       Extract uniquely mapped reads from the input alignment\n"
                    + "           gtf2struc  Transform GTF annotation into struc format\n"
                    + "           tojuncs    Combine and transform junction list(s) into a junction list in 'juncs' format\n"
                    + "           toevents   Combine the junction list(s) and generate alternative splicing events\n"
                    + "           countjunc  Given junction list in junc/bed/gtf format, or event list in events format, output table of read count of each junction for the input alignments\n"
                    + "           gui        Open HiTSeq GUI.\n");
    }
    
    /**
     * Run the program of "test" command
     * @param args the command line arguments
     */
    private static void runTest(String[] args){
        htsjdk.samtools.SamReader reader = htsjdk.samtools.SamReaderFactory.makeDefault().open(new File(args[1]));
        reader.getFileHeader().setSortOrder(SAMFileHeader.SortOrder.queryname);
        htsjdk.samtools.SAMFileWriter writer = new htsjdk.samtools.SAMFileWriterFactory().makeBAMWriter(reader.getFileHeader(), false, new File(args[2]));
        for (htsjdk.samtools.SAMRecord record : reader) {
            //writer.addAlignment(record);
            System.out.print(record.getSAMString());
        }
        htsjdk.samtools.util.CloserUtil.close(reader);
        htsjdk.samtools.util.CloserUtil.close(writer);
    }
    
    /**
     * Run the program of "info" command
     * @param args the command line arguments
     */
    private static void runInfo(String[] args){
        String pathMapping=args[1];
        File mappingFile = new File(pathMapping);
        MappingProcessor processor = new MappingProcessor(mappingFile);
        processor.collectMappingInformation();
        ArrayList<Integer> information = processor.getMappingInfo();

        System.out.println("\n" + args[1] + ":");
        System.out.printf("%45s|          %d\n", "Number of reads: ", information.get(0));
        System.out.printf("%45s|          %d\n", "Number of mapped reads: ", information.get(1));
        System.out.printf("%45s|          %d\n", "Number of uniquely mapped reads: ", information.get(2));
        System.out.println();
        if (information.size() > 3) {
            for (int i = 3; i < information.size(); i += 2) {
                System.out.printf("%45s|          %d\n", "Number of reads with " + (i - 3) / 2 + " mismatch(es)", information.get(i));
                System.out.printf("%46s          %d\n", "(UNIQUE READS)", information.get(i + 1));
            }
        }
        System.out.println();
    }
    
    /**
     * The program to run "uniq" command
     * @param args the command line arguments
     */
    private static void runUniq(String[] args){
        File inputSam = new File(args[1]);
        MappingProcessor processor = new MappingProcessor(inputSam);
        File outputSam;
        if (args.length > 2) {
            outputSam = new File(args[2]);
        } else {
            String parentDir = inputSam.getParent();
            String output = parentDir + "/unique." + inputSam.getName();
            outputSam = new File(output);
        }

        int numUniqueReads = processor.extractUniquelyMappedReads(outputSam, true);
        System.out.println("Total Number of Uniquely Mapped Reads of " + args[1] + ": " + numUniqueReads);
    }
    
    /**
     * The program to run "gtf2struc" command
     * @param args the command line arguments
     */
    private static void transformGTF2Struc(String[] args){
        String pathAnnotation = args[1];
        Annotation annotation = new Annotation(new File(pathAnnotation), "gtf");
        annotation.outputInStruc();
    }
    
    private static void transform4Junction(String[] args){
        String cmd=args[0];
        
        if (args.length == 1 || args[1].equals("-h")) {
            System.err.println("\nThis is the help of '" + cmd + "' command of HiTSeq.");
            System.err.println("Usage: java -jar HiTSeq.jar " + cmd + " [options] <in.file> [in2.file ...]\n"
                    + "   Or: HiTSeq.sh " + cmd + " [options] <in.file> [in2.file ...]");
            System.err.println("\n"
                    + "Options: -h        This help page\n"
                    + "         -a [str]  The file type of junction annotation file (default: juncs format; options: juncs/bed[tophat output]/gtf/bam/events)\n"
                    + "         -s [int]  Strandedness of BAM/SAM input (default: 0 [no strand information]; 1/-1 [same/opposite strandness]; only work with '-a bam')\n");
            System.exit(0);
        }

        // read the parameters
        ParameterSet parameters = new ParameterSet(cmd);
        parameters.readCommandLineArgs(args);

        String fileType = parameters.getAnnotFormat();
        int strandSpecific = parameters.getStrandedness();
        int firstInputIndex = parameters.getFirstSAMIdx();

        // read input
        JunctionSet junctions = new JunctionSet();
        for (int i = firstInputIndex; i < args.length; i++) {
            String pathBed = args[i];
            if (fileType.equalsIgnoreCase("bam")) {
                junctions.addJunctionSet(new File(pathBed), strandSpecific);
            } else if(fileType.equalsIgnoreCase("events")){
                ASEventSet eventSet=new ASEventSet(new File(pathBed));
                junctions.addJunctionSet(new JunctionSet(eventSet));
            } else {
                junctions.addJunctionSet(new File(pathBed), fileType);
            }
        }

        if (cmd.equalsIgnoreCase("tojuncs")) {
            junctions.outputInJuncs(true);
        } else {
            junctions.groupJuncSet();
            System.err.println("done junction grouping");
            //junctions.outputJuncGroups();
            ASEventSet juncEvents = new ASEventSet(junctions);
            System.err.println("done event generation");
            juncEvents.outputASEventSet();
        }
    }
    
    /**
     * The program to run "count" and "rpkm" command
     * @param args the command line arguments
     */
    private static void runReadCounting(String[] args){
        String cmd=args[0];
        if (args.length == 1 || args[1].equals("-h")) {
            System.err.println("\nThis is the help of '" + cmd.toLowerCase() + "' command of HiTSeq.");
            System.err.println("Usage: java -jar HiTSeq.jar " + cmd.toLowerCase() + " [options] <annotation.struc> <in.bam> [in2.bam ...]\n"
                    + "   Or: HiTSeq.sh " + cmd.toLowerCase() + " [options] <annotation.struc> <in.bam> [in2.bam ...]");
            System.err.println("\n"
                    + "Options: -h        This help page\n"
                    + "         -s [int]  Strandedness (default: 0 - no strand information; 1 - same strandness; -1 - opposite strandness)\n"
                    + "         -n        For reads mapped to n-loci, assign 1/n read to each hit\n"
                    + "         -c        Do read collapse to remove PCR duplicates\n"
                    + "         -m [int]  The mode to deal with multi-gene hits (default: mode 0 - abandon ambiguous reads; options: 0-3)\n"
                    + "         -t [int]  The maximum iteration time to assign ambiguous reads (default: 2). Only work with -m 3\n"
                    + "         -a [str]  The file type of annotation file (default: struc format; options: struc/gtf/bed)\n");
            System.exit(0);
        }

        // read the parameters
        ParameterSet parameters = new ParameterSet(cmd);
        parameters.readCommandLineArgs(args);

        int strandSpecific = parameters.getStrandedness();
        boolean considerNH = parameters.getConsiderNHTag();
        boolean readCollapse = parameters.getReadCollapseTag();
        int modeForMultiGenesOverlap = parameters.getModeForMultiGenesOverlap();
        int iterationLimit = parameters.getIterationLimit();
        String annotFormat = parameters.getAnnotFormat();

        int firstSAMIndex = parameters.getFirstSAMIdx();

        // read annotation
        String pathAnnotation = args[firstSAMIndex];
        Annotation annotation = new Annotation(new File(pathAnnotation), annotFormat);
        if (modeForMultiGenesOverlap == 0) {
            annotation.estimateAmbiguousGeneRegions();
        }
        HashMap<String, Double> totalNumMappedReads = new HashMap<>();
        HashMap<String, HashMap<String, Double>> readCount = new HashMap<>();
        HashMap<String, HashMap<String, Double>> fpkm = new HashMap<>();
        for (String gene : annotation.getGeneSet()) {
            readCount.put(gene, new HashMap<String, Double>());
            fpkm.put(gene, new HashMap<String, Double>());
        }
        firstSAMIndex++;

        for (int i = firstSAMIndex; i < args.length; i++) {
            String pathMapping = args[i];
            File mappingFile = new File(pathMapping);

            // Read counting
            annotation.resetPointer();
            ReadCounter counter = new ReadCounter(mappingFile, annotation, strandSpecific, modeForMultiGenesOverlap);
            counter.estimateCounts(considerNH, readCollapse, iterationLimit);
            HashMap<String, Double> count = counter.getCounts();
            for (String gene : count.keySet()) {
                readCount.get(gene).put(args[i], count.get(gene));
            }
            totalNumMappedReads.put(args[i], counter.getTotalNumReads());

            // Calculate RPKM if necessary
            if (cmd.equalsIgnoreCase("rpkm")) {
                counter.estimateRPKM();
                HashMap<String, Double> fpkmGene = counter.getRPKM();
                for (String gene : fpkmGene.keySet()) {
                    fpkm.get(gene).put(args[i], fpkmGene.get(gene));
                }
            }

            System.err.println("done " + args[i] + "\n");
        }

        String header = cmd.equalsIgnoreCase("count") ? "GENE_ID\tLENGTH" : "GENE_ID";
        for (int i = firstSAMIndex; i < args.length; i++) {
            header = header + "\t" + args[i];
        }
        System.out.println(header);

        java.util.TreeSet<String> sortedGeneNames = new java.util.TreeSet<>(new java.util.Comparator<String>() {

            @Override
            public int compare(String o1, String o2) {
                return o1.compareTo(o2);
            }
        });

        if (cmd.equalsIgnoreCase("count")) {
            String totalReads = String.valueOf(totalNumMappedReads.get(args[firstSAMIndex]).intValue());
            if (args.length > firstSAMIndex + 1) {
                for (int i = firstSAMIndex + 1; i < args.length; i++) {
                    totalReads = totalReads + "\t" + String.valueOf(totalNumMappedReads.get(args[i]).intValue());
                }
            }
            System.out.println("TOTAL_READS\tNA\t" + totalReads);


            for (String gene : readCount.keySet()) {
                sortedGeneNames.add(gene);
            }
            for (Iterator<String> it = sortedGeneNames.iterator(); it.hasNext();) {
                String gene = it.next();

                int geneLength;
                if (modeForMultiGenesOverlap == 0) { // abandom ambiguous reads, use exclusive length
                    if (strandSpecific == 0) // no strand information
                    {
                        geneLength = annotation.getExclusiveGeneLengthNoStrand(gene);
                    } else {
                        geneLength = annotation.getExclusiveGeneLength(gene);
                    }
                } else // not to abandom ambiguous reads, use total length
                {
                    geneLength = annotation.getGeneLength(gene);
                }

                String readNum = "";
                for (int i = firstSAMIndex; i < args.length; i++) {
                    if (modeForMultiGenesOverlap == 1) {
                        readNum = readNum + "\t" + String.valueOf(readCount.get(gene).get(args[i]));
                    } else {
                        readNum = readNum + "\t" + String.valueOf(readCount.get(gene).get(args[i]).intValue());
                    }
                }
                System.out.println(gene + "\t" + geneLength + readNum);
            }
        } else if (cmd.equalsIgnoreCase("rpkm")) {
            for (String gene : fpkm.keySet()) {
                sortedGeneNames.add(gene);
            }
            for (Iterator<String> it = sortedGeneNames.iterator(); it.hasNext();) {
                String gene = it.next();
                String fpkmGene = "";
                for (int i = firstSAMIndex; i < args.length; i++) {
                    fpkmGene = fpkmGene + "\t" + String.valueOf(fpkm.get(gene).get(args[i]));
                }
                System.out.println(gene + fpkmGene);
            }
        }
    }
    
    /**
     * The function to run read counting for junctions (command "countjunc")
     * @param args the command line arguments
     */
    private static void runJunctionCounting(String[] args){
        String cmd=args[0];
        
        if (args.length == 1 || args[1].equals("-h")) {
            System.err.println("\nThis is the help of '" + cmd + "' command of HiTSeq.");
            System.err.println("Usage: java -jar HiTSeq.jar " + cmd + " [options] <junc.file> <in.bam> [in2.bam ...]\n"
                    + "   Or: HiTSeq.sh " + cmd + " [options] <junc.file> <in.bam> [in2.bam ...]");
            System.err.println("\n"
                    + "Options: -h        This help page\n"
                    + "         -a [str]  The file type of junction annotation file (default: juncs format; options: juncs/bed[tophat output]/gtf/bam/events)\n"
                    + "         -n        For reads mapped to n-loci, only assign 1/n read to each hit\n"
                    + "         -c        Do read collapse to remove PCR duplicates\n"
                    + "         -s [int]  Strandedness of BAM/SAM input (default: 0 [no strand information]; 1/-1 [same/opposite strandness]; only work with '-a bam')\n"
                    + "         -e        Output the counting for events instead of junctions\n");
            System.exit(0);
        }

        // read the parameters
        ParameterSet parameters = new ParameterSet(cmd);
        parameters.readCommandLineArgs(args);

        int strandSpecific = parameters.getStrandedness();
        boolean considerNH = parameters.getConsiderNHTag();
        boolean readCollapse = parameters.getReadCollapseTag();
        String juncType = parameters.getAnnotFormat();
        boolean outputForEvents = parameters.getOutputForEventsTag();

        int firstSAMIndex = parameters.getFirstSAMIdx();

        // generate junction set and AS event set
        String pathJunctions = args[firstSAMIndex];
        JunctionSet junctions;
        ASEventSet events = null;
        if (!juncType.equals("event")) {
            junctions = new JunctionSet(new File(pathJunctions), juncType);
            System.err.println("done reading junction set.");
            if (outputForEvents) {
                events = new ASEventSet(junctions);
                System.err.println("done generating event set.");
            }
        } else {
            events = new ASEventSet(new File(pathJunctions));
            System.err.println("done reading event set.");
            junctions = new JunctionSet(events);
            System.err.println("done generating junction set.");
        }

        HashMap<String, Double> totalNumMappedReads = new HashMap<>();
        HashMap<Junction, HashMap<String, Double>> juncCount = new HashMap<>();
        HashMap<ASEvent, HashMap<String, ArrayList<Double>>> eventCount = new HashMap<>();
        for (String chrom : junctions.getJunctions().keySet()) {
            System.err.println("initializing counting data for chromosome: " + chrom);
            for (Junction junc : junctions.getJunctions().get(chrom)) {
                juncCount.put(junc, new HashMap<String, Double>());
            }
        }
        if (outputForEvents) {
            for (ASEvent event : events.getAllEvents()) {
                eventCount.put(event, new HashMap<String, ArrayList<Double>>());
            }
        }

        firstSAMIndex++;

        System.err.println("\nstart counting...");
        // start reading SAM/BAM files
        for (int i = firstSAMIndex; i < args.length; i++) {
            String pathMapping = args[i];
            File mappingFile = new File(pathMapping);

            // Read counting
            ReadCounter counter = new ReadCounter(mappingFile, junctions, strandSpecific);
            counter.estimateJunctionCounts(considerNH, readCollapse);
            HashMap<Junction, Double> count = counter.getJunctionCounts();
            for (Junction junc : count.keySet()) {
                juncCount.get(junc).put(args[i], count.get(junc));
            }
            totalNumMappedReads.put(args[i], counter.getTotalNumReads());

            if (outputForEvents) {
                HashMap<ASEvent, ArrayList<Double>> countEvents = events.quantifyInclusion(count);
                for (ASEvent event : countEvents.keySet()) {
                    eventCount.get(event).put(args[i], countEvents.get(event));
                }
            }

            System.err.println("done " + args[i] + "\n");
        }

        // output
        if (outputForEvents) {
            String header = "GROUP\tEVENT_TYPE\tINCLU_LEFT\tINCLU_RIGHT\tEXCLU_LEFT\tEXCLU_RIGHT";
            for (int i = firstSAMIndex; i < args.length; i++) {
                header = header + "\tINCLU:" + args[i] + "\tEXCLU:" + args[i];
            }
            System.out.println(header);

            java.util.TreeSet<ASEvent> sortedEvents = new java.util.TreeSet<>(new java.util.Comparator<ASEvent>() {

                @Override
                public int compare(ASEvent arg0, ASEvent arg1) {
                    return (arg0.compareTo(arg1));
                }
            });
            sortedEvents.addAll(eventCount.keySet());

            for (Iterator<ASEvent> it = sortedEvents.iterator(); it.hasNext();) {
                ASEvent event = it.next();
                String info = event.toString();
                String readNum = "";
                for (int i = firstSAMIndex; i < args.length; i++) {
                    ArrayList<Double> nums = eventCount.get(event).get(args[i]);
                    readNum += "\t" + nums.get(0).intValue() + "\t" + nums.get(1).intValue();
                }
                System.out.println(info + readNum);
            }
        } else {
            String header = "JUNC_CHROM\tJUNC_START\tJUNC_END\tJUNC_STRAND";
            for (int i = firstSAMIndex; i < args.length; i++) {
                header = header + "\t" + args[i];
            }
            System.out.println(header);

            java.util.TreeSet<Junction> sortedJunctions = new java.util.TreeSet<>(new java.util.Comparator<Junction>() {

                @Override
                public int compare(Junction o1, Junction o2) {
                    return o1.compareTo(o2);
                }
            });

            String totalReads = String.valueOf(totalNumMappedReads.get(args[firstSAMIndex]).intValue());
            if (args.length > firstSAMIndex + 1) {
                for (int i = firstSAMIndex + 1; i < args.length; i++) {
                    totalReads = totalReads + "\t" + String.valueOf(totalNumMappedReads.get(args[i]).intValue());
                }
            }
            System.out.println("TOTAL_READS\tNA\tNA\tNA\t" + totalReads);

            sortedJunctions.addAll(juncCount.keySet());
            for (Iterator<Junction> it = sortedJunctions.iterator(); it.hasNext();) {
                Junction junc = it.next();
                String readNum = "";
                for (int i = firstSAMIndex; i < args.length; i++) {
                    readNum = readNum + "\t" + String.valueOf(juncCount.get(junc).get(args[i]).intValue());
                }
                System.out.println(junc.getChrom() + "\t" + junc.getStartSite() + "\t" + junc.getEndSite() + "\t" + junc.getStrand() + readNum);
            }
        }
    }

    
    
    
    /**
     * The main function for command line running
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        if(args.length==0 || args[0].equals("help")){
            showHelp();
            System.exit(0);
        }
        
        String cmd=args[0];
        if(cmd.equalsIgnoreCase("info")){
            runInfo(args);
        }
        else if(cmd.equalsIgnoreCase("count") || cmd.equalsIgnoreCase("rpkm")){
            runReadCounting(args);
        }
        else if(cmd.equalsIgnoreCase("uniq")){
            runUniq(args);
        }
        else if(cmd.equalsIgnoreCase("gtf2struc")){
            transformGTF2Struc(args);
        }
        else if(cmd.equalsIgnoreCase("tojuncs") || cmd.equalsIgnoreCase("toevents")){
            transform4Junction(args);
        }
        else if(cmd.equalsIgnoreCase("countjunc")){
            runJunctionCounting(args);
        }
        
        else if(cmd.equalsIgnoreCase("test")){ // a command for test only
            runTest(args);
        }
        else if(cmd.equalsIgnoreCase("gui")){
            try{
                java.lang.reflect.Method method=hitseq.gui.GraphicUserInterface.class.getMethod("main", String[].class);
                method.invoke(null, (Object) args);
            } catch (Exception e){
                System.err.println(e);
            }
        }
    }
    
    
    
    
    /**
     * The inner class of parameter set, to initial, read and return the parameters for the commands
     */
    private static class ParameterSet {
        private int firstSAMIndex;
        private int strandSpecific;
        private boolean considerNH;
        private boolean readCollapse;
        private int modeForMultiGenesOverlap;
        private int iterationLimit;
        private String annotFormat;
        private boolean outputForEvents=false;
        
        ParameterSet(String cmd){
            firstSAMIndex=1;
            if(cmd.equalsIgnoreCase("count") || cmd.equalsIgnoreCase("rpkm")){
                strandSpecific = 0;
                considerNH = false;
                readCollapse = false;
                modeForMultiGenesOverlap = 0;
                iterationLimit = 2;
                annotFormat = "struc";
            } else if(cmd.equalsIgnoreCase("countjunc") || cmd.equalsIgnoreCase("tojuncs") || cmd.equalsIgnoreCase("toevents")){
                annotFormat = "juncs";
                strandSpecific = 0;
                considerNH = false;
                readCollapse = false;
                outputForEvents = false;
            }
        }
        
        int getFirstSAMIdx(){
            return firstSAMIndex;
        }
        
        int getStrandedness(){
            return strandSpecific;
        }
        
        boolean getConsiderNHTag(){
            return(considerNH);
        }
        
        boolean getReadCollapseTag(){
            return(readCollapse);
        }
        
        int getModeForMultiGenesOverlap(){
            return(modeForMultiGenesOverlap);
        }
        
        int getIterationLimit(){
            return(iterationLimit);
        }
        
        String getAnnotFormat(){
            return(annotFormat);
        }
        
        boolean getOutputForEventsTag(){
            return(outputForEvents);
        }
        
        void readCommandLineArgs(String[] args){
            String cmd=args[0];
            while(true){
                java.util.regex.Pattern pattern=java.util.regex.Pattern.compile("^-");
                java.util.regex.Matcher matcher=pattern.matcher(args[this.firstSAMIndex]);
                if(matcher.find()){
                    String optionsString=matcher.replaceAll("");
                    for(int i=0; i<optionsString.length(); i++){
                        String option=optionsString.substring(i, i+1);
                        switch (option) {
                            case "s":
                                int idx=optionsString.indexOf("s");
                                if(idx<optionsString.length()-1){
                                    System.err.println("\nParameter error. The strandness needs to be given.\n");
                                    System.exit(0);
                                }
                                this.firstSAMIndex++;
                                try{
                                   this.strandSpecific=Integer.parseInt(args[this.firstSAMIndex]);
                                   if(this.strandSpecific>1 || this.strandSpecific<-1){
                                       System.err.println("\nParameter error. The strandness should be int {-1,0,1}.\n");
                                       System.exit(0);
                                   }
                                } catch(java.lang.NumberFormatException e){
                                   System.err.println("\nParameter error. The strandness should be int {-1,0,1}.\n");
                                   System.exit(0);
                                }
                                break;
                            case "n":
                                this.considerNH=true;
                                break;
                            case "c":
                                this.readCollapse=true;
                                break;
                            case "m":
                                idx=optionsString.indexOf("m");
                                if(idx<optionsString.length()-1){
                                    System.err.println("\nParameter error. The mode needs to be given.\n");
                                    System.exit(0);
                                }
                                this.firstSAMIndex++;
                                try{
                                   this.modeForMultiGenesOverlap=Integer.parseInt(args[this.firstSAMIndex]);
                                   if(this.modeForMultiGenesOverlap>3 || this.modeForMultiGenesOverlap<0){
                                       System.err.println("\nParameter error. The mode should be int 0-3.\n");
                                       System.exit(0);
                                   }
                                } catch(java.lang.NumberFormatException e){
                                   System.err.println("\nParameter error. The mode should be int 0-3.\n");
                                   System.exit(0);
                                }
                                break;
                            case "t":
                                idx=optionsString.indexOf("t");
                                if(idx<optionsString.length()-1){
                                    System.err.println("\nParameter error. The iteration limit needs to be given.\n");
                                    System.exit(0);
                                }
                                this.firstSAMIndex++;
                                try{
                                   this.iterationLimit=Integer.parseInt(args[this.firstSAMIndex]);
                                   if(this.iterationLimit<=0){
                                       System.err.println("\nParameter error. The iteration limit should be positive integer.\n");
                                       System.exit(0);
                                   }
                                } catch(java.lang.NumberFormatException e){
                                   System.err.println("\nParameter error. The iteration limit should be positive integer.\n");
                                   System.exit(0);
                                }
                                break;
                            case "a":
                                idx=optionsString.indexOf("a");
                                if(idx<optionsString.length()-1){
                                    System.err.println("\nParameter error. The annotation file type needs to be given.\n");
                                    System.exit(0);
                                }
                                this.firstSAMIndex++;
                                this.annotFormat=args[this.firstSAMIndex];
                                if(cmd.equalsIgnoreCase("countjunc") || cmd.equalsIgnoreCase("tojuncs") || cmd.equalsIgnoreCase("toevents")){
                                    if ((!this.annotFormat.equalsIgnoreCase("gtf")) && (!this.annotFormat.equalsIgnoreCase("bed")) && (!this.annotFormat.equalsIgnoreCase("juncs")) && (!this.annotFormat.equalsIgnoreCase("bam")) && (!this.annotFormat.equalsIgnoreCase("events"))) {
                                        System.err.println("\nParameter error. The mode should be one of \"juncs\", \"gtf\", \"bed\", \"bam\" and \"events\".\n");
                                        System.exit(0);
                                    } 
                                } else if (cmd.equalsIgnoreCase("count") || cmd.equalsIgnoreCase("rpkm")) {
                                    if ((!this.annotFormat.equalsIgnoreCase("gtf")) && (!this.annotFormat.equalsIgnoreCase("bed")) && (!this.annotFormat.equalsIgnoreCase("struc"))) {
                                        System.err.println("\nParameter error. The mode should be one of \"struc\", \"gtf\" and \"bed\".\n");
                                        System.exit(0);
                                    }
                                }
                                break;
                            case "e":
                                outputForEvents=true;
                                break;
                            default:
                                System.err.println("\nParameter error. No parameter "+option+"\n");
                                break;
                        }
                        System.err.println("option added: '"+option+"'");
                    }
                    firstSAMIndex++;
                }
                else
                    break;
            }
        }
    }
}

