/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hitseq;

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
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        /**
         * The command of this package could be the following.
         * 1. info - to get the basic information of the mapping)
         * 2. count - count the number of reads covering each gene in annotation
         * 3. rpkm - estimate the RPKM for each gene in annotation
         * 4. uniq - extract the uniquely mapped reads in the mapping
         */
        
        // TODO code application logic here
        if(args.length==0 || args[0].equals("help")){
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
                    + "           countjunc  Given junction list in junc/bed/gtf format, output table of read count of each junction for the input alignments\n");
            System.exit(0);
        }
        
        String cmd=args[0];
        if(cmd.equalsIgnoreCase("info")){
            String pathMapping=args[1];
            File mappingFile=new File(pathMapping);
            MappingProcessor processor=new MappingProcessor(mappingFile);
            processor.collectMappingInformation();
            ArrayList<Integer> information=processor.getMappingInfo();
            
            System.out.println("\n"+args[1]+":");
            System.out.printf("%45s|          %d\n","Number of reads: ",information.get(0));
            System.out.printf("%45s|          %d\n","Number of mapped reads: ",information.get(1));
            System.out.printf("%45s|          %d\n","Number of uniquely mapped reads: ",information.get(2));
            System.out.println();
            if(information.size()>3)
                for(int i=3; i<information.size(); i+=2){
                    System.out.printf("%45s|          %d\n","Number of reads with "+(i-3)/2+" mismatch(es)",information.get(i));
                    System.out.printf("%46s          %d\n","(UNIQUE READS)",information.get(i+1));
                }
            System.out.println();
        }
        else if(cmd.equalsIgnoreCase("count") || cmd.equalsIgnoreCase("rpkm")){
            if(args.length==1 || args[1].equals("-h")){
                System.err.println("\nThis is the help of '"+cmd.toLowerCase()+"' command of HiTSeq.");
                System.err.println("Usage: java -jar HiTSeq.jar "+cmd.toLowerCase()+" [options] <annotation.struc> <in.bam> [in2.bam ...]\n"
                        + "   Or: HiTSeq.sh "+cmd.toLowerCase()+" [options] <annotation.struc> <in.bam> [in2.bam ...]");
                System.err.println("\n"
                        + "Options: -h        This help page\n"
                        + "         -s [int]  Strandedness (default: 0 - no strand information; 1 - same strandness; -1 - opposite strandness)\n"
                        + "         -n        For reads mapped to n-loci, assign 1/n read to each hit\n"
                        + "         -c        Do read collapse to remove PCR duplicates\n"
                        + "         -m [int]  The mode to deal with multi-gene hits (default: mode 0 - abandon ambiguous reads)\n"
                        + "         -t [int]  The maximum iteration time to assign ambiguous reads (default: 2). Only work with -m 3\n"
                        + "         -a [str]  The file type of annotation file (default: struc format)\n");
                System.exit(0);
            }
            
            int strandSpecific=0;
            boolean considerNH=false;
            boolean readCollapse=false;
            int modeForMultiGenesOverlap=0;
            int iterationLimit=2;
            String annotFormat="struc";
            
            int firstSAMIndex=1;
            while(true){
                java.util.regex.Pattern pattern=java.util.regex.Pattern.compile("^-");
                java.util.regex.Matcher matcher=pattern.matcher(args[firstSAMIndex]);
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
                                firstSAMIndex++;
                                try{
                                   strandSpecific=Integer.parseInt(args[firstSAMIndex]);
                                   if(strandSpecific>1 || strandSpecific<-1){
                                       System.err.println("\nParameter error. The strandness should be int {-1,0,1}.\n");
                                       System.exit(0);
                                   }
                                } catch(java.lang.NumberFormatException e){
                                   System.err.println("\nParameter error. The strandness should be int {-1,0,1}.\n");
                                   System.exit(0);
                                }
                                break;
                            case "n":
                                considerNH=true;
                                break;
                            case "c":
                                readCollapse=true;
                                break;
                            case "m":
                                idx=optionsString.indexOf("m");
                                if(idx<optionsString.length()-1){
                                    System.err.println("\nParameter error. The mode needs to be given.\n");
                                    System.exit(0);
                                }
                                firstSAMIndex++;
                                try{
                                   modeForMultiGenesOverlap=Integer.parseInt(args[firstSAMIndex]);
                                   if(modeForMultiGenesOverlap>3 || modeForMultiGenesOverlap<0){
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
                                firstSAMIndex++;
                                try{
                                   iterationLimit=Integer.parseInt(args[firstSAMIndex]);
                                   if(iterationLimit<=0){
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
                                firstSAMIndex++;
                                annotFormat=args[firstSAMIndex];
                                if((!annotFormat.equalsIgnoreCase("gtf")) && (!annotFormat.equalsIgnoreCase("bed")) && (!annotFormat.equalsIgnoreCase("struc"))){
                                    System.err.println("\nParameter error. The mode should be one of \"struc\", \"gtf\" and \"bed\".\n");
                                    System.exit(0);
                                }
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
            
            String pathAnnotation=args[firstSAMIndex];            
            Annotation annotation=new Annotation(new File(pathAnnotation), annotFormat);
            if(modeForMultiGenesOverlap==0)
                annotation.estimateExclusiveGeneLength();
            HashMap<String, Double> totalNumMappedReads=new HashMap<>();
            HashMap<String, HashMap<String, Double>> readCount=new HashMap<>();
            HashMap<String, HashMap<String, Double>> fpkm=new HashMap<>();
            for(String gene : annotation.getGeneSet()){
                readCount.put(gene, new HashMap<String, Double>());
                fpkm.put(gene, new HashMap<String, Double>());
            }
            firstSAMIndex++;
            
            for(int i=firstSAMIndex; i<args.length; i++){
                String pathMapping=args[i];
                File mappingFile=new File(pathMapping);
                
                // Read counting
                annotation.resetPointer();
                ReadCounter counter=new ReadCounter(mappingFile,annotation,strandSpecific,modeForMultiGenesOverlap);
                counter.estimateCounts(considerNH, readCollapse, iterationLimit);
                HashMap<String, Double> count=counter.getCounts();
                for(String gene : count.keySet())
                    readCount.get(gene).put(args[i], count.get(gene));
                totalNumMappedReads.put(args[i], counter.getTotalNumReads());
                
                // Calculate RPKM if necessary
                if(cmd.equalsIgnoreCase("rpkm")){
                    counter.estimateRPKM();
                    HashMap<String, Double> fpkmGene=counter.getRPKM();
                    for (String gene : fpkmGene.keySet()) {
                        fpkm.get(gene).put(args[i], fpkmGene.get(gene));
                    }
                }
                
                System.err.println("done "+args[i]+"\n");
            }
            
            String header=cmd.equalsIgnoreCase("count") ? "GENE_ID\tLENGTH" : "GENE_ID";
            for(int i=firstSAMIndex; i<args.length; i++)
                header=header+"\t"+args[i];
            System.out.println(header);
            
            java.util.TreeSet<String> sortedGeneNames=new java.util.TreeSet<>(new java.util.Comparator<String>() {
                @Override
                public int compare(String o1, String o2) {
                    return o1.compareTo(o2);
                }
            });
            
            if(cmd.equalsIgnoreCase("count")){
                String totalReads=String.valueOf(totalNumMappedReads.get(args[firstSAMIndex]).intValue());
                if(args.length>firstSAMIndex+1)
                    for(int i=firstSAMIndex+1; i<args.length; i++)
                        totalReads=totalReads+"\t"+String.valueOf(totalNumMappedReads.get(args[i]).intValue());
                System.out.println("TOTAL_READS\tNA\t"+totalReads);

                
                for(String gene : readCount.keySet())
                    sortedGeneNames.add(gene);
                for (Iterator<String> it = sortedGeneNames.iterator(); it.hasNext();) {
                    String gene = it.next();
                    
                    int geneLength;
                    if(modeForMultiGenesOverlap==0){ // abandom ambiguous reads, use exclusive length
                        if(strandSpecific==0) // no strand information
                            geneLength=annotation.getExclusiveGeneLengthNoStrand(gene);
                        else
                            geneLength=annotation.getExclusiveGeneLength(gene);
                    } else // not to abandom ambiguous reads, use total length
                        geneLength=annotation.getGeneLength(gene);
                    
                    String readNum="";
                    for(int i=firstSAMIndex; i<args.length; i++){
                        if(modeForMultiGenesOverlap==1)
                            readNum=readNum+"\t"+String.valueOf(readCount.get(gene).get(args[i]));
                        else
                            readNum=readNum+"\t"+String.valueOf(readCount.get(gene).get(args[i]).intValue());
                    }
                    System.out.println(gene+"\t"+geneLength+readNum);
                }
            }
            else if(cmd.equalsIgnoreCase("rpkm")){
                for(String gene : fpkm.keySet())
                    sortedGeneNames.add(gene);
                for (Iterator<String> it = sortedGeneNames.iterator(); it.hasNext();) {
                    String gene = it.next();
                    String fpkmGene="";
                    for(int i=firstSAMIndex; i<args.length; i++)
                        fpkmGene=fpkmGene+"\t"+String.valueOf(fpkm.get(gene).get(args[i]));
                    System.out.println(gene+fpkmGene);
                }
            }   
        }
        else if(cmd.equalsIgnoreCase("uniq")){
            File inputSam=new File(args[1]);
            MappingProcessor processor=new MappingProcessor(inputSam);
            File outputSam;
            if(args.length>2)
                outputSam=new File(args[2]);
            else
                outputSam=new File("unique."+args[1]);
            
            int numUniqueReads=processor.extractUniquelyMappedReads(outputSam);
            System.out.println("Total Number of Uniquely Mapped Reads of "+args[1]+": "+numUniqueReads);
        }
        else if(cmd.equalsIgnoreCase("gtf2struc")){
            String pathAnnotation=args[1];            
            Annotation annotation=new Annotation(new File(pathAnnotation), "gtf");
            annotation.outputInStruc();
        }
        else if(cmd.equalsIgnoreCase("tojuncs")){
            if(args.length==1 || args[1].equals("-h")){
                System.err.println("\nThis is the help of 'tojuncs' command of HiTSeq.");
                System.err.println("Usage: java -jar HiTSeq.jar tojuncs [options] <in.file> [in2.file ...]\n"
                        + "   Or: HiTSeq.sh tojunc [options] <in.file> [in2.file ...]");
                System.err.println("\n"
                        + "Options: -h        This help page\n"
                        + "         -a [str]  The file type of junction annotation file (default: juncs format)\n"
                        + "         -s [int]  Strandedness of BAM/SAM input (default: 0 [no strand information]; 1/-1 [same/opposite strandness]; only work with '-a bam')\n");
                System.exit(0);
            }
            
            String fileType="juncs";
            int strandSpecific=0;
            int firstInputIndex=1;
            while(true){
                java.util.regex.Pattern pattern=java.util.regex.Pattern.compile("^-");
                java.util.regex.Matcher matcher=pattern.matcher(args[firstInputIndex]);
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
                                firstInputIndex++;
                                try{
                                   strandSpecific=Integer.parseInt(args[firstInputIndex]);
                                   if(strandSpecific>1 || strandSpecific<-1){
                                       System.err.println("\nParameter error. The strandness should be int {-1,0,1}.\n");
                                       System.exit(0);
                                   }
                                } catch(java.lang.NumberFormatException e){
                                   System.err.println("\nParameter error. The strandness should be int {-1,0,1}.\n");
                                   System.exit(0);
                                }
                                break;
                            case "a":
                                idx=optionsString.indexOf("a");
                                if(idx<optionsString.length()-1){
                                    System.err.println("\nParameter error. The annotation file type needs to be given.\n");
                                    System.exit(0);
                                }
                                firstInputIndex++;
                                fileType=args[firstInputIndex];
                                if((!fileType.equalsIgnoreCase("gtf")) && (!fileType.equalsIgnoreCase("bed")) && (!fileType.equalsIgnoreCase("juncs"))){
                                    System.err.println("\nParameter error. The mode should be one of \"juncs\", \"gtf\" and \"bed\".\n");
                                    System.exit(0);
                                }
                                break;
                            default:
                                System.err.println("\nParameter error. No parameter "+option+"\n");
                                break;
                        }
                        System.err.println("option added: '"+option+"'");
                    }
                    firstInputIndex++;
                }
                else
                    break;
            }
            
            JunctionSet junctions=new JunctionSet();
            for(int i=firstInputIndex; i<args.length; i++){
                String pathBed=args[i];
                if(fileType.equalsIgnoreCase("bam"))
                    junctions.addJunctionSet(new File(pathBed), strandSpecific);
                else
                    junctions.addJunctionSet(new File(pathBed), fileType);
            }
            junctions.outputInJuncs(true); 
        }
        else if(cmd.equalsIgnoreCase("countjunc")){
            if(args.length==1 || args[1].equals("-h")){
                System.err.println("\nThis is the help of 'countjunc' command of HiTSeq.");
                System.err.println("Usage: java -jar HiTSeq.jar countjunc [options] <junc.file> <in.bam> [in2.bam ...]\n"
                        + "   Or: HiTSeq.sh countjunc [options] <junc.file> <in.bam> [in2.bam ...]");
                System.err.println("\n"
                        + "Options: -h        This help page\n"
                        + "         -a [str]  The file type of junction annotation file (default: juncs format)\n"
                        + "         -n        For reads mapped to n-loci, only assign 1/n read to each hit\n"
                        + "         -c        Do read collapse to remove PCR duplicates\n"
                        + "         -s [int]  Strandedness of BAM/SAM input (default: 0 [no strand information]; 1/-1 [same/opposite strandness]; only work with '-a bam')\n");
                System.exit(0);
            }
            
            String juncType="juncs";
            int strandSpecific=0;
            boolean considerNH=false;
            boolean readCollapse=false;
            
            int firstSAMIndex=1;
            while(true){
                java.util.regex.Pattern pattern=java.util.regex.Pattern.compile("^-");
                java.util.regex.Matcher matcher=pattern.matcher(args[firstSAMIndex]);
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
                                firstSAMIndex++;
                                try{
                                   strandSpecific=Integer.parseInt(args[firstSAMIndex]);
                                   if(strandSpecific>1 || strandSpecific<-1){
                                       System.err.println("\nParameter error. The strandness should be int {-1,0,1}.\n");
                                       System.exit(0);
                                   }
                                } catch(java.lang.NumberFormatException e){
                                   System.err.println("\nParameter error. The strandness should be int {-1,0,1}.\n");
                                   System.exit(0);
                                }
                                break;
                            case "n":
                                considerNH=true;
                                break;
                            case "c":
                                readCollapse=true;
                                break;
                            case "a":
                                idx=optionsString.indexOf("a");
                                if(idx<optionsString.length()-1){
                                    System.err.println("\nParameter error. The annotation file type needs to be given.\n");
                                    System.exit(0);
                                }
                                firstSAMIndex++;
                                juncType=args[firstSAMIndex];
                                if((!juncType.equalsIgnoreCase("gtf")) && (!juncType.equalsIgnoreCase("bed")) && (!juncType.equalsIgnoreCase("juncs"))){
                                    System.err.println("\nParameter error. The mode should be one of \"juncs\", \"gtf\" and \"bed\".\n");
                                    System.exit(0);
                                }
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
            
            String pathJunctions=args[firstSAMIndex];
            JunctionSet junctions=new JunctionSet(new File(pathJunctions), juncType);
            HashMap<String, Double> totalNumMappedReads=new HashMap<>();
            HashMap<Junction, HashMap<String, Double>> juncCount=new HashMap<>();
            for(String chrom : junctions.getJunctions().keySet())
                for(Junction junc : junctions.getJunctions().get(chrom))
                    juncCount.put(junc, new HashMap<String, Double>());
            firstSAMIndex++;
            System.err.println("done reading junction set...");
            
            for(int i=firstSAMIndex; i<args.length; i++){
                String pathMapping=args[i];
                File mappingFile=new File(pathMapping);
                
                // Read counting
                ReadCounter counter=new ReadCounter(mappingFile,junctions,strandSpecific);
                counter.estimateJunctionCounts(considerNH, readCollapse);
                HashMap<Junction, Double> count=counter.getJunctionCounts();
                for(Junction junc : count.keySet())
                    juncCount.get(junc).put(args[i], count.get(junc));
                totalNumMappedReads.put(args[i], counter.getTotalNumReads());
                
                System.err.println("done "+args[i]+"\n");
            }
            
            String header="JUNC_CHROM\tJUNC_START\tJUNC_END\tJUNC_STRAND";
            for(int i=firstSAMIndex; i<args.length; i++)
                header=header+"\t"+args[i];
            System.out.println(header);
            
            java.util.TreeSet<Junction> sortedJunctions=new java.util.TreeSet<>(new java.util.Comparator<Junction>() {
                @Override
                public int compare(Junction o1, Junction o2) {
                    return o1.compareTo(o2);
                }
            });
            
            String totalReads=String.valueOf(totalNumMappedReads.get(args[firstSAMIndex]).intValue());
            if(args.length>firstSAMIndex+1)
                for(int i=firstSAMIndex+1; i<args.length; i++)
                    totalReads=totalReads+"\t"+String.valueOf(totalNumMappedReads.get(args[i]).intValue());
            System.out.println("TOTAL_READS\tNA\tNA\tNA\t"+totalReads);

            sortedJunctions.addAll(juncCount.keySet());
            for (Iterator<Junction> it = sortedJunctions.iterator(); it.hasNext();) {
                Junction junc = it.next();
                String readNum="";
                for(int i=firstSAMIndex; i<args.length; i++)
                    readNum=readNum+"\t"+String.valueOf(juncCount.get(junc).get(args[i]).intValue());
                System.out.println(junc.getChrom()+"\t"+junc.getDonorSite()+"\t"+junc.getAcceptorSite()+"\t"+junc.getStrand()+readNum);
            }
        }
        
        else if(cmd.equalsIgnoreCase("annotoverlap")){ // a command for test only, to output overlapping gene pairs in the given struc annotation.
            Annotation annotation=new Annotation(new File(args[1]));
            annotation.estimateExclusiveGeneLength(true);
        }
        
        else if(cmd.equalsIgnoreCase("junctest")){ // a command for test only, to test junction set generation method
            /*
            String inputType=args[1];
            
            JunctionSet junctions=new JunctionSet();
            for(int i=2; i<args.length; i++){
                String pathBed=args[i];
                junctions.addJunctionSet(new File(pathBed), inputType);
            }
            junctions.outputInJunc(true); 
            */
            
            File juncFile=new File(args[1]);
            File inputFile=new File(args[2]);
            JunctionSet junctions=new JunctionSet(juncFile, "juncs");
            ReadCounter counter=new ReadCounter(inputFile, junctions, 0);
            counter.estimateJunctionCounts(false, true);
            HashMap<Junction, Double> juncCounts=counter.getJunctionCounts();
            
            java.util.TreeSet<Junction> sortedJunctions=new java.util.TreeSet<>(new java.util.Comparator<Junction>() {
                @Override
                public int compare(Junction o1, Junction o2) {
                    return o1.compareTo(o2);
                }
            });
            sortedJunctions.addAll(juncCounts.keySet());
            for(Iterator<Junction> it=sortedJunctions.iterator(); it.hasNext(); ){
                Junction junc=it.next();
                System.out.println(junc.getChrom()+"\t"+junc.getDonorSite()+"\t"+junc.getAcceptorSite()+"\t"+junc.getStrand()+"\t"+juncCounts.get(junc).intValue());
            }
        }
    }
}

