package varbin;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringReader;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;
import net.sourceforge.argparse4j.inf.Subparser;
import net.sourceforge.argparse4j.inf.Subparsers;

// TODO: Auto-generated Javadoc
//import org.apache.log4j.Logger;


/**
 * The Class Varbin.
 * Given variant alleles for a main sample (vcf and bam file)
 * Analyze same variant alleles in multiple background samples (bams)
 * Use this background bam info to predict if alleles are systematic errors or true variants
 * Generate output that includes "bin" values (1 thru 4) for each alt allele and likelihood histograms for each alt allele
 * @author jacobdurtschi
 */


public class Varbin {

    private static Namespace args = null;
    private static ArrayList<String> bkgdList = null;
    private static LeftAlignCall vcfLeftCall = null;
	private static VcfSplit mainVcfParts;
	private static BufferedWriter outWriter = null;

	private static ProcessedCalls processedSnv1 = null;
	private static ProcessedCalls processedSnv2 = null;
	private static ProcessedCalls processedIndel1 = null;
	private static ProcessedCalls processedIndel2 = null;
	
	private static BinData binDataSnv1 = null;
	private static BinData binDataSnv2 = null;
	private static BinData binDataIndel1 = null;
	private static BinData binDataIndel2 = null;

	/**
	 * The main method.
	 *
	 * @param args the arguments
	 */
	public static void main(String[] arguments) {


		
        ArgumentParser parser = ArgumentParsers.newArgumentParser("Jarbin")
                .description("JarBin (Java-Varbin)");
        Subparsers subparsers = parser.addSubparsers()
        		.title("JarBin tools")
        		.help("")
        		.dest("tool");
        Subparser parserBin = subparsers.addParser("bin")
        		.description("bin subcommands")
        		.help("calc bins and bin data for a vcf file using a list of bkgd sample files");
        Subparser parserBkgd = subparsers.addParser("bkgd")
        		.description("bkgd subcommands")
        		.help("precalculate bkgd data for one bkgd sample (i.e., make bkgd sample files for a bam file)");
        parserBkgd.addArgument("-bam")
                .metavar("BAM")
                .type(String.class)
                .required(true)
                .help("file path of bkgd sample bam file to process (required)");
        parserBkgd.addArgument("-p", "--output_prefix")
                .type(String.class)
                .required(true)
        		.dest("prefix")
                .help("prefix for bkgd output files for this sample");
        parserBkgd.addArgument("-f", "--folder")
                .type(String.class)
        		.dest("folder")
        		.setDefault("./")
                .required(false)
                .help("file path string for output (default = working folder)");
        parserBkgd.addArgument("--samtools")
                .type(String.class)
        		.dest("samtools")
        		.setDefault("")
                .help("path to samtools" +
                		" (default = assume in PATH)");
        parserBkgd.addArgument("--bcftools")
                .type(String.class)
        		.dest("bcftools")
        		.setDefault("")
                .help("path to samtools" +
                		" (default = assume in PATH)");
        parserBkgd.addArgument("--bgzip")
                .type(String.class)
        		.dest("bgzip")
        		.setDefault("")
                .help("path to bgzip" +
                		" (default = assume in PATH)");
        parserBkgd.addArgument("--tabix")
                .type(String.class)
        		.dest("tabix")
        		.setDefault("")
                .help("path to tabix" +
                		" (default = assume in PATH)");
        parserBkgd.addArgument("-r", "--reference")
                .type(String.class)
        		.dest("ref")
                .help("path to reference sequence");
        parserBkgd.addArgument("-gatk", "--gatk")
                .type(String.class)
        		.dest("gatk")
        		.setDefault("/Usr/local/bin/GenomeAnalysisTK.jar")
                .help("path to GenomeAnalysisTK.jar" +
                		" (default = /Usr/local/bin/GenomeAnalysisTK.jar)");
        parserBkgd.addArgument("--gatk_threads")
                .type(String.class)
        		.dest("threads")
        		.setDefault("4")
                .help("number of threads used by GATK UnifiedGenotyper" +
                		" (default = 4)");
        parserBkgd.addArgument("--gatk_mem_options")
                .type(String.class)
        		.dest("mem")
        		.setDefault("Xmx2g")
                .help("java memory option used by GATK UnifiedGenotyper." +
                		"Note that the dash (-) is not included to avoid" +
                		"shell parsing issues." +
                		" (default = Xmx1g)");
        parserBin.addArgument("-bam")
                .metavar("BAM")
                .type(String.class)
                .required(true)
                .help("file path of primary sample bam file (required)");
        parserBin.addArgument("-vcf")
                .metavar("VCF")
                .type(String.class)
                .required(true)
                .help("file path of variants of interest vcf file (required)");
        parserBin.addArgument("-bg", "--background_table")
        		.type(Arguments.fileType().acceptSystemIn().verifyCanRead())
        		.dest("bg")
        		.setDefault("-")
                .required(true)
                .help("file path to file containing list of background bams," +
                		" background files are identified based on order" +
                		" in this file. (required, - = stdin)");
        parserBin.addArgument("-p", "--output_prefix")
                .type(String.class)
                .required(true)
        		.dest("prefix")
                .help("prefix for bkgd output files for this sample");
        parserBin.addArgument("-f", "--folder")
                .type(String.class)
        		.dest("folder")
        		.setDefault("./varbin_results")
                .required(false)
                .help("file path string for output (default = ./varbin_results)");
        parserBin.addArgument("-o", "--output_table")
                .type(String.class)
        		.dest("out")
        		.setDefault("varbin_output.txt")
                .help("file name for output table" +
                		" (default = varbin_output.txt)");
        parserBin.addArgument("-gatk", "--gatk")
                .type(String.class)
        		.dest("gatk")
        		.setDefault("/Usr/local/bin/GenomeAnalysisTK.jar")
                .help("path to GenomeAnalysisTK.jar" +
                		" (default = /Usr/local/bin/GenomeAnalysisTK.jar)");
        parserBin.addArgument("-r", "--reference")
                .type(String.class)
        		.dest("ref")
                .required(true)
                .help("path to reference sequence");
        parserBin.addArgument("--gatk_threads")
                .type(String.class)
        		.dest("threads")
        		.setDefault("4")
                .help("number of threads used by GATK UnifiedGenotyper" +
                		" (default = 4)");
        parserBin.addArgument("--gatk_mem_options")
                .type(String.class)
        		.dest("mem")
        		.setDefault("Xmx2g")
                .help("java memory option used by GATK UnifiedGenotyper." +
                		"Note that the dash (-) is not included to avoid" +
                		"shell parsing issues." +
                		" (default = Xmx2g)");
		try {
			args= parser.parseArgs(arguments);
			System.out.println(args);
		} catch (ArgumentParserException e) {
				parser.handleError(e);
		}
		
		
		
		
		if (args.get("tool").toString().equals("bin")) {
			BinCalc();
		} else if (args.get("tool").toString().equals("bkgd")) {
			BkgdCalc();
		}
	} // end of main
	

	private static void BkgdCalc() {
		String[] alleleVcf;
		String refVcf;


		//use samtools and bcftools to get full set of variant/reference evidence
		System.out.println("Calling ");
		MpileupCall.makeCall(args);
		alleleVcf = MpileupCall.getVcfOut();
		refVcf = MpileupCall.getCovOut();
		

		//call gatk unifiedgenotyper to get final bkgd vcf files
		String[] ugOut = new String[3];
		ugOut[0] = args.getString("folder") + "/" + args.getString("prefix") + ".temp1.vcf";
		ugOut[1] = args.getString("folder") + "/" + args.getString("prefix") + ".temp2.vcf";
		ugOut[2] = args.getString("folder") + "/" + args.getString("prefix") + ".temp3.vcf";
		UnifiedGenotyperCall[] ugCall = new UnifiedGenotyperCall[3];
		ugCall[0] = new UnifiedGenotyperCall(alleleVcf[0], args.getString("bam"), VariantType.BOTH, ugOut[0], args);
		ugCall[1] = new UnifiedGenotyperCall(alleleVcf[1], args.getString("bam"), VariantType.BOTH, ugOut[1], args);
		ugCall[2] = new UnifiedGenotyperCall(alleleVcf[2], args.getString("bam"), VariantType.BOTH, ugOut[2], args);
		ugCall[0].callGatk();
		ugCall[1].callGatk();
		ugCall[2].callGatk();
		System.out.println("made ug calls");
		System.out.println(ugCall[0].vcfOut);
		System.out.println(ugCall[1].vcfOut);
		System.out.println(ugCall[2].vcfOut);
		System.out.println("completed ug calls");
		String[] ugLeftOut = new String[3];
		ugLeftOut[0] = args.getString("folder") + "/" + args.getString("prefix") + ".temp1Left.vcf";
		ugLeftOut[1] = args.getString("folder") + "/" + args.getString("prefix") + ".temp2Left.vcf";
		ugLeftOut[2] = args.getString("folder") + "/" + args.getString("prefix") + ".temp3Left.vcf";
		LeftAlignCall[] leftCall = new LeftAlignCall[3];
		leftCall[0] = new LeftAlignCall(ugOut[0], ugLeftOut[0], args);
		leftCall[1] = new LeftAlignCall(ugOut[1], ugLeftOut[1], args);
		leftCall[2] = new LeftAlignCall(ugOut[2], ugLeftOut[2], args);
		leftCall[0].makeCall();
		leftCall[1].makeCall();
		leftCall[2].makeCall();
		String comboVcf= args.getString("folder") + "/" + args.getString("prefix") + ".bkgd.temp";
		combineVcfs(args, comboVcf, refVcf, leftCall[0].vcfOut, leftCall[1].vcfOut, leftCall[2].vcfOut);
		String sortedVcf = args.getString("folder") + "/" + args.getString("prefix") + ".bkgd.vcf";
		sortVcf(args, comboVcf,  sortedVcf);
		String zipVcf = args.getString("folder") + "/" + args.getString("prefix") + ".bkgd.vcf.gz";
		bgzip(args, sortedVcf,  zipVcf);
		tabix(args, zipVcf);
		System.out.println("Completed Jarbin bkgd file generation");
	}
	
	private static void sortVcf(Namespace args, String vcf, String sortedVcf) {
		String sortComm = "cat <(grep \\^# " + vcf + ") <(grep -v \\^# " + vcf + " | sort -k1,1 -k2,2n) > " + sortedVcf;
		List<String> call = new ArrayList<String>(Arrays.asList("/bin/bash", "-c", sortComm));

		// launch the command and grab stdin/stdout and stderr
		ProcessBuilder pb = new ProcessBuilder(call);
		Process p = null;

		try {
			System.out.println(call.toString());
			System.out.println(pb.command().toString());
			p = pb.start();
		} catch (IOException e) {
			System.err.println("ERROR: starting system call process to sort, bgzip, and tabix");
			System.err.println(pb.command().toString());
			e.printStackTrace();
		}
		try {
			p.waitFor();
		} catch (InterruptedException e) {
			System.err.println("ERROR: executing system call process to sort, bgzip, and tabix");
			System.err.println(pb.command().toString());
			e.printStackTrace();
		}
		
	}

	private static void bgzip(Namespace args, String vcf, String zipVcf) {
		//TODO bgzip to final file
		//TODO tabix index
		String comm = args.getString("bgzip") + " -c " + vcf + " > " + zipVcf;
		List<String> call = new ArrayList<String>(Arrays.asList("/bin/bash", "-c", comm));

		// launch the command and grab stdin/stdout and stderr
		ProcessBuilder pb = new ProcessBuilder(call);
		Process p = null;

		try {
			System.out.println(call.toString());
			System.out.println(pb.command().toString());
			p = pb.start();
		} catch (IOException e) {
			System.err.println("ERROR: starting system call process to bgzip");
			System.err.println(pb.command().toString());
			e.printStackTrace();
		}
		try {
			p.waitFor();
		} catch (InterruptedException e) {
			System.err.println("ERROR: executing system call process to bgzip");
			System.err.println(pb.command().toString());
			e.printStackTrace();
		}
		
	}
	private static void tabix(Namespace args, String zipVcf) {
		//TODO tabix index
		String comm = args.getString("tabix") + " -p vcf " + zipVcf;
		List<String> call = new ArrayList<String>(Arrays.asList("/bin/bash", "-c", comm));

		// launch the command and grab stdin/stdout and stderr
		ProcessBuilder pb = new ProcessBuilder(call);
		Process p = null;

		try {
			System.out.println(call.toString());
			System.out.println(pb.command().toString());
			p = pb.start();
		} catch (IOException e) {
			System.err.println("ERROR: starting system call process to bgzip");
			System.err.println(pb.command().toString());
			e.printStackTrace();
		}
		try {
			p.waitFor();
		} catch (InterruptedException e) {
			System.err.println("ERROR: executing system call process to bgzip");
			System.err.println(pb.command().toString());
			e.printStackTrace();
		}
	}
	private static void combineVcfs(Namespace args, String comboVcfPath, String refVcf, String...vcf) {
		String line;
		
		BufferedWriter comboWtr = null;
		try {
			comboWtr = new BufferedWriter(new FileWriter(comboVcfPath));
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		BufferedReader[] vcfReader = new BufferedReader[vcf.length];
		for (int i = 0 ; i< vcf.length ; i++) {
			try {
				vcfReader[i] = new BufferedReader(new FileReader(vcf[i]));
			} catch (FileNotFoundException e1) {
				System.err.println("ERROR: temp vcf file not found");
				e1.printStackTrace();
			}

			try {
				while ((line = vcfReader[i].readLine()) != null) {
					if (i == 0) {
						comboWtr.write(line + "\n");
					} else if (! line.startsWith("#")) {
						comboWtr.write(line + "\n");
					}
				}
			} catch (IOException e) {
				System.err.println("ERROR: writing to file " + comboVcfPath);
				e.printStackTrace();
			}
			
			
			try {
				vcfReader[i].close();
			} catch (IOException e) {
				System.err.println("ERROR: closing file " + vcf[i]);
				e.printStackTrace();
			}
		}

		BufferedReader refVcfReader = null;
		try {
			refVcfReader = new BufferedReader(new FileReader(refVcf));
		} catch (FileNotFoundException e1) {
			System.err.println("ERROR: writing to file " + refVcf);
			e1.printStackTrace();
		}
		try {
			while ((line = refVcfReader.readLine()) != null) {
				if (! line.startsWith("#")) {
					comboWtr.write(line + "\n");
				}
			}
		} catch (IOException e) {
			System.err.println("ERROR: writing to file " + refVcf);
			e.printStackTrace();
		}
		try {
			refVcfReader.close();
		} catch (IOException e) {
			System.err.println("ERROR: closing file " + refVcf);
			e.printStackTrace();
		}
		try {
			comboWtr.flush();
			comboWtr.close();
		} catch (IOException e) {
			System.err.println("ERROR: closing file " + comboVcfPath);
			e.printStackTrace();
		}
	}
	private static void BinCalc() {

		//print some info
        System.out.println("--Input Parameters--");
        System.out.println("main sample bam = " + args.get("bam"));
        System.out.println("variants of iterest vcf = " + args.get("vcf"));
        System.out.println("background bam table = " + args.get("background_table"));
        System.out.println("output folder = " + args.get("folder"));
		System.out.println("");
		
		//create output directory if it does not exist
	    File dirOut = new File(args.getString("folder"));
	    if (dirOut.exists() && dirOut.isFile()){ //don't create if file with same name already exists
			System.err.println("ERROR: output directory not created, it exists as a file  " + (String) args.get("folder"));
			System.exit(1);
	    } else if (dirOut.exists() && dirOut.isDirectory()) { //continue if directory exists
			System.err.println("Warning: output directory already existed  " + (String) args.get("folder"));
	    } else if (!dirOut.mkdir()){ //otherwise try to make directory
			System.err.println("ERROR: output directory not created (perhaps permissions issue?)  " + (String) args.get("folder"));
			System.exit(1);
	    }

	    //run gatk leftalign and trim on the input vcf to propperly match background vcf variants
	    String vcfInName = Paths.get(args.getString("vcf")).getFileName().toString();
	    String vcfLeftOut = args.getString("folder") + "/" + vcfInName + ".leftaligned.vcf";
	    vcfLeftCall = new LeftAlignCall(args.getString("vcf"), vcfLeftOut, args);
	    vcfLeftCall.makeCall();

		//split main input vcf into 4 (snv alt allele 1, snv alt allele 2, indel alt allele 1, indel alt allele 2)
	    //this avoids any potential issues with how second alt alleles are handled or annotated
	    //each of these component vcf files will be used as allele targets for BAM variant calling
		mainVcfParts = new VcfSplit(vcfLeftCall.vcfOut , (String) args.get("folder"));
		mainVcfParts.makeCall();

		
		System.out.println("");
		System.out.println("--Component vcf variant counts--");
		System.out.println("snv1Count: " + mainVcfParts.snv1Count);
		System.out.println("snv2Count: " + mainVcfParts.snv2Count);
		System.out.println("indel1Count: " + mainVcfParts.indel1Count);
		System.out.println("indel2Count: " + mainVcfParts.indel2Count);

		//Read in List of background bam file paths
		bkgdList = ReadBkgdList((File) args.get("bg"));
		System.out.println("-------------Bkgd file list-------------");
		Integer counter = 1;
		for (String vcf : bkgdList) {
			System.out.println((counter++).toString() + ".   " + vcf);
		}

		//TODO fix Tabix search in VarsFromBkgd.java to handle no position entry, or wrong alt entry
		//			 There will be either a hit, a X alt, or nothing
		//			 If hit, use vcf info
		//			 If wrong alt, use (either X or wrong alt info?????)
		//			 If nothing, make a default bkgd VarObj entry that can be identified when calculating bin stats
		//			 One of these must go into the same collection of variant objects as before for bin calc
		//TODO fix ProcessedCalls.java so that the processing is not done in constructor
		//TODO fix BinData.java to account for any changes to VarObj.java or VarList.java
		//TODO Improve output (table?, other format?)

		//TODO polish an outline of process for presentation
		//break primary vcf into 4 component vcfs (snv alt1/alt2 and indel alt1/alt2)
		//for each of the 4 component vcfs
		//if the component vcf contains at least one variant:
		//1. recall component vcf variants in primary bam with GATK
		//2. convert the new vcf output to variant list object
		//3. Tabix search component vcf variants in bkgd vcf
		//4. generate variant list object from each bkgd vcf
		//5. compare coponent vcf variant list to bkgd variant lists --> JarBin stats and bin # for each variant
		
		
		
		//--for first alt allele SNVs--
		
		String outPathSnv1 = args.getString("folder") + "/" + args.get("prefix") + ".tempSnv1.vcf";
		if (mainVcfParts.snv1Count > 0) {
			System.out.println("Starting main and background BAM file GATK calls for first alt allele, snv variants...");
			processedSnv1 = new ProcessedCalls(mainVcfParts.snv1Path, VariantType.SNP, (String) args.get("bam"), bkgdList, outPathSnv1, args);
		} else {
			System.out.println("There were no first alt allele, snv variants in the input vcf file");
		}
		
		//--for second alt allele SNVs--
		String outPathSnv2 = args.getString("folder") + "/" + args.getString("prefix") + ".tempSnv1.vcf";
		
		if (mainVcfParts.snv2Count > 0) {
			System.out.println("Starting main and background BAM file GATK calls for second alt allele, snv variants ("
						+ args.get("bam") + ")...");
			processedSnv2 = new ProcessedCalls(mainVcfParts.snv2Path, VariantType.SNP, (String) args.get("bam"), bkgdList, outPathSnv2, args);
		} else {
			System.out.println("There were no second alt allele, snv variants in the input vcf file");
		}

		//--for first alt allele Indels--
		String outPathIndel1 = args.getString("folder") + "/" + args.getString("prefix") + ".tempSnv1.vcf";
		
		if (mainVcfParts.indel1Count > 0) {
			System.out.println("Starting main and background BAM file GATK calls for first alt allele, indel variants ("
						+ args.get("bam") + ")...");
			processedIndel1 = new ProcessedCalls(mainVcfParts.indel1Path, VariantType.INDEL, (String) args.get("bam"), bkgdList, outPathIndel1, args);
		} else {
			System.out.println("There were no first alt allele, indel variants in the input vcf file");
		}
		
		//--for second alt allele Indels--
		String outPathIndel2 = args.getString("folder") + "/" + args.getString("prefix") + ".tempSnv1.vcf";
		
		if (mainVcfParts.indel2Count > 0) {
			System.out.println("Starting main and background BAM file GATK calls for second alt allele, indel variants ("
						+ args.get("bam") + ")...");
			processedIndel2 = new ProcessedCalls(mainVcfParts.indel2Path, VariantType.INDEL, (String) args.get("bam"), bkgdList, outPathIndel2, args);
		} else {
			System.out.println("There were no second alt allele, indel variants in the input vcf file");
		}
		
		//process the main and background variant call data
		//to get bin results
		if (mainVcfParts.snv1Count > 0) {
			binDataSnv1 = new BinData(processedSnv1);
		}
		if (mainVcfParts.snv2Count > 0) {
			binDataSnv2 = new BinData(processedSnv2);
		}
		if (mainVcfParts.indel1Count > 0) {
			binDataIndel1 = new BinData(processedIndel1);
		}
		if (mainVcfParts.indel2Count > 0) {
			binDataIndel2 = new BinData(processedIndel2);
		}
		
		//print each set of results to file 
		String outPath = args.getString("folder") + "/" + args.getString("out");
		try{
			outWriter = new BufferedWriter(new FileWriter(outPath));
		} catch (IOException e) {
			System.err.println("ERROR: creating output table file");
			e.printStackTrace();
			System.exit(1);
		}
		
		// --for snv alt1--
		if (mainVcfParts.snv1Count > 0) {
			try {
				outWriter.write("Bin data for snv first alt alleles\n");
			} catch (IOException e) {
				System.err.println("ERROR: writing to output table file");
				e.printStackTrace();
				System.exit(1);
			}
			writeBinData(binDataSnv1, outWriter);
		} else {
			try {
				outWriter.write("\n");
				outWriter.write("No snv first alt alleles were processed\n");
				outWriter.write("\n");
			} catch (IOException e) {
				System.err.println("ERROR: writing to output table file");
				e.printStackTrace();
				System.exit(1);
			}
		}

		// --for snv alt2--
		if (mainVcfParts.snv2Count > 0) {
			try {
				outWriter.write("Bin data for snv second alt alleles\n");
			} catch (IOException e) {
				System.err.println("ERROR: writing to output table file");
				e.printStackTrace();
				System.exit(1);
			}
			writeBinData(binDataSnv2, outWriter);
		} else {
			try {
				outWriter.write("\n");
				outWriter.write("No snv second alt alleles were processed\n");
				outWriter.write("\n");
			} catch (IOException e) {
				System.err.println("ERROR: writing to output table file");
				e.printStackTrace();
				System.exit(1);
			}
		}

		// --for indel alt1--
		if (mainVcfParts.indel1Count > 0) {
			try {
				outWriter.write("Bin data for indel first alt alleles\n");
			} catch (IOException e) {
				System.err.println("ERROR: writing to output table file");
				e.printStackTrace();
				System.exit(1);
			}
			writeBinData(binDataIndel1, outWriter);
		} else {
			try {
				outWriter.write("\n");
				outWriter.write("No indel first alt alleles were processed\n");
				outWriter.write("\n");
			} catch (IOException e) {
				System.err.println("ERROR: writing to output table file");
				e.printStackTrace();
				System.exit(1);
			}
		}

		// --for indel alt2--
		if (mainVcfParts.indel2Count > 0) {
			try {
				outWriter.write("Bin data for indel second alt alleles\n");
			} catch (IOException e) {
				System.err.println("ERROR: writing to output table file");
				e.printStackTrace();
				System.exit(1);
			}
			writeBinData(binDataIndel2, outWriter);
		} else {
			try {
				outWriter.write("\n");
				outWriter.write("No indel second alt alleles were processed\n");
				outWriter.write("\n");
			} catch (IOException e) {
				System.err.println("ERROR: writing to output table file");
				e.printStackTrace();
				System.exit(1);
			}
		}
		
		try {
			outWriter.close();
		} catch (IOException e) {
			System.err.println("ERROR: closing output table file");
			e.printStackTrace();
			System.exit(1);
		}
		/*
		System.out.println("--Variant list object contents--");
		for (VarObj thisVar : processedIndel1.mainVars.vars) {
			System.out.println(thisVar.chrom + ":" + thisVar.position + ", "
					+ StringUtils.join(thisVar.refList, ',').toString() + "->"
					+ StringUtils.join(thisVar.altList, ',') + ", "
					+ thisVar.varType.toString() + ", Genotype: "
					+ thisVar.genotype.toString() + ","
					+ thisVar.formatMap.get("GT") + " PLRD: " + thisVar.plrd + " PassFilter: " + thisVar.passFilter);
		}
		*/
		
		//TODO VarBinTable
		//get local sequence context
		//make summaries of count and list of bkgd files that had problematic values (low qual, low cov, low MQ etc.)

		//TODO Output Data
		//VarBinTable makes output table file
		//VarBinPlots makes output plot pdf document
		
		//TODO LATER merge table with csv annotations (do this by hand for the time being)
		
	}
	
	private static void writeBinData(BinData binData, BufferedWriter writer) {
		try {
			writer.write("chr\tpos\tref\talt\tbin\tvarType\tGT\tfailedFilterCount\tmainPLRD\tbgMedianPLRD\tstdDevPLRD\tbg3sigmaPLRD\tbg4sigmaPLRD\tbg5sigmaPLRD\tbg6sigmaPLRD\n");
			for (int i = 0 ; i < binData.bin.size() ; i++ ) {
				//print chr, pos, ref, alt, bin, variantType, GT, mainPLRD, median, 3sigma, 4sigma, 5sigma, 6sigma, 
			writer.write(((binData.chrom.get(i) == null) ? "." : binData.chrom.get(i)) + "\t"
					+ ((binData.position.get(i) == null) ? "." : binData.position.get(i)) + "\t"
					+ ((binData.ref.get(i) == null) ? "." : binData.ref.get(i)) + "\t"
					+ ((binData.alt.get(i) == null) ? "." : binData.alt.get(i)) + "\t"
					+ ((binData.bin.get(i) == null) ? "." : binData.bin.get(i)) + "\t"
					+ ((binData.variantType.get(i) == null) ? "." : binData.variantType.get(i).toString()) + "\t"
					+ ((binData.genotype.get(i) == null) ? "." : binData.genotype.get(i).toString()) + "\t"
					+ ((binData.failCount.get(i) == null) ? "." : binData.failCount.get(i).toString()) + "\t"
					+ ((binData.mainPLRD.get(i) == null) ? "." : binData.mainPLRD.get(i).toString()) + "\t"
					+ ((binData.median.get(i) == null) ? "." : binData.median.get(i).toString()) + "\t"
					+ ((binData.sigma.get(i) == null) ? "." : binData.sigma.get(i).toString()) + "\t"
					+ ((binData.sigma3.get(i) == null) ? "." : binData.sigma3.get(i).toString()) + "\t"
					+ ((binData.sigma4.get(i) == null) ? "." : binData.sigma4.get(i).toString()) + "\t"
					+ ((binData.sigma5.get(i) == null) ? "." : binData.sigma5.get(i).toString()) + "\t"
					+ ((binData.sigma6.get(i) == null) ? "." : binData.sigma6.get(i).toString()) + "\n");
			}
		} catch (IOException e) {
			System.err.println("ERROR: writing to output table file");
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	private static ArrayList<String> ReadBkgdList(File bamListFile) {
		BufferedReader bamListReader = null;
		ArrayList<String> bamList = null;
		
		try {
			bamListReader = new BufferedReader(new FileReader(bamListFile));
		} catch (IOException e) {
			System.err.println("ERROR: accessing/reading background bam list file " + bamListFile.getPath());
			e.printStackTrace();
			System.exit(1);
		}

		try {
			bamList = new ArrayList<String>();
			for(String line; (line = bamListReader.readLine()) != null; ) {
				if(! line.startsWith("#")){
					bamList.add(line);
				}
			}
			bamListReader.close();
		} catch (IOException e) {
			System.err.println("ERROR: accessing/reading background bam list file " + bamListFile.getPath());
			e.printStackTrace();
			System.exit(1);
		}
		return bamList;
	}

}