package varbin;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

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

	/**
	 * The main method.
	 *
	 * @param args the arguments
	 */
	public static void main(String[] arguments) {

        Namespace args = null;
        ArrayList<String> bamList = null;
		VcfSplit mainVcfParts;
		BufferedWriter outWriter = null;

		ProcessedCalls processedSnv1 = null;
		ProcessedCalls processedSnv2 = null;
		ProcessedCalls processedIndel1 = null;
		ProcessedCalls processedIndel2 = null;
		
		BinData binDataSnv1 = null;
		BinData binDataSnv2 = null;
		BinData binDataIndel1 = null;
		BinData binDataIndel2 = null;
		
        ArgumentParser parser = ArgumentParsers.newArgumentParser("Varbin")
                .description("Run Java-Varbin program");
        parser.addArgument("-bam")
                .metavar("BAM")
                .type(String.class)
                .required(true)
                .help("file path of primary sample bam file (required)");
        parser.addArgument("-vcf")
                .metavar("VCF")
                .type(String.class)
                .required(true)
                .help("file path of variants of interest vcf file (required)");
        parser.addArgument("-bg", "--background_table")
        		.type(Arguments.fileType().acceptSystemIn().verifyCanRead())
        		.dest("bg")
        		.setDefault("-")
                .required(true)
                .help("file path to file containing list of background bams," +
                		" background files are identified based on order" +
                		" in this file. (required, - = stdin)");
        parser.addArgument("-f", "--folder")
                .type(String.class)
        		.dest("folder")
        		.setDefault("./varbin_results")
                .required(false)
                .help("file path string for output (default = ./varbin_results)");
        parser.addArgument("-o", "--output_table")
                .type(String.class)
        		.dest("out")
        		.setDefault("varbin_output.txt")
                .help("file name for output table" +
                		" (default = varbin_output.txt)");
        parser.addArgument("-gatk", "--gatk")
                .type(String.class)
        		.dest("gatk")
        		.setDefault("/Usr/local/bin/GenomeAnalysisTK.jar")
                .help("path to GenomeAnalysisTK.jar" +
                		" (default = /Usr/local/bin/GenomeAnalysisTK.jar)");
        parser.addArgument("-r", "--reference")
                .type(String.class)
        		.dest("ref")
                .help("path to reference sequence");
        parser.addArgument("--gatk_threads")
                .type(String.class)
        		.dest("threads")
        		.setDefault("4")
                .help("number of threads used by GATK UnifiedGenotyper" +
                		" (default = 4)");
        parser.addArgument("--gatk_mem_options")
                .type(String.class)
        		.dest("mem")
        		.setDefault("Xmx2g")
                .help("java memory option used by GATK UnifiedGenotyper." +
                		"Note that the dash (-) is not included to avoid" +
                		"shell parsing issues." +
                		" (default = Xmx2g)");
		try {
			args= parser.parseArgs(arguments);
			} catch (ArgumentParserException e) {
				parser.handleError(e);
			}
		
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

		//split main input vcf into 4 (snv alt allele 1, snv alt allele 2, indel alt allele 1, indel alt allele 2)
	    //this avoids any potential issues with how second alt alleles are handled or annotated
	    //each of these component vcf files will be used as allele targets for BAM variant calling
		mainVcfParts = new VcfSplit((String) args.get("vcf"), (String) args.get("folder"));
		
		System.out.println("");
		System.out.println("--Component vcf variant counts--");
		System.out.println("snv1Count: " + mainVcfParts.snv1Count);
		System.out.println("snv2Count: " + mainVcfParts.snv2Count);
		System.out.println("indel1Count: " + mainVcfParts.indel1Count);
		System.out.println("indel2Count: " + mainVcfParts.indel2Count);

		//Read in List of background bam file paths
		bamList = ReadBamList((File) args.get("bg"));

		//for each of the 4 component vcfs
		//if the component vcf contains at least one variant:
		//1. recall component vcf variants in primary bam with GATK
		//2. convert the new vcf output to variant list object
		//3. recall component vcf variants in background bams with GATK
		//2. convert the background vcf outputs to variant list objects
		
		//--for first alt allele SNVs--
		if (mainVcfParts.snv1Count > 0) {
			System.out.println("Starting main and background BAM file GATK calls for first alt allele, snv variants ("
						+ args.get("bam") + ")...");
			processedSnv1 = new ProcessedCalls(mainVcfParts.snv1Path, VariantType.SNP, (String) args.get("bam"), bamList, args);
		} else {
			System.out.println("There were no first alt allele, snv variants in the input vcf file");
		}
		
		//--for second alt allele SNVs--
		if (mainVcfParts.snv2Count > 0) {
			System.out.println("Starting main and background BAM file GATK calls for second alt allele, snv variants ("
						+ args.get("bam") + ")...");
			processedSnv2 = new ProcessedCalls(mainVcfParts.snv2Path, VariantType.SNP, (String) args.get("bam"), bamList, args);
		} else {
			System.out.println("There were no second alt allele, snv variants in the input vcf file");
		}

		//--for first alt allele Indels--
		if (mainVcfParts.indel1Count > 0) {
			System.out.println("Starting main and background BAM file GATK calls for first alt allele, indel variants ("
						+ args.get("bam") + ")...");
			processedIndel1 = new ProcessedCalls(mainVcfParts.indel1Path, VariantType.INDEL, (String) args.get("bam"), bamList, args);
		} else {
			System.out.println("There were no first alt allele, indel variants in the input vcf file");
		}
		
		//--for second alt allele Indels--
		if (mainVcfParts.indel2Count > 0) {
			System.out.println("Starting main and background BAM file GATK calls for second alt allele, indel variants ("
						+ args.get("bam") + ")...");
			processedIndel2 = new ProcessedCalls(mainVcfParts.indel2Path, VariantType.INDEL, (String) args.get("bam"), bamList, args);
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
				outWriter.write("");
				outWriter.write("No snv first alt alleles were processed\n");
				outWriter.write("");
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
				outWriter.write("");
				outWriter.write("No snv second alt alleles were processed\n");
				outWriter.write("");
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
				outWriter.write("");
				outWriter.write("No indel first alt alleles were processed\n");
				outWriter.write("");
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
				outWriter.write("");
				outWriter.write("No indel second alt alleles were processed\n");
				outWriter.write("");
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
	
	private static ArrayList<String> ReadBamList(File bamListFile) {
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