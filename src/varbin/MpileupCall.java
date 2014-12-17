package varbin;

//imports
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.lang3.StringUtils;

import net.sourceforge.argparse4j.inf.Namespace;

public class MpileupCall {
	// decs
	private static List<String> mpileupCall;
	private static List<String> bcftoolsCall;
	private static List<String> comboCall;
	private static String bamPath;
	private static String samtoolsPath; // assume on PATH;
	private static String bcftoolsPath; // assume on PATH;
	private static String refPath; // required
									// "/Users/jacobdurtschi/Reference/Genome/human_g1k_v37.fasta";
	private static File[] vcfFile;
	private static File refVcfFile;
	private static BufferedWriter[] vcfWriter;
	private static BufferedWriter refVcfWriter;
	private static ProcessBuilder pb;
	private static Process p;
	private static BufferedReader result;
	
	private static String[] vcfOut;
	private static String refVcfOut;
	private static String stderrOut;
	private static final String rawBaseQualityScore = "1";
	private static final String rawMapQualityScore = "5";
	private static final String minBaseQualityScore = "20";
	private static final String minMapQualityScore = "10";
	private static final int minDpQual = 3;

	// constructor private because utility class with only static methods.
	private MpileupCall() {
		throw new AssertionError();
	}

	public static String[] getVcfOut() {
		return vcfOut;
	}

	public static String getCovOut() {
		return refVcfOut;
	}

	public static void makeCall(Namespace args) {
		mpileupCall = null;
		bcftoolsCall = null;
		bamPath = args.getString("bam");
		refPath = args.getString("ref");
		vcfOut = new String[3];
		vcfOut[0] = args.getString("folder") + "/"
				+ args.getString("prefix") + ".alleles1.vcf";
		vcfOut[1] = args.getString("folder") + "/"
				+ args.getString("prefix") + ".alleles2.vcf";
		vcfOut[2] = args.getString("folder") + "/"
				+ args.getString("prefix") + ".alleles3.vcf";
		refVcfOut = args.getString("folder") + "/"
				+ args.getString("prefix") + ".ref.vcf";
		stderrOut = "";
		if (args.getString("samtools").equals("")) {
			ExecPath samPath = new ExecPath("samtools");
			samtoolsPath =  samPath.findPath();
		} else {
			samtoolsPath =  args.getString("samtools");
		}
		if (args.getString("bcftools").equals("")) {
			ExecPath bcfPath = new ExecPath("bcftools");
			bcftoolsPath =  bcfPath.findPath();
		} else {
			bcftoolsPath =  args.getString("bcftools");
		}

		vcfFile = new File[3];
		vcfFile[0] = new File(vcfOut[0]);
		vcfFile[1] = new File(vcfOut[1]);
		vcfFile[2] = new File(vcfOut[2]);
		refVcfFile = new File(refVcfOut);

		vcfWriter = new BufferedWriter[3];
		vcfWriter[0] = null;
		vcfWriter[1] = null;
		vcfWriter[2] = null;
		refVcfWriter = null;

		mpileupCall = buildMpileupCall();
		bcftoolsCall = buildBcftoolsCall();
		result = callSamtools(mpileupCall, bcftoolsCall);
		parse();
	}

	private static void parse() {
		try {
			vcfWriter[0] = new BufferedWriter(new FileWriter(vcfFile[0]));
			vcfWriter[1] = new BufferedWriter(new FileWriter(vcfFile[1]));
			vcfWriter[2] = new BufferedWriter(new FileWriter(vcfFile[2]));
		} catch (IOException e) {
			e.printStackTrace();
		}

		try {
			refVcfWriter = new BufferedWriter(new FileWriter(refVcfFile));
		} catch (IOException e) {
			e.printStackTrace();
		}

		String line;
		VarObj parsed = new VarObj();
		try {
			while ((line = result.readLine()) != null) {
				if (line.startsWith("#")) {
					vcfWriter[0].write(line + "\n");
					vcfWriter[1].write(line + "\n");
					vcfWriter[2].write(line + "\n");
					refVcfWriter.write(line + "\n");
					//System.out.println(line);
				} else {
					parsed.parse(line); //reuse same object (parsed)
					Integer dpQual = 0; //depth of quality filtered reads/bases
					if (parsed.infoMap.containsKey("I16")) {
						String[] i16Array = parsed.infoMap.get("I16").split(",");
						for (int i= 0 ; i < 4 ; i++ ) {
							dpQual += Integer.parseInt(i16Array[i]);
						}
					}
					

					if ((parsed.altList.size() == 1) && (parsed.altList.get(0).equals("X")) && (dpQual >= minDpQual)) {
						// this is only "X" for alt (means all reads are ref at this position)
						// write to reference only vcf file
						refVcfWriter.write(parsed.chrom + "\t"
								+ parsed.position + "\t" + "." + "\t"
								+ StringUtils.join(parsed.refList, ",")
								+ "\t"
								+ "X" // alt
								+ "\t" + parsed.qual + "\t" + parsed.filter
								+ "\t");
						List<String> refInfo = new ArrayList<String>();
						if (parsed.infoMap.containsKey("DP")) {
							refInfo.add("DP=" + parsed.infoMap.get("DP"));
						}
						if (parsed.infoMap.containsKey("I16")) {
							refInfo.add("DPQ=" + dpQual.toString());
							// refInfo.add("I16=" + parsed.infoMap.get("I16"));
						}
						if (refInfo.isEmpty()) {
							refVcfWriter.write(".");
						} else {
							refVcfWriter.write(StringUtils.join(refInfo, ";"));
						}
						refVcfWriter.write(parsed.chrom + "\t" + "PL" + "\t"
								+ parsed.formatMap.get("PL") + "\n");

					} else {
						if (parsed.altList.get(parsed.altList.size() - 1)
								.equals("X")) {
							parsed.altList
									.remove(parsed.altList.size() - 1);
						}
						
						// add new vcf line for each alt allele
						for (int i = 0 ; i < parsed.altList.size() ; i++) {
						vcfWriter[i]
									.write(parsed.chrom
											+ "\t"
											+ parsed.position
											+ "\t"
											+ "."
											+ "\t"
											+ StringUtils.join(
													parsed.refList, ",")
											+ "\t"
											// choose only first allele or all alleles or eachAlt
											+ parsed.altList.get(i).toString()
											//+ StringUtils.join(
											//		parsed.altList, ",")
											+ "\t"
											+ parsed.qual
											+ "\t"
											+ parsed.filter
											+ "\t"
											+ "."
											+ "\t"
											+ "PL"
											+ "\t"
											+ parsed.formatMap.get("PL")
											+ "\n"
											);
						}
						
					}
				}
					// all other lines, parse, get chr, pos, covRaw, covQual
					// alt X lines, to bucket
					// alt other, remove X send to vcfOut parsed.parse(line);

			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		try {
			vcfWriter[0].close();
			vcfWriter[1].close();
			vcfWriter[2].close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		try {
			refVcfWriter.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static List<String> buildMpileupCall() {

		List<String> call = new ArrayList<String>(Arrays.asList(samtoolsPath,
				"mpileup", "-B", "-C", "0", "-d", "1000", "-q", minMapQualityScore, "-Q",
				minBaseQualityScore, "-e", "10", "-gu", "-f", refPath, bamPath));

		return call;
	}

	private static List<String> buildBcftoolsCall() {
		// samtools mpileup -B -C 0 -d 1000 -f
		// ~/Reference/Genome/human_g1k_v37.fasta -q 10 -Q 20
		// jarbin_test_short_HHT12.bam -gu -L 1000 -e 10 \
		// | bcftools view -AN - \
		// | sed 's/,X / /' \
		// | sed 's/,[0-9]+:[0-9]+:[0-9]$//' \
		// | tee >(grep -v 'INDEL' > junk_snv.vcf) \
		// | grep 'INDEL\|^#' > junk_indel.vcf
		//
		// | grep -v "\tX\t" \
		List<String> call = new ArrayList<String>(Arrays.asList(bcftoolsPath,
				"view", "-A", "-N", "-"));

		return call;
	}

	private static BufferedReader callSamtools(List<String> mpileupCall,
			List<String> bcfCall) {
		BufferedReader stdout;
		String mpileupStr = StringUtils.join(mpileupCall, " ");
		String bcfStr = StringUtils.join(bcfCall, " ");
		List<String> fullCall = new ArrayList<String>(Arrays.asList("/bin/bash", "-c", mpileupStr + " | " + bcfStr));

		// launch the command and grab stdin/stdout and stderr
		pb = new ProcessBuilder(fullCall);
		p = null;

		try {
			System.out.println(fullCall.toString());
			System.out.println(pb.command().toString());
			p = pb.start();
		} catch (IOException e) {
			System.err.println("ERROR: starting system call process to samtools mpileup and bcftools");
			System.err.println(pb.command().toString());
			e.printStackTrace();
		}
		stdout = new BufferedReader(new InputStreamReader(
				p.getInputStream()));

		return stdout;
	}
}

/**
 * 
 * Standard GATK UnifiedGenotyper annotations Standard annotations in the list
 * below are marked with a '*'.
 * 
 * Available annotations for the VCF INFO field: AlleleBalance BaseCounts
 * BaseQualityRankSumTest ChromosomeCounts ClippingRankSumTest Coverage
 * FisherStrand GCContent HaplotypeScore HardyWeinberg HomopolymerRun
 * InbreedingCoeff LikelihoodRankSumTest LowMQ MVLikelihoodRatio
 * MappingQualityRankSumTest MappingQualityZero NBaseCount QualByDepth
 * RMSMappingQuality ReadPosRankSumTest SampleList SnpEff SpanningDeletions
 * StrandOddsRatio TandemRepeatAnnotator TransmissionDisequilibriumTest
 * VariantType
 * 
 * Available annotations for the VCF FORMAT field: AlleleBalanceBySample
 * DepthPerAlleleBySample DepthPerSampleHC MappingQualityZeroBySample
 * StrandBiasBySample
 * 
 * Available classes/groups of annotations: ActiveRegionBasedAnnotation
 * ExperimentalAnnotation RankSumTest RodRequiringAnnotation StandardAnnotation
 * WorkInProgressAnnotation
 **/
