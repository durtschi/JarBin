package varbin;

//imports
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.lang3.StringUtils;

import net.sourceforge.argparse4j.inf.Namespace;

public class MpileupCallPiped {
	// decs
	public static List<String> mpileupCall;
	public static List<String> bcftoolsCall;
	public static String bamPath;
	public static String samtoolsPath; // assume on PATH;
	public static String bcftoolsPath; // assume on PATH;
	public static String refPath; // required
									// "/Users/jacobdurtschi/Reference/Genome/human_g1k_v37.fasta";
	public static String vcfOut;
	public static String covOut;
	public static String stderrOut;
	public static final String rawBaseQualityScore = "1";
	public static final String rawMapQualityScore = "5";
	public static final String minBaseQualityScore = "1";
	public static final String minMapQualityScore = "5";

	/**
	 * Constructor MpileupCall
	 * 
	 * @param bamPath
	 *            file path to background bam
	 * @param args
	 *            argparse Namespase with user options
	 */
	private MpileupCallPiped() {
		throw new AssertionError();
	}

	public static void makeCall(Namespace args) {
		MpileupCallPiped.mpileupCall = null;
		MpileupCallPiped.bcftoolsCall = null;
		MpileupCallPiped.bamPath = args.getString("bam");
		MpileupCallPiped.refPath = args.getString("ref");
		MpileupCallPiped.vcfOut = args.getString("folder") + "/"
				+ args.getString("prefix") + ".alleles.vcf";
		MpileupCallPiped.covOut = args.getString("folder") + "/"
				+ args.getString("prefix") + ".cov.txt";
		MpileupCallPiped.stderrOut = "";
		if (args.getString("samtools").equals("")) {
			ExecPath samPath = new ExecPath("samtools");
			MpileupCallPiped.samtoolsPath =  samPath.findPath();
		} else {
			MpileupCallPiped.samtoolsPath =  args.getString("samtools");
		}
		if (args.getString("bcftools").equals("")) {
			ExecPath bcfPath = new ExecPath("bcftools");
			MpileupCallPiped.bcftoolsPath =  bcfPath.findPath();
		} else {
			MpileupCallPiped.bcftoolsPath =  args.getString("bcftools");
		}

		BufferedReader samOut;
		String line;
		VarObj parsed = new VarObj();
		File vcfFile = new File(vcfOut);
		File covFile = new File(covOut);
		BufferedWriter vcfWriter = null;
		BufferedWriter covWriter = null;

		mpileupCall = buildMpileupCall(samtoolsPath, refPath, bamPath);
		bcftoolsCall = buildBcftoolsCall(bcftoolsPath);
		samOut = callSamtools(mpileupCall, bcftoolsCall);

		try {
			vcfWriter = new BufferedWriter(new FileWriter(vcfFile));
		} catch (IOException e) {
			e.printStackTrace();
		}

		try {
			covWriter = new BufferedWriter(new FileWriter(covFile));
		} catch (IOException e) {
			e.printStackTrace();
		}

		try {
			while ((line = samOut.readLine()) != null) {
				if (line.startsWith("#")) {
					vcfWriter.write(line);
					System.out.println(line);
				} else {
					parsed.parse(line); //no "new" object, just keep using same one over and over (parsed)
					//write coverage to cov file
					covWriter.write(parsed.chrom
									+ "\t"
									+ parsed.position
									+ "\t"
									+ parsed.infoMap.get("DP")
									+ "\t"
									+ parsed.formatMap.get("DP")
									);

					if (parsed.altList.size() == 1) {
						if (!parsed.altList.get(0).equals("X")) {
							System.out.println(line);
						} else {
							if (parsed.altList.get(parsed.altList.size() - 1)
									.equals("X")) {
								parsed.altList
										.remove(parsed.altList.size() - 1);
								// write line to the intermediate, allele vcf
								vcfWriter
										.write(parsed.chrom
												+ "\t"
												+ parsed.position
												+ "\t"
												+ "."
												+ "\t"
												+ StringUtils.join(
														parsed.refList, ",")
												+ "\t"
												+ StringUtils.join(
														parsed.altList, ",")
												+ "\t"
												+ parsed.qual
												+ "\t"
												+ parsed.filter
												+ "\t"
												+ "."
												+ "\t"
												+ "GT"
												+ "\t"
												+ parsed.genotype
												);
							}
						}
					}
					// all other lines, parse, get chr, pos, covRaw, covQual
					// alt X lines, to bucket
					// alt other, remove X send to vcfOut parsed.parse(line);

				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		try {
			vcfWriter.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		try {
			covWriter.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static List<String> buildMpileupCall(String samtools, String ref,
			String bam) {

		List<String> call = new ArrayList<String>(Arrays.asList(samtools,
				"mpileup", "-B", "-C", "0", "-d", "1000", "-q", "10", "-Q",
				"20", "-e", "10", "-gu", "-f", ref, bam));

		return call;
	}

	private static List<String> buildBcftoolsCall(String bcftoolsPath) {
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

		// launch the command and grab stdin/stdout and stderr
		ProcessBuilder pbMpileup = new ProcessBuilder(mpileupCall);
		Process processMpileup = null;
		ProcessBuilder pbBcftools = new ProcessBuilder(bcfCall);
		Process processBcftools = null;

		try {
			System.out.println(mpileupCall.toString());
			System.out.println(pbMpileup.command().toString());
			processMpileup = pbMpileup.start();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		InputStream stdoutMpileup = processMpileup.getInputStream();
		System.out.println("I'm here");
		try {
			System.out.println(bcfCall.toString());
			System.out.println(pbBcftools.command().toString());
			processBcftools = pbBcftools.start();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.out.println("I'm here");
		OutputStream stdinBcftools = processBcftools.getOutputStream();
		// BufferedReader stderr = new BufferedReader(new InputStreamReader(
		// processBcftools.getErrorStream()));

		Pipe pipe = new Pipe(stdoutMpileup, stdinBcftools);
		System.out.println("I'm here");
		pipe.run();
		System.out.println("I'm here");
		//new Thread(pipe).start();
		System.out.println("I'm here");


		BufferedReader stdout = new BufferedReader(new InputStreamReader(
				processBcftools.getInputStream()));

		//TODO remove the following, just a test
		String line;
		try {
			while ((line = stdout.readLine()) != null) {
				System.out.println(line);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.out.println("I'm here");

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
