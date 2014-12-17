package varbin;
//imports
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import net.sourceforge.argparse4j.inf.Namespace;
 

public class UnifiedGenotyperCall {
	//decs
	public String alleleVcf;
	public List<String> gatkCall;
	public String bamPath;
	public VariantType genotypeLikelyhoodsModel;
	public String vcfOut;
	public String stderrOut;
	public String gatkPath; // = "/Users/jacobdurtschi/Tools/GenomeAnalysisTK.jar";
	public String refPath; // = "/Users/jacobdurtschi/Reference/Genome/human_g1k_v37.fasta";
	public String threads; // = "4";
	public String defaultMemOptions; // = "-Xmx2g";
	public static final String standCallConf = "0.0";
	public static final String standEmitConf = "0.0";
	public static final String minIndelFrac = "0.00001";
	public static final String downsampleToCoverage = "100000";
	public static final String output = "/dev/stdout"; //"/Users/jacobdurtschi/Data/jarbin_test/jarbin_test_output.vcf";
	public static final String maxDeletionFraction = "1.0";
	public static final String minBaseQualityScore = "17";
	//public static final String alleles = "/Users/jacobdurtschi/Data/jarbin_test/jarbin_test_short.vcf";
	public static final String outputMode = "EMIT_ALL_SITES";
	public static final String loggingLevel = "INFO";



	/**
	 * @param alleleVcf
	 * @param bamPath
	 * @param genotypeLikelyhoodsModel
	 * @param args
	 */
	public UnifiedGenotyperCall(String alleleVcf, String bamPath, VariantType genotypeLikelyhoodsModel, String vcfOut, Namespace args){
		
		this.gatkPath = args.getString("gatk"); // = "/Users/jacobdurtschi/Tools/GenomeAnalysisTK.jar";
		this.refPath = args.getString("ref"); // = "/Users/jacobdurtschi/Reference/Genome/human_g1k_v37.fasta";
		this.threads = args.getString("threads"); // = "4";
		this.defaultMemOptions = args.getString("mem"); // = "-Xmx2g";
		this.vcfOut = vcfOut;
		this.stderrOut = "";
		this.alleleVcf = alleleVcf;
		this.bamPath = bamPath;
		this.genotypeLikelyhoodsModel = genotypeLikelyhoodsModel;
		this.gatkCall = buildGatkCall(bamPath, genotypeLikelyhoodsModel);
	}
	
	private List<String> buildGatkCall(String bamPath, VariantType genotypeLikelyhoodsModel){
		
		List<String> call = new ArrayList<String>(Arrays.asList("java",
				"-jar", gatkPath,
				"-T", "UnifiedGenotyper",
				"-R", refPath,
				"-I", bamPath,
				"-o", output,
				"--genotype_likelihoods_model", genotypeLikelyhoodsModel.toString(),
				"--genotyping_mode", "GENOTYPE_GIVEN_ALLELES",
				"--alleles", alleleVcf,
				"-L:alleles,VCF", alleleVcf,
				"--output_mode", outputMode,
				"-stand_call_conf", standCallConf,
				"-stand_emit_conf", standEmitConf,
				"-rf", "BadCigar",
				"-nt", threads,
				"-minIndelFrac", minIndelFrac,
				"--downsample_to_coverage", downsampleToCoverage,
				"--max_deletion_fraction", maxDeletionFraction,
				"--min_base_quality_score", minBaseQualityScore,
				"--logging_level", loggingLevel,
				//"--annotation", "ReadDepthAndAllelicFractionBySample",
				"--annotation", "AlleleBalanceBySample",
				"--annotation", "MappingQualityZeroBySample",
				//"--annotation", "MappingQualityZeroFraction",
				"--annotation", "AlleleBalance" ,
				"--annotation", "BaseCounts" ,
				"--annotation", "FisherStrand" ,
				"--annotation", "QualByDepth" ,
				"--annotation", "Coverage"));
		return call;
	}
		

	public void callGatk(){
		
		// make temp output file
		BufferedWriter bw = null;
		try {
			bw = new BufferedWriter(new FileWriter(this.vcfOut));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			System.err.println("ERROR: creating file");
			e.printStackTrace();
		}

	    try {
		      String line;
		      //OutputStream stdin = null;
		      //InputStream stderr = null;
		      //InputStream stdout = null;

		      // launch the command and grab stdin/stdout and stderr
		      ProcessBuilder pb = new ProcessBuilder(this.gatkCall);
		      Process process = pb.start();
		      BufferedReader stdout = new BufferedReader(new  InputStreamReader(process.getInputStream()));
		      BufferedReader stderr = new BufferedReader(new InputStreamReader(process.getErrorStream()));

		      // get process stdout
		      while ((line = stdout.readLine()) != null) {
		        bw.write(line + "\n");
		      }
		      stdout.close();
		      bw.flush();
		      bw.close();

		      // get process stderr
		      while ((line = stderr.readLine()) != null) {
		        this.stderrOut = this.stderrOut + line + "\n";
		        System.out.println(line + "\n");
		      }
		      stderr.close();
		      process.destroy();
		    } catch (Exception err) {
		      err.printStackTrace();
		    }	
	}
}

/**
 * 
 * Standard GATK UnifiedGenotyper annotations
 * Standard annotations in the list below are marked with a '*'.

Available annotations for the VCF INFO field:
	AlleleBalance
	BaseCounts
	*BaseQualityRankSumTest
	*ChromosomeCounts
	ClippingRankSumTest
	*Coverage
	*FisherStrand
	GCContent
	*HaplotypeScore
	HardyWeinberg
	HomopolymerRun
	*InbreedingCoeff
	LikelihoodRankSumTest
	LowMQ
	MVLikelihoodRatio
	*MappingQualityRankSumTest
	*MappingQualityZero
	NBaseCount
	*QualByDepth
	*RMSMappingQuality
	*ReadPosRankSumTest
	SampleList
	SnpEff
	*SpanningDeletions
	StrandOddsRatio
	*TandemRepeatAnnotator
	TransmissionDisequilibriumTest
	VariantType

Available annotations for the VCF FORMAT field:
	AlleleBalanceBySample
	*DepthPerAlleleBySample
	DepthPerSampleHC
	MappingQualityZeroBySample
	StrandBiasBySample

Available classes/groups of annotations:
	ActiveRegionBasedAnnotation
	ExperimentalAnnotation
	RankSumTest
	RodRequiringAnnotation
	StandardAnnotation
	WorkInProgressAnnotation
**/