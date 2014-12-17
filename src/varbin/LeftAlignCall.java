package varbin;
//imports
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import net.sourceforge.argparse4j.inf.Namespace;

 

public class LeftAlignCall {
	//decs
	public String vcfIn;
	public String vcfOut;
	public List<String> call;
	public String stderrOut;
	public String gatkPath; // = "/Users/jacobdurtschi/Tools/GenomeAnalysisTK.jar";
	public String refPath; // = "/Users/jacobdurtschi/Reference/Genome/human_g1k_v37.fasta";
	public String defaultMemOptions; // = "-Xmx2g";
	public static final String loggingLevel = "INFO";
	public static final String strictness = "LENIENT";

	public LeftAlignCall(String vcfIn, String vcfOut, Namespace args) {
		this.vcfIn = vcfIn;
		this.vcfOut = vcfOut;
		this.gatkPath = args.getString("gatk");
		this.refPath = args.getString("ref");
		this.stderrOut = "";
		this.call = null;
	}

	public void makeCall(){
		this.call = buildCall();
		callGatk(this.call);
	}
	
	private List<String> buildCall(){
		String javaCall = "java -Xmx2g -jar " + this.gatkPath
				+ " -R " + this.refPath
				+ " -T LeftAlignAndTrimVariants "
				+ " --validation_strictness " + strictness
				+ " --unsafe LENIENT_VCF_PROCESSING "
				+ " --logging_level " + loggingLevel
				+ " --splitMultiallelics "
				+ " --trimAlleles "
				+ " -V " + this.vcfIn
				+ " -o " + this.vcfOut;
		
		List<String> call = new ArrayList<String>(Arrays.asList("/bin/bash", "-c", javaCall));
		
		return call;
	}

	private void callGatk(List<String> call) {
    	// launch the command and grab stdin/stdout and stderr
    	ProcessBuilder pb = new ProcessBuilder(call);
    	Process p = null;
		String line;
		BufferedReader stderr = null;

    	
    	
		try {
			p = pb.start();
			stderr = new BufferedReader(new InputStreamReader(p.getErrorStream()));
		} catch (IOException e) {
			System.err.println("ERROR: calling gatk leftalign");
			e.printStackTrace();
		}
		// get process stderr
		try {
		while ((line = stderr.readLine()) != null) {
		  this.stderrOut = this.stderrOut + line + "\n";
		  System.out.println(line + "\n");
		}
		  stderr.close();
		  p.destroy();
		} catch (Exception e1) {
			System.out.println("ERROR: reading call std error and closing");
			e1.printStackTrace();
		}	
	}
}	