package varbin;

import java.util.ArrayList;
import net.sourceforge.argparse4j.inf.Namespace;

public class ProcessedCalls {
	public String vcfPath;
	public String mainBamPath;
	public ArrayList<String> backgroundBams;
	public VariantType varType;

	public UnifiedGenotyperCaller mainCall = null;
	public VarList mainVars = null;
	public ArrayList<VarList> backgroundVars = null;
	
	
	
	public ProcessedCalls(String vcfPath, VariantType varType, String mainBamPath, ArrayList<String> backgroundBams, Namespace args) {
		UnifiedGenotyperCaller gatkCall = null;
		
		this.vcfPath = vcfPath;
		this.mainBamPath = mainBamPath;
		this.backgroundBams = backgroundBams;

		this.mainCall = new UnifiedGenotyperCaller(
				vcfPath, mainBamPath,
				varType, args);
		this.mainVars = new VarList(this.mainCall.vcfOut);
		
		this.backgroundVars = new ArrayList<VarList>();
		for (String bamPath : backgroundBams) {
			System.out.println("calling variant with GATK for background BAM file  "
					+ bamPath);
			gatkCall = new UnifiedGenotyperCaller(
				vcfPath, bamPath,
				varType, args);
				this.backgroundVars.add(new VarList(gatkCall.vcfOut));
		}
	}
}
