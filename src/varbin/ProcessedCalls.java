package varbin;

import java.util.ArrayList;
import net.sourceforge.argparse4j.inf.Namespace;

public class ProcessedCalls {
	public String alleleVcfPath;
	public String outVcfPath;
	public String mainBamPath;
	public VariantType varType;

	public UnifiedGenotyperCall mainCall = null;
	public VarList mainVars = null;

	public ArrayList<String> backgroundVcfs;
	public VarsFromBkgd backgroundVars;
	
	
	
	public ProcessedCalls(String alleleVcfPath, VariantType varType, String mainBamPath, ArrayList<String> bkgdFiles, String outVcfPath, Namespace args) {
		
		this.alleleVcfPath = alleleVcfPath;
		this.outVcfPath = outVcfPath;
		this.mainBamPath = mainBamPath;
		this.backgroundVcfs = bkgdFiles;

		this.mainCall = new UnifiedGenotyperCall(
				alleleVcfPath, mainBamPath,
				varType, outVcfPath, args);
		this.mainCall.callGatk();
		this.mainVars = new VarList();
		this.mainVars.makeList(this.mainCall.vcfOut);
		
		this.backgroundVars = new VarsFromBkgd(this.backgroundVcfs, this.mainVars);
		this.backgroundVars.getBkgdVars();
		

		}
	}
