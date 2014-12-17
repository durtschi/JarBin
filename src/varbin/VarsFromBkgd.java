package varbin;

import htsjdk.tribble.readers.TabixReader;

import java.io.IOException;
import java.util.ArrayList;

public class VarsFromBkgd {

		ArrayList<VarList> bkgdMatches;
		ArrayList<String> bkgdVcfs;
		VarList mainVars;
		String REF_ONLY = "X";
		String NO_DATA = ".\t0\t.\t.\t.\t.\t.\tDP=0\tGT:DP\t0/0:0";

	public VarsFromBkgd(ArrayList<String> bkgdVcfs, VarList mainVars) {
		this.bkgdMatches = new ArrayList<VarList>();
		this.bkgdVcfs = bkgdVcfs;
		this.mainVars = mainVars;
	}
	
	public void getBkgdVars(){
		
		for (String vcf : bkgdVcfs) {
			//Perform the lookup
			this.bkgdMatches.add(matchBkgd(vcf, mainVars));
		}
		
		
	

			//TODO
			//use tabix to grab right hit if exists
			//if does exist, grab in some way
			//if does not exist grab default in some way
			//complete list of these grabbed items
			//return in correct format (I think this will be a VarList)
	}
	
	private VarList matchBkgd(String vcf, VarList mainVars) {
		
		TabixReader reader = null;
		TabixReader.Iterator iter = null;
		ArrayList<VarObj> varArray = new ArrayList<VarObj>();
			try {
				reader = new TabixReader(vcf);
			} catch (IOException e) {
				System.err.println("ERROR: opening bkgdFile for Tabix reading");
				e.printStackTrace();
			}
			for (VarObj var : mainVars.vars) {
				String queryStr = var.chrom + ":" + var.position + "-" + Integer.toString(Integer.parseInt(var.position) + 0);
				VarObj bkgdVarMatch = null;

				iter = reader.query(queryStr);
				if(iter != null) {
					try {
						String bkgdVarStr = iter.next();
						VarObj bkgdVar = new VarObj();	
						while(bkgdVarStr != null) {
							System.out.println("next bkgd var string");
							System.out.println(bkgdVarStr);
							bkgdVar.parse(bkgdVarStr);
							// if pos, ref, alt match up, then use this VarObj
							if (bkgdVar.refList.equals(var.refList)
									&& bkgdVar.altList.equals(var.altList)) {
								bkgdVarMatch = bkgdVar;
								break;
							// else if alt = REF_ONLY ("X") then keep this (only ref evidence found in alignment)
							} else if (bkgdVar.altList.equals(REF_ONLY)) {
								bkgdVarMatch = bkgdVar;
								break;
							}
							bkgdVarStr = iter.next();
						}
						// else, pass back nothing in some recognizable way
						if (bkgdVarMatch == null) {
							bkgdVar.parse(NO_DATA);
							bkgdVarMatch = bkgdVar;
						}
						varArray.add(bkgdVarMatch);
						System.out.println("the bkgdVar line is:");
						System.out.println(bkgdVarMatch.line);
					} catch (IOException e) {
					System.err.println("ERROR: during Tabix query");
					e.printStackTrace();
					} // try/catch
				} // if	
			}	
			VarList bkgdVarList = new VarList();
			bkgdVarList.makeList(varArray);
			return bkgdVarList;
	}
}
