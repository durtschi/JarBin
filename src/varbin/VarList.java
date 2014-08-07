package varbin;

import java.util.ArrayList;
import java.util.Arrays;

public class VarList{

	public ArrayList<VarObj> vars;
	public VarObj oneVar;

	public VarList(String vcf){

		//initialize an empty ArrayList of variant objects
		vars = new ArrayList<VarObj>();
		//convert vcf string to list of variant objects
		ArrayList<String> vcfList = new ArrayList<String>(Arrays.asList(vcf.split("\n")));
		for(String line : vcfList) {
		    if( ! line.startsWith("#") && ! line.isEmpty() && ! line.equals("")){ //ignore header, metadata, or empty line
		    	vars.add(new VarObj(line));
		    }
		}
	}
}
