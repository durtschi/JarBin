package varbin;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

public class VarList {

	public ArrayList<VarObj> vars;
	public String vcfIn;
	public Boolean fromVcf;

	public VarList() {
	}
	
	public void makeList(String vcf) {
		this.fromVcf = true;
		this.vcfIn = vcf;
		this.vars = new ArrayList<VarObj>();
		privateMakeList();
	}
		
	public void makeList(ArrayList<VarObj> listOfVars) {
		this.fromVcf = false;
		this.vcfIn = "none";
		this.vars = listOfVars;
	}

	private void privateMakeList() {
		
		VarObj oneVar;
		String line;
		BufferedReader vcfReader = null;
		try {
			vcfReader = new BufferedReader(new FileReader(this.vcfIn));
		} catch (FileNotFoundException e1) {
			System.err.println("ERROR: temp vcf file not found");
			e1.printStackTrace();
		}

		//convert vcf to list of variant objects
		try {
			while ((line = vcfReader.readLine()) != null) {
				if( ! line.startsWith("#") && ! line.isEmpty() && ! line.equals("")){ //ignore header, metadata, or empty line
					//System.out.println("here is the variant line --->   " + line);
					oneVar = new VarObj();
					oneVar.parse(line);
					vars.add(oneVar);
				}
			}
		} catch (IOException e) {
				System.err.println("ERROR: reading from file " + this.vcfIn);
				e.printStackTrace();
		}
			
		try {
			vcfReader.close();
		} catch (IOException e) {
			System.err.println("ERROR: closing file " + this.vcfIn);
			e.printStackTrace();
		}
	}
}
