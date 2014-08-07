package varbin;

import java.awt.List;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

public class ProcessedCalls {
	public String vcfPath;
	public String mainBamPath;
	public ArrayList<String> backgroundBams;
	public VariantType varType;

	public UnifiedGenotyperCaller mainCall = null;
	public VarList mainVars = null;
	public ArrayList<VarList> backgroundVars = null;
	
	
	
	public ProcessedCalls(String vcfPath, VariantType varType, String mainBamPath, ArrayList<String> backgroundBams) {
		UnifiedGenotyperCaller gatkCall = null;
		
		this.vcfPath = vcfPath;
		this.mainBamPath = mainBamPath;
		this.backgroundBams = backgroundBams;

		this.mainCall = new UnifiedGenotyperCaller(
				vcfPath, mainBamPath,
				varType);
		this.mainVars = new VarList(this.mainCall.vcfOut);
		
		this.backgroundVars = new ArrayList<VarList>();
		for (String bamPath : backgroundBams) {
			System.out.println("calling variant with GATK for background BAM file  "
					+ bamPath);
			gatkCall = new UnifiedGenotyperCaller(
				vcfPath, bamPath,
				varType);
				this.backgroundVars.add(new VarList(gatkCall.vcfOut));
		}
	}
}
