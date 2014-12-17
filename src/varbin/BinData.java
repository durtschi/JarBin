package varbin;

import java.util.ArrayList;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

public class BinData {
				//print chr, pos, ref, alt, bin, variantType, GT, mainPLRD, median, 3sigma, 4sigma, 5sigma, 6sigma, 
	public ArrayList<String> chrom;
	public ArrayList<String> position;
	public ArrayList<String> ref;
	public ArrayList<String> alt;
	public ArrayList<VariantType> variantType;
	public ArrayList<Genotype> genotype;
	public ArrayList<Double> mainPLRD;
	public ArrayList<Integer> failCount;

	public ArrayList<Double> median;
	public ArrayList<Double> q1;
	public ArrayList<Double> q3;
	public ArrayList<Double> iqr;
	public ArrayList<Double> sigma;
	public ArrayList<Double> sigmas2main;
	public ArrayList<Double> sigma3;
	public ArrayList<Double> sigma4;
	public ArrayList<Double> sigma5;
	public ArrayList<Double> sigma6;
	public ArrayList<String> bin;
	
	public BinData(ProcessedCalls processedCalls) {
		this.chrom = new ArrayList<String>(); 
		this.position = new ArrayList<String>();
		this.ref = new ArrayList<String>();
		this.alt = new ArrayList<String>();
		this.variantType = new ArrayList<VariantType>();
		this.genotype = new ArrayList<Genotype>();
		this.mainPLRD = new ArrayList<Double>();
		this.failCount = new ArrayList<Integer>();
		
		this.median = new ArrayList<Double>();
		this.q1 = new ArrayList<Double>();
		this.q3 = new ArrayList<Double>();
		this.iqr = new ArrayList<Double>();
		this.sigma = new ArrayList<Double>();
		this.sigmas2main = new ArrayList<Double>();
		this.sigma3 = new ArrayList<Double>();
		this.sigma4 = new ArrayList<Double>();
		this.sigma5 = new ArrayList<Double>();
		this.sigma6 = new ArrayList<Double>();
		this.bin = new ArrayList<String>();

		ArrayList<Boolean>  passFilterList = new ArrayList<Boolean>();
		ArrayList<Double>   plrdList = new ArrayList<Double>();
		ArrayList<Double>   plrdFailList = new ArrayList<Double>();
		ArrayList<Integer>   whichFailList = new ArrayList<Integer>();
		ArrayList<Genotype> genotypeList = new ArrayList<Genotype>();

		for (int i = 0; i < processedCalls.mainVars.vars.size(); i++) { //loop through each allele
			//these are reset for each allele
			passFilterList.clear();
			plrdList.clear();
			plrdFailList.clear();
			whichFailList.clear();
			genotypeList.clear();

			// Gather useful info into lists
			this.chrom.add(processedCalls.mainVars.vars.get(i).chrom);
			this.position.add(processedCalls.mainVars.vars.get(i).position);
			this.ref.add(processedCalls.mainVars.vars.get(i).refList.get(0));
			this.alt.add(processedCalls.mainVars.vars.get(i).altList.get(0));
			this.variantType.add(processedCalls.mainVars.vars.get(i).varType);
			this.genotype.add(processedCalls.mainVars.vars.get(i).genotype);
			this.mainPLRD.add(processedCalls.mainVars.vars.get(i).plrd);
			
			//collect data for one allele in all background samples
			for (int j = 0; j < processedCalls.backgroundVars.size(); j++) {
				passFilterList
						.add(processedCalls.backgroundVars.get(j).vars.get(i).passFilter);
				plrdList.add(processedCalls.backgroundVars.get(j).vars.get(i).plrd);
				genotypeList
						.add(processedCalls.backgroundVars.get(j).vars.get(i).genotype);
				if (!processedCalls.backgroundVars.get(j).vars.get(i).passFilter
						&& !processedCalls.backgroundVars.get(j).vars.get(i).genotype
								.equals(Genotype.NOCALL)
						&& !processedCalls.backgroundVars.get(j).vars.get(i).genotype
								.equals(Genotype.HOM)
						&& processedCalls.backgroundVars.get(j).vars.get(i).plrd != null
						&& processedCalls.mainVars.vars.get(i).plrd != null) {
					plrdFailList
							.add(processedCalls.backgroundVars.get(j).vars.get(i).plrd);
					whichFailList.add(j);
				}
			}

			// make count of failed filter samples for this allele
			int failCountVal = plrdFailList.size();
			this.failCount.add(failCountVal);
			
			// Calc bin stats
			if (plrdFailList.size() >= 3) {
				DescriptiveStatistics plrdStats = new DescriptiveStatistics();
				for (double val : plrdFailList) {
					plrdStats.addValue(val);
				}
				// make first quartile list
				double q1val = plrdStats.getPercentile(25.0);
				this.q1.add(q1val);
				// make third quartile list
				double q3val = plrdStats.getPercentile(75.0);
				this.q3.add(q3val);
				// make iqr = q3 -q1 list
				double iqrVal = q3val - q1val;
				this.iqr.add(iqrVal);
				// make surogate sigma (iqr/1.35) list
				double sigmaVal = iqrVal / 1.35;
				this.sigma.add(sigmaVal);
				// make median list
				double medianVal = plrdStats.getPercentile(50.0);
				this.median.add(medianVal);
				// make 3sigma list 
				double sigma3Val = medianVal + 3*sigmaVal;
				this.sigma3.add(sigma3Val);
				// make 4sigma list 
				double sigma4Val = medianVal + 4*sigmaVal;
				this.sigma4.add(sigma4Val);
				// make 5sigma list 
				double sigma5Val = medianVal + 5*sigmaVal;
				this.sigma5.add(sigma5Val);
				// make 6sigma list 
				double sigma6Val = medianVal + 6*sigmaVal;
				this.sigma6.add(sigma6Val);
				// make sigmas to main sample ( (mainPLRD - median)/sigma )
				// list
				double sigmas2mainVal = (processedCalls.mainVars.vars.get(i).plrd - medianVal)
						/ sigmaVal;
				this.sigmas2main.add(sigmas2mainVal);

				// use above values for bin heuristic
				if (processedCalls.mainVars.vars.get(i).genotype.equals(Genotype.HOM) 
						|| (processedCalls.mainVars.vars.get(i).plrd >= 10.0
						&& processedCalls.mainVars.vars.get(i).plrd >= (medianVal + 6*sigmaVal))) {
					this.bin.add("1");
				} else if (processedCalls.mainVars.vars.get(i).plrd >= (medianVal + 6*sigmaVal)
						|| processedCalls.mainVars.vars.get(i).plrd >= 10.0
						&& processedCalls.mainVars.vars.get(i).plrd >= (medianVal + 3*sigmaVal)) {
					this.bin.add("2");
				} else if (processedCalls.mainVars.vars.get(i).plrd >= (medianVal + 3*sigmaVal)) {
					this.bin.add("3");
				} else {
					this.bin.add("4");
				};
				
			} else {
				this.median.add(null);
				this.q1.add(null);
				this.q3.add(null);
				this.iqr.add(null);
				this.sigma.add(null);
				this.sigma3.add(null);
				this.sigma4.add(null);
				this.sigma5.add(null);
				this.sigma6.add(null);
				this.sigmas2main.add(null);
				this.bin.add("?");
			}
		}
	}
}
