package varbin;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

// TODO: Auto-generated Javadoc
/**
 * The Class VcfSplit.
 * Takes input vcf and breaks into 4 vcfs
 * SNVs in first alt allele position
 * SNVs in second alt allele position
 * Indels in first alt allele position
 * Indels in second alt allele position
 *
 * Writes these split vcfs to files for later
 * forced allele calls to a variant caller
 * GATK UnifiedGenotyper
 * (Preferable would be to send these vcfs
 * to GATK as stdin but this is not supported
 * by GATK)
 * 
 * @author jacobdurtschi
 */
public class VcfSplit {
	
	/** vcf paths and counts for SNVs and Indels in first and second alt allele position. */
	public String snv1Path;
	public String snv2Path;
	public String indel1Path;
	public String indel2Path;
	public int snv1Count;
	public int snv2Count;
	public int indel1Count;
	public int indel2Count;
	public String inputVcf;
	
	/**
	 * Instantiates a new vcf split.
	 *
	 * @param vcfInPath the vcf in path
	 */
	public VcfSplit(String vcfInPath, String outFolder) {
		this.inputVcf = vcfInPath;
		this.snv1Count = 0;
		this.snv2Count = 0;
		this.indel1Count = 0;
		this.indel2Count = 0;
		
		// sort vcf lines into 4 buffers: (header into all), snv first alt, 
		//          snv second alt , indel first alt, indel second alt
		// Do this because gatk does not like multiple alt alleles for some processes.
		this.snv1Path = outFolder + "/" + "snvAlt1.vcf";
		this.snv2Path = outFolder + "/" + "snvAlt2.vcf";
		this.indel1Path = outFolder + "/" + "indelAlt1.vcf";
		this.indel2Path = outFolder + "/" + "indelAlt2.vcf";
	}
		
		
	public void makeCall() {
	
		BufferedReader br = null;
		BufferedWriter snv1Writer = null;
		BufferedWriter snv2Writer = null;
		BufferedWriter indel1Writer = null;
		BufferedWriter indel2Writer = null;
		String[] lineList;
		String[] altList;
		String[] refList;
		String newVcfLine;
		String chrom;
		String position;

		this.snv1Count = 0;
		this.snv2Count = 0;
		this.indel1Count = 0;
		this.indel2Count = 0;
		
		try{
			snv1Writer = new BufferedWriter( new FileWriter( this.snv1Path));
			snv2Writer = new BufferedWriter( new FileWriter( this.snv2Path));
			indel1Writer = new BufferedWriter( new FileWriter( this.indel1Path));
			indel2Writer = new BufferedWriter( new FileWriter( this.indel2Path));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			System.err.println("ERROR: creating component vcf files\n" + this.snv1Path + "\n" + this.snv2Path + "\n" + this.indel1Path + "\n" + this.indel2Path + "\n");
			e.printStackTrace();
			System.exit(1);
		}
		
		try{
			br = new BufferedReader(new FileReader(this.inputVcf));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			System.err.println("ERROR: accessing/reading input vcf file " + this.inputVcf);
			e.printStackTrace();
			System.exit(1);
		}

		try {
			// process each line of input vcf
			for(String line; (line = br.readLine()) != null; ) {
			    //if first char is # send to all files (private headerLine)
			    if(line.startsWith("##")){ //metadata line
			    	//System.out.println("This is a metadata line.");
			    	snv1Writer.write(line + "\n");
			    	snv2Writer.write(line + "\n");
			    	indel1Writer.write(line + "\n");
			    	indel2Writer.write(line + "\n");
			    //table header line
			    }else if(line.startsWith("#CHROM")){ 
			    	String newHeaderLine = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample";
			    	//System.out.println("This is the header line.");
			    	snv1Writer.write(newHeaderLine + "\n");
			    	snv2Writer.write(newHeaderLine + "\n");
			    	indel1Writer.write(newHeaderLine + "\n");
			    	indel2Writer.write(newHeaderLine + "\n");
			    //else variant line
			    }else{
			    	lineList = line.split("\t");
			    	chrom = lineList[0];
			    	position = lineList[1];
			    	refList = lineList[3].split(",");
			    	altList = lineList[4].split(",");

			    	if((refList.length  == 1) & (altList.length == 1)){ //single alt allele variant
			    		if(refList[0].length() == 1 & altList[0].length() == 1){ //single SNV variant
			    			//System.out.println("This is a single alt snv line.");
			    			newVcfLine = chrom + "\t" + position + "\t.\t" + refList[0] + "\t" + altList[0] + "\t.\tPASS\t.\tGT\t./.";
			    			snv1Writer.write(newVcfLine + "\n");
			    			this.snv1Count += 1;
			    		}else if(refList[0].length() > 1 | altList[0].length() > 1){ //single Indel variant
			    			//System.out.println("This is a single alt indel line.");
			    			newVcfLine = chrom + "\t" + position + "\t.\t" + refList[0] + "\t" + altList[0] + "\t.\tPASS\t.\tGT\t./.";
			    			indel1Writer.write(newVcfLine + "\n");
			    			this.indel1Count += 1;
			    		}
			    	}
			    	if(altList.length > 1){ //multi alt allele variant
			    		if(refList[0].length() == 1 & altList[0].length() == 1){ //first alt allele is SNV variant
			    			//System.out.println("First alt allele is an SNV.");
			    			newVcfLine = chrom + "\t" + position + "\t.\t" + refList[0] + "\t" + altList[0] + "\t.\tPASS\t.\tGT\t./.";
			    			snv1Writer.write(newVcfLine + "\n");
			    			this.snv1Count += 1;
			    		}else if(refList[0].length() > 1 | altList[0].length() > 1){ //first alt allele is Indel variant
			    			//System.out.println("First alt allele is an indel.");
			    			newVcfLine = chrom + "\t" + position + "\t.\t" + refList[0] + "\t" + altList[0] + "\t.\tPASS\t.\tGT\t./.";
			    			indel1Writer.write(newVcfLine + "\n");
			    			this.indel1Count += 1;
			    		}
			    		if(refList[0].length() == 1 & altList[1].length() == 1){ //second alt allele is SNV variant
			    			//System.out.println("Second alt allele is an SNV.");
			    			newVcfLine = chrom + "\t" + position + "\t.\t" + refList[0] + "\t" + altList[1] + "\t.\tPASS\t.\tGT\t./.";
			    			snv2Writer.write(newVcfLine + "\n");
			    			this.snv2Count += 1;
			    		}else if(refList[0].length() > 1 | altList[0].length() > 1){ //second alt allele is Indel variant
			    			//System.out.println("Second alt allele is an indel.");
			    			newVcfLine = chrom + "\t" + position + "\t.\t" + refList[0] + "\t" + altList[1] + "\t.\tPASS\t.\tGT\t./.";
			    			indel2Writer.write(newVcfLine + "\n");
			    			this.indel2Count += 1;
			    			//System.out.println();
			    		}
			    	}
			    }
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			System.err.println("ERROR: Problem writing component vcf files" + this.snv1Path);
			e.printStackTrace();
			System.exit(1);
		}

		try {
			snv1Writer.close();
			snv2Writer.close();
			indel1Writer.close();
			indel2Writer.close();
		} catch (IOException e) {
			System.err.println("ERROR: Problem closing component vcf files" + this.snv1Path);
			e.printStackTrace();
			System.exit(1);
		}
	}
}
