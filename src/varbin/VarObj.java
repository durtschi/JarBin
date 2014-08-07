package varbin;

import java.util.Map;

public class VarObj {
	
	//attributes directly from vcf
	private String[] varList;
	public String line;
	public String chrom;
	public String position;
	public String id;
	public String[] refList;
	public String[] altList;
	public String qual;
	public String filter;
    public DefaultHashMap<String,String> infoMap;
    public DefaultHashMap<String,String> formatMap;
    
    //attributes not pulled directly from vcf (calculated here)
    public Double plrd; //Likelihood ratio over depth value calculated from "PL" values from vcf
    public boolean passFilter; //Heuristic variant filter based on quality and bias values in the vcf
    public Genotype genotype; // enumerated type: WILDTYPE, HET, HOM, or NOCALL
    public VariantType varType; // enumerated type: SNP or INDEL
    
    //others
    public Map<String, Float> snvFilters;
    public Map<String, Float> indelFilters;
    
	public VarObj(String varString){

        // map snv filter parameters
        this.snvFilters = new DefaultHashMap<String, Float>();
        this.snvFilters.put("QualMin", new Float(30.0));
        this.snvFilters.put("QDmin", new Float(2.0));
        this.snvFilters.put("MQmin", new Float(40.0));
        this.snvFilters.put("FSmax", new Float(60.0));
        this.snvFilters.put("HaplotypeScoreMax", new Float(13.0));
        this.snvFilters.put("MQRankSumMin", new Float(-12.5));
        this.snvFilters.put("ReadPosRankSumMin", new Float(-8.0));
        // map indel filter parameters
        this.indelFilters = new DefaultHashMap<String, Float>();
        this.indelFilters.put("QualMin", new Float(30.0));
        this.indelFilters.put("QDmin", new Float(2.0));
        this.indelFilters.put("FSmax", new Float(200.0));
        this.indelFilters.put("ReadPosRankSumMin", new Float(-8.0));
        
        
		String[] infoList;
		String[] formatList;
		String[] sampleList;
	
		//parse vcf line
		line = varString.replaceAll("\n", "");
	    if( ! line.startsWith("#") && ! line.isEmpty() && ! line.equals("")){ //ignore blank, header, or metadata line
	    	this.varList = line.trim().split("\t");
	    	this.chrom = varList[0];
	    	this.position = varList[1];
	    	this.id = varList[2];
	    	this.refList = varList[3].split(",");
	    	this.altList = varList[4].split(",");
	    	this.qual = varList[5];
	    	this.filter = varList[6];
	    	infoList = varList[7].split(";");
	    	formatList = varList[8].split(":");
	    	sampleList = varList[9].split(":");

	    	//convert info field to map
	    	this.infoMap = new DefaultHashMap<String,String>();
	    	for (String infoItem : infoList) {
	    		String[] infoPair = infoItem.trim().split("=");
	    		if(infoPair.length == 2){
	    			this.infoMap.put(infoPair[0], infoPair[1]);
	    		}else if(infoPair.length == 1){
	    			this.infoMap.put(infoPair[0], "N/A");
	    		}else{
	    			System.err.println("ERROR: malformed info field in gatk output vcf file line\n" + line);
	    			System.exit(1);
	    		}
	    	}

	    	//convert format/sample fields to map
		    formatMap = new DefaultHashMap<String,String>();
	    	if(formatList.length == sampleList.length){
	    		for (int i=0; i<formatList.length; i++) {
	    			this.formatMap.put(formatList[i], sampleList[i]);
	    		}
	    	}else{
	    		System.err.println("ERROR: malformed format field in gatk output vcf file line\n" + line);
	    		System.exit(1);
	    	}
	    	
	    	//determine variant type
	    	if (refList[0].length() == 1 && altList[0].length() == 1){
	    		this.varType = VariantType.SNP;
	    	}else{
	    		this.varType = VariantType.INDEL;
	    	}

	    	//determine genotype (REF, HET, HOM, NOCALL)
	    	this.genotype = readGenotype(this.formatMap.get("GT"));

	    	//calculate likelihood ratio over depth values (plrd) for variant
	    	this.plrd = plrdCalc(this.formatMap.getDefault("PL", "0,0,0"), this.infoMap.getDefault("DP", "0"));

			}else{
				System.err.println("ERROR: malformed line in vcf \n" + line);
				System.exit(1);
			}
	    
	    	//determine passFilter 
	    	this.passFilter = PassFilter(this.snvFilters, this.indelFilters);
	}

	private Boolean PassFilter(Map<String, Float> snvFilters, Map<String, Float> indelFilters) {
		Boolean passFilter = null;
		
		// make sure we have all of the filter values
		if (!snvFilters.containsKey("QDmin")
				|| !snvFilters.containsKey("QualMin")
				|| !snvFilters.containsKey("MQmin")
				|| !snvFilters.containsKey("FSmax")
				|| !snvFilters.containsKey("HaplotypeScoreMax")
				|| !snvFilters.containsKey("MQRankSumMin")
				|| !snvFilters.containsKey("ReadPosRankSumMin")
				|| !indelFilters.containsKey("QualMin")
				|| !indelFilters.containsKey("QDmin")
				|| !indelFilters.containsKey("FSmax")
				|| !indelFilters.containsKey("ReadPosRankSumMin")
				) {
			System.err.println("ERROR: missing snv or indel filter values");
			System.exit(1);
		}
		if (this.qual.equals(".")) {
			this.qual = "0.0";
			
		}
		if (this.varType == VariantType.SNP) {
			if (Float.parseFloat(this.qual) < ((Float) snvFilters.get("QualMin")).floatValue() ||
				Float.parseFloat(this.infoMap.getDefault("QD", "0.0")) < ((Float) snvFilters.get("QDmin")).floatValue() ||
				Float.parseFloat(this.infoMap.getDefault("MQ", "0.0")) < ((Float) snvFilters.get("MQmin")).floatValue() ||
				Float.parseFloat(this.infoMap.getDefault("FS", "0.0")) > ((Float) snvFilters.get("FSmax")).floatValue() ||
				Float.parseFloat(this.infoMap.getDefault("HaplotypeScore", "0.0")) > ((Float) snvFilters.get("HaplotypeScoreMax")).floatValue() ||
				Float.parseFloat(this.infoMap.getDefault("MQRankSum", "0.0")) < ((Float) snvFilters.get("MQRankSumMin")).floatValue() ||
				Float.parseFloat(this.infoMap.getDefault("ReadPosRankSum", "0.0")) < ((Float) snvFilters.get("ReadPosRankSumMin")).floatValue()) {
				passFilter = false;
			} else {
				passFilter = true;
			}
		} else if (this.varType == VariantType.INDEL) {
			if (Float.parseFloat(this.qual) < ((Float) snvFilters.get("QualMin")).floatValue() ||
				Float.parseFloat(this.infoMap.getDefault("QD", "0.0")) < ((Float) indelFilters.get("QDmin")) ||
				Float.parseFloat(this.infoMap.getDefault("FS", "0.0")) > ((Float) indelFilters.get("FSmax")) ||
				Float.parseFloat(this.infoMap.getDefault("ReadPosRankSum", "0.0")) < ((Float) indelFilters.get("ReadPosRankSumMin"))) {
				passFilter = false;
			} else {
				passFilter = true;
			}
		}
		if (passFilter == null) {
    		System.err.println("ERROR: variant filter pass/fail not correctly determined\n");
    		System.exit(1);
		}
		return passFilter.booleanValue();
	}
	
	
    private Double plrdCalc(String plString, String depthString){
    	Double depth; //depth
    	Double plRef; // pl value for wildtype referene genotype (Phred scale)
    	Double plHet; // pl value for het genotype (Phred scale)
    	Double plHom; // pl value for hom genotype (Phred scale)
    	Double plHetLinear; // pl value for het genotype in linear scale rather than Phred scale
    	Double plHomLinear; // pl value for hom genotype in linear scale rather than Phred scale
    	Double plAlt; // pl value for hom and het genotype linearly combined (Phred scale)
    	Double plRatio; // pl ratio of ref/(het + hom) in Phred scale
    	Double plRatioOverDepth; // Phred scale pl ratio devided by depth
    	if (plString != null && depthString != null && Double.parseDouble(depthString) >= 3){
    		String[] plStringArray = plString.split(",");
    		depth = Double.parseDouble(depthString);
    		if (plStringArray.length != 3){
    			System.err.println("ERROR: malformed format field in gatk output vcf file line\n" + line);
				System.exit(1);
    		}
    		plRef = Double.parseDouble(plStringArray[0]);
    		plHet = Double.parseDouble(plStringArray[1]);
    		plHom = Double.parseDouble(plStringArray[2]);
    		depth = Double.parseDouble(depthString);
    		plHetLinear = Math.pow(10.0, -plHet/10.0);
    		plHomLinear = Math.pow(10.0, -plHom/10.0);
    		//need to make sure we don't attempt log(0) error
    		//this might happen if we have numbers smaller than double precision
    		//instead we will estimate an answer for plAlt based on log(really small number)
    		if (plHetLinear == 0.0 && plHomLinear == 0.0){
    			plAlt = 4000.0; //corresponds to log10(10e-400)
    		}else{
    			plAlt = -10.0 * Math.log10(plHetLinear + plHomLinear);
    		}
    		plRatio = plRef - plAlt; //same as PLR
    		plRatioOverDepth = plRatio / depth;
    	}else{
    		plRatioOverDepth = (Double) null ;
    	}

		return plRatioOverDepth;
    }
    
    private Genotype readGenotype(String gtString){
    	Genotype genotypeResult = null;
    	if (gtString.equals("0/0")){
    		genotypeResult = Genotype.WILDTYPE;
    	}else if(gtString.equals("0/1")){
    		genotypeResult = Genotype.HET;
    	}else if(gtString.equals("1/1")){
    		genotypeResult = Genotype.HOM;
    	}else if(gtString.equals("./.")){
    		genotypeResult = Genotype.NOCALL;
    	}
    		
    	return genotypeResult;	
    }
	
}
