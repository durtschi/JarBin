package varbin;

import java.util.HashMap;

public class FilterLimits {

    public final static HashMap<String, Float> snvFilters = new HashMap<String, Float>();
	static { 
        // map snv filter parameters
        snvFilters.put("QualMin", new Float(30.0));
        snvFilters.put("QDmin", new Float(2.0));
        snvFilters.put("MQmin", new Float(40.0));
        snvFilters.put("FSmax", new Float(60.0));
        snvFilters.put("HaplotypeScoreMax", new Float(13.0));
        snvFilters.put("MQRankSumMin", new Float(-12.5));
        snvFilters.put("ReadPosRankSumMin", new Float(-8.0));
	}

    public final static HashMap<String, Float> indelFilters = new HashMap<String, Float>();
	static { 
        // map indel filter parameters
        indelFilters.put("QualMin", new Float(30.0));
        indelFilters.put("QDmin", new Float(2.0));
        indelFilters.put("FSmax", new Float(200.0));
        indelFilters.put("ReadPosRankSumMin", new Float(-8.0));
	}
	
}
