package phylo.ref;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.Map;

public class RegionUtils {

	public static Map<Integer, String> readBedFileIntoRegionString(Path bedFile, int regionLen) throws IOException{
		Map<Integer,String> ret = new LinkedHashMap<>();
		
		BufferedReader lineReader = new BufferedReader(new FileReader(bedFile.toFile()));
	    String line = null;
	    
	    int index=1;
	    
	    while ((line = lineReader.readLine()) != null) {
	    	String[] splits = line.split("\\s+");
	    	String chrom=splits[0];
	    	int start=Integer.parseInt(splits[1]);
	    	int end=Integer.parseInt(splits[2]);
	    	
	    	int ri=0;
	    	while(start+(ri+1)*regionLen -1 <=end) {
	    		Region r = new Region(chrom, start+ri*regionLen, start+(ri+1)*regionLen-1);
	    		ret.put(index, r.toString());
	    		index++;
	    		ri++;
	    	}
	    	//process trailing region with length < regionLen
	    	if(start+ri*regionLen<end) {
	    		Region r = new Region(chrom, start+ri*regionLen, end);
	    		ret.put(index, r.toString());
	    		index++;
	    	}
	    	
	    }
	    
	    lineReader.close();
		
		
		return ret;
	}
	
	public static void main(String[] args) throws IOException {
		Path bedFile=Path.of("C:\\Users\\tanxu\\Desktop\\scratch\\bed\\Sbicolor_454_v3.0.1.1-10.chrom.bed");
		int len=1000000;
		
		Map<Integer, String> ret = readBedFileIntoRegionString(bedFile, len);
		
		ret.forEach((index, rs)->{
			System.out.println(index+"==="+rs);
		});
	}
	
}
