package phylo.tree.dist.diff.merge;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.function.Predicate;

import phylo.ref.Region;

public class MergerUtils {
	
	
	/**
	 * build and return a map from chrom name to the region filter to the chrom 
	 * @return
	 * @throws IOException 
	 */
	public static Map<String, Predicate<Region>> makeChromNamePredicateMap(Path chromContainingFile, int chromColIndex) throws IOException{
		Map<String, Predicate<Region>> ret = new LinkedHashMap<>();
		
		/////////////////
	    BufferedReader lineReader = new BufferedReader(new FileReader(chromContainingFile.toFile()));
	    String line = null;
	 
	    while ((line = lineReader.readLine()) != null) {
	    	String[] splits = line.split("\\s+");
	    	
	    	if(!ret.containsKey(splits[chromColIndex])) {
	    		ret.put(splits[chromColIndex], e->{return e.getReferenceName().equals(splits[chromColIndex]);});
	    	}
	    }
	 
	    lineReader.close();
		
	    ////////////////////
	    return ret;
	}
	
	
	
}
