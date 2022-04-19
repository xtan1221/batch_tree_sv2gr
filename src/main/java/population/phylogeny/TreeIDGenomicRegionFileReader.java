package population.phylogeny;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;
import genomics.utils.SimpleGenomicRegion;

/**
 * reader for a file containing the regional tree id and the genomic region 
 * 
 * 1	Chr01	1	1000000
 * 
 * @author tanxu
 *
 */
public class TreeIDGenomicRegionFileReader {
	private final Path treeIDRegionFile;
	
	////////////////////////
	private Map<Integer, SimpleGenomicRegion> treeIDRegionMap;
	
	/**
	 * 
	 * @param treeIDRegionFile
	 */
	public TreeIDGenomicRegionFileReader(Path treeIDRegionFile) {
		super();
		this.treeIDRegionFile = treeIDRegionFile;
		
		////////////////////////////////////////////
		this.run();
	}
	
	
	void run() {
		this.treeIDRegionMap=new HashMap<>();
		
		try {
			////////////////////
			BufferedReader lineReader = new BufferedReader(new FileReader(this.treeIDRegionFile.toFile()));
			String line = null;
			
			while ((line = lineReader.readLine()) != null) {
				if(line.isEmpty() || line.startsWith("\"\"")) { //skip first line: "","V1","V2"
					continue;
				}
				
				//10	Chr01	9000001	10000000
				String[] splits=line.split("\\s+");
				
				
				int treeID=Integer.parseInt(splits[0]);
				String chrom=splits[1];
				int start=Integer.parseInt(splits[2]);
				int end=Integer.parseInt(splits[3]);
				
				this.treeIDRegionMap.put(treeID, new SimpleGenomicRegion(chrom, start, end));
			}
			
			lineReader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	///////////////////////////////////
	/**
	 * @return the treeIDRegionMap
	 */
	public Map<Integer, SimpleGenomicRegion> getTreeIDRegionMap() {
		return treeIDRegionMap;
	}
	
	
}
