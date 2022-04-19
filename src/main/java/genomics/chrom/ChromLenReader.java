package genomics.chrom;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ChromLenReader {
	
	/**
	 * column 1 is the chrom name, column2 is the length
	 */
	private final Path chromLengthFile;
	
	/**
	 * comparator to sort chrom names, if null , use the default one
	 */
	private final Comparator<String> chromNameComparator;
	
	///////////////////////////////////
	private List<String> sortedChromNames;
	private Map<String, Integer> chromNameLengthMap;
	private int totalChromLen;
	
	/**
	 * 
	 * @param chromLenFile
	 */
	public ChromLenReader(Path chromLenFile) {
		super();
		this.chromLengthFile = chromLenFile;
		this.chromNameComparator=null;
		
		this.readFile();
	}
	
	public ChromLenReader(Path chromLengthFile, Comparator<String> chromNameComparator) {
		super();
		this.chromLengthFile = chromLengthFile;
		this.chromNameComparator=chromNameComparator;
		this.readFile();
	}
	


	void readFile() {
		this.sortedChromNames=new ArrayList<>();
		this.chromNameLengthMap=new HashMap<>();
		
		
		try {
			////////////////////
			BufferedReader lineReader = new BufferedReader(new FileReader(this.chromLengthFile.toFile()));
			String line = null;
			
			while ((line = lineReader.readLine()) != null) {
				if(line.isEmpty()) { //first line is id	family	order
					continue;
				}
				
				String[] splits = line.split("\\s+");
				//Sobic.005G019600.5.p ===> Sobic.005G019600.5
				
				String chromName = splits[0];
				int len = Integer.parseInt(splits[1]);
				
				this.sortedChromNames.add(chromName);
				this.chromNameLengthMap.put(chromName, len);
				
				this.totalChromLen+=len;
			}

			lineReader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		if(this.chromNameComparator!=null) {
			Collections.sort(this.sortedChromNames, this.chromNameComparator);
		}else {
			Collections.sort(this.sortedChromNames);
		}
	}



	/**
	 * @return the sortedChromNames
	 */
	public List<String> getSortedChromNames() {
		return sortedChromNames;
	}

	

	/**
	 * @return the chromNameLengthMap
	 */
	public Map<String, Integer> getChromNameLengthMap() {
		return chromNameLengthMap;
	}
	
	public int getTotalChromLen() {
		return this.totalChromLen;
	}
}
