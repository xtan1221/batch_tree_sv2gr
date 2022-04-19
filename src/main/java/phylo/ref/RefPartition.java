package phylo.ref;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;


/**
 * output a bed file containing partition of a genome into windows of a specific size
 * 
 * @author tanxu
 *
 */
public class RefPartition {
	/**
	 * each line contains two columns delimited by space
	 * first column is the name of the chromosome/reference
	 * second column is the length;
	 */
	private final Path chromLengthFile;
	
	/**
	 * the minimal length of a chromosome to be included in the resulted partition
	 */
	private final int minChromLen;
	
	/**
	 * window size
	 */
	private final int windowSize;
	/**
	 * whether or not keep the last window with size less than the {@link #windowSize}
	 */
	private final boolean toKeepTrailingShortWindow;
	
	///////////////////////////
	private List<Region> windowList;
	
	/**
	 * 
	 * @param chromLengthFile
	 * @param minChromLen
	 * @param windowSize
	 */
	RefPartition(Path chromLengthFile, int minChromLen, int windowSize, boolean toKeepTrailingShortWindow){
		
		
		this.chromLengthFile = chromLengthFile;
		this.minChromLen = minChromLen;
		this.windowSize = windowSize;
		this.toKeepTrailingShortWindow = toKeepTrailingShortWindow;
		
		this.run();
	}
	

	void run() {
		this.windowList = new ArrayList<>();
		
		try {
		    BufferedReader lineReader = new BufferedReader(new FileReader(this.chromLengthFile.toFile()));
		    String line = null;
		 
		    while ((line = lineReader.readLine()) != null) {
//		        System.out.println(line);
		        if(line.trim().isEmpty()) {
					continue;
				}	
				String[] splits = line.split("\\s+");
				
				String chrom = splits[0];
				int len = Integer.parseInt(splits[1]);
				
				if(len<this.minChromLen) {
					continue;
				}
				
				int index=0;
				
				while((index+1)*this.windowSize<len) {
					this.windowList.add(new Region(chrom, index*this.windowSize+1, (index+1)*this.windowSize));
					index++;
				}
				
				if(index*this.windowSize==len) {
					//the end of last full window is exactly the same with the length of the chromosome
				}else {
					if(this.toKeepTrailingShortWindow)
						this.windowList.add(new Region(chrom, index*this.windowSize+1, len));
				}
		    }
		    
		    lineReader.close();
		} catch (IOException ex) {
		    System.err.println(ex);
		}
	}
	
	
	/**
	 * @return the windowList
	 */
	public List<Region> getWindowList() {
		return windowList;
	}
	
	
}
