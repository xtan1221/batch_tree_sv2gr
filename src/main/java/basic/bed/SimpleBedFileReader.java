package basic.bed;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import genomics.utils.SimpleGenomicRegion;


/**
 * reader for the first three columns of a bed file into a set of {@link SimpleGenomicRegion};
 * 		chrom, start, end
 * 
 * note that bed file's 'start' column is exclusive while 'end' column is inclusive
 * 
 * @author tanxu
 *
 */
public class SimpleBedFileReader {
	private final Path bedFile;
	
	//////////////////////////////////
	/**
	 * 
	 */
	private List<SimpleGenomicRegion> regions;
	
	
	public SimpleBedFileReader(Path refMaskedRegionBed) {
		super();
		this.bedFile = refMaskedRegionBed;
		this.read();
	}

	
	private void read() {
		this.regions=new ArrayList<>();
		
		try {
			////////////////////
			BufferedReader lineReader = new BufferedReader(new FileReader(this.bedFile.toFile()));
			String line = null;
			
			while ((line = lineReader.readLine()) != null) {
				if(line.isEmpty()) { //first line is id	family	order
					continue;
				}
				
				//id	family	order
				//Sobic.004G180200.1_1_Pkinase	Sobic.004G180200.1_1_Pkinase	Rgene
				String[] splits=line.split("\\s+");
				
				String chrom=splits[0];
				int start = Integer.parseInt(splits[1])+1; //start column is exclusive in bed file, thus need to fix it
				int end=Integer.parseInt(splits[2]); //end column is inclusive
				
				SimpleGenomicRegion gr = new SimpleGenomicRegion(chrom, start, end);
				this.regions.add(gr);
			}
			
			lineReader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	

	/**
	 * @return the maskedRegions
	 */
	public List<SimpleGenomicRegion> getRegions() {
		return regions;
	}
	
	
}
