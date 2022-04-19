package genomics.geneAnnotation;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;

import genomics.utils.Strand;
import genomics.utils.Transcript;

/**
 * reader for a bed file containing a set of transcripts
 * 
 * the 'name' column is the transcript id and should be unique
 * 
 * @author tanxu
 *
 */
public class TranscriptGenomicRegionBedReader {
	private final Path transcriptGenomicRegionBedFile;
	
	//////////////////
	private Map<String, Transcript> transcriptIDGenomicRegionMap;
	
	public TranscriptGenomicRegionBedReader(Path transcriptGenomicRegionBedFile) {
		super();
		this.transcriptGenomicRegionBedFile = transcriptGenomicRegionBedFile;
		
		this.run();
	}


	private void run() {
		this.transcriptIDGenomicRegionMap=new HashMap<>();
		
		try {
			////////////////////
			BufferedReader lineReader = new BufferedReader(new FileReader(transcriptGenomicRegionBedFile.toFile()));
			String line = null;
			
			while ((line = lineReader.readLine()) != null) {
				if(line.isEmpty()) {
					continue;
				}
				
				//Chr01	22390	37967	Sobic.001G000400.6	.	-
				String[] splits=line.split("\t");
				String chrom=splits[0];
				int start = Integer.parseInt(splits[1])+1; //bed file format's start column is 0-based
				int end=Integer.parseInt(splits[2]);
				//
				String transcriptID = splits[3];
				if(this.transcriptIDGenomicRegionMap.containsKey(transcriptID)) {
					lineReader.close();
					throw new IllegalArgumentException("duplicate transcript id "+transcriptID+" is found in the bed file, which is not expected!");
				}
				Strand strand=Strand.find(splits[5]);
				
				this.transcriptIDGenomicRegionMap.put(
						transcriptID, 
						new Transcript(chrom, start, end, strand, transcriptID)
						);
			}
			
			lineReader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}


	/**
	 * @return the transcriptIDGenomicRegionMap
	 */
	public Map<String, Transcript> getTranscriptIDGenomicRegionMap() {
		return transcriptIDGenomicRegionMap;
	}
	
	
}
