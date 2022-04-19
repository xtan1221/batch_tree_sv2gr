package genomics.utils;

import java.util.ArrayList;
import java.util.List;

/**
 * a gene on a chromosome that may contains one or more {@link Transcript}s
 * 
 * @author tanxu
 *
 */
public class Gene extends SimpleGenomicRegion{
	private final String name;
	private final List<Transcript> transcripts;
	
	private Transcript longestTranscript;
	
	
	public Gene(String name, String chrom, int start, int end, Strand strand) {
		super(chrom, start, end, strand);
		if(name==null||name.isEmpty()) {
			throw new IllegalArgumentException("given gene name cannot be null or empty!");
		}
		
		this.name=name;
		this.transcripts=new ArrayList<>();
	}
	
	
	public void addTranscript(Transcript transcript) {
		this.transcripts.add(transcript);
	}


	/**
	 * @return the transcripts
	 */
	public List<Transcript> getTranscripts() {
		return transcripts;
	}


	/**
	 * @return the name
	 */
	public String getName() {
		return name;
	}

	
	/**
	 * find out and return the longest transcript
	 * @return
	 */
	public Transcript getLongestTranscript() {
		if(this.longestTranscript==null) {
			for(Transcript t:this.transcripts) {
				if(this.longestTranscript==null||this.longestTranscript.getLen()<t.getLen())
					this.longestTranscript=t;
			}
		}
		
		return this.longestTranscript;
	}
	
}
