package phylogeny;

import java.util.List;

/**
 * calculate the mean distance between all pairs of individuals between two given group of individuals
 * 
 * to calculate within group mean distance, simply set both the value of {@link #group1IDs} and {@link #group2IDs} to that group
 * @author tanxu
 *
 */
public class BetweenGroupMeanDistanceCalculator {
	private final PhylipDistanceMatrixFileReader phylipDistanceMatrixFileReader;
	private final List<String> group1IDs;
	private final List<String> group2IDs;
	
	/////////////////////////////////
	private double meanDist;
	
	public BetweenGroupMeanDistanceCalculator(PhylipDistanceMatrixFileReader phylipDistanceMatrixFileReader,
			List<String> group1iDs, List<String> group2iDs) {
		super();
		this.phylipDistanceMatrixFileReader = phylipDistanceMatrixFileReader;
		group1IDs = group1iDs;
		group2IDs = group2iDs;
		
		/////////////////////////
		this.calculate();
	}


	void calculate() {
		double sum=0;
		int pair=0;
		
		for(int i=0;i<this.group1IDs.size();i++) {
			String id1=this.group1IDs.get(i);
			for(int j=0;j<this.group2IDs.size();j++) {
				String id2=this.group2IDs.get(j);
				
				sum+=this.phylipDistanceMatrixFileReader.lookupDistance(id1, id2);
				pair++;
			}
		}
		
		this.meanDist=sum/pair;
	}

	
	////////////////////////////////
	/**
	 * @return the meanDist
	 */
	public double getMeanDist() {
		return meanDist;
	}
	
}
