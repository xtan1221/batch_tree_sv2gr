package population.phylogeny.adjacentTreeDist;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;

public class AdjacentRegionalTreeDistStatisticalSummaryFileReader {
	private final Path adjacentRegionalTreeDistStatisticalSummaryFile;
	
	///////////////////////////
	private ChromSummary fullGenomeChromSummary;
	private Map<String, ChromSummary> chromNameChromSummaryMap;
	
	
	public AdjacentRegionalTreeDistStatisticalSummaryFileReader(Path adjacentRegionalTreeDistStatisticalSummaryFile) {
		super();
		this.adjacentRegionalTreeDistStatisticalSummaryFile = adjacentRegionalTreeDistStatisticalSummaryFile;
		
		
		//////////////////////
		this.run();
	}



	void run() {
		chromNameChromSummaryMap = new HashMap<>();
		
		try {
			////////////////////
			BufferedReader lineReader = new BufferedReader(new FileReader(this.adjacentRegionalTreeDistStatisticalSummaryFile.toFile()));
			String line = null;
			
			while ((line = lineReader.readLine()) != null) {
				if(line.isEmpty() || line.startsWith("\"\"")) { //skip first line: "","trees","total_num","mean_dist","sd_dist","variance_dist"
					continue;
				}
				
				//"411",-0.376731077051497,0.149388208125863
				//"all_chrom_missing_data_0",-0.00715884945308318,0.00533199453320135
				String[] splits=line.split(",");
				
				String regionString=splits[1].substring(1, splits[1].length()-1);
				int totalNum=Integer.parseInt(splits[2]);
				double meanDist=Double.parseDouble(splits[3]);
				double sdDist=Double.parseDouble(splits[4]);
				double varianceDist=Double.parseDouble(splits[5]);
				
				if(regionString.equals("all_chrom")) {//pair of trees across full genome
					if(this.fullGenomeChromSummary==null) {
						this.fullGenomeChromSummary=new ChromSummary();
						this.fullGenomeChromSummary.setChrom("full_genome");
					}
					this.fullGenomeChromSummary.setTotalPairNum(totalNum);
					this.fullGenomeChromSummary.setMeanDist(meanDist);
					this.fullGenomeChromSummary.setSd(sdDist);
					this.fullGenomeChromSummary.setVariance(varianceDist);
				}else if(regionString.equals("adjacent_all_chrom")) {//pairs of adjacent regional trees across full genome
					if(this.fullGenomeChromSummary==null) {
						this.fullGenomeChromSummary=new ChromSummary();
						this.fullGenomeChromSummary.setChrom("full_genome");
					}
					this.fullGenomeChromSummary.setTotalAdjacentPairNum(totalNum);
					this.fullGenomeChromSummary.setMeanDistOfAdjacentPairs(meanDist);
					this.fullGenomeChromSummary.setSdOfAdjacentPairs(sdDist);
					this.fullGenomeChromSummary.setVarianceOfAdjacentPairs(varianceDist);
				}else {//chromosome: all_Chr1 or adjacent_Chr2
					String[] splits2=regionString.split("_");
					String chromName=splits2[1];
					if(!this.chromNameChromSummaryMap.containsKey(chromName)) {
						ChromSummary cs=new ChromSummary();
						cs.setChrom(chromName);
						this.chromNameChromSummaryMap.put(chromName, cs);
					}
					
					if(regionString.startsWith("adjacent")) {//pairs of adjacent regional trees on a chrom
						this.chromNameChromSummaryMap.get(chromName).setTotalAdjacentPairNum(totalNum);
						this.chromNameChromSummaryMap.get(chromName).setMeanDistOfAdjacentPairs(meanDist);
						this.chromNameChromSummaryMap.get(chromName).setSdOfAdjacentPairs(sdDist);
						this.chromNameChromSummaryMap.get(chromName).setVarianceOfAdjacentPairs(varianceDist);
					}else {////pairs of all regional trees on a chrom
						this.chromNameChromSummaryMap.get(chromName).setTotalPairNum(totalNum);
						this.chromNameChromSummaryMap.get(chromName).setMeanDist(meanDist);
						this.chromNameChromSummaryMap.get(chromName).setSd(sdDist);
						this.chromNameChromSummaryMap.get(chromName).setVariance(varianceDist);
					}
				}
				
			}
			
			lineReader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	
	/**
	 * @return the fullGenomeChromSummary
	 */
	public ChromSummary getFullGenomeChromSummary() {
		return fullGenomeChromSummary;
	}



	/**
	 * @return the chromNameChromSummaryMap
	 */
	public Map<String, ChromSummary> getChromNameChromSummaryMap() {
		return chromNameChromSummaryMap;
	}


	
	static class ChromSummary{
		private String chrom;
		private int totalPairNum;
		private int totalAdjacentPairNum;
		/**
		 * mean distance between all pairs of trees on the chrom?
		 */
		private double meanDist;
		private double meanDistOfAdjacentPairs;
		/**
		 * sd of distance between all pairs of trees on the chrom?
		 */
		private double sd;
		private double sdOfAdjacentPairs;
		/**
		 * variance of distance between all pairs of trees on the chrom?
		 */
		private double variance;
		private double varianceOfAdjacentPairs;
		
		public ChromSummary() {
			super();
		}

		/**
		 * @return the chrom
		 */
		public String getChrom() {
			return chrom;
		}

		/**
		 * @param chrom the chrom to set
		 */
		public void setChrom(String chrom) {
			this.chrom = chrom;
		}

		/**
		 * @return the totalPairNum
		 */
		public int getTotalPairNum() {
			return totalPairNum;
		}

		/**
		 * @param totalPairNum the totalPairNum to set
		 */
		public void setTotalPairNum(int totalPairNum) {
			this.totalPairNum = totalPairNum;
		}

		/**
		 * @return the totalAdjacentPairNum
		 */
		public int getTotalAdjacentPairNum() {
			return totalAdjacentPairNum;
		}

		/**
		 * @param totalAdjacentPairNum the totalAdjacentPairNum to set
		 */
		public void setTotalAdjacentPairNum(int totalAdjacentPairNum) {
			this.totalAdjacentPairNum = totalAdjacentPairNum;
		}

		/**
		 * @return the meanDist
		 */
		public double getMeanDist() {
			return meanDist;
		}

		/**
		 * @param meanDist the meanDist to set
		 */
		public void setMeanDist(double meanDist) {
			this.meanDist = meanDist;
		}

		/**
		 * @return the meanDistOfAdjacentPairs
		 */
		public double getMeanDistOfAdjacentPairs() {
			return meanDistOfAdjacentPairs;
		}

		/**
		 * @param meanDistOfAdjacentPairs the meanDistOfAdjacentPairs to set
		 */
		public void setMeanDistOfAdjacentPairs(double meanDistOfAdjacentPairs) {
			this.meanDistOfAdjacentPairs = meanDistOfAdjacentPairs;
		}

		/**
		 * @return the sd
		 */
		public double getSd() {
			return sd;
		}

		/**
		 * @param sd the sd to set
		 */
		public void setSd(double sd) {
			this.sd = sd;
		}

		/**
		 * @return the sdOfAdjacentPairs
		 */
		public double getSdOfAdjacentPairs() {
			return sdOfAdjacentPairs;
		}

		/**
		 * @param sdOfAdjacentPairs the sdOfAdjacentPairs to set
		 */
		public void setSdOfAdjacentPairs(double sdOfAdjacentPairs) {
			this.sdOfAdjacentPairs = sdOfAdjacentPairs;
		}

		/**
		 * @return the variance
		 */
		public double getVariance() {
			return variance;
		}

		/**
		 * @param variance the variance to set
		 */
		public void setVariance(double variance) {
			this.variance = variance;
		}

		/**
		 * @return the varianceOfAdjacentPairs
		 */
		public double getVarianceOfAdjacentPairs() {
			return varianceOfAdjacentPairs;
		}

		/**
		 * @param varianceOfAdjacentPairs the varianceOfAdjacentPairs to set
		 */
		public void setVarianceOfAdjacentPairs(double varianceOfAdjacentPairs) {
			this.varianceOfAdjacentPairs = varianceOfAdjacentPairs;
		}
		
		
	}
}
