package population.phylogeny.adjacentTreeDist;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * reader for a file containing the distance between adjacent tree pairs
 * 
 * 
 * @author tanxu
 *
 */
public class AdjacentRegionalTreeDistFileReader {
	private final Path adjacentRegionalTreeDistFile;
	
	////////////////////////////////
	private Map<String, List<AdjacentTreePairDist>> chromAdjacentTreePairDistsMap;
	
	
	public AdjacentRegionalTreeDistFileReader(Path adjacentRegionalTreeDistFile) {
		super();
		this.adjacentRegionalTreeDistFile = adjacentRegionalTreeDistFile;
		
		
		///////////////////////////
		this.run();
	}
	
	void run() {
		this.chromAdjacentTreePairDistsMap=new HashMap<>();
		
		try {
			////////////////////
			BufferedReader lineReader = new BufferedReader(new FileReader(this.adjacentRegionalTreeDistFile.toFile()));
			String line = null;
			
			while ((line = lineReader.readLine()) != null) {
				if(line.isEmpty() || line.startsWith("\"\"")) { //skip first line: "","dist","Chr","tree1ID","tree2ID","tree1Start","tree1End","tree2Start","tree2End"
					continue;
				}
				
				//"411",-0.376731077051497,0.149388208125863
				//"all_chrom_missing_data_0",-0.00715884945308318,0.00533199453320135
				String[] splits=line.split(",");
				
				double distance=Double.parseDouble(splits[1]);
				String chrom=splits[2].substring(1, splits[2].length()-1);
				int tree1ID=Integer.parseInt(splits[3]);
				int tree2ID=Integer.parseInt(splits[4]);
				int tree1Start=Integer.parseInt(splits[5]);
				int tree1End=Integer.parseInt(splits[6]);
				int tree2Start=Integer.parseInt(splits[7]);
				int tree2End=Integer.parseInt(splits[8]);
				
				AdjacentTreePairDist adjacentTreePairDist = 
						new AdjacentTreePairDist(
								chrom, distance, 
								tree1ID, tree2ID, 
								tree1Start, tree1End, 
								tree2Start, tree2End);
				
				if(!this.chromAdjacentTreePairDistsMap.containsKey(chrom)) {
					this.chromAdjacentTreePairDistsMap.put(chrom, new ArrayList<>());
				}
				this.chromAdjacentTreePairDistsMap.get(chrom).add(adjacentTreePairDist);
			}
			
			lineReader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public List<AdjacentTreePairDist> getAdjacentTreePairDists(){
		List<AdjacentTreePairDist> ret = new ArrayList<>();
		this.getChromAdjacentTreePairDistsMap().forEach((chrom, pairs)->{
			ret.addAll(pairs);
		});
		return ret;
	}
	
	/**
	 * @return the chromAdjacentTreePairDistsMap
	 */
	public Map<String, List<AdjacentTreePairDist>> getChromAdjacentTreePairDistsMap() {
		return chromAdjacentTreePairDistsMap;
	}
	

	static class AdjacentTreePairDist{
		private final String chrom;
		private final double distance;
		private final int regionalTree1ID;
		private final int regionalTree2ID;
		private final int regionalTree1Start;
		private final int regionalTree1End;
		private final int regionalTree2Start;
		private final int regionalTree2End;
		
		public AdjacentTreePairDist(String chrom, double distance, int regionalTree1ID, int regionalTree2ID, int regionalTree1Start,
				int regionalTree1End, int regionalTree2Start, int regionalTree2End) {
			super();
			this.chrom = chrom;
			this.distance=distance;
			this.regionalTree1ID = regionalTree1ID;
			this.regionalTree2ID = regionalTree2ID;
			this.regionalTree1Start = regionalTree1Start;
			this.regionalTree1End = regionalTree1End;
			this.regionalTree2Start = regionalTree2Start;
			this.regionalTree2End = regionalTree2End;
		}

		/**
		 * @return the chrom
		 */
		public String getChrom() {
			return chrom;
		}

		/**
		 * @return the distance
		 */
		public double getDistance() {
			return distance;
		}

		/**
		 * @return the regionalTree1ID
		 */
		public int getRegionalTree1ID() {
			return regionalTree1ID;
		}

		/**
		 * @return the regionalTree2ID
		 */
		public int getRegionalTree2ID() {
			return regionalTree2ID;
		}

		/**
		 * @return the regionalTree1Start
		 */
		public int getRegionalTree1Start() {
			return regionalTree1Start;
		}

		/**
		 * @return the regionalTree1End
		 */
		public int getRegionalTree1End() {
			return regionalTree1End;
		}

		/**
		 * @return the regionalTree2Start
		 */
		public int getRegionalTree2Start() {
			return regionalTree2Start;
		}

		/**
		 * @return the regionalTree2End
		 */
		public int getRegionalTree2End() {
			return regionalTree2End;
		}
		
		
	}
}
