package population.phylogeny;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;


/**
 * reader for CSV file containing the calculated coordinates of all regional trees, chromosomal trees and the full genome tree
 * 
 * 
 * regional tree
 * 		"33",0.0487103762400359,-0.0171330380442861 //regional trees' id are only numbers
 * chrom tree
 * 		"Chr01_missing_data_27",0.105232444209899,-0.0309108382945484 //missing data 27 indicates at most 27 samples can have missing data at any specific nucleotide site
 * full genome tree
 * 		"all_chrom_missing_data_10",0.0142333098233633,0.0403515132890494 //missing data 10 indicates at most 10 samples can have missing data at any specific nucleotide site
 * 
 * @author tanxu
 *
 */
public class AllTreeMDSCoordinateFileReader {
	private final Path treeMDSCoordinateFile;
	/**
	 * check if a tree id string is of a chrom tree
	 */
	private final Predicate<String> chromTreeIDPredicate;
	/**
	 * check if a tree id string is of a full genome tree
	 */
	private final Predicate<String> fullGenomeTreeIDPredicate; 
	
	//////////////////////////////////
	private Map<Integer, Coordinate> regionalTreeIDCoordinateMap;
	private Map<String, Map<Integer, Coordinate>> chromNameMissingDataNumCoordinateMapMap;
	private Map<Integer, Coordinate> fullGenomeMissingDataNumCoordinateMap;
	private Set<Integer> missingDataNums;
	
	
	public AllTreeMDSCoordinateFileReader(
			Path treeMDSCoordinateFile, 
			Predicate<String> chromTreeIDPredicate,
			Predicate<String> fullGenomeTreeIDPredicate) {
		super();
		this.treeMDSCoordinateFile = treeMDSCoordinateFile;
		this.chromTreeIDPredicate = chromTreeIDPredicate;
		this.fullGenomeTreeIDPredicate = fullGenomeTreeIDPredicate;
		
		//////////////////////////
		this.run();
	}
	
	void run() {
		this.regionalTreeIDCoordinateMap=new HashMap<>();
		this.chromNameMissingDataNumCoordinateMapMap = new HashMap<>();
		this.fullGenomeMissingDataNumCoordinateMap=new HashMap<>();
		this.missingDataNums=new HashSet<>();
		
		
		
		try {
			////////////////////
			BufferedReader lineReader = new BufferedReader(new FileReader(this.treeMDSCoordinateFile.toFile()));
			String line = null;
			
			while ((line = lineReader.readLine()) != null) {
				if(line.isEmpty() || line.startsWith("\"\"")) { //skip first line: "","V1","V2"
					continue;
				}
				
				//"411",-0.376731077051497,0.149388208125863
				//"all_chrom_missing_data_0",-0.00715884945308318,0.00533199453320135
				String[] splits=line.split(",");
				
				
				String treeIDString=splits[0].substring(1, splits[0].length()-1);//retrieve the tree id string
				double x=Double.parseDouble(splits[1]);
				double y=Double.parseDouble(splits[2]);
				Coordinate coord=new Coordinate(x,y);
				
				try {
					int regionalTreeID=Integer.parseInt(treeIDString);
					this.regionalTreeIDCoordinateMap.put(regionalTreeID, coord);
				}catch(NumberFormatException e) {//not a regional tree id
					if(this.chromTreeIDPredicate.test(treeIDString)) {//Chr01_missing_data_27
						String[] splits2=treeIDString.split("_missing_data_");
						String chromName=splits2[0];
						int missingDataNum=Integer.parseInt(splits2[1]);
						
						//
						if(!this.chromNameMissingDataNumCoordinateMapMap.containsKey(chromName)) {
							this.chromNameMissingDataNumCoordinateMapMap.put(chromName, new HashMap<>());
						}
						this.chromNameMissingDataNumCoordinateMapMap.get(chromName).put(missingDataNum, coord);
						
						//
						this.missingDataNums.add(missingDataNum);
					}else if(this.fullGenomeTreeIDPredicate.test(treeIDString)){
						String[] splits2=treeIDString.split("_missing_data_");
						int missingDataNum=Integer.parseInt(splits2[1]);
						
						//
						this.fullGenomeMissingDataNumCoordinateMap.put(missingDataNum, coord);
						
						//
						this.missingDataNums.add(missingDataNum);
					}else {
						lineReader.close();
						throw new UnsupportedOperationException("unrecognized tree id string:"+treeIDString);
					}
				}
			}
			
			lineReader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	/////////////////////////////////////////
	/**
	 * @return the regionalTreeIDCoordinateMap
	 */
	public Map<Integer, Coordinate> getRegionalTreeIDCoordinateMap() {
		return regionalTreeIDCoordinateMap;
	}
	
	/**
	 * @return the chromNameMissingDataNumCoordinateMapMap
	 */
	public Map<String, Map<Integer, Coordinate>> getChromNameMissingDataNumCoordinateMapMap() {
		return chromNameMissingDataNumCoordinateMapMap;
	}

	/**
	 * @return the fullGenomeMissingDataNumCoordinateMap
	 */
	public Map<Integer, Coordinate> getFullGenomeMissingDataNumCoordinateMap() {
		return fullGenomeMissingDataNumCoordinateMap;
	}

	/**
	 * @return the missingDataNums
	 */
	public Set<Integer> getMissingDataNums() {
		return missingDataNums;
	}


	/////////////////////////////////////////
	public static class Coordinate{
		private final double x;
		private final double y;
		
		public Coordinate(double x, double y) {
			super();
			this.x = x;
			this.y = y;
		}

		/**
		 * @return the x
		 */
		public double getX() {
			return x;
		}

		/**
		 * @return the y
		 */
		public double getY() {
			return y;
		}
		
		
		/**
		 * calculate the Euclidean distance between the two given coordinate
		 * @param a
		 * @param b
		 * @return
		 */
		public static double calculateEuclideanDistance(Coordinate a, Coordinate b) {
			return Math.sqrt(Math.pow(a.x-b.x, 2) + Math.pow(a.y-b.y, 2));
		}
	}
}
