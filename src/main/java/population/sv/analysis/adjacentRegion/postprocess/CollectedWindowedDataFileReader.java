package population.sv.analysis.adjacentRegion.postprocess;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import population.sv.analysis.adjacentRegion.piAndThetaW.CollectedPiAndThetaWDataWriter;
import population.sv.utils.SimpleSVType;

/**
 * see the {@link CollectedPiAndThetaWDataWriter} for details about the data file
 * 
 * @author tanxu
 *
 */
public class CollectedWindowedDataFileReader {
	/**
	 * 
	 */
	private final Path collectedDataFile;
	/**
	 * number of window on each side of the SV
	 */
	private final int windowNum;
	
	/////////////////////
	private Map<SimpleSVType, List<SVRecord>> svTypeSVRecordsMap;
	
	public CollectedWindowedDataFileReader(Path collectedDataFile, int windowNum) {
		super();
		this.collectedDataFile = collectedDataFile;
		this.windowNum = windowNum;
		
		this.readFile();
	}
	
	void readFile() {
		this.svTypeSVRecordsMap=new HashMap<>();
		
		try {
			////////////////////
			BufferedReader lineReader = new BufferedReader(new FileReader(this.collectedDataFile.toFile()));
			String line = null;
			
			while ((line = lineReader.readLine()) != null) {
				if(line.startsWith("sv_type")||line.isEmpty()) { //header line starts with #
					continue;
				}
				
				//#sample_index	population_level	sample_name
				//1		1	propinquum_4
				String[] splits=line.split("\\s+");
				
				SimpleSVType svType=SimpleSVType.valueOf(splits[0]);
				String chrom=splits[1];
				int start=Integer.parseInt(splits[2]);
				int end=Integer.parseInt(splits[3]);
				int size=Integer.parseInt(splits[4]);
				Double derivedAlleleProportion=splits[5].equals("N/A")?null:Double.parseDouble(splits[5]);
				
				List<Double> leftSideWindowProportion = new ArrayList<>();
				for(int i=0;i<this.windowNum;i++) {
					leftSideWindowProportion.add(splits[i+6].equals("N/A")?null:Double.parseDouble(splits[i+6]));
				}
				
				List<Double> rightSideWindowProportion = new ArrayList<>();
				for(int i=0;i<this.windowNum;i++) {
					rightSideWindowProportion.add(splits[i+this.windowNum+6].equals("N/A")?null:Double.parseDouble(splits[i+this.windowNum+6]));
				}
				
				SVRecord record = 
						new SVRecord(
							svType, chrom, start ,end, size, derivedAlleleProportion, 
							leftSideWindowProportion, rightSideWindowProportion);
				
				if(!this.svTypeSVRecordsMap.containsKey(svType)) {
					this.svTypeSVRecordsMap.put(svType, new ArrayList<>());
				}
				this.svTypeSVRecordsMap.get(svType).add(record);
			}
			
			lineReader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	
	/**
	 * @return the svTypeSVRecordsMap
	 */
	public Map<SimpleSVType, List<SVRecord>> getSvTypeSVRecordsMap() {
		return svTypeSVRecordsMap;
	}



	/////////////////////
	public static class SVRecord{
		private final SimpleSVType type;
		private final String chrom;
		private final int start;
		private final int end;
		private final int size;
		private final Double derivedAlleleProportion;
		/**
		 * start from the right-most window
		 */
		private final List<Double> leftSideWindowValues;
		/**
		 * start from the left-most window
		 */
		private final List<Double> rightSideWindowValues;
		
		
		public SVRecord(
				SimpleSVType type, 
				String chrom, int start, int end, int size,
				Double derivedAlleleProportion,
				List<Double> leftSideWindowProportion,
				List<Double> rightSideWindowProportion) {
			super();
			this.type = type;
			this.chrom=chrom;
			this.start = start;
			this.end = end;
			this.size=size;
			this.derivedAlleleProportion=derivedAlleleProportion;
			this.leftSideWindowValues = leftSideWindowProportion;
			this.rightSideWindowValues = rightSideWindowProportion;
		}
		
		public int getSize() {
			return this.size;
		}
		/**
		 * @return the type
		 */
		public SimpleSVType getType() {
			return type;
		}


		/**
		 * @return the start
		 */
		public int getStart() {
			return start;
		}


		/**
		 * @return the end
		 */
		public int getEnd() {
			return end;
		}

		
		/**
		 * @return the chrom
		 */
		public final String getChrom() {
			return chrom;
		}

		/**
		 * @return the derivedAlleleProportion
		 */
		public final Double getDerivedAlleleProportion() {
			return derivedAlleleProportion;
		}

		/**
		 * @return the leftSideWindowProportion
		 */
		public List<Double> getLeftSideWindowProportion() {
			return leftSideWindowValues;
		}


		/**
		 * @return the rightSideWindowProportion
		 */
		public List<Double> getRightSideWindowProportion() {
			return rightSideWindowValues;
		}
		
		
		
	}
	
	
}
