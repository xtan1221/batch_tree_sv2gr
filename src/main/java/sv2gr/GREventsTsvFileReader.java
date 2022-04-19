package sv2gr;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import population.sv.utils.SimpleSVType;
import sv2gr.gr.GRType;

/**
 * reader of the tsv file written by {@link GREventsWriter}
 * @author tanxu
 * 
 */
public class GREventsTsvFileReader {
	private final Path inferredGREventsTsvFile;
	
	////////////////////////
	private List<SimpleGREvent> grEvents;
	
	public GREventsTsvFileReader(Path inferredGREventsTsvFile) {
		super();
		this.inferredGREventsTsvFile = inferredGREventsTsvFile;
		
		///////////////////////////////
		this.run();
	}

	void run() {
		this.grEvents=new ArrayList<>();
		
		try {
			////////////////////
			BufferedReader lineReader = new BufferedReader(new FileReader(this.inferredGREventsTsvFile.toFile()));
			String line = null;
			
			while ((line = lineReader.readLine()) != null) {
				if(line.isEmpty() || line.startsWith("chrom")) { //skip header line
					continue;
				}
				
				//id	family	order
				//Sobic.004G180200.1_1_Pkinase	Sobic.004G180200.1_1_Pkinase	Rgene
				String[] splits=line.split("\t");
				
//				.append("chrom").append("\t")
				String chrom=splits[0];
//				.append("SV_type").append("\t")
				SimpleSVType svType=SimpleSVType.valueOf(splits[1]);
//				.append("SV_start").append("\t")
				int svStart=Integer.parseInt(splits[2]);
//				.append("SV_end").append("\t")
				int svEnd=Integer.parseInt(splits[3]);
//				.append("size").append("\t") //size of the sv and GR?
				int svSize=Integer.parseInt(splits[4]);
//				//gr related
//				.append("GR_type").append("\t")
				GRType grEventType=GRType.valueOf(splits[5]);
//				.append("reverse_of_SV").append("\t")
				boolean reverseOfSV=Boolean.parseBoolean(splits[6]);
//				.append("anc_dist").append("\t") //average distance from ancestral node to all leaf nodes
				double ancDist=Double.parseDouble(splits[7]);
//				.append("desc_dist").append("\t")//average distance from descendant node to all leaf nodes
				double descDist=Double.parseDouble(splits[8]);
//				.append("branch_num").append("\t")//number of branches on the path where the GR event occured
				int branchNum=Integer.parseInt(splits[9]);
//				.append("max_bootstrap").append("\t")//max bootstrap value of branch on the path where the GR event occured
				double maxBootstrap=Double.parseDouble(splits[10]);
//				//regional tree related
//				.append("window_size").append("\t")//window size of the regional tree on which the GR event was inferred
				int treeWindowSize=Integer.parseInt(splits[11]);
//				.append("tree_id").append("\t")//regional tree id
				int treeID=Integer.parseInt(splits[12]);
//				.append("tree_start").append("\t")
				int treeStart=Integer.parseInt(splits[13]);
//				.append("tree_end").append("\t")
				int treeEnd=Integer.parseInt(splits[14]);
//				.append("anc_leaf_samples").append("\t")//the list of samples corresponding to the leaf nodes of the clade rooted with the SECOND node to ancestor node of the path where the GR event occured
				String ancLeafSamplesString = splits[15];
				String[] ancSamplesSplits=ancLeafSamplesString.split(":");
				List<String> ancLeafSamples=new ArrayList<>();
				for(String sample:ancSamplesSplits) {
					if(!sample.isEmpty())
						ancLeafSamples.add(sample);
				}
//				.append("desc_leaf_samples");//the list of samples corresponding to the leaf nodes of the clade rooted with the descendant node of the path where the GR event occured
				String descLeafSamplesString = splits[16];
				String[] descSamplesSplits=descLeafSamplesString.split(":");
				List<String> descLeafSamples=new ArrayList<>();
				for(String sample:descSamplesSplits) {
					if(!sample.isEmpty())
						descLeafSamples.add(sample);
				}
				
				SimpleSV sv = new SimpleSV(chrom, svType, svStart, svEnd, svSize);
				SimpleGREvent gr=new SimpleGREvent(
						sv, grEventType, reverseOfSV, ancDist, descDist,
						branchNum, maxBootstrap, treeWindowSize, treeID, treeStart, treeEnd,
						ancLeafSamples, descLeafSamples);
				
				this.grEvents.add(gr);
			}
			
			lineReader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	
	/**
	 * @return the grEvents
	 */
	public List<SimpleGREvent> getGrEvents() {
		return grEvents;
	}


	////////////////////////////////////////////
	public static class SimpleSV {
//		.append("chrom").append("\t")
		private final String chrom;
//		.append("SV_type").append("\t")
		private final SimpleSVType svType;
//		.append("SV_start").append("\t")
		private final int svStart;
//		.append("SV_end").append("\t")
		private final int svEnd;
//		.append("size").append("\t") //size of the sv and GR?
		private final int svSize;
		
		public SimpleSV(String chrom, SimpleSVType svType, int svStart, int svEnd, int svSize) {
			super();
			this.chrom = chrom;
			this.svType = svType;
			this.svStart = svStart;
			this.svEnd = svEnd;
			this.svSize = svSize;
		}

		/**
		 * @return the chrom
		 */
		public String getChrom() {
			return chrom;
		}

		/**
		 * @return the svType
		 */
		public SimpleSVType getSvType() {
			return svType;
		}

		/**
		 * @return the svStart
		 */
		public int getSvStart() {
			return svStart;
		}

		/**
		 * @return the svEnd
		 */
		public int getSvEnd() {
			return svEnd;
		}

		/**
		 * @return the svSize
		 */
		public int getSvSize() {
			return svSize;
		}
		
		
	}
	
	public static class SimpleGREvent{
		private final SimpleSV sv;
//		//gr related
//		.append("GR_type").append("\t")
		private final GRType grEventType; 
//		.append("reverse_of_SV").append("\t")
		private final boolean reverseOfSV;
//		.append("anc_dist").append("\t") //average distance from ancestral node to all leaf nodes
		private final double ancDist;
//		.append("desc_dist").append("\t")//average distance from descendant node to all leaf nodes
		private final double descDist;
//		.append("branch_num").append("\t")//number of branches on the path where the GR event occured
		private final int branchNum;
//		.append("max_bootstrap").append("\t")//max bootstrap value of branch on the path where the GR event occured
		private final double maxBootstrap;
//		//regional tree related
//		.append("window_size").append("\t")//window size of the regional tree on which the GR event was inferred
		private final int treeWindowSize;
//		.append("tree_id").append("\t")//regional tree id
		private final int treeID;
//		.append("tree_start").append("\t")
		private final int treeStart;
//		.append("tree_end").append("\t")
		private final int treeEnd;
//		.append("anc_leaf_samples").append("\t")//the list of samples corresponding to the leaf nodes of the clade rooted with the SECOND node to ancestor node of the path where the GR event occured
		private final List<String> ancLeafSamples;
//		.append("desc_leaf_samples");//the list of samples corresponding to the leaf nodes of the clade rooted with the descendant node of the path where the GR event occured
		private final List<String> descLeafSamples;
		
		public SimpleGREvent(SimpleSV sv, GRType grEventType, boolean reverseOfSV, double ancDist, double descDist,
				int branchNum, double maxBootstrap, int treeWindowSize, int treeID, int treeStart, int treeEnd,
				List<String> ancLeafSamples, List<String> descLeafSamples) {
			super();
			this.sv = sv;
			this.grEventType = grEventType;
			this.reverseOfSV = reverseOfSV;
			this.ancDist = ancDist;
			this.descDist = descDist;
			this.branchNum = branchNum;
			this.maxBootstrap = maxBootstrap;
			this.treeWindowSize = treeWindowSize;
			this.treeID = treeID;
			this.treeStart = treeStart;
			this.treeEnd = treeEnd;
			this.ancLeafSamples = ancLeafSamples;
			this.descLeafSamples = descLeafSamples;
		}

		/**
		 * @return the sv
		 */
		public SimpleSV getSv() {
			return sv;
		}

		/**
		 * @return the grEventType
		 */
		public GRType getGrEventType() {
			return grEventType;
		}

		/**
		 * @return the reverseOfSV
		 */
		public boolean isReverseOfSV() {
			return reverseOfSV;
		}

		/**
		 * @return the ancDist
		 */
		public double getAncDist() {
			return ancDist;
		}

		/**
		 * @return the descDist
		 */
		public double getDescDist() {
			return descDist;
		}

		/**
		 * @return the branchNum
		 */
		public int getBranchNum() {
			return branchNum;
		}

		/**
		 * @return the maxBootstrap
		 */
		public double getMaxBootstrap() {
			return maxBootstrap;
		}

		/**
		 * @return the treeWindowSize
		 */
		public int getTreeWindowSize() {
			return treeWindowSize;
		}

		/**
		 * @return the treeID
		 */
		public int getTreeID() {
			return treeID;
		}

		/**
		 * @return the treeStart
		 */
		public int getTreeStart() {
			return treeStart;
		}

		/**
		 * @return the treeEnd
		 */
		public int getTreeEnd() {
			return treeEnd;
		}

		/**
		 * @return the ancLeafSamples
		 */
		public List<String> getAncLeafSamples() {
			return ancLeafSamples;
		}

		/**
		 * @return the descLeafSamples
		 */
		public List<String> getDescLeafSamples() {
			return descLeafSamples;
		}
	}
}
