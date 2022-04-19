package sv2gr;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import population.sv.utils.SimpleSVType;
import sv2gr.GREventsTsvFileReader.SimpleGREvent;
import sv2gr.GREventsTsvFileReader.SimpleSV;
import sv2gr.gr.GRType;

/**
 * reader of the tsv file written by {@link GREventsWriterForRicePop}
 * @author tanxu
 * 
 */
public class GREventsTsvFileReaderForRicePop {
	private final Path inferredGREventsTsvFile;
	
	////////////////////////
	private List<SimpleGREvent> allGREvents;
	private List<SimpleGREvent> indicaGREvents;
	private List<SimpleGREvent> japonicaGREvents;
	
	
	
	public GREventsTsvFileReaderForRicePop(Path inferredGREventsTsvFile) {
		super();
		this.inferredGREventsTsvFile = inferredGREventsTsvFile;
		
		///////////////////////////////
		this.run();
	}

	void run() {
		this.allGREvents=new ArrayList<>();
		this.indicaGREvents=new ArrayList<>();
		this.japonicaGREvents=new ArrayList<>();
		
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
				
				String variety=splits[17];
				this.allGREvents.add(gr);
				if(variety.equals("indica")) {
					this.indicaGREvents.add(gr);
				}else if(variety.equals("japonica")) {
					this.japonicaGREvents.add(gr);
				}
			}
			
			lineReader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	

	/**
	 * @return the allGREvents
	 */
	public List<SimpleGREvent> getAllGREvents() {
		return allGREvents;
	}

	/**
	 * @return the indicaGREvents
	 */
	public List<SimpleGREvent> getIndicaGREvents() {
		return indicaGREvents;
	}

	/**
	 * @return the japonicaGREvents
	 */
	public List<SimpleGREvent> getJaponicaGREvents() {
		return japonicaGREvents;
	}


}
