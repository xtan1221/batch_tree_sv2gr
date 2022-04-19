package sv2gr;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Collections;
import java.util.List;

import sv2gr.gr.GREvent;
import sv2gr.gr.GREventUtils;

/**
 * 
 * @author tanxu
 *
 */
public class GREventsWriter {
	private final Path outDir;
	private final List<GREvent> grEvents;
	
	////////////////////////////////////
	private Path outFile;
	/**
	 * 
	 * @param outDir
	 * @param grEvents
	 */
	public GREventsWriter(Path outDir, List<GREvent> grEvents) {
		super();
		this.outDir = outDir;
		this.grEvents = grEvents;
		
		/////////////////////////////
		this.prepare();
		this.run();
	}
	
	void prepare() {
		this.outFile = Path.of(this.outDir.toString(), "GR.events.tsv");

		if(this.outFile.toFile().exists()) {
			System.out.println("outfile exists, delete it...");
			this.outFile.toFile().delete();
		}
		
		//////////////////////
		//sort GR events by their SV
		Collections.sort(
				this.grEvents, 
				(a,b)->{
					if(!a.getOriginalSVLocus().getChrom().equals(b.getOriginalSVLocus().getChrom())) {
						return a.getOriginalSVLocus().getChrom().compareTo(b.getOriginalSVLocus().getChrom());
					}else {
						return a.getOriginalSVLocus().getStart()-b.getOriginalSVLocus().getStart();
					}
				}
				);
	}
	
	
	void run() {
	    try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(outFile.toString()));
			writer.append(this.getHeaderLine());
			writer.newLine();
			
			
			for(GREvent gr:this.grEvents) {
				StringBuilder sb=new StringBuilder();
				
				sb
				//sv related
				.append(gr.getOriginalSVLocus().getChrom()).append("\t")
				.append(gr.getOriginalSVLocus().getType()).append("\t")
				.append(gr.getOriginalSVLocus().getStart()).append("\t")
				.append(gr.getOriginalSVLocus().getEnd()).append("\t")
				.append(gr.getOriginalSVLocus().getSize()).append("\t")
				//gr related
				.append(GREventUtils.inferGRType(gr.getOriginalSVLocus().getType(), gr.isReversalOfSV()).getStringValue()).append("\t")//GR type
				.append(gr.isReversalOfSV()).append("\t")
				.append(gr.getAncestralNodeAverageDistToLeafNodes()).append("\t")
				.append(gr.getDescendantNodeAverageDistToLeafNodes()).append("\t")
				.append(gr.getBranchNumOnDescendantToAncestralTreePath()).append("\t")
				.append(gr.getAssignedTreePath().getMaxBoostrap()==null?"N/A":gr.getAssignedTreePath().getMaxBoostrap()).append("\t") //max_bootstrap on the path
				//regional tree related
				.append(gr.getRegionalTree().getWindowSize()).append("\t")
				.append(gr.getRegionalTree().getID()).append("\t")
				.append(gr.getRegionalTree().getStart()).append("\t")
				.append(gr.getRegionalTree().getEnd()).append("\t")
				//mis
				.append(gr.getAssignedTreePath().getLeafSampleNameListStringOfAncestorNode()).append("\t") //anc_leaf_samples; the list of samples corresponding to the leaf nodes of the clade rooted with the SECOND node to ancestor node of the path where the GR event occured
				.append(gr.getAssignedTreePath().getLeafSampleNameListStringOfDescendantNode());//desc_leaf_samples; the list of samples corresponding to the leaf nodes of the clade rooted with the descendant node of the path where the GR event occured
				
				
				writer.append(sb.toString());
				writer.newLine();
			}
			
			writer.flush();
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	String getHeaderLine() {
		StringBuilder sb=new StringBuilder();
		
		sb
		//sv related
		.append("chrom").append("\t")
		.append("SV_type").append("\t")
		.append("SV_start").append("\t")
		.append("SV_end").append("\t")
		.append("size").append("\t") //size of the sv and GR?
		//gr related
		.append("GR_type").append("\t")
		.append("reverse_of_SV").append("\t")
		.append("anc_dist").append("\t") //average distance from ancestral node to all leaf nodes
		.append("desc_dist").append("\t")//average distance from descendant node to all leaf nodes
		.append("branch_num").append("\t")//number of branches on the path where the GR event occured
		.append("max_bootstrap").append("\t")//max bootstrap value of branch on the path where the GR event occured
		//regional tree related
		.append("window_size").append("\t")//window size of the regional tree on which the GR event was inferred
		.append("tree_id").append("\t")//regional tree id
		.append("tree_start").append("\t")
		.append("tree_end").append("\t")
		.append("anc_leaf_samples").append("\t")//the list of samples corresponding to the leaf nodes of the clade rooted with the SECOND node to ancestor node of the path where the GR event occured
		.append("desc_leaf_samples");//the list of samples corresponding to the leaf nodes of the clade rooted with the descendant node of the path where the GR event occured
		
		return sb.toString();
	}
}
