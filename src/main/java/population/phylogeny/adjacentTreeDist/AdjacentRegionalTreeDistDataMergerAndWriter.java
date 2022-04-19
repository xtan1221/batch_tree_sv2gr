package population.phylogeny.adjacentTreeDist;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;

import population.phylogeny.adjacentTreeDist.AdjacentRegionalTreeDistFileReader.AdjacentTreePairDist;
import population.phylogeny.adjacentTreeDist.AdjacentRegionalTreeDistStatisticalSummaryFileReader.ChromSummary;

/**
 * 
 * @author tanxu
 *
 */
public class AdjacentRegionalTreeDistDataMergerAndWriter {
	private final Path outFile;
	private final AdjacentRegionalTreeDistFileReader adjacentRegionalTreeDistFileReader;
	private final AdjacentRegionalTreeDistStatisticalSummaryFileReader adjacentRegionalTreeDistStatisticalSummaryFileReader;
	
	
	public AdjacentRegionalTreeDistDataMergerAndWriter(
			Path outFile,
			AdjacentRegionalTreeDistFileReader adjacentRegionalTreeDistFileReader,
			AdjacentRegionalTreeDistStatisticalSummaryFileReader adjacentRegionalTreeDistStatisticalSummaryFileReader) {
		super();
		this.outFile = outFile;
		this.adjacentRegionalTreeDistFileReader = adjacentRegionalTreeDistFileReader;
		this.adjacentRegionalTreeDistStatisticalSummaryFileReader = adjacentRegionalTreeDistStatisticalSummaryFileReader;
		
		/////////////////////////////////
		if(this.outFile.toFile().exists()) {
			System.out.println("outfile exists, delete it...");
			this.outFile.toFile().delete();
		}
		
		this.run();
	}
	
	void run() {
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(outFile.toFile()));
			writer.append(this.getHeaderLine());
		    writer.newLine();
		    
		    for(String chrom:this.adjacentRegionalTreeDistFileReader.getChromAdjacentTreePairDistsMap().keySet()) {
		    	
		    	for(AdjacentTreePairDist pair:this.adjacentRegionalTreeDistFileReader.getChromAdjacentTreePairDistsMap().get(chrom)) {
		    		StringBuilder sb = new StringBuilder();
		    		ChromSummary chromSummary=this.adjacentRegionalTreeDistStatisticalSummaryFileReader.getChromNameChromSummaryMap().get(chrom);
		    		
		    		sb.append(chrom).append("\t")
			    	.append(pair.getRegionalTree1ID()).append("\t")
			    	.append(pair.getRegionalTree2ID()).append("\t")
			    	.append(pair.getDistance()).append("\t")
			    	.append(pair.getRegionalTree1Start()).append("\t")
			    	.append(pair.getRegionalTree1End()).append("\t")
			    	.append(pair.getRegionalTree2Start()).append("\t")
			    	.append(pair.getRegionalTree2End()).append("\t")
			    	.append(chromSummary.getTotalPairNum()).append("\t")
			    	.append(chromSummary.getMeanDist()).append("\t")
			    	.append(chromSummary.getSd()).append("\t")
			    	.append(chromSummary.getVariance()).append("\t")
			    	.append(chromSummary.getTotalAdjacentPairNum()).append("\t")
			    	.append(chromSummary.getMeanDistOfAdjacentPairs()).append("\t")
			    	.append(chromSummary.getSdOfAdjacentPairs()).append("\t")
			    	.append(chromSummary.getVarianceOfAdjacentPairs());
		    		
		    		writer.append(sb.toString());
		    		writer.newLine();
		    	}
		    	
		    }
		    
		    writer.flush();
		    writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	private String getHeaderLine() {
		StringBuilder sb = new StringBuilder();
		
		sb.append("chrom").append("\t")
		.append("tree1_ID").append("\t")
		.append("tree2_ID").append("\t")
		.append("dist").append("\t")
		.append("tree1_start").append("\t")
		.append("tree1_end").append("\t")
		.append("tree2_start").append("\t")
		.append("tree2_end").append("\t")
		.append("total_tree_pair_num").append("\t")//total number of tree pairs on the chrom of the tree 1 and 2
		.append("total_tree_pair_mean_dist").append("\t")//
		.append("total_tree_pair_sd_dist").append("\t")
		.append("total_tree_pair_var_dist").append("\t")
		.append("total_adjacent_tree_pair_num").append("\t")//total number of adjacent tree pairs on the chrom of the tree1 and 2
		.append("total_adjacent_tree_pair_mean_dist").append("\t")
		.append("total_adjacent_tree_pair_sd_dist").append("\t")
		.append("total_adjacent_tree_pair_var_dist");
		
		return sb.toString();
	}
}
