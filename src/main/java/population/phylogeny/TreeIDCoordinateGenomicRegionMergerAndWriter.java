package population.phylogeny;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import genomics.utils.SimpleGenomicRegion;
import population.phylogeny.AllTreeMDSCoordinateFileReader.Coordinate;


/**
 * 
 * combine the tree MDS coordinate and the regional tree's genomic region into a single output file for figure plot with R
 * 
 * the output file contains the following columns:
 * 
 * 1. tree_id 
 * 		for regional trees, tree id is the index number of the tree
 * 		for chromosome trees, tree id is the name of the chromsomoe
 * 		for full genome trees, tree id is the 'full_genome'
 * 2. type
 * 		type of the tree, 'full_genome', 'chrom', 'regional'
 * 3. x
 * 		the calculated x coordinate of the tree by MDS
 * 4. y
 * 		the calculated y coordinate of the tree by MDS
 * 5. chrom
 * 		the name of the chrom of the tree
 * 		for full genome tree, chrom is 'full_genome'
 * 6. start
 * 		start position of the regional tree on the corresponding chromosome
 * 		'N/A' for chromosmome trees and full genome trees
 * 7. end
 * 		end position of the regional tree on the corresponding chromosome
 * 		'N/A' for chromosmome trees and full genome trees
 * 8. missing data
 * 		number of missing data when building the full genome and chromosome trees
 * 			multiple missing data may be present for chromsome and full genome trees in the same dataset
 * 		'N/A' for regional trees (their missing data are all predefined and fixed)
 * 
 * 9. dist_to_chrom_tree
 * 		Euclidean distance from the current regional tree to the chrom tree on the 2D space (the one with the missing data num equal to {@link #selectedFullGenomeChromosomeTreeMissingDataNum})
 * 		'N/A' for chrom trees and full genome tree
 * 10. dist_to_full_genome_tree
 * 		Euclidean distance from the current regional/chrom tree to the full genome tree on the 2D space (the one with the missing data num equal to {@link #selectedFullGenomeChromosomeTreeMissingDataNum})
 * 		'N/A' for the full genome tree
 * 
 * 
 * @author tanxu
 *
 */
public class TreeIDCoordinateGenomicRegionMergerAndWriter {
	private final Path outFile;
	
	private final AllTreeMDSCoordinateFileReader allTreeMDSCoordinateFileReader;
	private final TreeIDGenomicRegionFileReader treeIDGenomicRegionFileReader;
	
	/**
	 * the number of allowed missing data of the full genome tree and chrom tree to calculate the columns 'dist_to_chrom_tree' and 'dist_to_full_genome_tree'
	 */
	private final int selectedFullGenomeChromosomeTreeMissingDataNum;
	/////////////////////////////////////
	
	//////////////////////
	private List<Integer> regionalTreeIDs;
	
	/**
	 * 
	 * @param outFile
	 * @param allTreeMDSCoordinateFileReader
	 * @param treeIDGenomicRegionFileReader
	 */
	public TreeIDCoordinateGenomicRegionMergerAndWriter(
			Path outFile, 
			AllTreeMDSCoordinateFileReader allTreeMDSCoordinateFileReader,
			TreeIDGenomicRegionFileReader treeIDGenomicRegionFileReader,
			int selectedFullGenomeChromosomeTreeMissingDataNum
			) {
		super();
		this.outFile = outFile;
		this.allTreeMDSCoordinateFileReader = allTreeMDSCoordinateFileReader;
		this.treeIDGenomicRegionFileReader = treeIDGenomicRegionFileReader;
		this.selectedFullGenomeChromosomeTreeMissingDataNum=selectedFullGenomeChromosomeTreeMissingDataNum;
		
		
		if(this.outFile.toFile().exists()) {
			System.out.println("outfile exists, delete it...");
			this.outFile.toFile().delete();
		}
		
		this.prepare();
		this.run();
	}

	void prepare() {
		this.regionalTreeIDs=new ArrayList<>();
		this.regionalTreeIDs.addAll(this.allTreeMDSCoordinateFileReader.getRegionalTreeIDCoordinateMap().keySet());
		
		Collections.sort(this.regionalTreeIDs);
	}
	
	void run() {
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(outFile.toFile()));
			writer.append(this.getHeaderLine());
		    writer.newLine();
		    //first process full genome trees
		    for(int missingDataNum:this.allTreeMDSCoordinateFileReader.getFullGenomeMissingDataNumCoordinateMap().keySet()) {
		    	Coordinate coord=this.allTreeMDSCoordinateFileReader.getFullGenomeMissingDataNumCoordinateMap().get(missingDataNum);
		    	
		    	StringBuilder sb=new StringBuilder();
		    	sb.append("full_genome").append("\t")//tree id
		    	.append("full_genome").append("\t")//type
		    	.append(coord.getX()).append("\t")
		    	.append(coord.getY()).append("\t")
		    	.append("full_genome").append("\t")//chrom
				.append("N/A").append("\t") //start
				.append("N/A").append("\t") //end
				.append(missingDataNum).append("\t")
				.append("N/A").append("\t")//dist_to_chrom_tree
				.append("N/A");//dist_to_full_genome_tree
		    	
		    	writer.append(sb.toString());
		    	writer.newLine();
		    }
		    
		    //process chrom trees
		    for(String chromName:this.allTreeMDSCoordinateFileReader.getChromNameMissingDataNumCoordinateMapMap().keySet()) {
		    	for(int missingDataNum:this.allTreeMDSCoordinateFileReader.getChromNameMissingDataNumCoordinateMapMap().get(chromName).keySet()) {
		    		Coordinate coord=this.allTreeMDSCoordinateFileReader.getChromNameMissingDataNumCoordinateMapMap().get(chromName).get(missingDataNum);
			    	double distance = 
			    			Coordinate.calculateEuclideanDistance(
			    					coord, 
			    					this.allTreeMDSCoordinateFileReader.getFullGenomeMissingDataNumCoordinateMap().get(this.selectedFullGenomeChromosomeTreeMissingDataNum)
			    					);
		    		
			    	StringBuilder sb=new StringBuilder();
			    	sb.append(chromName).append("\t")//tree id
			    	.append("chrom").append("\t")//type
			    	.append(coord.getX()).append("\t")
			    	.append(coord.getY()).append("\t")
			    	.append(chromName).append("\t")//chrom
					.append("N/A").append("\t") //start
					.append("N/A").append("\t") //end
					.append(missingDataNum).append("\t")
					.append("N/A").append("\t")//dist_to_chrom_tree
					.append(distance);//dist_to_full_genome_tree
			    	
			    	writer.append(sb.toString());
			    	writer.newLine();
		    	}
		    }
		    
		    //process regional trees
		    for(int treeID:this.regionalTreeIDs) {
		    	Coordinate coord = this.allTreeMDSCoordinateFileReader.getRegionalTreeIDCoordinateMap().get(treeID);
		    	SimpleGenomicRegion region=this.treeIDGenomicRegionFileReader.getTreeIDRegionMap().get(treeID);
		    	double distanceToFullGenomeTree = 
		    			Coordinate.calculateEuclideanDistance(
		    					coord, 
		    					this.allTreeMDSCoordinateFileReader.getFullGenomeMissingDataNumCoordinateMap().get(this.selectedFullGenomeChromosomeTreeMissingDataNum)
		    					);
		    	double distanceToChromTree = 
		    			Coordinate.calculateEuclideanDistance(
		    					coord, 
		    					this.allTreeMDSCoordinateFileReader.getChromNameMissingDataNumCoordinateMapMap().get(region.getChrom()).get(this.selectedFullGenomeChromosomeTreeMissingDataNum)
		    					);
	    		
		    	StringBuilder sb=new StringBuilder();
		    	sb.append(treeID).append("\t")//tree id
		    	.append("regional").append("\t")//type
		    	.append(coord.getX()).append("\t")
		    	.append(coord.getY()).append("\t")
		    	.append(region.getChrom()).append("\t")//chrom
				.append(region.getStart()).append("\t") //start
				.append(region.getEnd()).append("\t") //end
				.append("N/A").append("\t")//missing_data
				.append(distanceToChromTree).append("\t")//dist_to_chrom_tree
				.append(distanceToFullGenomeTree);//dist_to_full_genome_tree
		    	
		    	
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
	
	
	private String getHeaderLine() {
		StringBuilder sb=new StringBuilder();
		sb.append("tree_id").append("\t")
		.append("type").append("\t")
		.append("x").append("\t")
		.append("y").append("\t")
		.append("chrom").append("\t")
		.append("start").append("\t")
		.append("end").append("\t")
		.append("missing_data").append("\t")
		.append("dist_to_chrom_tree").append("\t")
		.append("dist_to_full_genome_tree");
		
		
		return sb.toString();
	}
}
