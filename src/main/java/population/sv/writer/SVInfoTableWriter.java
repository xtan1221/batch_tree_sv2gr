package population.sv.writer;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.List;

import population.sv.utils.SimpleSVLocus;

/**
 * write information of a set of given SVs to a table file with each line containing one SV
 * 
 * the table file can be used as input to perform statistical analysis with R
 * 
 * @author tanxu
 *
 */
public class SVInfoTableWriter {
	/**
	 * 
	 */
	private final List<SimpleSVLocus> svLocusList;
	
	private final Path outTableFile;
	
	/**
	 * 
	 * @param svLocusList
	 * @param outTableFile
	 */
	public SVInfoTableWriter(List<SimpleSVLocus> svLocusList, Path outTableFile) {
		super();
		this.svLocusList = svLocusList;
		this.outTableFile = outTableFile;
		
		if(this.outTableFile.toFile().exists()) {
			System.out.println("given outTableFile exists, delete it...");
			this.outTableFile.toFile().delete();
		}
		
		this.run();
	}
	
	
	void run() {
	    
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(this.outTableFile.toString()));
			
			writer.append(this.getHeaderLine());
			writer.newLine();
//			int i=0;
//			String outputString="";
			for(SimpleSVLocus locus:this.svLocusList) {
				StringBuilder sb = new StringBuilder();
				sb.append(locus.getChrom()).append("\t")
				.append(locus.getStart()).append("\t")
				.append(locus.getEnd()).append("\t")
				.append(locus.getEnd()-locus.getStart()+1).append("\t")
				.append(locus.getType());
				
//				i++;
				
//				System.out.println(sb.toString()+"\t"+i);
//				outputString=outputString.concat(sb.toString()).concat("\n");
				
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
		//chrom
		//start
		//end
		//type
		//
		StringBuilder sb = new StringBuilder();
		
		sb.append("chrom").append("\t")
		.append("start").append("\t")
		.append("end").append("\t")
		.append("size").append("\t")
		.append("type");
		
		return sb.toString();
	}
	
}
