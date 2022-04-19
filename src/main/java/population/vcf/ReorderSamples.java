package population.vcf;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * 
 * 
 * 1. filter out all breakend type SVs
 * 
 * @author tanxu
 * 
 */
public class ReorderSamples {
	private final Path inVCF;
	private final List<String> reorderedSampleNameList;
	private final Path outVCF;
	
	/**
	 * map from sample name to the col index in the {@link #inVCF}
	 */
	private Map<String, Integer> sampleNameColIndexMap;
	
	public ReorderSamples(Path inVCF, List<String> reorderedSampleNameList, Path outVCF) {
		super();
		this.inVCF = inVCF;
		this.reorderedSampleNameList = reorderedSampleNameList;
		this.outVCF = outVCF;
		
		
		if(this.outVCF.toFile().exists()) {
			System.out.println("outVCF file exists, delete it...");
			this.outVCF.toFile().delete();
		}
		
		//////////////////////
		this.run();
	}


	

	void run() {
		this.sampleNameColIndexMap=new HashMap<>();
		try {
		    BufferedWriter writer = new BufferedWriter(new FileWriter(outVCF.toFile()));
			////////////////////
			BufferedReader lineReader = new BufferedReader(new FileReader(this.inVCF.toFile()));
			String line = null;
			
			while ((line = lineReader.readLine()) != null) {
				if(line.isEmpty()) { //first line is id	family	order
					continue;
				}
				
				if(line.startsWith("#")&& !line.startsWith("#CHROM")) {
					writer.append(line);
					writer.newLine();
					continue;
				}
				
				
				String[] splits = line.split("\\s+");
				
				if(line.startsWith("#CHROM")) {//header line for column names
					//#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2	...
					//keep first 9 columns
					//reorder columns starting from 10
					for(int i=9;i<splits.length;i++) {
						this.sampleNameColIndexMap.put(splits[i], i);
					}
				}
				
				//add the first 9 columns without reordering
				StringBuilder sb=new StringBuilder();
				for(int i=0;i<9;i++) {
					sb.append(splits[i]).append("\t");
				}
				//add the columns for samples with reordering
				for(int i=9;i<splits.length;i++) {
					String newSample=this.reorderedSampleNameList.get(i-9);
					int index=this.sampleNameColIndexMap.get(newSample);
					
					sb.append(splits[index]).append("\t");
				}
				
				writer.append(sb.toString().trim());
				writer.newLine();
				
			}

			lineReader.close();
			writer.flush();
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
