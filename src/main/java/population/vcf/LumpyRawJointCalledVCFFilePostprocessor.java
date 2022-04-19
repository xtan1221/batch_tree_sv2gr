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
public class LumpyRawJointCalledVCFFilePostprocessor {
	private final Path inVCF;
	private final Path outVCF;
	
	
	public LumpyRawJointCalledVCFFilePostprocessor(Path inVCF, Path outVCF) {
		super();
		this.inVCF = inVCF;
		this.outVCF = outVCF;
		
		
		if(this.outVCF.toFile().exists()) {
			System.out.println("outVCF file exists, delete it...");
			this.outVCF.toFile().delete();
		}
		
		//////////////////////
		this.run();
	}


	void run() {
		try {
		    BufferedWriter writer = new BufferedWriter(new FileWriter(outVCF.toFile()));
			////////////////////
			BufferedReader lineReader = new BufferedReader(new FileReader(this.inVCF.toFile()));
			String line = null;
			
			while ((line = lineReader.readLine()) != null) {
				if(line.isEmpty()) { //first line is id	family	order
					continue;
				}
				
				if(line.startsWith("#")) {
					writer.append(line);
					writer.newLine();
					continue;
				}
				
				//data line
				String[] splits = line.split("\\s+");
				
				String type=splits[4]; //<DEL>, <DUP>, <INV>, <INS>, <CNV>
				if(!type.contains("<")) {//breakpoint type SV, skip
					continue;
				}
				
				writer.append(line);
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
