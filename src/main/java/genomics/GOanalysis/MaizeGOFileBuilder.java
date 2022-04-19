package genomics.GOanalysis;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * build gene id GO term file to perform GO enrichment analysis required by http://systemsbiology.cau.edu.cn/agriGOv2/c_SEA.php
 * 
 * @author tanxu
 * 
 */
public class MaizeGOFileBuilder {
	private final Path subsetGeneListFile;
	private final int geneIDColIndex;
	private final Path outputSubsetGeneGOTermsFile;
	/////////////////////////////
	private Path inFile=Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\data processing and analysis\\input data\\maize\\reference\\RefGen_v5(Zm-B73-REFERENCE-NAM-5.0) - used in this study\\annotation\\Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.interproscan.tsv");
	private Path outFile=Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\data processing and analysis\\input data\\maize\\reference\\RefGen_v5(Zm-B73-REFERENCE-NAM-5.0) - used in this study\\annotation\\gene.id.GO.terms.txt");
	
	
	private Map<String, Set<String>> idGOTermsMap;
	/**
	 * 
	 */
	public MaizeGOFileBuilder(Path subsetGeneListFile, int geneIDColIndex, Path outputSubsetGeneGOTermsFile){
		
		this.subsetGeneListFile=subsetGeneListFile;
		this.geneIDColIndex=geneIDColIndex;
		this.outputSubsetGeneGOTermsFile=outputSubsetGeneGOTermsFile;
		
		if(this.outFile.toFile().exists()) {
			this.outFile.toFile().delete();
		}
		
		if(this.outputSubsetGeneGOTermsFile.toFile().exists()) {
			this.outputSubsetGeneGOTermsFile.toFile().delete();
		}
		
		this.readGOTerms();
		this.writeToFile();
	}
	
	
	
	void readGOTerms() {
		this.idGOTermsMap=new HashMap<>();
		
		try {
			BufferedReader lineReader = new BufferedReader(new FileReader(this.inFile.toFile()));
			String line = null;
			
			while ((line = lineReader.readLine()) != null) {
				if(line.isEmpty()) { //first line is id	family	order
					continue;
				}
				
				String[] splits=line.split("\t");
				
				//transform from 'Zm00001eb019340_T003' to 'Zm00001eb019340'
				String id=splits[0].split("_")[0];
				
				if(splits.length<=13) {
					continue;
				}
				
				String[] goTerms=splits[13].split("\\|");
				

				if(!this.idGOTermsMap.containsKey(id)) {
					this.idGOTermsMap.put(id, new HashSet<>());
				}
				this.idGOTermsMap.get(id).addAll(Arrays.asList(goTerms));
			}

			lineReader.close();
			
			////////////////////
			BufferedWriter writer = new BufferedWriter(new FileWriter(outFile.toString()));
			for(String id:this.idGOTermsMap.keySet()) {
				for(String go:this.idGOTermsMap.get(id)) {
					writer.append(id.concat("\t").concat(go));
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
	
	
	void writeToFile() {
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(this.outputSubsetGeneGOTermsFile.toString()));
			
			BufferedReader lineReader = new BufferedReader(new FileReader(this.subsetGeneListFile.toString()));
			String line = null;
			
			while ((line = lineReader.readLine()) != null) {
				if(line.isEmpty() || line.startsWith("#")) { //first line is id	family	order
					continue;
				}
				
				String[] splits=line.split("\\s+");
				
				String id=splits[this.geneIDColIndex];
				
				if(this.idGOTermsMap.containsKey(id)) {
					for(String go:this.idGOTermsMap.get(id)) {
						writer.append(id.concat("\t").concat(go));
						writer.newLine();
					}
				}
			}
			lineReader.close();
			
			writer.flush();
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		Path subsetGeneListFile = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv\\maize\\GO_analysis\\output\\gene.list.overlapping.ingroup.SVs.tsv");
		int geneIDColIndex=0;
		Path outputSubsetGeneGOTermsFile=Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv\\maize\\GO_analysis\\output\\test.go.list.txt");
		MaizeGOFileBuilder builder = new MaizeGOFileBuilder(subsetGeneListFile, geneIDColIndex, outputSubsetGeneGOTermsFile);
	}
}
