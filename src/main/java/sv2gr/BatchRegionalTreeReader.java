package sv2gr;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import phylo.tree.reader.NewickFileFormatType;
import sv2gr.tree.RegionalTree;
import sv2gr.tree.Tree;

public class BatchRegionalTreeReader {
	/**
	 * 
	 */
	private final Path batchRegionalTreeFile;
	/**
	 * 
	 */
	private final int windowSize;
	
	////////////////////////////////
	private List<RegionalTree> regionalTrees;
	
	
	public BatchRegionalTreeReader(Path batchRegionalTreeFile, int windowSize) {
		super();
		this.batchRegionalTreeFile = batchRegionalTreeFile;
		this.windowSize=windowSize;
		
		
		this.read();
	}
	
	
	void read() {
		this.regionalTrees = new ArrayList<>();
		
		try {
			BufferedReader lineReader = new BufferedReader(new FileReader(this.batchRegionalTreeFile.toFile()));
			String line = null;
			
			while ((line = lineReader.readLine()) != null) {
				if(line.isEmpty() || line.startsWith("#")) { //header line starts with #
					continue;
				}
				//tree_id	chrom_regions	newick_string
				//29	Chr01:28000001-29000000	(3237:0.011546625,(((((((((11657:1.816E-5,(((13153:6.06E-6,Gp0018195:0.0)0.1:0.0,Gp0009948:0.0)0.1:0.0,13158:0.0)0.25:1.0E-8)0.24:4.0E-8,Gp0018209:1.207E-5)0.27:1.2E-7,(11652:3.036E-5,(13128:1.858E-5,Gp0018189:3.593E-5)0.88:1.203E-5)0.85:1.196E-5)0.32:3.5E-7,11647:3.588E-5)0.75:2.73E-6,(13068:0.0,Gp0018212:0.0)1.0:3.327E-5)0.81:1.079E-5,(((13148:0.0,Gp0009978:0.0)0.65:6.06E-6,Gp0018183:6.06E-6)0.35:8.0E-8,Gp0018197:5.443E-5)1.0:8.374E-5)1.0:6.383E-5,(((11757:9.089E-5,Gp0018186:1.812E-5)1.0:0.00176273,Gp0018192:0.00175222)1.0:2.976E-4,((13133:8.495E-5,Gp0018187:7.859E-5)1.0:6.923E-4,13138:6.0496E-4)1.0:8.1208E-4)0.93:9.817E-5)0.57:2.936E-5,13123:0.00275905)0.45:1.505E-5,13143:0.00234498):0.011546625);
				String[] splits=line.split("\\s+");
				int ID=Integer.parseInt(splits[0]);
				String[] chromRegionSplits=splits[1].split(":");
				String[] startEndSplits=chromRegionSplits[1].split("-");
				
				String chrom=chromRegionSplits[0];
				int start=Integer.parseInt(startEndSplits[0]);
				int end=Integer.parseInt(startEndSplits[1]);
				
				//////////
				String newickTreeString=splits[2];
				Tree tree = Tree.fromNewickString(newickTreeString, NewickFileFormatType.SIMPLE_NEWICK_2);
				
				
				RegionalTree regionalTree = new RegionalTree(ID, this.windowSize, chrom, start, end, tree);
				
//				System.out.println("read tree for "+chrom+" "+start+" "+end+" == "+newickTreeString);
				
				this.regionalTrees.add(regionalTree);
			}
			
			lineReader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}


	/**
	 * @return the regionalTrees
	 */
	public List<RegionalTree> getRegionalTrees() {
		return regionalTrees;
	}
	
	
}
