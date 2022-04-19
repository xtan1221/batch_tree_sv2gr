package phylo.tree.distance_ditched;

import java.nio.file.Path;

import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProcessExecutor;

public class TreeDistanceRscriptRunner {
	public static final Log log = Log.getInstance(TreeDistanceRscriptRunner.class);
	
	static final Path TREE_DIST_CALCULATOR_R_SCRIPT_PATH=Path.of("/home/tanxu/phylogeny/tree_dist/Generalized_Robinson_Foulds_distances_2.r"); 
	
	
	/**
	 * run the R script to calculate the distance between the two given trees;
	 * then output the calculated distance and the tree ids to the given output file;
	 * @param tree1ID
	 * @param tree1NewickString
	 * @param tree2ID
	 * @param tree2NewickString
	 * @param outputFilePathString
	 * @return
	 */
	static int run(String tree1ID, String tree1NewickString, String tree2ID, String tree2NewickString, String outputFilePathString) {
		String commandLineString="Rscript ".concat(TREE_DIST_CALCULATOR_R_SCRIPT_PATH.toString()).concat(" ")
				.concat(tree1ID).concat(" ")
				.concat(tree1NewickString).concat(" ")
				.concat(tree2ID).concat(" ")
				.concat(tree2NewickString).concat(" ")
				.concat(outputFilePathString).concat(" ");
		
		log.info("run command to calcualte pairwise tree distance: '".concat(commandLineString).concat("'"));
		
		//
		return ProcessExecutor.execute(commandLineString);
	}
	
	public static void main(String[] args) {
		String tree1ID = "tree1";
		String tree1NewickString="(3237:0.005875445,((((((((((((((11647:1.107E-5,Gp0018183:1.015E-5)0.164:2.2E-7,Gp0018209:0.0)0.394:1.9E-7,13133:2.182E-5)0.13:2.6E-7,Gp0018212:0.0)0.046:2.6E-7,13068:0.0)0.136:5.2E-7,Gp0018195:0.0)0.258:1.22E-6,13148:0.0)0.296:3.62E-6,13158:0.0)0.348:3.69E-6,Gp0009948:0.0)0.558:1.326E-5,Gp0018186:2.105E-5)1.0:2.7384E-4,(13123:4.7793E-4,13153:9.17E-5)0.762:1.4888E-4)1.0:3.6125E-4,((((((11652:2.883E-5,(11657:0.0,11757:3.04E-6)0.98:1.857E-5)0.998:1.87E-5,13128:1.532E-5)0.712:1.7E-5,Gp0009978:0.0)0.586:1.141E-5,Gp0018189:1.015E-5)0.698:1.257E-5,Gp0018192:2.079E-5)0.286:4.5E-7,Gp0018197:1.01E-5)1.0:1.6131E-4)1.0:1.9884E-4,(13138:3.6014E-4,Gp0018187:5.35E-6)1.0:4.2105E-4)1.0:0.00221341,13143:0.00374712):0.005875445);"; 
		String tree2ID = "tree2";
		String tree2NewickString = "(3237:0.004332975,((((((((((((Gp0018186:9.7E-7,Gp0018209:1.399E-5)0.946:1.51E-6,Gp0018212:0.0)0.468:1.12E-6,((13068:0.0,(13158:0.0,(Gp0009948:1.4757E-4,Gp0018195:0.0)1.0:4.222E-5)0.878:8.12E-6)0.618:1.53E-6,13148:0.0)0.2:7.1E-7)0.632:7.1E-7,13133:0.0)0.558:1.75E-6,13123:1.229E-5)0.796:4.3E-6,13128:2.95E-6)0.776:3.06E-6,((11647:1.457E-5,((11657:1.2E-7,11757:0.0)1.0:9.46E-5,((Gp0009978:6.87E-6,Gp0018197:3.68E-6)0.996:1.265E-5,Gp0018192:2.298E-5)0.778:2.148E-5)0.858:8.42E-6)0.958:5.55E-6,Gp0018183:7.87E-6)0.834:5.24E-6)0.756:6.74E-6,Gp0018189:1.052E-5)0.84:1.174E-5,(13138:0.0,Gp0018187:5.65E-6)0.63:2.455E-5)0.644:3.864E-5,13143:1.507E-5)1.0:5.3047E-4,13153:8.1894E-4)1.0:0.00134629,11652:0.00211289):0.004332975);"; 
		String outputFilePathString = "/home/tanxu/phylogeny/tree_dist/outfile";
		
		
		TreeDistanceRscriptRunner.run(tree1ID, tree1NewickString, tree2ID, tree2NewickString, outputFilePathString);
		System.exit(0);
	}
}
