package phylo.tree.dist;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang3.StringUtils;

import basic.Pair;

public class DistMatrixUtils {
	public static final int SEQ_NAME_LENGTH=10;
	
	/**
	 * 
	 * @param matrix
	 */
	public static void printMatrix(String[][] matrix) {
		for(int i=0;i<matrix.length;i++) {
			for(int j=0;j<matrix[i].length;j++) {
				System.out.print(matrix[i][j]+"	");
			}
			System.out.println();
		}
	}
	
	
	public static String[][] toStringMatrix(int[][] mat){
		String[][] ret = new String[mat.length][mat[0].length];
		
		for(int i=0;i<mat.length;i++) {
			for(int j=0;j<mat[i].length;j++) {
				ret[i][j] = Integer.toString(mat[i][j]);
			}
		}
		
		return ret;
	}
	
	public static String[][] toStringMatrix(double[][] mat){
		String[][] ret = new String[mat.length][mat[0].length];
		
		for(int i=0;i<mat.length;i++) {
			for(int j=0;j<mat[i].length;j++) {
				ret[i][j] = Double.toString(mat[i][j]);
			}
		}
		
		return ret;
	}
	
	/**
	 * write the given matrix to the given output file
	 * 		each lines contains the data from matrix with first column being the sequence name;
	 * 
	 * for example:
	 * seq1		4	5	2
	 * seq2		3	2	9
	 * seq3		0	18	2
	 * 
	 * @param outFile
	 * @param orderedNameList
	 * @param matrix
	 * @throws IOException 
	 */
	public static void writeMatrixToFile(Path outFile, List<String> orderedNameList, int[][] matrix) throws IOException {
		BufferedWriter writer = new BufferedWriter(new FileWriter(outFile.toFile(), true));
	    
	    for(int i=0;i<orderedNameList.size();i++) {
			String line = orderedNameList.get(i).concat("\t");
			
			String row = "";
			for(int j=0;j<matrix[i].length;j++) {
				if(!row.isEmpty()) {
					row = row.concat("\t");
				}
				row = row.concat(Integer.toString(matrix[i][j]));
			}
			
			line = line.concat(row);
			
			writer.append(line);
			writer.newLine();
		}
		
		//////////
		writer.flush();
		writer.close();
	}
	
	/**
	 * write the given matrix to the given output file
	 * 		each lines contains the data from matrix with first column being the sequence name;
	 * 
	 * for example:
	 * seq1		4	5	2
	 * seq2		3	2	9
	 * seq3		0	18	2
	 * @param outFile
	 * @param orderedNameList
	 * @param matrix
	 * @throws IOException
	 */
	public static void writeMatrixToFile(Path outFile, List<String> orderedNameList, long[][] matrix) throws IOException {
		BufferedWriter writer = new BufferedWriter(new FileWriter(outFile.toFile(), true));
	    
	    for(int i=0;i<orderedNameList.size();i++) {
			String line = orderedNameList.get(i).concat("\t");
			
			String row = "";
			for(int j=0;j<matrix[i].length;j++) {
				if(!row.isEmpty()) {
					row = row.concat("\t");
				}
				row = row.concat(Long.toString(matrix[i][j]));
			}
			
			line = line.concat(row);
			
			writer.append(line);
			writer.newLine();
		}
		
		//////////
		writer.flush();
		writer.close();
	}
	
	/**
	 * read from the file written by {@link #writeMatrixToFile} method
	 * 
	 * note that there is no header line for total sites and there is no padding for seq id
	 * 
	 * @throws IOException 
	 * 
	 */
	public static Pair<List<String>, int[][]> readFromMatrixFile(Path diffMatrixAndSeqLenFile) throws IOException{
		
	    BufferedReader lineReader = new BufferedReader(new FileReader(diffMatrixAndSeqLenFile.toFile()));
	    String line = null;
	    
//	    //read the seq length from first line
//	    line = lineReader.readLine();
//	    int seqLen=Integer.parseInt(line.trim());
	    
	    //read the diff matrix from following lines
	    List<String> seqNameList = new ArrayList<>();
	    int[][] diffMatrix=null;
	    int index = 0;
	    while ((line = lineReader.readLine()) != null) {
	    	if(line.isEmpty()) {
	    		continue;
	    	}
	    	String[] splits = line.split("\\s+");
	    	if(diffMatrix==null) {
	    		diffMatrix = new int[splits.length-1][splits.length-1];
	    	}
	    	
	    	seqNameList.add(splits[0]);
	    	
	    	for(int i=1;i<splits.length;i++) {
	    		diffMatrix[index][i-1] = Integer.parseInt(splits[i]);
	    	}
	    	index++;
	    }
	    lineReader.close();
	    
	    
	    /////
	    return new Pair<>(seqNameList, diffMatrix);
		
	}
	
	/**
	 * write the given distance matrix to a output file with the format of phylip distance matrix
	 * 	see https://evolution.gs.washington.edu/phylip/doc/distance.html
	 * The input format for distance data is straightforward. 
	 * 1. The first line of the input file contains the number of species. 
	 * 2. There follows species data, starting, as with all other programs, with a species name. 
	 * 		1. The species name is ten characters long, and must be padded out with blanks if shorter. 
	 * 		2. For each species there then follows a set of distances to all the other species (options selected in the programs' menus allow the distance matrix to be upper or lower triangular or square). The distances can continue to a new line after any of them. If the matrix is lower-triangular, the diagonal entries (the distances from a species to itself) will not be read by the programs. If they are included anyway, they will be ignored by the programs, except for the case where one of them starts a new line, in which case the program will mistake it for a species name and get very confused.
	 * @param outFile
	 * @param num
	 * @param matrix
	 * @throws IOException 
	 */
	public static void writeToPhylipDistMatrixFile(Path outFile, List<String> orderedNameList, String[][] matrix) throws IOException {
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(outFile.toFile(), true));
	    
		writer.append("\t").append(Integer.toString(orderedNameList.size()));
	    writer.newLine();
	    
	    ////////
		for(int i=0;i<orderedNameList.size();i++) {
			String paddedSeqName = StringUtils.rightPad(orderedNameList.get(i), SEQ_NAME_LENGTH);
			String line = paddedSeqName.concat("\t");
			
			String row = "";
			for(int j=0;j<matrix[i].length;j++) {
				if(!row.isEmpty()) {
					row = row.concat("\t");
				}
				row = row.concat(matrix[i][j]);
			}
			
			line = line.concat(row);
			
			writer.append(line);
			writer.newLine();
		}
		
		//////////
		writer.flush();
		writer.close();
	}
	
	
	/**
	 * calculate p-distance matrix from pairwise difference matrix and the total sequence length;
	 * 
	 * @param totalLen
	 * @param diffMatrix
	 * @return
	 */
	public static double[][] calculatePDistaFromDiffMatrix(int totalLen, int[][] diffMatrix){
		double[][] pDistMatrix = new double[diffMatrix.length][diffMatrix.length];
		
		for(int i=0;i<diffMatrix.length;i++) {
			for(int j=0;j<diffMatrix[i].length;j++) {
				pDistMatrix[i][j] = (double)diffMatrix[i][j]/totalLen;
			}
		}
		
		return pDistMatrix;
	}
	
	/**
	 * calculate p-distance matrix from pairwise difference matrix and the total sequence length;
	 * 
	 * @param totalLen
	 * @param diffMatrix
	 * @return
	 */
	public static double[][] calculatePDistaFromDiffMatrix(long[][] totalSitesNumMatrix, int[][] diffMatrix){
		double[][] pDistMatrix = new double[diffMatrix.length][diffMatrix.length];
		
		for(int i=0;i<diffMatrix.length;i++) {
			for(int j=0;j<diffMatrix[i].length;j++) {
				pDistMatrix[i][j] = (double)diffMatrix[i][j]/totalSitesNumMatrix[i][j];
			}
		}
		
		return pDistMatrix;
	}
	
	/**
	 * calculate J-C model distance matrix from p-distance matrix
	 * 
	 * see https://www.megasoftware.net/mega1_manual/Distance.html
	 * 		d = -3loge(1 - 4p/3)/4
	 * @param totalLen
	 * @param pDistMatrix
	 * @return
	 */
	public static double[][] calcualteJCDistMatrixFromPDistMatrix(double[][] pDistMatrix){
		double[][] jcDistMatrix=new double[pDistMatrix.length][pDistMatrix.length];
		for(int i=0;i<pDistMatrix.length;i++) {
			for(int j=0;j<pDistMatrix[i].length;j++) {
				jcDistMatrix[i][j] = -0.75*Math.log(1-(double)4*pDistMatrix[i][j]/3);
			}
		}
		return jcDistMatrix;
	}
	
	
	//////////////////////////////////////////
	public static void main(String[] args) {
		double pdist=1.0267572951105817E-4;
		
		double jcdist=-0.75*Math.log(1-(double)4/3*pdist);
		
		System.out.println(jcdist);
	}
	
}
