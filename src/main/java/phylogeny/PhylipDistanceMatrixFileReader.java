package phylogeny;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/**
 * read in a phylip distance matrix file
 * 
 * first line is the number of individuals
 * 
 * the following lines have the first column containing the name/ID of the individual
 * 
 * 
 * @author tanxu
 *
 */
public class PhylipDistanceMatrixFileReader {
	private final Path distanceMatrixFile;
	
	/////////////////////////
	private int individualNum;
	private List<String> IDs;
	
	private Double[][] matrix;
	
	
	public PhylipDistanceMatrixFileReader(Path distanceMatrixFile) {
		super();
		this.distanceMatrixFile = distanceMatrixFile;
		
		this.read();
	}

	void read() {
		this.IDs=new ArrayList<>();
		
		try {
			////////////////////
			BufferedReader lineReader = new BufferedReader(new FileReader(this.distanceMatrixFile.toFile()));
			String line = null;
			
			boolean headerLineProcessed=false;
			
			while ((line = lineReader.readLine()) != null) {
				if(line.isEmpty()) { //first line is id	family	order
					continue;
				}
				if(!headerLineProcessed) {
					this.individualNum=Integer.parseInt(line.trim());
					headerLineProcessed=true;
					this.initializeDistMatrix();
				}else {
					String[] splits=line.split("\\s+");
					
					this.IDs.add(splits[0]);
					int individualIndex=this.IDs.size()-1;
					
					for(int i=1;i<splits.length;i++) {
						this.matrix[individualIndex][i-1]=Double.parseDouble(splits[i]);
					}
					
					
				}
				
			}

			lineReader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		/////////////////check
		if(this.IDs.size()!=this.individualNum) {
			throw new IllegalArgumentException("unequal individual number found!");
		}
		
		for(int i=0;i<this.matrix.length;i++) {
			for(int j=0;j<this.matrix[i].length;j++) {
				if(this.matrix[i][j]==null)
					throw new IllegalArgumentException("null matrix cell value found!");
			}
		}
	}
	
	/**
	 * initialize the matrix and set every cell value to null
	 */
	private void initializeDistMatrix() {
		this.matrix=new Double[this.individualNum][];
		
		for(int i=0;i<this.matrix.length;i++) {
			this.matrix[i]=new Double[this.individualNum];
			for(int j=0;j<this.matrix[i].length;j++) {
				this.matrix[i][j]=null;
			}
		}
	}

	/**
	 * @return the individualNum
	 */
	public int getIndividualNum() {
		return individualNum;
	}

	/**
	 * @return the iDs
	 */
	public List<String> getIDs() {
		return IDs;
	}
	
	/**
	 * @return the matrix
	 */
	public Double[][] getMatrix() {
		return matrix;
	}
	
	/**
	 * lookup and return the distance between the two given individual IDs
	 * @param id1
	 * @param id2
	 * @return
	 */
	public double lookupDistance(String id1, String id2) {
		int index1=this.IDs.indexOf(id1);
		int index2=this.IDs.indexOf(id2);
		
		return this.matrix[index1][index2];
	}
	//////////////////////
	
}
