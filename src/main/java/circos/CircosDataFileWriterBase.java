package circos;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;


/**
 * 
 * @author tanxu
 *
 * @param <T>
 */
public abstract class CircosDataFileWriterBase<T> {
	/**
	 * 
	 */
	protected final Path outFile;
	
	/////////////////////
	protected BufferedWriter writer;
	
	/**
	 * 
	 * @param outFile
	 * @throws IOException
	 */
	public CircosDataFileWriterBase(Path outFile) throws IOException {
		super();
		this.outFile = outFile;
		
		this.prepare();
	}
	
	private void prepare() throws IOException {
		if(this.outFile.toFile().exists()) {
			System.out.println("outFile already exists, delete it...");
			this.outFile.toFile().delete();
		}
		
		this.writer = new BufferedWriter(new FileWriter(this.outFile.toString()));
	}
	
	/**
	 * 
	 * @param record
	 * @return
	 */
	protected abstract String getDataLineString(T record);
	
	/**
	 * 
	 * @param record
	 * @throws IOException 
	 */
	public void add(T record) throws IOException {
		this.writer.append(this.getDataLineString(record));
		this.writer.newLine();
	}
	
	/**
	 * flush and close the writer;
	 * 
	 * @throws IOException 
	 * 
	 */
	public void close() throws IOException {
		this.writer.flush();
		this.writer.close();
	}
}
