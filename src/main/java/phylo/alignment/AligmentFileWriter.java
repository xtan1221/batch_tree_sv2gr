package phylo.alignment;

import java.io.File;
import java.nio.file.Path;

public abstract class AligmentFileWriter {
	private final MultipleAlignment multipleAlignment;
	private final Path outputFileDir;
	private final String outFileBaseName;
	
	/**
	 * constructor
	 * @param multipleAlignment
	 * @param outputFileDir directory to put the generated alignment file
	 * @param outFileBaseName file base name of the generated alignment file, together with the file extension to build full file name
	 */
	AligmentFileWriter(MultipleAlignment multipleAlignment, Path outputFileDir, String outFileBaseName){
		
		if(multipleAlignment==null)
			throw new IllegalArgumentException("given multipleAlignment cannot be null or empty!");
		
		this.multipleAlignment=multipleAlignment;
		this.outputFileDir=outputFileDir;
		this.outFileBaseName=outFileBaseName;
		
		if(this.getOutputFile().exists()) {
			throw new IllegalArgumentException("given output file already exists!");
		}
	}
	
	abstract String getFileExtension();
	
	public abstract void write();
	
	/**
	 * return the output multiple alignment File
	 * @return
	 */
	public File getOutputFile() {
		return Path.of(this.getOutputFileDir().toString(), this.getOutFileBaseName().concat(this.getFileExtension())).toFile();
	}
	
	

	/**
	 * @return the outputFileDir
	 */
	public Path getOutputFileDir() {
		return outputFileDir;
	}

	/**
	 * @return the outFileBaseName
	 */
	public String getOutFileBaseName() {
		return outFileBaseName;
	}

	/**
	 * @return the multipleAlignment
	 */
	public MultipleAlignment getMultipleAlignment() {
		return multipleAlignment;
	}
}
