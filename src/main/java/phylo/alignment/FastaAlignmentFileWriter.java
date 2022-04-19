package phylo.alignment;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;

public class FastaAlignmentFileWriter extends AligmentFileWriter{
	public static final String FILE_EXTENSION=".fasta";
	
	/**
	 * 
	 * @param seqNameAlignmentSeqMap
	 * @param outputFilePath
	 * @param outFileBaseName
	 */
	public FastaAlignmentFileWriter(MultipleAlignment multipleAlignment, Path outputFilePath, String outFileBaseName) {
		super(multipleAlignment, outputFilePath, outFileBaseName);
		// TODO Auto-generated constructor stub
	}
	
	@Override
	String getFileExtension() {
		return FILE_EXTENSION;
	}
	
	
	@Override
	public void write() {
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(getOutputFile(), false));
			
			for(String seqName:this.getMultipleAlignment().getSeqNameAlignmentSeqMap().keySet()) {
				String seq=this.getMultipleAlignment().getSeqNameAlignmentSeqMap().get(seqName);
				writer.append(">".concat(seqName));
				writer.newLine();
				writer.append(seq);
				writer.newLine();
			}
			
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
}
