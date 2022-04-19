package phylo.alignment;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;


/**
 * write multiple sequence alignment to a sequential phylip format file;
 * 
 * @author tanxu
 *
 */
public class SequentialPhylipAlignmentFileWriter extends AligmentFileWriter{
	public static final int SEQ_NAME_LENGTH=10;
	public static final String FILE_EXTENSION=".phy";
	
	
	/**
	 * map from the original seq names to the recoded seq names that are equal to or less than 10 characters long
	 * note that the recoded seq names should all be unique
	 */
	private final Map<String,String> seqNameRecodedSeqNameMap;
	
	
	/**
	 * 
	 * @param seqNameAlignmentSeqMap
	 * @param outputFilePath
	 * @param outFileBaseName
	 * @param seqNameRecodedSeqNameMap
	 */
	SequentialPhylipAlignmentFileWriter(MultipleAlignment multipleAlignment, Path outputFilePath, String outFileBaseName, Map<String,String> seqNameRecodedSeqNameMap) {
		super(multipleAlignment, outputFilePath, outFileBaseName);
		// TODO Auto-generated constructor stub
		
		for(String seqName:seqNameRecodedSeqNameMap.keySet()) {
			if(seqNameRecodedSeqNameMap.get(seqName).length()>SEQ_NAME_LENGTH)
				throw new IllegalArgumentException("given recoded seq name length cannot exceed "+ SEQ_NAME_LENGTH + " characters!");
		}
		
		this.seqNameRecodedSeqNameMap=seqNameRecodedSeqNameMap;
	}
	
	/**
	 * make a SequentialPhylipAlignmentFileWriter with original seq names;
	 * if seq name is longer than {@link SEQ_NAME_LENGTH}, trim it;
	 * no seq name uniqueness is checked!
	 * @param multipleAlignment
	 * @param outputFileDir
	 * @param outFileBaseName
	 * @return
	 */
	public static SequentialPhylipAlignmentFileWriter makeWithOriginalSeqName(MultipleAlignment multipleAlignment, Path outputFileDir, String outFileBaseName) {
		Map<String,String> seqNameRecodedSeqNameMap=new HashMap<>();
		
		for(String seqName:multipleAlignment.getSeqNameAlignmentSeqMap().keySet()) {
			if(seqName.length()>SEQ_NAME_LENGTH) {
				seqNameRecodedSeqNameMap.put(seqName, seqName.substring(0, 10)); //trim the first ten characters
			}else {
				seqNameRecodedSeqNameMap.put(seqName, seqName);
			}
		}
		
		return new SequentialPhylipAlignmentFileWriter(multipleAlignment, outputFileDir, outFileBaseName, seqNameRecodedSeqNameMap);
	}
	
	/**
	 * make a SequentialPhylipAlignmentFileWriter with default recoded seq names
	 * seq_index
	 * 
	 * @param multipleAlignment
	 * @param outputFileDir
	 * @param outFileBaseName
	 * @return
	 */
	static SequentialPhylipAlignmentFileWriter makeWithDefaultSeqNameRecoding(MultipleAlignment multipleAlignment, Path outputFileDir, String outFileBaseName) {
		Map<String,String> seqNameRecodedSeqNameMap=new HashMap<>();
		
		int index=1;
		for(String seqName:multipleAlignment.getSeqNameAlignmentSeqMap().keySet()) {
			String recodedName="seq".concat(Integer.toString(index));
			seqNameRecodedSeqNameMap.put(seqName, recodedName);
			index++;
		}
		
		return new SequentialPhylipAlignmentFileWriter(multipleAlignment, outputFileDir, outFileBaseName, seqNameRecodedSeqNameMap);
	}
	
	@Override
	String getFileExtension() {
		return FILE_EXTENSION;
	}
	
	@Override
	public void write() {
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(getOutputFile(), false));
			
			//write header line
			//
			writer.append("\t".concat(Integer.toString(this.getMultipleAlignment().getSeqNum())).concat("\t").concat(Integer.toString(this.getMultipleAlignment().getAlignmentLen())));
			writer.newLine();
			//write sequences  
			for(String seqName:this.getMultipleAlignment().getSeqNameAlignmentSeqMap().keySet()) {
				String paddedRecodedSeqName = StringUtils.rightPad(this.seqNameRecodedSeqNameMap.get(seqName), SEQ_NAME_LENGTH);
				
				String line = paddedRecodedSeqName.concat(this.getMultipleAlignment().getSeqNameAlignmentSeqMap().get(seqName));
				
				writer.append(line);
				writer.newLine();
			}
			
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
}
