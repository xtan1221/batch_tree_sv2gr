package circos;

import java.io.IOException;
import java.nio.file.Path;


/**
 * writer for a data file containing numeric value for windows with following columns
 * 
 * 		col 1: chrom name
 * 		col 2: start
 * 		col 3: end
 * 		col 4: numeric value
 * 	
 * the data file can be used to create LINE plot, SCATTER plot, HISTOGRAM and HEAT MAP with circos
 * 
 * @author tanxu
 *
 */
public class WindowedNumericValueDataFileWriter extends CircosDataFileWriterBase<WindowedNumericValueRecord>{
	
	public WindowedNumericValueDataFileWriter(Path outFile) throws IOException {
		super(outFile);
		// TODO Auto-generated constructor stub
	}
	
	
	@Override
	protected String getDataLineString(WindowedNumericValueRecord record) {
		StringBuilder sb=new StringBuilder();
		sb.append(record.getChr()).append("\t")
		.append(record.getStart()).append("\t")
		.append(record.getEnd()).append("\t")
		.append(record.getValue());
		
		return sb.toString();
	}
}
