package phylo.assembly;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import picard.sam.SamToFastq;

/**
 * a simple wrapper of {@link SamToFastq} for pair end reads;
 * 
 * all default settings of {@link SamToFastq} will be used;
 * 
 * note that based on the source code, {@link SamToFastq} does not require input bam file to be sorted, it will store all SAMRecords of the first encountered mate of pair end reads in a HashMap until the second mate is encountered
 * 
 * 
 * @author tanxu
 *
 */
public class SimpleSamToFastqWrapper extends SamToFastq{
	private final Path inputBamFile;
	private final Path firstEndFastqFile; 
	private final Path secondEndFastqFile;
	
	///////////
	private String[] samToFastqArgs;
	private SamToFastq samToFastq;
	
	/**
	 * 
	 * @param inputBamFile
	 * @param firstEndFastqFile
	 * @param secondEndFastqFile
	 */
	SimpleSamToFastqWrapper(Path inputBamFile, Path firstEndFastqFile, Path secondEndFastqFile){
		this.inputBamFile = inputBamFile;
		this.firstEndFastqFile = firstEndFastqFile;
		this.secondEndFastqFile=secondEndFastqFile;
		
		this.buildSamToFastqArgs();
		this.runSamToFastq();
	}
	
	/**
	 * build the arguments like the following:
	 * java -jar $EBROOTPICARD/picard.jar SamToFastq I=input.bam FASTQ=out_1.fastq SECOND_END_FASTQ=out_2.fastq
	 */
	private void buildSamToFastqArgs() {
		List<String> argList = new ArrayList<>();
		
		//
		argList.add("I=".concat(this.inputBamFile.toString())); //I=input.bam
		argList.add("FASTQ=".concat(this.firstEndFastqFile.toString())); //FASTQ=out_1.fastq
		argList.add("SECOND_END_FASTQ=".concat(this.secondEndFastqFile.toString())); //SECOND_END_FASTQ=out_2.fastq
		
		this.samToFastqArgs = argList.toArray(new String[argList.size()]);
		
		//for testing
//		for(String a:CommandLineSyntaxTranslater.convertPicardStyleToPosixStyle(this.samToFastqArgs)) {
//			System.out.println(a);
//		}
	}
	
	private void runSamToFastq() {
		this.samToFastq = new SamToFastq();
		this.samToFastq.instanceMainWithExit(this.samToFastqArgs);
	}
	
	
	public static void main(String[] args) {
		Path inputBamFile = Path.of("a");
		Path firstEndFastqFile = Path.of("b");
		Path secondEndFastqFile = Path.of("c");
		
		SimpleSamToFastqWrapper b2f=new SimpleSamToFastqWrapper(inputBamFile, firstEndFastqFile, secondEndFastqFile);
	}
}
