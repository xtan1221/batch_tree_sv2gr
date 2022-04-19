package phylo.assembly;

import java.nio.file.Path;

public class SimpleSamToFastqWrapperTest {
	
	

	public static void main(String[] args) {
		Path inputBamFile = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\assembly\\reads extraction\\SRR486615.Chr01.1972000.1983000.sorted.bam");
		Path firstEndFastqFile = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\assembly\\reads extraction\\SRR486615.Chr01.1972000.1983000.1.fastq");
		Path secondEndFastqFile = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\assembly\\reads extraction\\SRR486615.Chr01.1972000.1983000.2.fastq");
		
		SimpleSamToFastqWrapper b2f=new SimpleSamToFastqWrapper(inputBamFile, firstEndFastqFile, secondEndFastqFile);
	}
}
