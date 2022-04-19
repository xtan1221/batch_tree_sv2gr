package genomics.chrom;

import java.nio.file.Path;

public class ChromLenReaderTest {

	
	public static void main(String[] args) {
		Path chromNameLenFile = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\gene_family_presence_absence_detection\\sorghum\\PoPoolationTE2_new_pipe\\output_data\\chrom1-10.length.txt");
		
		
		ChromLenReader reader = new ChromLenReader(chromNameLenFile,null);
		
		
		for(String chrName:reader.getSortedChromNames()) {
			System.out.println(chrName+"\t"+reader.getChromNameLengthMap().get(chrName));
		}
	}
}
