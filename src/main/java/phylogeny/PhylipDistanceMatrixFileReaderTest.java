package phylogeny;

import java.nio.file.Path;

public class PhylipDistanceMatrixFileReaderTest {
	
	
	
	
	public static void main(String[] args) {
		Path file=Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\phylogeny\\test\\phylip.dist.matrix.test.txt");
		
		PhylipDistanceMatrixFileReader reader = new PhylipDistanceMatrixFileReader(file);
		
		System.out.println(reader.getIndividualNum());
		reader.getIDs().forEach(id->{
			System.out.println(id);
		});
		
		System.out.println(reader.lookupDistance("11652","11652"));
		
	}
}
