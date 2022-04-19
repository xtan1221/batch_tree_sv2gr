package phylogeny;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class BetweenGroupMeanDistanceCalculatorTest {

	public static void main(String[] args) {
		Path file=Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\phylogeny\\test\\phylip.dist.matrix.test.txt");
		
		PhylipDistanceMatrixFileReader reader = new PhylipDistanceMatrixFileReader(file);
		
		List<String> group1=new ArrayList<>();
		group1.add("10946");
		group1.add("11282");
		List<String> group2=new ArrayList<>();
		group2.add("11282");
		
		BetweenGroupMeanDistanceCalculator calculator=new BetweenGroupMeanDistanceCalculator(reader, group1, group1);
		
		System.out.println(calculator.getMeanDist());
	}
}
