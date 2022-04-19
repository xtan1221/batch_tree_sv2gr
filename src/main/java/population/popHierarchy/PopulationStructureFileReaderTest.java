package population.popHierarchy;

import java.nio.file.Path;

public class PopulationStructureFileReaderTest {

	
	public static void main(String[] args) {
		Path samplePopulationInfoTableFile = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\gene_family_presence_absence_detection\\sorghum\\PoPoolationTE2_new_pipe\\output_data\\sorghum.27.sample.sample_name.sample_index.table.txt");		
		
		PopulationStructureFileReader ps=new PopulationStructureFileReader(samplePopulationInfoTableFile,false);
		
		System.out.println(ps.getOrderedSampleIndices());
	}
}

