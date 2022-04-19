package misc;

import java.nio.file.Path;

public class MaizeAccessionDownSamplingTest {

	public static void main(String[] args) {
		Path metadataFile=Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\data processing and analysis\\input data\\maize\\reseq\\HapMap 3.2.1 wgs reseq dataset\\HapMap3.2.1_NCBI_accession_metadata(merged by Xu).tsv");
		Path outputTSVFile =Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\data processing and analysis\\input data\\maize\\reseq\\HapMap 3.2.1 wgs reseq dataset\\HapMap3.2.1_NCBI_accession_metadata(merged by Xu)_highest_base_num_accession_per_genotype.tsv");
		int dataSetColIndex = 2;
		int accessionColIndex = 3;
		int baseNumColIndex = 4;
		int bioSampleColIndex = 7;
		int genotypeColIndex = 11;
		int libraryNameColIndex = 12;
		//////
		MaizeAccessionDownSampling ds = new MaizeAccessionDownSampling(
				metadataFile,outputTSVFile,
				dataSetColIndex, accessionColIndex, baseNumColIndex, bioSampleColIndex, genotypeColIndex, libraryNameColIndex
				);
		
		
		ds.readFile();
		ds.process();
		ds.print();
		ds.outputToFile();
	}
}
