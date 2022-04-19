package bed;

import java.nio.file.Path;

public class BedIntersectTest {

	public static void main(String[] args) {
		System.out.println("hello");
		Path bedFile1=Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\gene_family_presence_absence_detection\\sorghum\\PoPoolationTE2\\testing_data\\Sbicolor_454_v3.1.1.gene.bed");
		Path bedFile2=Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\gene_family_presence_absence_detection\\sorghum\\PoPoolationTE2\\testing_data\\nhmmer_against_full_ref_hits.bed");
		Path realPepNameList=Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\gene_family_presence_absence_detection\\sorghum\\PoPoolationTE2\\testing_data\\true_domain_containing_pep_id_list.txt");
		
		
		BedIntersect bi=new BedIntersect(bedFile1, bedFile2, realPepNameList);
		
		bi.readFiles();
		
		bi.analysis();
		
	}
}
