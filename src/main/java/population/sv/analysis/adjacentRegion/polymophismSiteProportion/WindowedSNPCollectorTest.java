package population.sv.analysis.adjacentRegion.polymophismSiteProportion;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class WindowedSNPCollectorTest {

	
	public static void main(String[] args) {
		Path allSitesVCFFile = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv\\sorghum\\SV_SNP_relatedness\\testing_data\\sorghum.27.sample.first1M.lines.vcf");
		////////////
		List<Integer> includedSampleIndexList = new ArrayList<>();
		includedSampleIndexList.add(1);
		includedSampleIndexList.add(2);
		includedSampleIndexList.add(3);
		includedSampleIndexList.add(4);
		
		Map<String, List<Window>> chromBuiltWindowsMap = new HashMap<>();
		chromBuiltWindowsMap.put("Chr01", new ArrayList<>());
		chromBuiltWindowsMap.put("Chr02", new ArrayList<>());
		chromBuiltWindowsMap.get("Chr01").add(new Window("Chr01", 154949, 154949));
//		chromBuiltWindowsMap.get("Chr01").add(new Window("Chr01", 154463, 154480));
//		chromBuiltWindowsMap.get("Chr01").add(new Window("Chr01", 154470, 154488));
//		chromBuiltWindowsMap.get("Chr02").add(new Window("Chr02", 1, 100));
		
		
		int minSampleNumWithNonMissingDataToInclude=3;
		
		
//		AllSiteVCFFileWindowedSNPCollector collector = new AllSiteVCFFileWindowedSNPCollector(allSitesVCFFile, includedSampleIndexList, chromBuiltWindowsMap, minSampleNumWithNonMissingDataToInclude);
//		
//		
//		System.out.println(collector.getTotalSiteNumWithNonMissingData());
//		System.out.println(collector.getTotalPolymorphismSiteNum());
//	
//		for(String chrom:collector.getChromBuiltWindowsMap().keySet()) {
//			for(Window w:collector.getChromBuiltWindowsMap().get(chrom)) {
//				System.out.println(w.toString());
//			}
//		}
	}
}
