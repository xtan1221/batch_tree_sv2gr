package phylo.tree.dist.diff.merge;

import java.io.IOException;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.function.Predicate;

import phylo.ref.Region;
import phylo.tree.dist.DistMatrixUtils;

public class RegionDiffMatrixAndNonMissingSiteNumMatrixFilesMergerTest {
	
	public static void main(String[] args) {
		
//		Path diffMatrixOutDir = Path.of("/scratch/tanxu/reseq/test/sv2gr/batch_matrix/sb.chrom1-10.100k.region.diff.matrix/region_diff_matrix_out/");
		Path diffMatrixOutDir = Path.of("C:\\Users\\tanxu\\Desktop\\scratch\\jc\\region_diff_matrix_out\\");
		
		String diffMatrixFileSuffix ="_totalSites_diffMatrix.txt";
		
//		Path regionIndexFile = Path.of("/scratch/tanxu/reseq/test/sv2gr/batch_matrix/sb.chrom1-10.100k.region.diff.matrix/region_index.txt");
		Path regionIndexFile = Path.of("C:\\Users\\tanxu\\Desktop\\scratch\\jc\\region_index.txt");
		
		try {
			Map<String,Predicate<Region>> outputDataNameRegionFilterMap = new LinkedHashMap<>();
			
			outputDataNameRegionFilterMap.put("all_chrom", e->{return true;});
			outputDataNameRegionFilterMap.putAll(MergerUtils.makeChromNamePredicateMap(regionIndexFile, 0));
			
			int seqNum = 25;
			
			RegionDiffMatrixAndNonMissingSiteNumMatrixFilesMerger merger = 
					new RegionDiffMatrixAndNonMissingSiteNumMatrixFilesMerger(
						diffMatrixOutDir, regionIndexFile, 
						outputDataNameRegionFilterMap,
						seqNum);
			
			merger.run();
			
			merger.getDataNameSummedPairwiseNonMissingSitesNumMatrixMap().forEach((dn,len)->{
				System.out.println(dn + "=="+len);
				
				int[][] diffMatrix=merger.getDataNameSummedPairwiseDiffMatrixMap().get(dn);
				
				DistMatrixUtils.printMatrix(DistMatrixUtils.toStringMatrix(diffMatrix));
				
			});
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		
		
		System.out.println("done");
	}
}
