package phylo.tree.dist.nj;

import java.io.IOException;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.function.Predicate;

import phylo.ref.Region;
import phylo.tree.dist.diff.merge.MergerUtils;

public class DiffMatrixToNJTreeTest {

	
	public static void main(String[] args) throws IOException {
//		Path diffMatrixOutDir = Path.of("C:\\Users\\tanxu\\Desktop\\scratch\\jc\\region_diff_matrix_out\\");
//		String diffMatrixFileSuffix ="_totalSites_diffMatrix.txt";
//		Path regionIndexFile = Path.of("C:\\Users\\tanxu\\Desktop\\scratch\\jc\\region_index.txt");
//		Path outputDir = Path.of("C:\\Users\\tanxu\\Desktop\\scratch\\jc\\jc_dist_matrix");
		
		Path diffMatrixOutDir = Path.of("/scratch/tanxu/reseq/test/sv2gr/batch_matrix/sb.chrom1-10.100k.region.diff.matrix/region_diff_matrix_out/");
		Path regionIndexFile = Path.of("/scratch/tanxu/reseq/test/sv2gr/batch_matrix/sb.chrom1-10.100k.region.diff.matrix/region_index.txt");
		Path outputDir = Path.of("/scratch/tanxu/reseq/test/sv2gr/batch_matrix/sb.chrom1-10.100k.region.diff.matrix/region_jc_dist_based_nj_tree/");
		
		int seqNum = 25;
		String outgroupName = "3237";
		
		
		///////////////////////
		Map<String,Predicate<Region>> outputDataNameRegionFilterMap = new LinkedHashMap<>();
		outputDataNameRegionFilterMap.put("all_chrom", e->{return true;});
		outputDataNameRegionFilterMap.putAll(MergerUtils.makeChromNamePredicateMap(regionIndexFile, 0));
		
		
		DiffMatrixToNJTree diffMatrixToNJTree = new DiffMatrixToNJTree(
				diffMatrixOutDir, regionIndexFile, 
				outputDataNameRegionFilterMap,
				seqNum,
				outgroupName, outputDir);
		
		
		//////////////////
//		long startTime=System.nanoTime();
		diffMatrixToNJTree.run();
		//
//		System.exit(0);
		
	}
}
