package genomics.utils;

import java.util.ArrayList;
import java.util.List;

public class GenomicRegionUtilsTest {
	
	
	
	static void testMergeOverlappingRegions() {
		List<SimpleGenomicRegion> regions=new ArrayList<>();
		regions.add(new SimpleGenomicRegion("1", 1, 100));
		regions.add(new SimpleGenomicRegion("1", 90, 105));
		regions.add(new SimpleGenomicRegion("2", 1, 100));
		regions.add(new SimpleGenomicRegion("1", 10, 100));
		regions.add(new SimpleGenomicRegion("1", 900, 1025));
		
		List<SimpleGenomicRegion> merged = GenomicRegionUtils.mergeOverlappingRegions(regions);
		
		merged.forEach(r->{
			System.out.println(r.toString());
		});
	}
	
	public static void main(String[] args) {
		testMergeOverlappingRegions();
	}
}
