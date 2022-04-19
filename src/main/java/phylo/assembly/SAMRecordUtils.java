package phylo.assembly;

import htsjdk.samtools.SAMRecord;

public class SAMRecordUtils {
	
	/**
	 * representative: no 2048 flag
	 * @param rec
	 * @return
	 */
	public static boolean isPrimaryAndRepresentative(SAMRecord rec) {
		
		return !rec.isSecondaryOrSupplementary();
		
	}
}
