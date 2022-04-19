package phylo.alignment;

import java.util.Map;

/**
 * wrapper class for a set of multiple aligned sequences
 * @author tanxu
 *
 */
public class MultipleAlignment {
	private final Map<String, String> seqNameAlignmentSeqMap;
	
	private final int seqNum;
	private final int alignmentLen;
	
	/**
	 * 
	 * @param seqNameAlignmentSeqMap
	 */
	public MultipleAlignment(Map<String, String> seqNameAlignmentSeqMap){
		if(seqNameAlignmentSeqMap==null||seqNameAlignmentSeqMap.isEmpty())
			throw new IllegalArgumentException("given seqNameAlignmentSeqMap cannot be null or empty!");
		
		Integer alignLen=null;
		for(String seqName:seqNameAlignmentSeqMap.keySet()) {
			int seqLen=seqNameAlignmentSeqMap.get(seqName).length();
			if(seqNameAlignmentSeqMap.get(seqName).length()==0)
				throw new IllegalArgumentException("sequence length cannot be 0!");
			
			if(alignLen==null)
				alignLen=seqLen;
			else
				if(seqLen!=alignLen)
					throw new IllegalArgumentException("sequence are not of the same length!"+alignLen+" "+seqLen+" "+seqName);
		}
		
		this.seqNameAlignmentSeqMap=seqNameAlignmentSeqMap;
		this.seqNum=this.seqNameAlignmentSeqMap.size();
		this.alignmentLen=alignLen;
	}

	/**
	 * @return the seqNameAlignmentSeqMap
	 */
	public Map<String, String> getSeqNameAlignmentSeqMap() {
		return seqNameAlignmentSeqMap;
	}

	/**
	 * @return the seqNum
	 */
	public int getSeqNum() {
		return seqNum;
	}

	/**
	 * @return the alignmentLen
	 */
	public int getAlignmentLen() {
		return alignmentLen;
	}
	
	
}
