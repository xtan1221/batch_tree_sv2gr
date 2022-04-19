package phylo.tree.phylo;

import phylo.tree.phylo.PhylipBashRunner.DistanceModel;

public class PhylipBasedTreeBuilderFactory extends MultipleAlignment2TreeFactory<PhylipBasedTreeBuilder>{
	private final int outgroupIndex;
	private final DistanceModel distanceModel;
	
	/**
	 * constructor
	 * @param uniqueID
	 * @param msa
	 * @param tmpDir
	 * @param outgroupIndex
	 */
	public PhylipBasedTreeBuilderFactory(int outgroupIndex, DistanceModel distanceModel) {
		
		this.outgroupIndex=outgroupIndex;
		this.distanceModel=distanceModel;
	}
	
	/**
	 * do not specify the outgroup index, use the default one (1)
	 * @param dataName
	 * @param msa
	 * @param tmpDir
	 * @param distanceModel
	 */
	public PhylipBasedTreeBuilderFactory(DistanceModel distanceModel) {
		this(1, distanceModel);
	}

	
	@Override
	public PhylipBasedTreeBuilder make() {
		return new PhylipBasedTreeBuilder(this.outgroupIndex, this.distanceModel);
	}
	
	
}
