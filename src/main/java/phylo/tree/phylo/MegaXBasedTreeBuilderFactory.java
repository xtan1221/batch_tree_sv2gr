package phylo.tree.phylo;

import java.nio.file.Path;

public class MegaXBasedTreeBuilderFactory extends MultipleAlignment2TreeFactory<MegaXBasedTreeBuilder>{
	/**
	 * the file containing all settings/configurations of running MegaX from command line;
	 */
	private final Path maoFile;
	
	public MegaXBasedTreeBuilderFactory(Path maoFile){
		this.maoFile=maoFile;
	}
	
	@Override
	public MegaXBasedTreeBuilder make() {
		return new MegaXBasedTreeBuilder(this.maoFile);
	}
	
}
