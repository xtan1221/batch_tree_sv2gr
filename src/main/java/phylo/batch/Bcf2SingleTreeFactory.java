package phylo.batch;

import java.io.FileWriter;
import java.nio.file.Path;

import phylo.tree.phylo.MultipleAlignment2TreeFactory;
import phylo.vcf.Vcf2AlignmentSeqBuilderFactory;



public class Bcf2SingleTreeFactory {
	/////////////////////////////
	/**
	 * bcfFile from which the regional VCF is queried
	 * must be indexed;
	 */
	private final Path bcfFile;
	
	/**
	 * contains the specific strategy to filter sites in bcf file and recode genotype to sequence
	 */
	private final Vcf2AlignmentSeqBuilderFactory vcf2AlignmentSeqBuilderFactory;
	
	/**
	 * construct phylogenetic tree from a MultipleAlignment
	 */
	private final MultipleAlignment2TreeFactory<?> msa2treeFactory;
	
	/**
	 * sample/individual name of the outgroup to reroot the phylogenetic tree;
	 * the name should be consistent with the name in the bcf file
	 */
	private final String outgroupName;
	
	/**
	 * the minimal alignment length so that the multiple alignment will be used to build tree
	 */
	private final int minAlignmentLen;
	/**
	 * 
	 * @param bcfFile
	 * @param vcf2AlignmentSeqBuilder
	 * @param multipleAlignment2Tree
	 * @param outgroupName
	 * @param tmpDir
	 * @param outputTreeFileWriter
	 */
	public Bcf2SingleTreeFactory(
			Path bcfFile, 
			Vcf2AlignmentSeqBuilderFactory vcf2AlignmentSeqBuilderFactory, 
			MultipleAlignment2TreeFactory<?> msa2treeFactory,
			String outgroupName,
			int minAlignmentLen){
		//TODO validation
		this.bcfFile = bcfFile;
		this.vcf2AlignmentSeqBuilderFactory=vcf2AlignmentSeqBuilderFactory;
		this.msa2treeFactory = msa2treeFactory;
		this.outgroupName = outgroupName;
		this.minAlignmentLen=minAlignmentLen;
	}
	

	
	/**
	 * 
	 * @param regionListString
	 * @param uniqueRegionIdentifier
	 * @param outputTreeFileWriter the file writer to the tree output file; shared by all Bcf2SingleTree2 running in the same batch thus synchronized
	 * @param tmpDir the temporary directory for all temporary folders and files generated during the process of the tree construction from bcf file
	 * @return
	 */
	public Bcf2SingleTree makeNew(String regionListString, String uniqueRegionIdentifier, FileWriter outputTreeFileWriter, Path tmpDir) {
		return new Bcf2SingleTree(
				this.bcfFile, this.vcf2AlignmentSeqBuilderFactory.make(), this.msa2treeFactory.make(),
				this.outgroupName, this.minAlignmentLen,
				regionListString, uniqueRegionIdentifier, outputTreeFileWriter, tmpDir
				);
	}
	
}
