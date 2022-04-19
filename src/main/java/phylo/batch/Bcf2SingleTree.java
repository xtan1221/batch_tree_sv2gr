package phylo.batch;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;

import org.apache.commons.io.FileUtils;

import htsjdk.samtools.util.Log;
import htsjdk.tribble.TribbleException;
import htsjdk.variant.vcf.VCFFileReader;
import phylo.alignment.MultipleAlignment;
import phylo.exception.MultipleAlignmentLengthIsTooShortException;
import phylo.ref.Region;
import phylo.tree.phylo.MegaXBasedTreeBuilder;
import phylo.tree.phylo.MultipleAlignment2TreeBase;
import phylo.tree.reader.NewickFileFormatType;
import phylo.tree.reader.Tree;
import phylo.vcf.Vcf2AlignmentSeqBuilder;
import population.vcf.utils.GenotypeRecoderFactory;
import population.vcf.utils.VariantContextFilterFactory;

/**
 * this class will do
 * 1. build a phylogenetic tree based on data in a list of regions from a bcf file
 * 
 * 2. reroot the phylogenetic tree with a specific outgroup
 * 
 * 3. append the newick tree string of the rerooted tree and the unique identifier of the list of regions to a output file
 * 
 * @author tanxu
 * 
 */
public class Bcf2SingleTree implements Callable<Integer>{
	public static final Log log=Log.getInstance(Bcf2SingleTree.class);
	
	/////////////////////////////
	/**
	 * bcfFile from which the regional VCF is queried
	 * must be indexed;
	 */
	private final Path bcfFile;
	
	/**
	 * contains the specific strategy to filter sites in bcf file and recode genotype to sequence
	 */
	private final Vcf2AlignmentSeqBuilder vcf2AlignmentSeqBuilder;
	
	/**
	 * construct phylogenetic tree from a MultipleAlignment
	 */
	private final MultipleAlignment2TreeBase multipleAlignment2Tree;
	
	/**
	 * sample/individual name of the outgroup to reroot the phylogenetic tree;
	 * the name should be consistent with the name in the bcf file
	 */
	private final String outgroupName;
	
	/**
	 * the minimal alignment length so that the multiple alignment will be used to build tree
	 */
	private final int minAlignmentLen;
	
	//////////////////
	/**
	 * a list of disjoint regions SNP within which are to be used to build regional tree;
	 * see -r of bcftools
	 */
	private final String regionListString;
	
	/**
	 * a unique string identifier of the regionList that distinguish it from other BcfRegionalSNP2SingleTree running in the same batch;
	 * help to create temporary folder and name output files
	 */
	private final String uniqueRegionIdentifier;
	
	/**
	 * the file writer that append to file where to store the built phylgenetic tree's newick string;
	 * each line contains information for a tree delimited by tab
	 * 		first column is the uniqueRegionIdentifier
	 * 		second column is the regionListString
	 * 		third column is the newick tree string
	 */
	private final FileWriter outputTreeFileWriter;
	
	/**
	 * the directory where to put the temporary files and folders
	 */
	private final Path tmpDir;
	
	
	///////////////////
	/**
	 * the queried VCF file from the {@link #bcfFile} with {@link #disjointRegionList};
	 */
	private Path queriedVcfFile;
	/**
	 * 
	 */
	private MultipleAlignment msa;
	/**
	 * 
	 */
	private Tree tree;
	
	//////////////////////////////////////job status
	private Long startTime=null;
	private Long endTime=null;
	private boolean successfullyFinished = false;
	
	/**
	 * 
	 * @param bcfFile
	 * @param vcf2AlignmentSeqBuilder
	 * @param multipleAlignment2Tree
	 * @param outgroupName
	 * @param regionListString
	 * @param uniqueRegionIdentifier
	 * @param tmpDir
	 * @param outputTreeFile
	 */
	Bcf2SingleTree(
			Path bcfFile, 
			Vcf2AlignmentSeqBuilder vcf2AlignmentSeqBuilder, MultipleAlignment2TreeBase multipleAlignment2Tree,
			String outgroupName, int minAlignmentLen,
			String regionListString, String uniqueRegionIdentifier, FileWriter outputTreeFileWriter, Path tmpDir){
		
		//TODO validation
		if(outgroupName!=null&&outgroupName.isEmpty())
			throw new IllegalArgumentException("given outgroupName cannot be empty string!");
		
		
		/////////////////////////
		this.bcfFile = bcfFile;
		this.vcf2AlignmentSeqBuilder=vcf2AlignmentSeqBuilder;
		this.multipleAlignment2Tree = multipleAlignment2Tree;
		this.outgroupName = outgroupName;
		this.minAlignmentLen=minAlignmentLen;
		//////////
		this.regionListString = regionListString;
		this.uniqueRegionIdentifier= uniqueRegionIdentifier;
		this.outputTreeFileWriter=outputTreeFileWriter;
		this.tmpDir=tmpDir;
	}
	
	/**
	 * @return the uniqueRegionIdentifier
	 */
	public String getUniqueRegionIdentifier() {
		return uniqueRegionIdentifier;
	}
	
	/**
	 * return whether this job has started or not
	 * @return
	 */
	public boolean started() {
		return this.startTime!=null;
	}
	
	/**
	 * return the start time of this job in nanoseconds
	 * @return the startTime
	 */
	public Long getStartTime() {
		return startTime;
	}
	
	/**
	 * return the time length of how long this job has being running since it started in nanoseconds
	 * @return
	 */
	public long hasRunFor() {
		if(!this.started()) {
			throw new UnsupportedOperationException("this job is not started yet!");
		}
		return System.nanoTime()-this.startTime;
	}
	
	/**
	 * @return the successfullyFinished
	 */
	public boolean isSuccessfullyFinished() {
		return successfullyFinished;
	}
	
	/**
	 * 
	 * @return
	 */
	public long getFullSuccessfulRunTime() {
		if(!this.isSuccessfullyFinished()) {
			throw new UnsupportedOperationException("this job is not successfully done!");
		}
		return this.endTime-this.startTime;
	}
	
	/**
	 * invoked after this job is canceled and need to be re-submitted;
	 * reset this job for correctly tracking the status of this job and re-submit the job
	 * @throws IOException 
	 */
	public void reset() throws IOException {
		this.startTime=null;
		this.successfullyFinished=false;
		this.endTime=null;
		//kill the process that runs bash code to construct tree
		//note that 
		log.info("process pid=("+this.multipleAlignment2Tree.getProcess().pid()+") to build tree is alive? for tree id=="+this.uniqueRegionIdentifier+" "+this.multipleAlignment2Tree.getProcess().isAlive());
//				+";process to build tree exit code:"+this.multipleAlignment2Tree.getProcess().exitValue());
		
		Runtime.getRuntime().exec("kill -9 ".concat(Long.toString(this.multipleAlignment2Tree.getProcess().pid())));
		
//		this.multipleAlignment2Tree.getProcess().destroyForcibly();
		try {
			Thread.sleep(3000);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
//		this.multipleAlignment2Tree.getProcess().destroy();
		if(this.multipleAlignment2Tree.getProcess()!=null)
			log.info("process pid=("+this.multipleAlignment2Tree.getProcess().pid()+") to build tree after destroyForcibly is alive? for tree id=="+this.uniqueRegionIdentifier+" "+this.multipleAlignment2Tree.getProcess().isAlive());
		
		FileUtils.cleanDirectory(this.tmpDir.toFile());
	}
	
	
	
	//////////////////////////////////
	/**
	 * return 1 if terminated during step 1;
	 * return 2 if terminated during step 2;
	 * return 3 if terminated during step 3
	 * return 0 if all steps are finished
	 */
	@Override
	public Integer call() throws Exception {
		this.startTime=System.nanoTime();
//		this.successfullyFinished=false;
//		this.endTime=null;
//		//clear 
//		FileUtils.cleanDirectory(this.tmpDir.toFile());
		
//		Thread.sleep(10000);
		if(this.queryVCF()!=0) { //step 1
			return 1;
		}
		
		if(this.buildMSA()!=0) {//step 2
			return 2;
		}
		
		if(this.inferTree()!=0) {//step 3
			return 3;
		}
		
		if(this.rerootTree()!=0) {
			return 4;
		}
		
		if(this.appendToFile()!=0) {
			return 5;
		}
		
		////////
		this.endTime=System.nanoTime();
		this.successfullyFinished=true;
		log.info("job is successfully done for tree id=="+this.uniqueRegionIdentifier+"  running time = "+this.getFullSuccessfulRunTime());
		return 0;
	}
	
	/**
	 * return 0 if succeed; return 1 otherwise
	 */
	private int queryVCF() {
		log.info("start querying regional VCF from bcf file"+" for tree id=="+this.uniqueRegionIdentifier);
		this.queriedVcfFile = Path.of(this.tmpDir.toString(), this.uniqueRegionIdentifier.concat(".vcf"));
		
		BcfUtils.queryRegionIntoVcfFileFromIndexedBcfFile(this.bcfFile.toString(), this.regionListString, this.queriedVcfFile.toString());
		
		if(this.queriedVcfFile.toFile().exists()) {
			log.info("querying regional VCF from bcf file is successfully finished "+" for tree id=="+this.uniqueRegionIdentifier);
			return 0;
		}else {
			log.error("queried regional vcf file is not found"+" for tree id=="+this.uniqueRegionIdentifier);
			return 1;
		}
	}
	
	/**
	 * build MultipleAlignment
	 * return 0 if succeed;
	 * return 1 otherwise
	 * @throws MultipleAlignmentLengthIsTooShortException 
	 */
	private int buildMSA() {
		log.info("start building multiple sequence alignment"+" for tree id=="+this.uniqueRegionIdentifier);
		////
		try {
			VCFFileReader reader = new VCFFileReader(this.queriedVcfFile, false);
			this.vcf2AlignmentSeqBuilder.run(reader);
			
			try {
				this.msa=this.vcf2AlignmentSeqBuilder.getMultipleAlignment();
			}catch(IllegalArgumentException e) {
//				e.printStackTrace();
				log.error(e.getMessage()+" for tree id=="+this.uniqueRegionIdentifier);
				return 1;
			}
			
			if(this.msa.getAlignmentLen()<this.minAlignmentLen){
				log.error("built alignment is too short (built="+this.msa.getAlignmentLen()+"; min len="+this.minAlignmentLen+") for tree id=="+this.uniqueRegionIdentifier);
				return 1;
			}
			
			log.info("building multiple sequence alignment is successfully finished "+" for tree id=="+this.uniqueRegionIdentifier);
			return 0;
			
		}catch(TribbleException e) {
			log.error("Exception is thrown when when building multiple sequence alignment for tree id=="+this.uniqueRegionIdentifier+" :"+e.getMessage());
			return 1;
		}
	}
	
	/**
	 * return 0 if succeed; any non-0 value will be considered abnormal termination
	 */
	private int inferTree() {
		log.info("start inferring tree"+" for tree id=="+this.uniqueRegionIdentifier);
		try {
			this.multipleAlignment2Tree.build(this.uniqueRegionIdentifier, this.msa, this.tmpDir);
		} catch (Exception e) { //tree newick file is not successfully created
			e.printStackTrace();
			log.error(e.getMessage()+" for tree id=="+this.uniqueRegionIdentifier);
			return 1;
		}
		
		log.info("tree inferring is successfully finished "+" for tree id=="+this.uniqueRegionIdentifier);
		
		
		this.tree=this.multipleAlignment2Tree.getTree();
		
		return 0;
	}
	
	
	/**
	 * 
	 * @return
	 */
	private int rerootTree() {
		log.info("start reroot tree"+" for tree id=="+this.uniqueRegionIdentifier);
		try {
			this.tree = this.tree.reroot(this.outgroupName);
		}catch(Exception e) {
			log.error("Exception found when rerooting tree "+" for tree id=="+this.uniqueRegionIdentifier);
			e.printStackTrace();
			return 1;
		}
		
		log.info("tree rerooting is successfully finished "+" for tree id=="+this.uniqueRegionIdentifier);
		return 0;
	}
	
	/**
	 * write the full set of information of the tree to the output file by using the {@link #outputTreeFileWriter}
	 * @return
	 */
	private int appendToFile() {
		log.info("start append to final output tree file"+" for tree id=="+this.uniqueRegionIdentifier);
		StringBuilder lineSB = new StringBuilder();
		lineSB.append(this.uniqueRegionIdentifier).append("\t").append(this.regionListString).append("\t").append(this.tree.toFullNewickString(NewickFileFormatType.SIMPLE_NEWICK_2));
		
		synchronized(this.outputTreeFileWriter)
		{
			log.info("write newick tree string to output tree file "+" for tree id=="+this.uniqueRegionIdentifier);
			BufferedWriter writer = new BufferedWriter(this.outputTreeFileWriter);
            // synchronizing the outputTreeFileWriter object
            try {
				writer.append(lineSB.toString());
				writer.newLine();
				writer.flush();
			} catch (IOException e) {
				e.printStackTrace();
				return 1;
			}
		}
		
		log.info("appending to final output file is successfully finished "+" for tree id=="+this.uniqueRegionIdentifier);
		return 0;
	}
	
	public static void main(String[] args) throws Exception {
		Path bcfFile = Path.of("/scratch/tanxu/reseq/test/batch_local_tree/sorghum.1.versi.1.timo.1.prop.all.24.bicolor.jointly.called.raw.first1000000lines.bcf");
		String outgroupName = "3237";
		int minAlignmentLen = 100;
		List<Region> disjointRegionList = new ArrayList<>(); 
		disjointRegionList.add(new Region("Chr01", 1, 100000));
		String uniqueRegionIdentifier = "test";
//		Path workingDir = Path.of("/scratch/tanxu/reseq/test/batch_local_tree/working");
		Path outputDir=Path.of("/scratch/tanxu/reseq/test/batch_local_tree/output");
		
		Path outputTreeFile = Path.of("");//TODO
		
		int maxMissingGenotype = 0;
		Vcf2AlignmentSeqBuilder builder = new Vcf2AlignmentSeqBuilder(
//				VariantContextFilterFactory.onlySNPSites().and(VariantContextFilterFactory.allIndividualsAreGenericHomozygous()).and(VariantContextFilterFactory.maxMissingCount(maxMissingGenotype)), 
				VariantContextFilterFactory.nonIndelSite().and(VariantContextFilterFactory.allIndividualsAreGenericHomozygous()).and(VariantContextFilterFactory.maxMissingCount(maxMissingGenotype)),
//				GenotypeRecoderFactory.biBaseSeqRecoder()
				GenotypeRecoderFactory.singleBaseSeqRecoder()
				);
		
		Path maoFile = Path.of("/home/tanxu/phylogeny/megaX/infer_NJ_nucleotide_pairwise_deletion.mao");
		
		MegaXBasedTreeBuilder treeBuilder = new MegaXBasedTreeBuilder(maoFile);
		
		FileWriter outputTreeFileWriter = new FileWriter(outputTreeFile.toFile(), true);
		
		Bcf2SingleTree tree=
				new Bcf2SingleTree(
						bcfFile, 
						builder, treeBuilder, 
						outgroupName, minAlignmentLen,
						BcfUtils.buildBcftoolsRegionString(disjointRegionList), uniqueRegionIdentifier, outputTreeFileWriter, outputDir);
		
		tree.call();
	}

}
