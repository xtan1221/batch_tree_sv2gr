package phylo.vcf;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.stream.Stream;

import basic.Pair;
import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import phylo.alignment.MultipleAlignment;
import phylo.batch.BatchRegionalTreeManager;
import population.vcf.utils.GenotypeRecoder;
import population.vcf.utils.VariantContextFilter;


/**
 * sequence builder that transform records in a vcf file for each individual/sample into a sequence that can be used for phylogenetic tree inference for the group of individuals/samples
 * 
 * whether a specific site should be included or not is fully decided by the given {@link @VariantContextFilter};
 * 
 * note that the VCF file can contain non-variant sites where individuals' genotype are either homozygous reference or missing;
 * 
 * @author tanxu
 * 
 */
public class Vcf2AlignmentSeqBuilder {
	public static Log log=Log.getInstance(BatchRegionalTreeManager.class);
	
	///////////////////////
	/**
	 * only process until the first N sites (before applying site filterings)
	 */
	public static Integer firstN=null; //if non-null, only the firstN SNP sites will be processed; non-null only for testing
	
	//////////////
	private final GenotypeRecoder genotypeRecoder;
	
	private final VariantContextFilter variantFilter;
	
	/////////////
	private Map<String, String> sampleNameSequenceStringMap;
	
	private String ref;
	
	/**
	 * 
	 * @param genotypeRecoder
	 * @param variantFilter
	 */
	public Vcf2AlignmentSeqBuilder(VariantContextFilter variantFilter, GenotypeRecoder genotypeRecoder){
		//
		
		
		this.variantFilter = variantFilter;
		this.genotypeRecoder = genotypeRecoder;
		
	}
	
	
	/**
	 * transform the records in the given VCFFileReader into sequences with this Vcf2AlignmentSeqBuilder
	 * @param reader
	 */
	public void run(VCFFileReader reader) {
		this.sampleNameSequenceStringMap=new LinkedHashMap<>();
		reader.getHeader().getGenotypeSamples().forEach(s->{
			sampleNameSequenceStringMap.put(s, "");
		});
		this.ref="";
		
		////////
		Stream<VariantContext> stream=reader.iterator().stream();
		if(firstN!=null) {
			stream=stream.limit(firstN);
		}
		
		stream=stream.filter(this.variantFilter); //pass the given filter
		
		stream.forEach(v->{
			if(v.getStart()%100000==0)
				log.info("vcf file progress:"+v.getContig()+" "+v.getStart());
			
			v.getGenotypesOrderedByName().forEach(g->{
				this.sampleNameSequenceStringMap.put(g.getSampleName(), this.sampleNameSequenceStringMap.get(g.getSampleName()).concat(this.genotypeRecoder.apply(new Pair<>(v,g))));
			});
			
		});
	}
	
	
	void print() {
		sampleNameSequenceStringMap.forEach((k,v)->{
			
			System.out.println(v+"\t"+k);
		});
	}
	
	/**
	 * @return the sampleNameSequenceStringMap
	 */
	public Map<String, String> getSampleNameSequenceStringMap() {
		return sampleNameSequenceStringMap;
	}
	
	
	/**
	 * 
	 * @return
	 */
	public MultipleAlignment getMultipleAlignment() {
		return new MultipleAlignment(this.getSampleNameSequenceStringMap());
	}
}
