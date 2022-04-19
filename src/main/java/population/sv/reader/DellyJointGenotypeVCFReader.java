package population.sv.reader;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Predicate;
import java.util.stream.Stream;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import population.sv.utils.SimpleSVLocus;
import population.sv.utils.SimpleSVType;
import population.utils.Genotype;
import population.vcf.utils.VariantContextFilter;


/**
 * reader for a vcf file containing jointly called SVs for a population by Delly
 * 
 * sv types of delly
 * 		INS: length of inserted sequence is not accurate??????
 * 		DEL:
 * 		INV:
 * 		DUP: all are found to be tandem duplications!?
 * 		BND: no record find to be of this type!!!
 * 
 * only sv with FILTER column being 'PASS' are included
 * 
 * @author tanxu
 * 
 */
public class DellyJointGenotypeVCFReader {
	private final Path vcfFile;
	/**
	 * filter using the {@link VariantContext}
	 * for example
	 * 		only include biallelic locus
	 * 
	 */
	private final VariantContextFilter variantContextFilter;
	
	/**
	 * filter by the built SimpleSVLocus
	 * for example
	 * 		filter these SV with length > threshold 
	 * 
	 */
	private final Predicate<SimpleSVLocus> simpleSVLocusFileter;
	///////////////////////////////////////////
	private VCFFileReader VCFFileReader;

	private Map<SimpleSVType, Map<String,List<SimpleSVLocus>>> svTypeChromLocusListMapMap;
	
	
	public DellyJointGenotypeVCFReader(
			Path vcfFile, 
			VariantContextFilter variantContextFilter,
			Predicate<SimpleSVLocus> simpleSVLocusFileter) {
		super();
		this.vcfFile = vcfFile;
		this.variantContextFilter = variantContextFilter;
		this.simpleSVLocusFileter = simpleSVLocusFileter;
		
		this.preprocess();
		this.run();
		this.sortByGenomicPosition();
	}
	

	void preprocess() {
		this.VCFFileReader = new VCFFileReader(this.vcfFile, false);
		
		this.svTypeChromLocusListMapMap=new HashMap<>();
		
	}
	
	
	void run() {
		Stream<VariantContext> stream=VCFFileReader.iterator().stream();
		
		stream.filter(this.variantContextFilter).forEach(vc->{
			//check if passed filters (FILTER column being 'PASS')
			boolean passedFilter=vc.isNotFiltered();
			if(!passedFilter) return;
			
			String CHROM = vc.getContig();
			int POS = vc.getStart();
			String ID=vc.getID();
			Double QUAL = vc.getPhredScaledQual();
			String ALT=vc.getAlternateAllele(0).getDisplayString();//
			
			//##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="PE confidence interval around POS">
			@SuppressWarnings("unchecked")
			List<Integer> CIPOS=(ArrayList<Integer>)vc.getAttribute("CIPOS");
//			ConfidenceInterval cipos = new ConfidenceInterval(CIPOS==null?null:CIPOS.get(0), CIPOS==null?null:CIPOS.get(1));
			//##INFO=<ID=CIEND,Number=2,Type=Integer,Description="PE confidence interval around END">
			@SuppressWarnings("unchecked")
			List<Integer> CIEND=(ArrayList<Integer>)vc.getAttribute("CIEND");
			//##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for POS2 coordinate in case of an inter-chromosomal translocation">
			String CHR2=vc.getAttribute("CHR2")==null?null:(String)vc.getAttribute("CHR2");
			//##INFO=<ID=POS2,Number=1,Type=Integer,Description="Genomic position for CHR2 in case of an inter-chromosomal translocation">
			Integer POS2=vc.getAttribute("POS2")==null?null:Integer.parseInt((String)vc.getAttribute("POS2"));
			//##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
			Integer END=vc.getAttribute("END")==null?null:Integer.parseInt((String)vc.getAttribute("END"));
			//##INFO=<ID=PE,Number=1,Type=Integer,Description="Paired-end support of the structural variant">
			Integer PE=vc.getAttribute("PE")==null?null:Integer.parseInt((String)vc.getAttribute("PE"));
			//##INFO=<ID=MAPQ,Number=1,Type=Integer,Description="Median mapping quality of paired-ends">
			Integer MAPQ=vc.getAttribute("MAPQ")==null?null:Integer.parseInt((String)vc.getAttribute("MAPQ"));
			//##INFO=<ID=SRMAPQ,Number=1,Type=Integer,Description="Median mapping quality of split-reads">
			Integer SRMAPQ=vc.getAttribute("SRMAPQ")==null?null:Integer.parseInt((String)vc.getAttribute("SRMAPQ"));
			//##INFO=<ID=SR,Number=1,Type=Integer,Description="Split-read support">
			Integer SR=vc.getAttribute("SR")==null?null:Integer.parseInt((String)vc.getAttribute("SR"));
			//##INFO=<ID=SRQ,Number=1,Type=Float,Description="Split-read consensus alignment quality">
			Double SRQ=vc.getAttribute("SRQ")==null?null:Double.parseDouble((String)vc.getAttribute("SRQ"));
			//##INFO=<ID=CONSENSUS,Number=1,Type=String,Description="Split-read consensus sequence">
			String CONSENSUS=(String)vc.getAttribute("CONSENSUS");
			//##INFO=<ID=CE,Number=1,Type=Float,Description="Consensus sequence entropy">
			Double CE=vc.getAttribute("CE")==null?null:Double.parseDouble((String)vc.getAttribute("CE"));
			//##INFO=<ID=CT,Number=1,Type=String,Description="Paired-end signature induced connection type">
			String CT=(String)vc.getAttribute("CT");
			//##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Insertion length for SVTYPE=INS.">
			Integer SVLEN=vc.getAttribute("SVLEN")==null?null:Integer.parseInt((String)vc.getAttribute("SVLEN"));
			//##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
			Boolean IMPRECISE=vc.getAttribute("IMPRECISE")==null?null:(Boolean)vc.getAttribute("IMPRECISE");
			//##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variation">
			Boolean PRECISE=vc.getAttribute("PRECISE")==null?null:(Boolean)vc.getAttribute("PRECISE");
			//##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
			String SVTYPE=(String)vc.getAttribute("SVTYPE");
			//##INFO=<ID=INSLEN,Number=1,Type=Integer,Description="Predicted length of the insertion">
			Integer INSLEN=vc.getAttribute("INSLEN")==null?null:Integer.parseInt((String)vc.getAttribute("INSLEN"));
			//##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description="Predicted microhomology length using a max. edit distance of 2">
			Integer HOMLEN=vc.getAttribute("HOMLEN")==null?null:Integer.parseInt((String)vc.getAttribute("HOMLEN"));
			
			
			///////////////////////////////////////////
			
			SimpleSVType svType=null;
			Integer start=null;
			Integer end=null;
			String chrom=CHROM;
			if(SVTYPE.equals("DEL")) {
				svType=SimpleSVType.DEL;
				start=POS;
				end=END;
				
			}else if(SVTYPE.equals("INS")) {
				//====================important notes on delly's INS type SV record in the output VCF file=====================
				//delly's VCF file has INS type SV's END attribute equal to the POS + inserted_seq_lenght, which is incorrect
				//for example, for a INS sv with CHROM=Chr02 and POS=64285230 and inserted seq length=29, its END attribute = 64285258
				//note that, normally for insertion type SV, the two breakends on the reference genome should have independent evidence; however, for those insertion SV of size << read length, the two breakends can be detected by the same set of reads, thus the evidence of them will be the same
				//=============================================================================================================
				svType=SimpleSVType.INS;
				start=POS;
				end=POS; //set the end the same as the start, length of INS is not accurate
				
			}else if(SVTYPE.equals("INV")) {
				svType=SimpleSVType.INV;
				start=POS;
				end=END;
				
			}else if(SVTYPE.equals("DUP")) {
				//DELLY's DUP sv are in a single line with only END attribute
				//thus, by reasoning, delly's DUP must be tandem duplication with the duplicated region between the POS and END
				svType=SimpleSVType.DUP;
				start=POS;
				end=END; //
				
			}else {
//				System.out.println("skip sv type:"+SVTYPE);
				return; 
			}
			
			///////////////////build the genotype map for all samples
			Map<Integer, Genotype> sampleIndexGenotypeMap = new HashMap<>();
			
			for(int i=0;i<vc.getNSamples();i++) {
				if(vc.getGenotype(i).isNoCall()) {
					sampleIndexGenotypeMap.put(i+1, Genotype.MISSING);
				}else {
					if(vc.getGenotype(i).isHomRef()) {
						sampleIndexGenotypeMap.put(i+1, Genotype.ABSENCE);
					}else if(vc.getGenotype(i).isHet()) {
						sampleIndexGenotypeMap.put(i+1, Genotype.HETER);
					}else if(vc.getGenotype(i).isHomVar()) {
						sampleIndexGenotypeMap.put(i+1, Genotype.PRESENCE);
					}else {
						htsjdk.variant.variantcontext.Genotype gt=vc.getGenotype(i);
						throw new UnsupportedOperationException("check genotype!");
					}
				}
			}
			
			
			//////////////////////////
			SimpleSVLocus svlocus = new SimpleSVLocus(
					svType, chrom, start, end, sampleIndexGenotypeMap
					);
			
			if(this.simpleSVLocusFileter.test(svlocus)) {
				
				if(!this.svTypeChromLocusListMapMap.containsKey(svType)) {
					this.svTypeChromLocusListMapMap.put(svType, new HashMap<>());
				}
				
				if(!this.svTypeChromLocusListMapMap.get(svType).containsKey(chrom)) {
					this.svTypeChromLocusListMapMap.get(svType).put(chrom, new ArrayList<>());
				}
				
				this.svTypeChromLocusListMapMap.get(svType).get(chrom).add(svlocus);
				
//				svlocus.setVariantContext(vc);
			}
		});
	}


	void sortByGenomicPosition() {
		for(SimpleSVType type:this.svTypeChromLocusListMapMap.keySet()) {
			
			for(String chrom:this.svTypeChromLocusListMapMap.get(type).keySet()) {
				Collections.sort(this.svTypeChromLocusListMapMap.get(type).get(chrom), 
						(a,b)->{
							return a.getStart()-b.getStart();
						});
			}
			
		}
	}
	
	
	/**
	 * @return the svTypeChromLocusListMapMap
	 */
	public Map<SimpleSVType, Map<String, List<SimpleSVLocus>>> getSvTypeChromLocusListMapMap() {
		return svTypeChromLocusListMapMap;
	}
	
	
	public List<SimpleSVLocus> getSimpleSVLocusList(){
		List<SimpleSVLocus> ret = new ArrayList<>();
		
		for(SimpleSVType type:this.getSvTypeChromLocusListMapMap().keySet()) {
			for(String chrom:this.getSvTypeChromLocusListMapMap().get(type).keySet()) {
				ret.addAll(this.getSvTypeChromLocusListMapMap().get(type).get(chrom));
			}
		}
		
		return ret;
	}
	
	/**
	 * 
	 */
	public Map<SimpleSVType, Integer> getSvTypeTotalLocusNumMap(){
		Map<SimpleSVType, Integer>  ret = new HashMap<>();
		
		for(SimpleSVType type:this.svTypeChromLocusListMapMap.keySet()) {
			ret.put(type, 0);
			for(String chrom:this.svTypeChromLocusListMapMap.get(type).keySet()) {
				ret.put(type, ret.get(type)+this.svTypeChromLocusListMapMap.get(type).get(chrom).size());
			}
		}
		
		return ret;
	}
	
	/**
	 * @return the vCFFileReader
	 */
	public VCFFileReader getVCFFileReader() {
		return VCFFileReader;
	}
	
	
}
