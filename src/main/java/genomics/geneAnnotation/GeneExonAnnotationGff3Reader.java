package genomics.geneAnnotation;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import basic.gff3.Gff3Hit;
import basic.gff3.Gff3Utils;
import genomics.utils.CDS;
import genomics.utils.Exon;
import genomics.utils.FivePrimeUTR;
import genomics.utils.Gene;
import genomics.utils.GenomicRegionUtils;
import genomics.utils.Intron;
import genomics.utils.ThreePrimeUTR;
import genomics.utils.Transcript;
import genomics.utils.UTR;

/**
 * reader class for gene annotations from a gff3 file;
 * 
 * @author tanxu
 *
 */
public class GeneExonAnnotationGff3Reader {
	private final Path gff3File;
	
	/**
	 * the attribute in the attribute column that contains the name of the gene for a data line of 'gene' type
	 * for example, for rice, the attribute name is 'Name' as shown below
	 * 		Chr1	phytozomev11	gene	2903	10817	.	+	.	ID=LOC_Os01g01010.MSUv7.0;Name=LOC_Os01g01010
	 * 
	 * note that different gff3 file of different genome annotation may have different attribute name! 
	 */
	private final String geneNameAttributeName;
	
	/**
	 * the attribute in the attribute column that contains the name of the transcript for a data line of 'mRNA' type
	 * for example, for rice, the attribute name is 'Name' as shown below
	 * 		Chr1	phytozomev11	mRNA	2903	10817	.	+	.	ID=LOC_Os01g01010.1.MSUv7.0;Name=LOC_Os01g01010.1;pacid=33123311;longest=1;Parent=LOC_Os01g01010.MSUv7.0
	 * 
	 * note that different gff3 file of different genome annotation may have different attribute name! 
	 */
	private final String transcriptIDAttributeName;
	
	///////////////
	private List<Gene> genes;
	
	private Map<String, List<Gene>> chromNameSortedGeneListMap;
	
	private List<UTR> utrs;
	private List<UTR> utrs3;
	private List<UTR> utrs5;
	
	private List<Intron> introns;

	private List<CDS> CDSs;
	
	private List<Exon> exons;
	/**
	 * 
	 * @param gff3File
	 * @param geneNameAttributeName
	 * @param transcriptIDAttributeName
	 */
	public GeneExonAnnotationGff3Reader(Path gff3File, String geneNameAttributeName, String transcriptIDAttributeName) {
		super();
		this.gff3File = gff3File;
		this.geneNameAttributeName = geneNameAttributeName;
		this.transcriptIDAttributeName = transcriptIDAttributeName;
		
		
		this.read();
	}

	/**
	 * 
	 */
	private void read() {
		System.out.println("Start reading the gff3 file...");
		this.genes = new ArrayList<>();
		
		try {
			BufferedReader lineReader = new BufferedReader(new FileReader(gff3File.toFile()));
			String line = null;
			
			///////////////
			Gene currentGene=null;
			Transcript currentTranscript=null;
			while ((line = lineReader.readLine()) != null) {
				if(line.startsWith("#") ||line.isEmpty()) {
					continue;
				}
				///
				Gff3Hit hit = Gff3Hit.build(line);
				
				
				///////////////
				if(hit.getType().equals(Gff3Utils.GENE_TYPE)) {
					currentGene = new Gene(
							hit.getAttributeNameStringValueMap().get(this.geneNameAttributeName),
							hit.getChrom(), hit.getStart(),hit.getEnd(),hit.getStrand()
							);
					this.genes.add(currentGene);
//					System.out.println("gene:"+currentGene.getName());
				}else if(hit.getType().equals(Gff3Utils.TRANSCRIPT_TYPE)) {
					currentTranscript = new Transcript(
							hit.getChrom(), hit.getStart(),hit.getEnd(), hit.getStrand(),
							hit.getAttributeNameStringValueMap().get(this.transcriptIDAttributeName)
							);
					currentGene.addTranscript(currentTranscript);
//					System.out.println("transcript:"+currentTranscript.getId());
				}else if(hit.getType().equals(Gff3Utils.EXON_TYPE)) {
					Exon exon=new Exon(hit.getChrom(), hit.getStart(), hit.getEnd(), hit.getStrand());
					currentTranscript.addExon(exon);
					
				}else if(hit.getType().equals(Gff3Utils.CDS_TYPE)) {
					CDS cds=new CDS(hit.getChrom(), hit.getStart(), hit.getEnd(), hit.getStrand(), hit.getPhase());
					
					currentTranscript.addCDS(cds);
					
				}else if(hit.getType().equals(Gff3Utils.FIVE_PRIME_UTR_TYPE)) {
					FivePrimeUTR fivePrimeUTR=new FivePrimeUTR(hit.getChrom(), hit.getStart(), hit.getEnd(), hit.getStrand());
					currentTranscript.addFivePrimeUTR(fivePrimeUTR);
				}else if(hit.getType().equals(Gff3Utils.THREE_PRIME_UTR_TYPE)) {
					ThreePrimeUTR threePrimeUTR=new ThreePrimeUTR(hit.getChrom(), hit.getStart(), hit.getEnd(), hit.getStrand());
					currentTranscript.addThreePrimeUTR(threePrimeUTR);
				}
			}
			
			lineReader.close();
		} catch (IOException ex) {
			System.err.println(ex);
		}
		
		System.out.println("Reading the gff3 file is done...");
	}
	
	private void buildChromSortedGeneMap() {
		this.chromNameSortedGeneListMap=new HashMap<>();
		
		for(Gene gene:this.genes) {
			if(!this.chromNameSortedGeneListMap.containsKey(gene.getChrom())) {
				this.chromNameSortedGeneListMap.put(gene.getChrom(), new ArrayList<>());
			}
			this.chromNameSortedGeneListMap.get(gene.getChrom()).add(gene);
		}
		
		for(String chrom:this.chromNameSortedGeneListMap.keySet()) {
			Collections.sort(this.chromNameSortedGeneListMap.get(chrom), GenomicRegionUtils.sorterByChromAndStartPos());
		}
	}
	
	
	/**
	 * lookup and return the list of {@link Gene}s that is overlapping with the given genomic region
	 * @param chrom
	 * @param start
	 * @param end
	 * @param contained if true, only return genes that are fully contained in the given region; otherwise, return all genes that are overlapping with the given region
	 */
	public List<Gene> queryGene(String chrom, int start, int end, boolean contained) {
		if(this.chromNameSortedGeneListMap==null) {
			this.buildChromSortedGeneMap();
		}
		
		List<Gene> ret = new ArrayList<>();
		if(!this.chromNameSortedGeneListMap.containsKey(chrom)) {
			return ret;
		}
		
		////////////////////
		for(Gene gene:this.chromNameSortedGeneListMap.get(chrom)) {
			if(gene.getEnd()<start) {
				continue;
			}else if(gene.getStart()>end) {
				break;
			}else {
				if(contained) {
					if(gene.getStart()>=start && gene.getEnd()<=end) {
						ret.add(gene);
					}
				}else {
					ret.add(gene);
				}
			}
		}
		
		return ret;
	}
	
	/**
	 * @return the genes
	 */
	public List<Gene> getGenes() {
		return genes;
	}
	
	
	/**
	 * 
	 * @return
	 */
	public List<Transcript> getLongestTranscripts(){
		List<Transcript> longestTranscripts = new ArrayList<>();
		this.getGenes().forEach(g->{
			longestTranscripts.add(g.getLongestTranscript());
		});
		return longestTranscripts;
	}
	/**
	 * return the {@link UTR}s of longest transcripts of all genes
	 * @return
	 */
	public List<UTR> getUTRs(){
		if(this.utrs==null) {
			this.utrs=new ArrayList<>();
			
			for(Gene gene:this.genes) {
				this.utrs.addAll(gene.getLongestTranscript().getFivePrimeUTRs());
				this.utrs.addAll(gene.getLongestTranscript().getThreePrimeUTRs());
			}
		}
		
		return this.utrs;
	}
	
	/**
	 * return the {@link UTR}s of longest transcripts of all genes
	 * @return
	 */
	public List<UTR> get3PrimeUTRs(){
		if(this.utrs3==null) {
			this.utrs3=new ArrayList<>();
			
			for(Gene gene:this.genes) {
				this.utrs3.addAll(gene.getLongestTranscript().getThreePrimeUTRs());
			}
		}
		
		return this.utrs3;
	}
	
	/**
	 * return the {@link UTR}s of longest transcripts of all genes
	 * @return
	 */
	public List<UTR> get5PrimeUTRs(){
		if(this.utrs5==null) {
			this.utrs5=new ArrayList<>();
			
			for(Gene gene:this.genes) {
				this.utrs5.addAll(gene.getLongestTranscript().getFivePrimeUTRs());
			}
		}
		
		return this.utrs5;
	}
	
	
	/**
	 * return the {@link UTR}s of longest transcripts of all genes
	 * @return
	 */
	public List<Intron> getIntrons(){
		if(this.introns==null) {
			this.introns=new ArrayList<>();
			
			for(Gene gene:this.genes) {
				this.introns.addAll(gene.getLongestTranscript().getIntrons());
			}
		}
		
		return this.introns;
	}
	
	/**
	 * return the list of CDS of longest transcripts of all genes
	 * @return
	 */
	public List<CDS> getCDSs(){
		if(this.CDSs==null) {
			this.CDSs=new ArrayList<>();
			
			for(Gene gene:this.genes) {
				this.CDSs.addAll(gene.getLongestTranscript().getCDSList());
			}
		}
		
		return this.CDSs;
	}
	
	/**
	 * return the list of Exons of longest transcripts of all genes
	 * @return
	 */
	public List<Exon> getExons(){
		if(this.exons==null) {
			this.exons=new ArrayList<>();
			
			for(Gene gene:this.genes) {
				this.exons.addAll(gene.getLongestTranscript().getExons());
			}
		}
		
		return this.exons;
	}
}
