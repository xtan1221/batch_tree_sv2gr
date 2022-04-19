package genomics.utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * a transcript on genome of a {@link Gene}, including all UTRs, introns and CDSs
 * 
 * note that exon is composed of UTR and CDS
 * 
 * @author tanxu
 * 
 */
public class Transcript extends SimpleGenomicRegion{
	/**
	 * transcript id
	 */
	private final String id;
	
	private final List<Exon> exonList;
	
	private final List<CDS> CDSList;
	
	/**
	 * there may be one or multiple {@link ThreePrimeUTR}s for a transcript
	 * 		for example, 3' UTR divided by introns
	 */
	private final List<ThreePrimeUTR> threePrimeUTRs;
	/**
	 * there may be one or multiple {@link FivePrimeUTR}s for a transcript
	 * 		for example, 5' UTR divided by introns
	 */
	private final List<FivePrimeUTR> fivePrimeUTRs;
	
	/////////////////////////////
	private List<Intron> introns;
	
	public Transcript(String chrom, int start, int end, Strand strand, String id) {
		super(chrom, start, end, strand);
		
		if(id==null||id.isEmpty()) {
			throw new IllegalArgumentException("given transcript id cannot be null or empty!");
		}
		
		this.id = id;
		this.exonList=new ArrayList<>();
		this.CDSList = new ArrayList<>();
		this.threePrimeUTRs=new ArrayList<>();
		this.fivePrimeUTRs=new ArrayList<>();
	}
	


	public void addExon(Exon exon) {
		this.exonList.add(exon);
	}
	
	
	
	/**
	 * @return the exonList
	 */
	public List<Exon> getExons() {
		return exonList;
	}



	/**
	 * @return the cDSList
	 */
	public List<CDS> getCDSList() {
		return CDSList;
	}



	/**
	 * @return the id
	 */
	public String getId() {
		return id;
	}



	public void addCDS(CDS cds) {
		this.CDSList.add(cds);
	}
	
	public void addThreePrimeUTR(ThreePrimeUTR utr) {
		this.threePrimeUTRs.add(utr);
	}
	
	public void addFivePrimeUTR(FivePrimeUTR utr) {
		this.fivePrimeUTRs.add(utr);
	}

	/**
	 * @return the threePrimeUTRs
	 */
	public List<ThreePrimeUTR> getThreePrimeUTRs() {
		return threePrimeUTRs;
	}



	/**
	 * @return the fivePrimeUTRs
	 */
	public List<FivePrimeUTR> getFivePrimeUTRs() {
		return fivePrimeUTRs;
	}
	
	/**
	 * 
	 * @return
	 */
	public List<Intron> getIntrons(){
		if(this.introns==null) {
			this.introns=new ArrayList<>();
			Collections.sort(this.exonList, GenomicRegionUtils.sorterByChromAndStartPos());
			Exon previous=null;
			for(Exon exon:this.exonList) {
				if(previous==null) {
					//
				}else {
					Intron intron=new Intron(this.getChrom(), previous.getEnd()+1, exon.getStart()-1, this.getStrand());
					this.introns.add(intron);
				}
				previous=exon;
			}
		}
		
		return this.introns;
	}
	
}
