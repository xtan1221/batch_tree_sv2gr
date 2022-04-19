package phylo.assembly;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import phylo.ref.Region;


/**
 * extract reads in a set of target region R so that at least one segment of the two ends of the read pair is uniquely mapped to R with mapq score >= threshold
 * 
 * then output all records to a bam file
 * 
 * ===================
 * this is to extract read sequences for regional de novo assembly, thus the full sequence of extracted reads should be present in the extracted bam file;
 * also, sequence of each included read should only be present once!
 * 
 * note that for records with soft clipping in CIGAR string, the full read sequence will be put in the column 10 of the bam file
 * for records that is flagged as supplementary alignment (with flag 2048) will have hard clipped CIGAR string and only the matched segment will be put in the bam file
 * 
 * thus, for a read end that has split alignments (chimeric alignment), only need to extract the one without flag 2048 (representative alignment)
 * 
 * the strategy to extract the representative alignment of a supplementary alignment r is as following:
 * 1. query the mate of r, m with {@link SAMRecordQueryUtils#queryNonSupplementaryMate(SamReader, SAMRecord)} rather than {@link SamReader#queryMate(SAMRecord)}
 * 2. query the mate of m, r', {@link SAMRecordQueryUtils#queryNonSupplementaryMate(SamReader, SAMRecord)} rather than {@link SamReader#queryMate(SAMRecord)}
 * 		r' will be the representative alignment of r!
 * 
 * ==============================
 * note that for read with both primary alignment and secondary alignments, the read should not be included;
 * to exclude the primary alignment of such read, simply filter out all alignments with mapq = 0, since the primary alignment of such reads will always be 0!
 * 
 * @author tanxu
 * 
 */
public class RegionalReadExtractor {
	private final Path inputSortedAndIndexedBamFile;
	
	private final List<Region> regionList;
	
	private final Path outputBamFile;
	
	/**
	 * must be positive integer;
	 * cannot be 0!
	 */
	private final int minMapq;
	
	/**
	 * 
	 * @param inputSortedAndIndexedBamFile
	 * @param regionList
	 * @param outputBamFile
	 * @param minMapq
	 */
	RegionalReadExtractor(Path inputSortedAndIndexedBamFile, List<Region> regionList, Path outputBamFile, int minMapq){
		//TODO
		if(minMapq<=0)
			throw new IllegalArgumentException("given minMapq must be positive integer (so that read with multiple alignments can be filtered out)!");
		
		
		this.inputSortedAndIndexedBamFile = inputSortedAndIndexedBamFile;
		this.regionList =regionList;
		this.outputBamFile=outputBamFile;
		this.minMapq = minMapq;
		
		
	}
	
	/**
	 * TODO need further testing and debug
	 * 
	 * basic strategy to ensure each qualified record is included in output bam file only once
	 * 
	 * 1. if a record r passes all filters, need to find its representative alignment because the supplementary alignment record's column 10 (SEQ) only contain the alignment segment (hard clipping)
	 * 		if r is a supplementary alignment record, need to 
	 * 			1. query the mate record m of r with {@link SAMRecordQueryUtils#queryNonSupplementaryMate(SamReader, SAMRecord)}
	 * 			2. then query the mate r' of m with {@link SAMRecordQueryUtils#queryNonSupplementaryMate(SamReader, SAMRecord)}
	 * 				r' is the representative alignment record of r
	 * 			3. 
	 * 			if r' passes all filter (including the region filter)
	 * 				skip it because it will be included by itself
	 * 			else //the r' does not pass the filters
	 * 				add the mate record to output bam file as well
	 * 		else
	 * 			add it to output bam file
	 * 		
	 * 		============
	 * 		process the mate end of the read of r
	 * 		1. query the mate record m of r with {@link SAMRecordQueryUtils#queryNonSupplementaryMate(SamReader, SAMRecord)}
	 * 		2. 
	 * 		if m also passes all filter (including the region filter)
	 * 			skip it because it will be included by itself
	 * 		else //the mate does not pass the filters
	 * 			add m to output bam file as well
	 * 		
	 * 2. else, skip it;
	 * 		it will be taken care of by its mate if it should be included!
	 * @throws IOException
	 */
	public void run() {
		SamReader reader = SamReaderFactory.makeDefault().open(inputSortedAndIndexedBamFile.toFile());
		if(!reader.hasIndex())
			throw new UnsupportedOperationException("input bam file is not indexed!");
		
		SamReader reader2 = SamReaderFactory.makeDefault().open(inputSortedAndIndexedBamFile.toFile());//for mate record querying
		
//		reader.getFileHeader().getSequenceDictionary().getSequences().forEach(seq->{
//			System.out.println(seq.getSequenceIndex()+" "+seq.getSequenceName());
//		});
		
		SAMFileWriterFactory SAMFileWriterFactory = new SAMFileWriterFactory();
//		SAMFileWriterFactory.setTempDirectory(tmpDir); 
//		SAMFileWriterFactory.setBufferSize(bufferSize)
		
		SAMFileWriter writer = SAMFileWriterFactory.makeBAMWriter(reader.getFileHeader(), false, outputBamFile);
		
		QueryInterval[] queryIntervals = new QueryInterval[this.regionList.size()];
		for(int i=0;i<this.regionList.size();i++){
			Region r = this.regionList.get(i);
			queryIntervals[i]=new QueryInterval(
					reader.getFileHeader().getSequenceIndex(r.getReferenceName()),
					r.getStart(),
					r.getEnd()
					);
		}
		
		//contained 
		//If true, each SAMRecord returned is will have its alignment completely contained in one of theintervals of interest. If false, the alignment of the returned SAMRecords need only overlap one ofthe intervals of interest.
		SAMRecordIterator recordIterator = reader.query(queryIntervals, false);
		
		recordIterator.forEachRemaining(r->{
//			System.out.println(r.getSAMString());
			SAMRecord mate = SAMRecordQueryUtils.queryNonSupplementaryMate(reader2, r);
			
			if(recordPassingAllNonRegionalFilter(r)) {
				//first process the end of the read of r
				if(r.getSupplementaryAlignmentFlag()) {//the record is supplementary, need to find its representative alignment because the supplementary alignment record's column 10 (SEQ) only contain the alignment segment (hard clipping)
					if(mate==null) {
						System.out.println("warning: mate record is not found!");
					}else {
						SAMRecord representativeRec = SAMRecordQueryUtils.queryNonSupplementaryMate(reader2, mate);
						if(representativeRec==null) {
							System.out.println("warning: representativeRec record is not found!");
						}else {
							if(this.recordPassingAllNonRegionalFilter(representativeRec) && this.recordInTargetRegion(representativeRec)) {
								//mate record passes all filters, it will be included by itself
							}else {//add the mate record by current record
								if(representativeRec.getReadString().equals("*"))
									System.out.println("Error: extracted alignment record has '*' sequence!");
								writer.addAlignment(representativeRec);
							}
						}
					}
				}else {//the record is representative, simply include it
					if(r.getReadString().equals("*"))
						System.out.println("Error: extracted alignment record has '*' sequence!");
					writer.addAlignment(r);
				}
				
				//then process the mate end of the read of r
				if(mate==null) {
					System.out.println("warning: mate record is not found!");
				}else {
					if(this.recordPassingAllNonRegionalFilter(mate) && this.recordInTargetRegion(mate)) {
						//mate record passes all filters, it will be included by itself
					}else {//add the mate record by current record
						if(mate.getReadString().equals("*"))
							System.out.println("Error: extracted alignment record has '*' sequence!");
						
						writer.addAlignment(mate);
					}
				}
			}
			
		});
		
		try {
			reader.close();
			reader2.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		writer.close();
		
		System.out.println("extraction is done!");
	}
	
	
	/**
	 * check if the given record passes all filters so that it can be included in the output BAM file
	 * note that whether the record is mapped onto the target regions or not is not checked by this method
	 * @param record
	 * @return
	 */
	private boolean recordPassingAllNonRegionalFilter(SAMRecord record) {
		return !record.getReadUnmappedFlag()&& 
				!record.isSecondaryAlignment()&& 
				record.getMappingQuality()>=minMapq&&
				!record.getDuplicateReadFlag();
	}
	
	/**
	 * 
	 * @param record
	 * @return
	 */
	private boolean recordInTargetRegion(SAMRecord record) {
		if(record.getReadUnmappedFlag()) 
			//this is needed because for unmapped read, it may still be assigned reference and position if its mate is mapped; for example:
			//SRR486615.97563620	69	Chr01	22	0	*	=	22	0	... //this is unmapped but has reference and position and is put together with its mapped mate in the sorted bam file!
			//SRR486615.97563620	137	Chr01	22	13	4S52M37S	=	22	0	...
			throw new IllegalArgumentException("given record is not mapped!");
		
		for(Region region:this.regionList) {
			if(region.getReferenceName().equals(record.getReferenceName())) {//on the same reference
				if(region.getStart()>record.getEnd() || region.getEnd()<record.getStart()) {
					//record and the region is non-overlapping; 
				}else {//
					return true;
				}
			}
		}
		
		return false;
	}
	
	/**
	 * return true if record1 is sorted at the upstream of record 2
	 * in detail, return true if the record 1's alignment start is at the upstream of record 2's
	 * 
	 * @param record1
	 * @param record2
	 * @return
	 */
	private boolean record1IsStrictlyMappedBeforeRecord2(SAMRecord record1, SAMRecord record2) {
		if(record1.getReadUnmappedFlag()||record2.getReadUnmappedFlag())
			throw new IllegalArgumentException("given record 1 and/or record 2 is not mapped!");
		
		return record1.getReferenceIndex()<record2.getReferenceIndex() || 
				record1.getReferenceIndex()==record2.getReferenceIndex() && record1.getAlignmentStart()<record2.getAlignmentStart();
	}
	
	
	public static void main(String[] args) {
		Path inputSortedBamFile = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\data processing and analysis\\SV detection\\result\\analysis\\SRR486615.trimmed.coord.sorted.dup.marked.rg.added.chr01.1-1000000.bam");
		List<Region> regionList = new ArrayList<>(); 
		regionList.add(new Region("Chr01", 619500, 619530));
		Path outputBamFile = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv2gr_java\\output\\out.bam");
		int minMapq = 10;
		
		RegionalReadExtractor extractor = new RegionalReadExtractor(inputSortedBamFile, regionList, outputBamFile, minMapq);
		
		extractor.run();
	}
}
