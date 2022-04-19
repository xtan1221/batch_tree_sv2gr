package phylo.assembly;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import phylo.ref.Region;


/**
 * extract reads in a set of target region R so that at least one segment of the two ends of the read pair is uniquely mapped to R
 * 
 * then output the read names in a list file
 * 
 * @author tanxu
 *
 */
public class RegionalReadExtractor2 {
	private final Path inputSortedAndIndexedBamFile;
	
	private final List<Region> regionList;
	
	private final Path outputReadNameListFile;
	
	private final int minMapq;
	
	/////////////
	private Set<String> readNameSet;
	/**
	 * 
	 * @param inputSortedAndIndexedBamFile
	 * @param regionList
	 * @param outputBamFile
	 * @param minMapq
	 */
	RegionalReadExtractor2(Path inputSortedAndIndexedBamFile, List<Region> regionList, Path outputReadNameListFile, int minMapq){
		//TODO
		this.inputSortedAndIndexedBamFile=inputSortedAndIndexedBamFile;
		this.regionList=regionList;
		this.outputReadNameListFile=outputReadNameListFile;
		this.minMapq = minMapq;
		
	}
	
	/**
	 * basic strategy to ensure each qualified record is included in output bam file only once
	 * 
	 * 1. if a record passes all filters, add it to output bam file
	 * 		if its mate also passes all filter (including the region filter)
	 * 			skip it because it will be included by itself
	 * 		else //the mate does not pass the filters
	 * 			add the mate record to output bam file as well
	 * 		
	 * 2. else, skip it;
	 * 		it will be taken care of by its mate if it should be included!
	 * @throws IOException
	 */
	public void run() {
		this.readNameSet = new HashSet<>();
		SamReader reader = SamReaderFactory.makeDefault().open(inputSortedAndIndexedBamFile.toFile());
		if(!reader.hasIndex())
			throw new UnsupportedOperationException("input bam file is not indexed!");
		
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
			System.out.println(r.getSAMString());
			if(recordPassingAllNonRegionalFilter(r)) {
				this.readNameSet.add(r.getReadName());
			}
		});
		
		try {
			reader.close();
			
		    BufferedWriter writer = new BufferedWriter(new FileWriter(this.outputReadNameListFile.toString()));

		    for(String rn:this.readNameSet) {
		    	writer.write(rn);
		    	writer.newLine();
		    }
			
			writer.flush();
			
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
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
		regionList.add(new Region("Chr01", 58, 100));
		Path outputBamFile = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv2gr_java\\output\\read.list.txt");
		int minMapq = 1;
		
		RegionalReadExtractor2 extractor = new RegionalReadExtractor2(inputSortedBamFile, regionList, outputBamFile, minMapq);
		
		extractor.run();
	}
}
