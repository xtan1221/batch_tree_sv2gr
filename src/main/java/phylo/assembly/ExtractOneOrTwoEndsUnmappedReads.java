package phylo.assembly;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SAMValidationError;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import phylo.ref.Region;

/**
 * extract all primary records of reads with both ends unmapped or only one end uniquely mapped while the other unmapped
 * 
 * for reads with both ends unmapped, this is trivial;
 * 
 * for reads with one end mapped, other end unmapped
 * 		if the mapped end is uniquely mapped (<==> mapq score > 0) and not marked as duplicate
 * 			retrieve the single representative alignment (trivial if not split aligned)
 * 		else //multiply mapped; the mapq == 0
 * 			do not include!
 * 
 * ======================
 * note that for reads with one end mapped, the other not mapped, if it is marked as duplicate, only the mapped end will be flagged 1024, while the unmapped end will not be flagged 1024
 * need to deal with above issue to avoid unexpected result
 * for both ends unmapped reads, they will never be marked as duplicate (validated with real data processed by picard.MarkDuplicates)
 * 
 * @author tanxu
 *
 */
public class ExtractOneOrTwoEndsUnmappedReads {
	/**
	 * the input bam file does not need to be sorted and indexed
	 */
	private final Path inputBamFile;
	
	
	private final Path outputBamFile;
	
	/**
	 * the minimal mapq for the mapped end to be included!
	 * must be positive integer;
	 * cannot be 0!
	 */
	private final int minMapq;
	
	/////////////////////////
	private final Log log = Log.getInstance(ExtractOneOrTwoEndsUnmappedReads.class);
	
	/**
	 * 
	 * @param inputBamFile
	 * @param regionList
	 * @param outputBamFile
	 * @param minMapq
	 */
	ExtractOneOrTwoEndsUnmappedReads(Path inputBamFile, Path outputBamFile, int minMapq){
		//TODO
		if(minMapq<=0)
			throw new IllegalArgumentException("given minMapq must be positive integer (so that read with multiple alignments can be filtered out)!");
		
		
		this.inputBamFile = inputBamFile;
		this.outputBamFile=outputBamFile;
		this.minMapq = minMapq;
		
	}
	
	/**
	 * strategy:
	 * 
	 * for a record r, 
	 * 		if it is unmapped, add it to the output bam
	 * 		else
	 * 			if its mate is unmapped
	 * 				if r is not marked as duplicate and uniquely mapped (mapq>0) and representative alignment (no flag 2048)
	 * 					add r to the output bam
	 * 
	 * @throws IOException
	 */
	public void run() {
		SamReader reader = SamReaderFactory.makeDefault().open(inputBamFile.toFile());
		
		SAMFileWriterFactory SAMFileWriterFactory = new SAMFileWriterFactory();
//		SAMFileWriterFactory.setTempDirectory(tmpDir); 
//		SAMFileWriterFactory.setBufferSize(bufferSize)
		
		SAMFileWriter writer = SAMFileWriterFactory.makeBAMWriter(reader.getFileHeader(), false, outputBamFile);
		
		final ProgressLogger progress = new ProgressLogger(log);
		
		Map<String,SAMRecord> readNameUnmappedMateRecordMap = new HashMap<>();
		Set<String> oneEndMappedReadNameSetWithMappedEndProcessed=new HashSet<>(); //the set of read names with one end mapped (meet all the requirements) and the mapped end is processed (either write to output bam file or skipped if marked as duplicate)
		Set<String> oneEndMappedReadNameSetWithMappedEndOutToFile=new HashSet<>();
		//note that for reads with one end mapped, the other not mapped, if it is marked as duplicate, only the mapped end will be flagged 1024, while the unmapped end will not be flagged 1024
		//need to deal with above issue to avoid unexpected result
		//for both ends unmapped reads, they will never be marked as duplicate (validated with real data)
		reader.iterator().forEachRemaining(currentRecord->{
			if(currentRecord.getReadUnmappedFlag()) {//end unmapped
				if(currentRecord.getMateUnmappedFlag()){//mate also unmapped, skip
					writer.addAlignment(currentRecord);
				}else{//mate is mapped
					if(oneEndMappedReadNameSetWithMappedEndProcessed.contains(currentRecord.getReadName())){//mapped mate has already been reached
						if(oneEndMappedReadNameSetWithMappedEndOutToFile.contains(currentRecord.getReadName())) {//mapped mate has been written to output file
							//mapped mate is already written out to file (thus meet all requirements)
							writer.addAlignment(currentRecord);
						}
					}else{//mapped mate is not reached yet, not sure if it meets all requirements, store currentRecord in readNameUnmappedMateRecordMap
						readNameUnmappedMateRecordMap.put(currentRecord.getReadName(), currentRecord);
					}
				}
			}else {//current record is mapped, 
				if(currentRecord.getMateUnmappedFlag()) {//mate is unmapped
					//only process for the primary and representative record of the mapped end
					if(SAMRecordUtils.isPrimaryAndRepresentative(currentRecord)) {
						if(currentRecord.getDuplicateReadFlag()||currentRecord.getMappingQuality()<this.minMapq) {
							//do not write to file
						}else {
							//write to file
							writer.addAlignment(currentRecord);
							
							//
							oneEndMappedReadNameSetWithMappedEndOutToFile.add(currentRecord.getReadName());
							
							if(readNameUnmappedMateRecordMap.containsKey(currentRecord.getReadName())) {
								writer.addAlignment(readNameUnmappedMateRecordMap.get(currentRecord.getReadName()));
							}
						}
						
						//add the read name to the processedReadNameSetWithOneEndMapped
						oneEndMappedReadNameSetWithMappedEndProcessed.add(currentRecord.getReadName());
						
						//always check if unmapped mate record is in readNameUnmappedMateRecordMap or not if current mapped record is primary and not supplementary; 
						//if yes, remove it
						if(readNameUnmappedMateRecordMap.containsKey(currentRecord.getReadName())) {
							readNameUnmappedMateRecordMap.remove(currentRecord.getReadName());
						}
						
					}else {
						//do nothing
					}
					
				}
			}
			//
			progress.record(currentRecord);
		});
		
		//there should be no unmapped mate left in the readNameUnmappedMateRecordMap
		if(!readNameUnmappedMateRecordMap.isEmpty())
			SAMUtils.processValidationError(new SAMValidationError(SAMValidationError.Type.MATE_NOT_FOUND,
                    "Found " + readNameUnmappedMateRecordMap.size() + " unpaired mates", null), ValidationStringency.DEFAULT_STRINGENCY);
		
		try {
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		writer.close();
		
		log.info("record extraction from bam is done!");
	}
	
	
	public static void main(String[] args) {
		Path inputSortedBamFile = Path.of("/scratch/tanxu/reseq/sb/single_sample_hc_pipeline_all_32/result/bam_files/SRR486615/SRR486615.trimmed.coord.sorted.dup.marked.rg.added.bam");
		List<Region> regionList = new ArrayList<>(); 
		regionList.add(new Region("Chr01", 619500, 619530));
		Path outputBamFile = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv2gr_java\\output\\out.bam");
		int minMapq = 10;
		
		RegionalReadExtractor extractor = new RegionalReadExtractor(inputSortedBamFile, regionList, outputBamFile, minMapq);
		
		extractor.run();
	}
	
	
	
	
	
	
}
