package phylo.assembly;

import htsjdk.samtools.SAMFormatException;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloseableIterator;

/**
 * 
 * @author tanxu
 *
 */
public class SAMRecordQueryUtils {
	
	/**
	 * this method is the same with {@link SamReader#queryMate(SAMRecord)} except for the following
	 * 
	 * if there are multiple records of mate read mapped at the same coordinate (with one being the representative other being supplementary with flag 2048), (which should be very rare but not impossible)
	 * 		this method will select the representative one (without flag 2048)
	 * in the original {@link SamReader#queryMate(SAMRecord)}, SAMFormatException will be thrown for such cases;
	 * 
	 * for example, if query the mate of the third record, since both representative and supplementary alignment are mapped at the same position (Chr01:906292), need to find out the one without flag 2048
	 * SRR486615.17723780	163	Chr01	906292	60	38S55M	=	906354	174	CTTACATCATTGACTGCAAAAGAGAGTTCAAAGAAGAACTTGGCTTTGTGCACTAAGCCAAGTTCTTCTTTGAACTCTCTTTTGCAGTCAATN	DFHHDFFIIGJBCD@BC@<HICE@CDHBE?@DG;FD;D9??@@?B?FGBD<8=<B=37C8CG=D=D@@=CHA=??C>?BE<C>(6;@######	SA:Z:Chr01,906292,-,31S62M,60,3;	MC:Z:7M5D100M	MD:Z:54G0	PG:Z:MarkDuplicates	RG:Z:SRR486615	NM:i:1	AS:i:54	XS:i:0	
	 * SRR486615.17723780	2227	Chr01	906292	60	31H62M	=	906354	113	CTTGGCTTAGTGCACAAAGCCAAGTTCTTCTTTGAACTCTCTTTTGCAGTCAATGATGTAAG	GC8C73=B<=8<DBGF?B?@@??9D;DF;GD@?EBHDC@ECIH<@CB@DCBJGIIFFDHHFD	SA:Z:Chr01,906292,+,38S55M,60,1;	MC:Z:7M5D100M	MD:Z:8T6T43T2	PG:Z:MarkDuplicates	RG:Z:SRR486615	NM:i:3	AS:i:49	XS:i:0	
	 * SRR486615.17723780	83	Chr01	906354	60	7M5D100M	=	906292	-174	TTTGACCCAAAGATGACTCTGTTTTTTGTGGTCCAAATCACCCAACAAGCAGTAATGAGGATTTCACTCAACAAAGGATTCCCAAAGCGACTTCGTGCTTCTATCAC	DBDB@@<A@:CCC>>:9395?35552B?;@@@?7;;7;:-==E@E=;IF@AHBBEIIHF??8B9BGB4>GF???<?F?;C?;??BIG=HGC@:DBD?:=:ABDA=??	MC:Z:38S55M	MD:Z:7^ATGAT100PG:Z:MarkDuplicates	RG:Z:SRR486615	NM:i:5	AS:i:100	XS:i:0
	 * 
	 * @param reader
	 * @param rec
	 * @return
	 */
	public static SAMRecord queryNonSupplementaryMate(final SamReader reader, final SAMRecord rec) {
        if (!rec.getReadPairedFlag()) {
            throw new IllegalArgumentException("queryMate called for unpaired read.");
        }
        if (rec.getFirstOfPairFlag() == rec.getSecondOfPairFlag()) {
            throw new IllegalArgumentException("SAMRecord must be either first and second of pair, but not both.");
        }
        final boolean firstOfPair = rec.getFirstOfPairFlag();
        final CloseableIterator<SAMRecord> it;
        if (rec.getMateReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
            it = reader.queryUnmapped();
        } else {
            it = reader.queryAlignmentStart(rec.getMateReferenceName(), rec.getMateAlignmentStart());
        }
        try {
            SAMRecord mateRec = null;
            while (it.hasNext()) {
                final SAMRecord next = it.next();
                if (!next.getReadPairedFlag()) {
                    if (rec.getReadName().equals(next.getReadName())) {
                        throw new SAMFormatException("Paired and unpaired reads with same name: " + rec.getReadName());
                    }
                    continue;
                }
                /////==============================XT added below
                if(next.getSupplementaryAlignmentFlag()) //skip the supplementary alignment record that happens to be mapped at exactly the same position with the representative alignment that is the real target of the query (which should be very rare but not impossible)
                	continue;
                /////==============================XT added above
                
                
                if (firstOfPair) {
                    if (next.getFirstOfPairFlag()) continue;
                } else {
                    if (next.getSecondOfPairFlag()) continue;
                }
                if (rec.getReadName().equals(next.getReadName())) {
                    if (mateRec != null) {
                        throw new SAMFormatException("Multiple SAMRecord with read name " + rec.getReadName() +
                                " for " + (firstOfPair ? "second" : "first") + " end.");
                    }
                    mateRec = next;
                }
            }
            return mateRec;
        } finally {
            it.close();
        }
    }
}
