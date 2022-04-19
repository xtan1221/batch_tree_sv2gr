package genomics.geneAnnotation.analysis;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.List;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import genomics.geneAnnotation.GeneExonAnnotationGff3Reader;
import genomics.utils.Gene;
import genomics.utils.SimpleGenomicRegion;

public class GeneAnnotationStatSummary {
	private final GeneExonAnnotationGff3Reader geneExonAnnotationGff3Reader;
	private final Path outDir;
	//////////////////////////////////
	private Path geneLenOutFile;
	private Path intronLenOutFile;
	private Path exonLenOutFile;
	private Path CDSLenOutFile;
	private Path UTRLenOutFile;
	private Path UTR3LenOutFile;
	private Path UTR5LenOutFile;
	
	private SummaryStatistics geneLenSumStat;
	private SummaryStatistics intronLenSumStat;
	private SummaryStatistics exonLenSumStat;
	private SummaryStatistics CDSLenSumStat;
	private SummaryStatistics UTRLenSumStat;
	private SummaryStatistics UTR3LenSumStat;
	private SummaryStatistics UTR5LenSumStat;
	
	private SummaryStatistics exonNumPerGeneSumStat;
	private SummaryStatistics intronNumPerGeneSumStat;
	
	
	
	
	public GeneAnnotationStatSummary(GeneExonAnnotationGff3Reader geneExonAnnotationGff3Reader, Path outDir) {
		super();
		this.geneExonAnnotationGff3Reader = geneExonAnnotationGff3Reader;
		this.outDir = outDir;
		
		this.prepare();
		this.run();
		this.report();
	}
	
	
	void prepare() {
		this.geneLenSumStat=new SummaryStatistics();
		this.intronLenSumStat=new SummaryStatistics();
		this.exonLenSumStat=new SummaryStatistics();
		this.CDSLenSumStat=new SummaryStatistics();
		this.UTRLenSumStat=new SummaryStatistics();
		this.UTR3LenSumStat=new SummaryStatistics();
		this.UTR5LenSumStat=new SummaryStatistics();
		
		this.exonNumPerGeneSumStat=new SummaryStatistics();
		this.intronNumPerGeneSumStat=new SummaryStatistics();
		
		this.geneLenOutFile=Path.of(this.outDir.toString(), "gene.len.txt");
		this.intronLenOutFile=Path.of(this.outDir.toString(), "intron.len.txt");
		this.exonLenOutFile=Path.of(this.outDir.toString(), "exon.len.txt");
		this.CDSLenOutFile=Path.of(this.outDir.toString(), "CDS.len.txt");
		this.UTRLenOutFile=Path.of(this.outDir.toString(), "UTR.len.txt");
		this.UTR3LenOutFile=Path.of(this.outDir.toString(), "UTR3.len.txt");
		this.UTR5LenOutFile=Path.of(this.outDir.toString(), "UTR5.len.txt");
		
		this.checkOutFile(this.geneLenOutFile);
		this.checkOutFile(this.intronLenOutFile);
		this.checkOutFile(this.exonLenOutFile);
		this.checkOutFile(this.CDSLenOutFile);
		this.checkOutFile(this.UTRLenOutFile);
		this.checkOutFile(this.UTR3LenOutFile);
		this.checkOutFile(this.UTR5LenOutFile);
	}
	
	private void checkOutFile(Path outFile) {
		if(outFile.toFile().exists()) {
			System.out.println("outFile exists, delete it ...");
			outFile.toFile().delete();
		}
	}
	
	
	void run() {
		for(Gene gene:this.geneExonAnnotationGff3Reader.getGenes()) {
			this.intronNumPerGeneSumStat.addValue(gene.getLongestTranscript().getIntrons().size());
			this.exonNumPerGeneSumStat.addValue(gene.getLongestTranscript().getExons().size());
		}
		
		
		///////////////////////
		this.processLen(this.geneExonAnnotationGff3Reader.getLongestTranscripts(), this.geneLenSumStat, this.geneLenOutFile);
		this.processLen(this.geneExonAnnotationGff3Reader.getIntrons(), this.intronLenSumStat, this.intronLenOutFile);
		this.processLen(this.geneExonAnnotationGff3Reader.getExons(), this.exonLenSumStat, this.exonLenOutFile);
		this.processLen(this.geneExonAnnotationGff3Reader.getCDSs(), this.CDSLenSumStat, this.CDSLenOutFile);
		this.processLen(this.geneExonAnnotationGff3Reader.getUTRs(), this.UTRLenSumStat, this.UTRLenOutFile);
		////////////////
		this.processLen(this.geneExonAnnotationGff3Reader.get3PrimeUTRs(), this.UTR3LenSumStat, this.UTR3LenOutFile);
		this.processLen(this.geneExonAnnotationGff3Reader.get5PrimeUTRs(), this.UTR5LenSumStat, this.UTR5LenOutFile);
	}
	
	private <R extends SimpleGenomicRegion> void processLen(List<R> regions, SummaryStatistics stat, Path outFile) {
		try {
		    BufferedWriter writer = new BufferedWriter(new FileWriter(outFile.toFile()));

			for(R region:regions) {
				stat.addValue(region.getLen());
				writer.append(Integer.toString(region.getLen()));
				writer.newLine();
			}
		    
			writer.close();
		} catch (IOException ex) {
			System.err.println(ex);
		}
		
	}
	
	void report() {
		System.out.println("average intron per transcript:"+this.intronNumPerGeneSumStat.getMean());
		System.out.println("average exon per transcript:"+this.exonNumPerGeneSumStat.getMean());
		/////////////////////////
		this.print(this.geneLenSumStat, "gene");
		this.print(this.intronLenSumStat, "intron");
		this.print(this.exonLenSumStat, "exon");
		this.print(this.CDSLenSumStat, "CDS");
		this.print(this.UTRLenSumStat, "UTR");
		this.print(this.UTR3LenSumStat, "UTR3");
		this.print(this.UTR5LenSumStat, "UTR5");
	}
	
	private void print(SummaryStatistics stat, String type) {
		System.out.println("==========================");
		System.out.println(type);
		System.out.println("total num:" + stat.getN());
		System.out.println("total len:" + stat.getSum());
		System.out.println("mean:" + stat.getMean());
		System.out.println("sd:" + stat.getStandardDeviation());
	}
}
