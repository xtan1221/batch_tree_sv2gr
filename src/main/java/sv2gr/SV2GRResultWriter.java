package sv2gr;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import basic.Pair;
import population.sv.utils.SimpleSVLocus;
import population.sv.utils.SimpleSVType;
import sv2gr.gr.GREvent;
import sv2gr.tree.RegionalTree;

public class SV2GRResultWriter {
	private final GREventsInfererAndTemporalOrderResolver GREventsInfererAndTemporalOrderResolver;
	
	private final Path outDir;
	
	///////////////////////////////
	private Path svToGREventTableFile;
	private Path temporallyOrderedGREventPairTableFile;
	private Path temporarllyOrderedGREventMatrixFile;
	
	public SV2GRResultWriter(
			sv2gr.GREventsInfererAndTemporalOrderResolver gREventsInfererAndTemporalOrderResolver,
			Path outDir) {
		super();
		GREventsInfererAndTemporalOrderResolver = gREventsInfererAndTemporalOrderResolver;
		this.outDir = outDir;
		
		this.prepare();
		this.writeToSvToGREventTableFile();
		this.writeToTemporallyOrderedGREventPairTableFile();
		this.writeToTemporarllyOrderedGREventMatrixFile();
	}

	
	

	void prepare() {
		this.svToGREventTableFile=Path.of(this.outDir.toString(), "SV2GR.recode.table.txt");
		
		
		this.temporallyOrderedGREventPairTableFile = Path.of(this.outDir.toString(), "temporally.ordered.GR.event.pairs.table.txt");
		
		this.temporarllyOrderedGREventMatrixFile=Path.of(this.outDir.toString(), "temporally.ordered.GR.event.pairs.matrix.txt");
		
		
		if(this.svToGREventTableFile.toFile().exists()) {
			System.out.println("svToGREventTableFile already exists, delete it...");
			this.svToGREventTableFile.toFile().delete();
		}
		
		if(this.temporallyOrderedGREventPairTableFile.toFile().exists()) {
			System.out.println("temporallyOrderedGREventPairTableFile already exists, delete it...");
			this.temporallyOrderedGREventPairTableFile.toFile().delete();
		}
		

		if(this.temporarllyOrderedGREventMatrixFile.toFile().exists()) {
			System.out.println("temporarllyOrderedGREventMatrixFile already exists, delete it...");
			this.temporarllyOrderedGREventMatrixFile.toFile().delete();
		}
		
	}
	
	
	void writeToSvToGREventTableFile() {
		
	    try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(this.svToGREventTableFile.toString()));
			writer.append(this.getSvToGREventTableFileHeaderLine());
			writer.newLine();
			
			for(RegionalTree tree:this.GREventsInfererAndTemporalOrderResolver.getRegionalTrees()) {
				for(SimpleSVLocus sv:tree.getCoveredSVLocusList()) {
					StringBuilder sb = new StringBuilder();
					
					sb.append(sv.getType()).append("\t")
					.append(sv.getChrom()).append("\t")
					.append(sv.getStart()).append("\t")
					.append(sv.getEnd()).append("\t")
					.append(tree.getSvLocusGREventMap().containsKey(sv)).append("\t")
					.append(tree.getSvLocusGREventMap().containsKey(sv)?
							tree.getSvLocusGREventMap().get(sv).isReversalOfSV()?"-":"+"
							:"N/A");
					
					writer.append(sb.toString());
					writer.newLine();
				}
			}
			
			writer.flush();
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	void writeToTemporallyOrderedGREventPairTableFile() {
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(this.temporallyOrderedGREventPairTableFile.toString()));
			writer.append(this.getTemporallyOrderedGREventPairTableFileHeaderLine());
			writer.newLine();
			
			for(RegionalTree tree:this.GREventsInfererAndTemporalOrderResolver.getTreeTemporallyOrderedGREventPairsMap().keySet()) {
				for(Pair<GREvent,GREvent> pair:this.GREventsInfererAndTemporalOrderResolver.getTreeTemporallyOrderedGREventPairsMap().get(tree)) {
					StringBuilder sb = new StringBuilder();
					
					sb.append(pair.getFirst().getOriginalSVLocus().getType()).append("\t")
					.append(pair.getFirst().getOriginalSVLocus().getChrom()).append("\t")
					.append(pair.getFirst().getOriginalSVLocus().getStart()).append("\t")
					.append(pair.getFirst().getOriginalSVLocus().getEnd()).append("\t")
					.append(tree.getSvLocusGREventMap().get(pair.getFirst().getOriginalSVLocus()).isReversalOfSV()?"-":"+").append("\t")
					
					.append(pair.getSecond().getOriginalSVLocus().getType()).append("\t")
					.append(pair.getSecond().getOriginalSVLocus().getChrom()).append("\t")
					.append(pair.getSecond().getOriginalSVLocus().getStart()).append("\t")
					.append(pair.getSecond().getOriginalSVLocus().getEnd()).append("\t")
					.append(tree.getSvLocusGREventMap().get(pair.getSecond().getOriginalSVLocus()).isReversalOfSV()?"-":"+");
					
					writer.append(sb.toString());
					writer.newLine();
				}
			}
			
			writer.flush();
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private void writeToTemporarllyOrderedGREventMatrixFile() {
		List<Pair<SimpleSVType, Boolean>> allIdentifiedGREventTypeList = new ArrayList<>();
		
		for(RegionalTree tree:this.GREventsInfererAndTemporalOrderResolver.getTreeTemporallyOrderedGREventPairsMap().keySet()) {
			for(Pair<GREvent, GREvent> pair:this.GREventsInfererAndTemporalOrderResolver.getTreeTemporallyOrderedGREventPairsMap().get(tree)) {
				Pair<SimpleSVType, Boolean> p1=new Pair<>(pair.getFirst().getOriginalSVLocus().getType(), pair.getFirst().isReversalOfSV());
				if(!allIdentifiedGREventTypeList.contains(p1))
					allIdentifiedGREventTypeList.add(p1);
				
				Pair<SimpleSVType, Boolean> p2=new Pair<>(pair.getSecond().getOriginalSVLocus().getType(), pair.getSecond().isReversalOfSV());
				if(!allIdentifiedGREventTypeList.contains(p2))
					allIdentifiedGREventTypeList.add(p2);
			}
		}
		
		Collections.sort(allIdentifiedGREventTypeList, 
				(a,b)->{
					if(a.getFirst().equals(b.getFirst())) {
						return a.getSecond().compareTo(b.getSecond());
					}else {
						return a.getFirst().compareTo(b.getFirst());
					}
				});
		
		/////////////////////////////
		List<List<Integer>> matrix=new ArrayList<>();
		
		for(int i=0;i<allIdentifiedGREventTypeList.size();i++) {
			List<Integer> list=new ArrayList<>();
			
			for(int j=0;j<allIdentifiedGREventTypeList.size();j++) {
				list.add(0);
			}
			
			matrix.add(list);
		}
		
		////////////////////////////////
		for(RegionalTree tree:this.GREventsInfererAndTemporalOrderResolver.getTreeTemporallyOrderedGREventPairsMap().keySet()) {
			for(Pair<GREvent, GREvent> pair:this.GREventsInfererAndTemporalOrderResolver.getTreeTemporallyOrderedGREventPairsMap().get(tree)) {
				Pair<SimpleSVType, Boolean> p1=new Pair<>(pair.getFirst().getOriginalSVLocus().getType(), pair.getFirst().isReversalOfSV());
				int p1Index=allIdentifiedGREventTypeList.indexOf(p1);
				
				Pair<SimpleSVType, Boolean> p2=new Pair<>(pair.getSecond().getOriginalSVLocus().getType(), pair.getSecond().isReversalOfSV());
				int p2Index=allIdentifiedGREventTypeList.indexOf(p2);
				
				matrix.get(p1Index).set(p2Index, matrix.get(p1Index).get(p2Index)+1);
			}
		}
		
		////////////////////////////////
		StringBuilder sb=new StringBuilder();
		sb.append("GR_type");
		for(Pair<SimpleSVType, Boolean> pair:allIdentifiedGREventTypeList) {
			sb.append("\t").append(pair.getFirst()).append(pair.getSecond()?"-":"+");
		}
		
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(this.temporarllyOrderedGREventMatrixFile.toString()));
			writer.append(sb.toString());
			writer.newLine();
			
			for(int i=0;i<matrix.size();i++) {
				Pair<SimpleSVType, Boolean> pair=allIdentifiedGREventTypeList.get(i);
				StringBuilder lineSb=new StringBuilder();
				lineSb.append(pair.getFirst()).append(pair.getSecond()?"-":"+");
				
				for(int j=0;j<matrix.get(i).size();j++) {
					lineSb.append("\t").append(matrix.get(i).get(j));
				}
				
				writer.append(lineSb.toString());
				writer.newLine();
			}
			
			writer.flush();
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	
	
	
	private String getSvToGREventTableFileHeaderLine() {
		StringBuilder sb = new StringBuilder();
		
		sb.append("sv_type").append("\t")
		.append("chrom").append("\t")
		.append("start").append("\t")
		.append("end").append("\t")
		.append("gr_event_inferred").append("\t") //whether the gr event inferred or not
		.append("gr_event_reversed"); //whether the inferred gr event is reversal of sv or not
		
		return sb.toString();
	}
	

	private String getTemporallyOrderedGREventPairTableFileHeaderLine() {
		StringBuilder sb = new StringBuilder();
		
		sb.append("sv_1_type").append("\t")
		.append("chrom1").append("\t")
		.append("start1").append("\t")
		.append("end1").append("\t")
		.append("gr_event_reversed_1").append("\t") //whether the inferred gr event is reversal of sv or not
		.append("sv_2_type").append("\t")
		.append("chrom2").append("\t")
		.append("start2").append("\t")
		.append("end2").append("\t")
		.append("gr_event_reversed_2");
		
		
		return sb.toString();
		
		
		
	}
	
	
}
