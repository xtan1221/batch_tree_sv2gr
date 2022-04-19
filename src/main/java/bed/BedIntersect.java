package bed;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * note that bed file format is 0-based coordinate system
 * 
 * a bed file record with start = 0, end =100 is the region [1,100] or (0,100]
 * 		see https://genome.ucsc.edu/FAQ/FAQformat.html#format1
 * @author tanxu
 *
 */
public class BedIntersect {
	private final Path bedFile1;
	private final Path bedFile2;
	private final Path realPepNameList;
	
	
	private List<BedRecord> recordSet1;
	private List<BedRecord> recordSet2;
	private Set<String> realPepNameSet;
	
	private Map<BedRecord, List<BedRecord>> recordFromSet1OverlappingRecordsFromSet2ListMap;
	
	private Set<String> geneNameSetFromBedFile1OverlappingWithRecordFromSet2;
	
	
	public BedIntersect(Path bedFile1, Path bedFile2,Path realPepNameList) {
		super();
		this.bedFile1 = bedFile1;
		this.bedFile2 = bedFile2;
		this.realPepNameList=realPepNameList;
	}
	
	
	void readFiles() {
		this.recordSet1=new ArrayList<>();
		this.recordSet2=new ArrayList<>();
		this.realPepNameSet=new HashSet<>();
		
		try {
			////////////////////
			BufferedReader lineReader = new BufferedReader(new FileReader(this.bedFile1.toFile()));
			String line = null;
			
			while ((line = lineReader.readLine()) != null) {
				if(line.isEmpty()) {
					continue;
				}
				
				BedRecord bedRecord = new BedRecord(line);
				
				this.recordSet1.add(bedRecord);
			}
			
			lineReader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		try {
			////////////////////
			BufferedReader lineReader = new BufferedReader(new FileReader(this.bedFile2.toFile()));
			String line = null;
			
			while ((line = lineReader.readLine()) != null) {
				if(line.isEmpty()) {
					continue;
				}
				
				BedRecord bedRecord = new BedRecord(line);
				
				this.recordSet2.add(bedRecord);
			}
			
			lineReader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		try {
			////////////////////
			BufferedReader lineReader = new BufferedReader(new FileReader(this.realPepNameList.toFile()));
			String line = null;
			
			while ((line = lineReader.readLine()) != null) {
				if(line.isEmpty()) {
					continue;
				}
				
				//convert from pep name Sobic.004G007800.1.p to gene name Sobic.004G007800
				Pattern p=Pattern.compile("^(.+)\\.[\\d]+\\.p$");
				Matcher m=p.matcher(line);
				if(m.matches()) {
					this.realPepNameSet.add(m.group(1));
//					System.out.println(m.group(1));
				}else {
					System.out.println("invalid peptide name:"+line);
				}
				
			}
			
			lineReader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	void analysis() {
		this.geneNameSetFromBedFile1OverlappingWithRecordFromSet2=new HashSet<>();
		this.recordFromSet1OverlappingRecordsFromSet2ListMap = new LinkedHashMap<>();
		//identify any bed record from set 1 that is overlapping with any bed record from set 2
		
		//
		Collections.sort(this.recordSet1);
		Collections.sort(this.recordSet2);
		
		
		for(BedRecord rec2:this.recordSet2) {
			
			for(BedRecord rec1:this.recordSet1) {
				if(rec1.ref.equals(rec2.ref)) {
					if(rec1.overlappingWith(rec2)) {
						System.out.println(rec1.toString()+"   ====   "+rec2.toString());
						
						if(!this.recordFromSet1OverlappingRecordsFromSet2ListMap.containsKey(rec1)) {
							this.recordFromSet1OverlappingRecordsFromSet2ListMap.put(rec1, new ArrayList<>());
						}
						
						this.recordFromSet1OverlappingRecordsFromSet2ListMap.get(rec1).add(rec2);
						this.geneNameSetFromBedFile1OverlappingWithRecordFromSet2.add(rec1.name);
					}
				}
			}
		}
		
		
		System.out.println(this.recordFromSet1OverlappingRecordsFromSet2ListMap.size());
		
		//////////
		for(String s:this.geneNameSetFromBedFile1OverlappingWithRecordFromSet2) {
			if(!this.realPepNameSet.contains(s)) {
//				System.out.println(s);
			}
		}
	}
	
	
	
	static class BedRecord implements Comparable<BedRecord>{
		
		static final String NULL_VALUE_STRING=".";
		private final String ref;
		private final int start;
		private final int end;
		private final String name;
		private final Double score;
		private final String strand;
		
		BedRecord(String dataLine){
			String[] splits = dataLine.split("\t");
			
			this.ref=splits[0];
			int start1=Integer.parseInt(splits[1]);
			int end1=Integer.parseInt(splits[2]);
			if(start1<end1) {
				this.start=start1;
				this.end=end1;
			}else {
				this.start=end1;
				this.end=start1;
			}
			
			
			this.name=splits[3].equals(NULL_VALUE_STRING)?null:splits[3];
			
			this.score=splits[4].equals(NULL_VALUE_STRING)?null:Double.parseDouble(splits[4]);
			
			this.strand=splits[5].equals(NULL_VALUE_STRING)?null:splits[5];
		}

		boolean overlappingWith(BedRecord r2) {
			if(!r2.ref.equals(this.ref)) {
				return false;
			}
			//r2 fully cover this
			if(r2.start<=this.start && r2.end>=this.end) {
				return true;
			}
			//this fully cover r2
			if(r2.start>=this.start && r2.end<=this.end) {
				return true;
			}
			
			//partially overlapping
			if(r2.end>=this.end) { 
				return r2.start<this.end;
			}else {
				return r2.end>this.start;
			}
			
		}
		
		@Override
		public int compareTo(BedRecord r2) {
			if(this.ref.compareTo(r2.ref)!=0) {
				return this.ref.compareTo(r2.ref);
			}
			
			if(this.start!=r2.start) {
				return this.start-r2.start;
			}
			
			return this.end-r2.end;
		}



		@Override
		public String toString() {
			StringBuilder sb=new StringBuilder();
			sb.append(this.ref).append("\t")
			.append(this.start).append("\t")
			.append(this.end).append("\t")
			.append(this.name==null?NULL_VALUE_STRING:this.name).append("\t")
			.append(this.score==null?NULL_VALUE_STRING:this.score).append("\t")
			.append(this.strand==null?NULL_VALUE_STRING:this.strand);
			
			return sb.toString();
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + end;
			result = prime * result + ((name == null) ? 0 : name.hashCode());
			result = prime * result + ((ref == null) ? 0 : ref.hashCode());
			result = prime * result + ((score == null) ? 0 : score.hashCode());
			result = prime * result + start;
			result = prime * result + ((strand == null) ? 0 : strand.hashCode());
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (!(obj instanceof BedRecord))
				return false;
			BedRecord other = (BedRecord) obj;
			if (end != other.end)
				return false;
			if (name == null) {
				if (other.name != null)
					return false;
			} else if (!name.equals(other.name))
				return false;
			if (ref == null) {
				if (other.ref != null)
					return false;
			} else if (!ref.equals(other.ref))
				return false;
			if (score == null) {
				if (other.score != null)
					return false;
			} else if (!score.equals(other.score))
				return false;
			if (start != other.start)
				return false;
			if (strand == null) {
				if (other.strand != null)
					return false;
			} else if (!strand.equals(other.strand))
				return false;
			return true;
		}
		
		
	}
}
