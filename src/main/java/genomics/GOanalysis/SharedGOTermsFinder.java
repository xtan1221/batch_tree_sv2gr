package genomics.GOanalysis;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.HashSet;
import java.util.Set;

public class SharedGOTermsFinder {
	/**
	 * tab-delimited file with the following columns
	 * col1: population	
	 * col2: gene set	
	 * col3: GO term	
	 * col4: Ontology	
	 * col5: Description	
	 * col6: Number in input list	
	 * col7: Number in BG/Ref	
	 * col8: p-value	
	 * col9: FDR
	 */
	private final Path file = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv\\GO_analysis_all\\all.GO.terms.in.sb.zm.os.tsv");
	
	/////////////////////////////////
	private Set<String> sbTerms;
	private Set<String> osTerms;
	private Set<String> zmTerms;
	private Set<String> allTerms;
	
	private Set<String> termsSharedBySbAndZm;
	private Set<String> termsSharedBySbAndOs;
	private Set<String> termsSharedByOsAndZm;
	private Set<String> termsSharedBySbOsZm;
	
	
	public SharedGOTermsFinder() {
		super();
		
		this.readFile();
		this.findoutSharedGOTerms();
		this.outputResult();
	}


	void readFile() {
		this.sbTerms=new HashSet<>();
		this.osTerms=new HashSet<>();
		this.zmTerms=new HashSet<>();
		this.allTerms=new HashSet<>();
		
		try {
			BufferedReader lineReader = new BufferedReader(new FileReader(file.toFile()));
			String line = null;
			
			while ((line = lineReader.readLine()) != null) {
				if(line.startsWith("population") ||line.isEmpty()) {
					continue;
				}
				///
				String[] splits=line.split("\t");
				
				String population=splits[0];
				String GOTerm=splits[2];
				
				this.allTerms.add(GOTerm);
				
				if(population.equals("zm")) {
					this.zmTerms.add(GOTerm);
				}else if(population.equals("sb")) {
					this.sbTerms.add(GOTerm);
				}else if(population.equals("os")) {
					this.osTerms.add(GOTerm);
				}
				
			}
			
			lineReader.close();
		} catch (IOException ex) {
			System.err.println(ex);
		}
	}
	
	
	void findoutSharedGOTerms() {
		this.termsSharedByOsAndZm=new HashSet<>();
		this.termsSharedBySbAndOs=new HashSet<>();
		this.termsSharedBySbAndZm=new HashSet<>();
		this.termsSharedBySbOsZm=new HashSet<>();
		
		for(String term:this.allTerms) {
			if(this.sbTerms.contains(term) && this.osTerms.contains(term) && this.zmTerms.contains(term)) {
				this.termsSharedBySbOsZm.add(term);
			}else if(this.sbTerms.contains(term) && this.osTerms.contains(term) && !this.zmTerms.contains(term)) {
				this.termsSharedBySbAndOs.add(term);
			}else if(!this.sbTerms.contains(term) && this.osTerms.contains(term) && this.zmTerms.contains(term)) {
				this.termsSharedByOsAndZm.add(term);
			}else if(this.sbTerms.contains(term) && !this.osTerms.contains(term) && this.zmTerms.contains(term)) {
				this.termsSharedBySbAndZm.add(term);
			}
		}
	}
	
	void outputResult() {
		System.out.println("GO terms shared by all three populations:"+this.termsSharedBySbOsZm.size());
		this.termsSharedBySbOsZm.forEach(t->{System.out.println(t);});
		System.out.println("======================");
		System.out.println("GO terms shared by sb and os:"+this.termsSharedBySbAndOs.size());
		this.termsSharedBySbAndOs.forEach(t->{System.out.println(t);});
		System.out.println("======================");
		System.out.println("GO terms shared by sb and zm:"+this.termsSharedBySbAndZm.size());
		this.termsSharedBySbAndZm.forEach(t->{System.out.println(t);});
		System.out.println("======================");
		System.out.println("GO terms shared by os and zm:"+this.termsSharedByOsAndZm.size());
		this.termsSharedByOsAndZm.forEach(t->{System.out.println(t);});
		System.out.println("======================");
		
	}
	
	public static void main(String[] args) {
		SharedGOTermsFinder finder = new SharedGOTermsFinder();
	}
}
