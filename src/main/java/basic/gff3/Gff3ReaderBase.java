package basic.gff3;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/**
 * base class for reading in data from a generic gff3 file;
 * 
 * @author tanxu
 * 
 */
public class Gff3ReaderBase {
	private final Path gff3File;
	
	///////////////
	
	private List<Gff3Hit> gff3Hits;
	
	
	public Gff3ReaderBase(Path gff3File) {
		super();
		this.gff3File = gff3File;
		
		this.read();
	}


	/**
	 * 
	 */
	private void read() {
		System.out.println("Start reading the gff3 file...");
		this.gff3Hits = new ArrayList<>();
		
		try {
			BufferedReader lineReader = new BufferedReader(new FileReader(gff3File.toFile()));
			String line = null;
			
			while ((line = lineReader.readLine()) != null) {
				if(line.startsWith("#")||line.isEmpty()) {
					continue;
				}
				
				///
				Gff3Hit hit = Gff3Hit.build(line);
				
				this.gff3Hits.add(hit);
			}
			
			lineReader.close();
		} catch (IOException ex) {
			System.err.println(ex);
		}
		
		System.out.println("Reading the gff3 file is done...");
	}


	/**
	 * @return the gff3Hits
	 */
	public List<Gff3Hit> getGff3Hits() {
		return gff3Hits;
	}
	
}
