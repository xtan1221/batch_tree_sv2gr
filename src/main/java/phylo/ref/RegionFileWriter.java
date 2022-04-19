package phylo.ref;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Collection;
import java.util.List;

public class RegionFileWriter {
	private final Path outputFile;
	
	///////////////
	private FileWriter fileWriter;
	private BufferedWriter bufferedWriter;
	
	RegionFileWriter(Path outputFile) throws IOException{
		this.outputFile=outputFile;
		
		if(this.outputFile.toFile().exists()) {
			if(this.outputFile.toFile().isFile()) {
				this.outputFile.toFile().delete();
			}
		}
		
		this.fileWriter = new FileWriter(this.outputFile.toFile(), true);
		this.bufferedWriter = new BufferedWriter(this.fileWriter);
	}
	
	
	void run(List<Region> simpleRegionList) throws IOException {
		int index = 1;
		for(Region r:simpleRegionList) {
			String line = Integer.toString(index).concat("\t").concat(r.toString());
			this.bufferedWriter.append(line);
			this.bufferedWriter.newLine();
			index++;
		}
		this.bufferedWriter.close();
		this.fileWriter.close();
	}
	
	void run(Collection<Collection<Region>> regions) {
		//TODO
		throw new UnsupportedOperationException("not implemented yet!");
		
	}
	
	
}
