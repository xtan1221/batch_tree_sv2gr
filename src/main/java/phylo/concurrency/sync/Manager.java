package phylo.concurrency.sync;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class Manager {
	
	public static void main(String[] args) {
		ExecutorService es = Executors.newFixedThreadPool(1000); 
		
		Collector c = new Collector();
		
		for(int i=0;i<100;i++) {
			Test1 t = new Test1(c);
			es.submit(t);
		}
		
		
	}
	
}
