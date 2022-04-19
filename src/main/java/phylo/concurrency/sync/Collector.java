package phylo.concurrency.sync;

public class Collector {
	private int num=0;
	
	synchronized void collect() {
		this.num++;
		
		System.out.println(this.num);
	}
}
