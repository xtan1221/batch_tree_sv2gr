package phylo.concurrency.sync;

public class Test1 implements Runnable{
	private final Collector collector;
	
	Test1(Collector collector){
		this.collector = collector;
	}
	
	
	@Override
	public void run() {
		try {
			Thread.sleep(3000);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		this.collector.collect();
	}
	
}
