package circos;

public class WindowedNumericValueRecord {
	private final String chr;
	private final int start;
	private final int end;
	private final double value;
	
	public WindowedNumericValueRecord(String chr, int start, int end, double value) {
		super();
		this.chr = chr;
		this.start = start;
		this.end = end;
		this.value = value;
	}

	/**
	 * @return the chr
	 */
	public final String getChr() {
		return chr;
	}

	/**
	 * @return the start
	 */
	public final int getStart() {
		return start;
	}

	/**
	 * @return the end
	 */
	public final int getEnd() {
		return end;
	}

	/**
	 * @return the value
	 */
	public final double getValue() {
		return value;
	}
}
