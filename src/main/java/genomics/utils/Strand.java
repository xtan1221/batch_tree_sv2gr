package genomics.utils;

/**
 * 
 * @author tanxu
 *
 */
public enum Strand {
	POSITIVE("+"),
	NEGATIVE("-"),
	UNKNOWN(".");
	
	/**
	 * string value consistent with the commonly used file format including bed, gff3 and nhmmer's output file formats
	 */
	private final String stringValue;
	
	Strand(String stringValue){
		this.stringValue=stringValue;
	}

	/**
	 * @return the stringValue
	 */
	public String getStringValue() {
		return stringValue;
	}
	
	/**
	 * return the corresponding {@link Strand} type based on given string value;
	 * @param stringValue
	 * @return
	 */
	public static Strand find(String stringValue) {
		if(stringValue == null) {
			throw new IllegalArgumentException("given stringValue cannot be null!");
		}
		if(stringValue.equals(POSITIVE.getStringValue())) {
			return POSITIVE;
		}else if(stringValue.equals(NEGATIVE.getStringValue())) {
			return NEGATIVE;
		}else if(stringValue.equals(UNKNOWN.getStringValue())) {
			return UNKNOWN;
		}else {
			throw new IllegalArgumentException("given stringValue cannot be recognized!");
		}
	}
}
