package sv2gr.gr;


/**
 * 
 * @author tanxu
 *
 */
public enum GRType {
	DEL ("DEL"),
	INS ("INS"),
//	DUP ("DUP"), DUP type GR events is not included in current study
	INV ("INV");
	
	private final String stringValue;
	
	
	GRType(String stringValue){
		this.stringValue=stringValue;
	}

	/**
	 * @return the stringValue
	 */
	public String getStringValue() {
		return stringValue;
	}
	
	
}
