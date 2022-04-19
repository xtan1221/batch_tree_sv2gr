package sv2gr.gr;

import population.sv.utils.SimpleSVType;

public class GREventUtils {
	
	
	/**
	 * 
	 * @param svType
	 * @param reverse
	 * @return
	 */
	public static GRType inferGRType(SimpleSVType svType, boolean reverse) {
		if(svType.equals(SimpleSVType.DEL)) {
			if(reverse) {
				return GRType.INS;
			}else {
				return GRType.DEL;
			}
		}else if(svType.equals(SimpleSVType.INS)) {
			if(reverse) {
				return GRType.DEL;
			}else {
				return GRType.INS;
			}
		}else if(svType.equals(SimpleSVType.INV)) {
			return GRType.INV;
		}else {
			throw new UnsupportedOperationException("given SV type is not supported for GR event inference:"+svType);
		}
		
		
	}
}
