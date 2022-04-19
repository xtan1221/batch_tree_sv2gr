package population.vcf.utils;

import java.util.function.Predicate;

import htsjdk.variant.variantcontext.VariantContext;

public interface VariantContextFilter extends Predicate<VariantContext>{
	
	
	default VariantContextFilter and(VariantContextFilter f) {
		return e->{
			return this.test(e) && f.test(e);
		};
	}
}
