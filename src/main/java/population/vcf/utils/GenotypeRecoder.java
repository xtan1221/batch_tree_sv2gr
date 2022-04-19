package population.vcf.utils;

import java.util.function.Function;

import basic.Pair;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public interface GenotypeRecoder extends Function<Pair<VariantContext, Genotype>, String>{
	
}
