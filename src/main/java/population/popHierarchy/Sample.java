package population.popHierarchy;

public class Sample {
	/**
	 * order index in the population
	 */
	private final int index;
	/**
	 * population level;
	 * see {@link PopulationStructureFileReader}
	 */
	private final int populationLevel;
	/**
	 * name of the sample
	 */
	private final String name;
	/**
	 * taxon
	 */
	private final String taxon;
	
	
	public Sample(int index, int populationLevel, String name, String taxon) {
		super();
		this.index = index;
		this.populationLevel = populationLevel;
		this.name = name;
		this.taxon = taxon;
	}


	/**
	 * @return the index
	 */
	public final int getIndex() {
		return index;
	}


	/**
	 * @return the populationLevel
	 */
	public final int getPopulationLevel() {
		return populationLevel;
	}


	/**
	 * @return the name
	 */
	public final String getName() {
		return name;
	}


	/**
	 * @return the taxon
	 */
	public final String getTaxon() {
		return taxon;
	}
	
	
}
