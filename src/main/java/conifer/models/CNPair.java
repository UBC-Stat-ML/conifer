package conifer.models;

import conifer.io.CNObservationFactory;

// TODO: add checks on the valid value of copy numbers
public class CNPair {

	// the wildtype read count
	private int rA;

	// the mutant read count
	private int ra;

	public CNPair(int rA, int ra) {

		if (ra < 0 || rA < 0)
			throw new RuntimeException("Counts can't be negative.");
		if (ra > CNObservationFactory.maximumNormalCopyNumber || rA > CNObservationFactory.maximumMutantCopyNumber) {
			// throw new ValueException("Counts can't be over " +
			// CNObservationFactory.maximumMutantCopyNumber);
			//System.err.println("Value over CTMC bound!");
		}

		this.rA = rA;
		this.ra = ra;
	}

	// wild type
	public int getrA() {
		return rA;
	}
	
	// mutatnt type
	public int getRa() {
		return ra;
	}

	public int getN() {
		return (rA + ra);
	}

	public static CNPair withCounts(int rA, int ra) {
		return new CNPair(rA, ra);
	}

	@Override
	public String toString() {
		return "(" + rA + ", " + ra + ")";
	}
}
