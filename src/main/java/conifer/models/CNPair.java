package conifer.models;

import com.sun.tools.internal.ws.wsdl.document.jaxws.Exception;

// TODO: add checks on the valid value of copy numbers
public class CNPair {

	// the normal read count
	private int ra;
	
	// the mutant read count
	private int rA;
	
	public CNPair(int ra, int rA) {
		
		if (ra < 0 || rA < 0) throw new RuntimeException("Counts can't be negative.");
		
		this.ra = ra;
		this.rA = rA;
	}
	
	public int getRa() {
		return ra;
	}
	
	public int getrA() {
		return rA;
	}
	
	public static CNPair withCounts(int ra, int rA) {
		return new CNPair(ra, rA);
	}
	
	@Override
	public String toString() {
		return "(" + ra + ", " + rA + ")";
	}
}
