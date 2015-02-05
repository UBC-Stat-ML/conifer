package conifer.ctmc.cnv;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import org.ejml.simple.SimpleMatrix;

import com.google.common.collect.Maps;
import com.sun.org.apache.xalan.internal.xsltc.util.IntegerArray;

/**
 * Class to generate the simple copy number + point mutation process The size
 * parameter is the ploidy number of wildtype (A), mutant (a) respectively State
 * space is (A, a, b) where A, a \in {0,1,..., Size} and b \in {0,1}
 * 
 * For biological reasons, state space is restrict in b = 0 state to:
 * 
 * (0,0,0), (1,0,0), ... , (Size, 0, 0)
 * 
 * For b = 1 no restrictions are imposed, that is space is of cardinality (Size
 * + 1) * (Size + 1)
 * 
 * 
 * Total size of matrix is thus, Size + 1 + (Size + 1) * (Size + 1), which is
 * significantly less than naive implementation of 2 * (Size + 1) * (Size + 1)
 * 
 * @author jewellsean
 *
 */

public class GenerateCNMatrices {
	private final int size;
	private double[][] mutationQ;
	private double[][] increaseQ;
	private double[][] decreaseQ;

	/**
	 * 
	 * @param size
	 *            where Z_{0-size} * Z_{0-size}
	 */
	public GenerateCNMatrices(int size) {
		this.size = size + 1;
		this.mutationQ = null;
		this.increaseQ = null;
		this.decreaseQ = null;
	}

	private Integer[] ind2Space(int index) {
		Integer[] ind2Space = new Integer[3];

		// constrained state space
		if (index < size) {
			ind2Space[0] = index;
			ind2Space[1] = 0;
			ind2Space[2] = 0;
			return ind2Space;
		}

		index -= size;
		ind2Space[0] = index / size; // int division
		ind2Space[1] = index % size;
		ind2Space[2] = 1;
		return ind2Space;
	}

	private int space2Ind(Integer[] state) {
		if (state[2] == 0)
			return state[0];

		return state[0] * size + state[1] + size;
	}

	private <T extends MutationSpecification> double[][] generateSpecification(T mutationSpecification) {
		int len = size * (size + 1);
		double[][] rateMatrix = new double[len][len];
		for (int row = 0; row < len; row++) {
			Integer[] state = ind2Space(row);
			// System.out.println("----- Start state:( " + state[0] + ", " +
			// state[1] + ", " + state[2] + ")-----");
			List<Integer[]> transitions = mutationSpecification.getValues(state, size);

			for (Integer[] transState : transitions) {
				if (transState != null) {
					// System.out.println("transition state:( " + transState[0]
					// + ", " + transState[1] + ", " + transState[2] + ")");
					int col = space2Ind(transState);
					rateMatrix[row][col] = 1;
				}
			}
			// System.out.println("   ");

		}

		return rateMatrix;
	}

	private void generateIncrease() {
		double[][] increaseQ = generateSpecification(new increaseSpecification());
		this.increaseQ = increaseQ;
	}

	private void generateDecrease() {
		double[][] decreaseQ = generateSpecification(new decreaseSpecification());
		this.decreaseQ = decreaseQ;
	}

	private void generateMutation() {
		double[][] mutationQ = generateSpecification(new simpleMutationSpecification());
		this.mutationQ = mutationQ;
	}

	public void ensureInitalize() {
		generateDecrease();
		generateIncrease();
		generateMutation();
	}

	/**
	 * Assumes that the state is (A,a, b)
	 * 
	 * @author Sean Jewell (jewellsean@gmail.com)
	 *
	 */
	private interface MutationSpecification {
		List<Integer[]> getValues(Integer[] state, int size);
	}

	private static class simpleMutationSpecification implements MutationSpecification {
		@Override
		public List<Integer[]> getValues(Integer[] state, int size) {
			List<Integer[]> transition = new ArrayList<Integer[]>();
		    if (state[2].intValue() == 0)
		        transition.add(modState(-1, 1, state, size));
			return transition;
		}
	}

	private static class increaseSpecification implements MutationSpecification {
		@Override
		public List<Integer[]> getValues(Integer[] state, int size) {
			List<Integer[]> transition = new ArrayList<Integer[]>();
			transition.add(modState(0, 1, state, size));
			transition.add(modState(1, 0, state, size));
			return transition;
		}
	}

	private static class decreaseSpecification implements MutationSpecification {
		@Override
		public List<Integer[]> getValues(Integer[] state, int size) {
			List<Integer[]> transition = new ArrayList<Integer[]>();
			transition.add(modState(-1, 0, state, size));
			transition.add(modState(0, -1, state, size));
			return transition;
		}
	}

	private static Integer[] modState(int modWild, int modMutant, Integer[] state, int size) {
		boolean unrestrict = false;
		// cannot get anywhere from this state
		if (state[0] == 0 && state[1] == 0)
			return null;

		// Point mutation
		if (modWild + modMutant == 0) {
			state[2] = 1;
			unrestrict = true;
		}

		int s0 = state[0] + modWild;
		int s1 = state[1] + modMutant;

		// test for superfluous cn changes
		if (!unrestrict && ((state[0] == 0 && modWild != 0) || (state[1] == 0 && modMutant != 0))) {
			return null;
		}

		// boundaries in state space
		if (s0 < 0 || s0 > (size - 1) || s1 < 0 || s1 > (size - 1))
			return null;

		Integer[] modState = new Integer[3];
		modState[0] = s0;
		modState[1] = s1;
		modState[2] = state[2];
		return modState;
	}

	public double[][] getMutationQ() {
		return mutationQ;
	}

	public double[][] getIncreaseQ() {
		return increaseQ;
	}

	public double[][] getDecreaseQ() {
		return decreaseQ;
	}

	private String getJSONString(double[][] matrix) {
		String a = Arrays.deepToString(matrix);
		a = a.replace(", [", ",\n[");
		return ("{\"rateMatrix\":" + a + "}");
	}

	public void help() {
		
	}
	
	private String nameForType(String type) {
		return "cn-" + String.valueOf(size) + "-" + type + ".txt";
	}
	
	public void saveMatricesAsJSON() {
		Map<String, double[][]> matrices = Maps.newHashMap();

		matrices.put(nameForType("decreaseQ"), getDecreaseQ());
		matrices.put(nameForType("increaseQ"), getIncreaseQ());
		matrices.put(nameForType("mutationQ"), getMutationQ());

		for (String m : matrices.keySet()) {
			String jsonString = getJSONString(matrices.get(m));
			try {
				PrintWriter out = new PrintWriter("src/main/resources/conifer/ctmc/" + m);
				out.write(jsonString);
				out.close();
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}
	}

	public static void main(String args[]) {
		int size = 5;

		GenerateCNMatrices cn = new GenerateCNMatrices(size);
		cn.ensureInitalize();
		cn.saveMatricesAsJSON();
		
		SimpleMatrix decrease = new SimpleMatrix(cn.getDecreaseQ());
		System.out.println("Decrease cn:");
		System.out.println(decrease.toString());
		
		SimpleMatrix increase = new SimpleMatrix(cn.getIncreaseQ());
		System.out.println("Increase cn:");
		System.out.println(increase.toString());

		SimpleMatrix mutation = new SimpleMatrix(cn.getMutationQ());
		System.out.println("Mutation cn:");
		System.out.println(mutation.toString());
		 
	}

	/**
	 * 
	 * @param matrixType should be one of mutationQ, increaseQ, or decreaseQ
	 * @param size
	 * @return
	 */
	public static String getPathForMatrix(String matrixType, int size) 
	{
		GenerateCNMatrices cn = new GenerateCNMatrices(size);
		String filePath = "src/main/resources/conifer/ctmc/" + cn.nameForType(matrixType);
		File f = new File(filePath);
		if(!f.exists() || f.isDirectory()) {
			cn.ensureInitalize();
			cn.saveMatricesAsJSON();
		}
		
		//return ("src/main/resources/conifer/ctmc/" + cn.nameForType(matrixType));
		return ("/conifer/ctmc/" + cn.nameForType(matrixType));
//		return ("/Users/sohrab/project/conifercp/src/main/resources/conifer/ctmc/" + cn.nameForType(matrixType));
	}
}
