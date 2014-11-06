package conifer.ctmc.cnv;

import org.ejml.simple.SimpleMatrix;

public class GenerateCNMatrices
{
  private final int size;
  private double [][] mutationQ;
  private double [][] increaseQ; 
  private double [][] decreaseQ; 

  /**
   * 
   * @param size where Z_{0-size} * Z_{0-size}
   */
  public GenerateCNMatrices(int size)
  {
    this.size = size;  
    this.mutationQ = null; 
    this.increaseQ = null;
    this.decreaseQ = null;
  }

  private int[] ind2Space(int index)
  {
    int [] ind2Space = new int[2];
    ind2Space[0] = index / size; // int division
    ind2Space[1] = index % size;
    return ind2Space;
  }

  private int space2Ind(int [] state)
  {
    return state[0] * size + state[1];
  }

  private <T extends MutationSpecification> double [][] generateSpecification(T mutationSpecification) 
  {
    double [][] rateMatrix = new double[size * size][size * size];
    for(int row = 0; row < size * size; row++)
    {
      int [] state = ind2Space(row); 
      state  = mutationSpecification.getValue(state, size);
      int col = space2Ind(state);
      rateMatrix[row][col] = 1;
    }

    return rateMatrix;
  }

  private void generateIncrease()
  {
    double [][] increaseQ = generateSpecification(new increaseSpecification());
    this.increaseQ = increaseQ;
  }

  private void generateDecrease()
  {
    double [][] decreaseQ = generateSpecification(new decreaseSpecification());
    this.decreaseQ = decreaseQ;
  }

  private void generateMutation()
  {
    double [][] mutationQ = generateSpecification(new simpleMutationSpecification());
    this.mutationQ  = mutationQ;
  }

  public void ensureInitalize()
  {
    generateDecrease();
    generateIncrease();
    generateMutation();
  }

  /**
   * Assumes that the state is (m,n)
   *   
   * @author Sean Jewell (jewellsean@gmail.com)
   *
   */
  private interface MutationSpecification
  {
    int[] getValue(int[] state, int size);
  }

  private static class simpleMutationSpecification implements MutationSpecification
  {
    @Override
    public int[] getValue(int[] state, int size)
    {
      if(state[1] > 0 & state[0] < (size - 1))
      {
        state[0] += 1;
        state[1] -= 1;
        return state;
      }
      return state;
    }
  }

  private static class decreaseSpecification implements MutationSpecification
  {
    @Override
    public int[] getValue(int[] state, int size)
    {
      if(state[0] > 0 && state[0] < (size - 1))
        state[0] += 1; 
      if(state[1] > 0 && state[1] < (size - 1))
        state[1] += 1;
      return state;
    }
  }

  private static class increaseSpecification implements MutationSpecification
  {
    @Override
    public int[] getValue(int[] state, int size)
    {
      if(state[0] > 0)
        state[0]-= 1; 
      if(state[1] > 0) 
        state[1] -= 1;
      return state;
    } 
  }

  public double [][] getMutationQ()
  {
    return mutationQ;
  }

  public double [][] getIncreaseQ()
  {
    return increaseQ;
  }


  public double [][] getDecreaseQ()
  {
    return decreaseQ;
  }

  public static void main(String args[])
  {
    int size = 3;
    GenerateCNMatrices cn = new GenerateCNMatrices(size);
    cn.ensureInitalize(); 
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
}
