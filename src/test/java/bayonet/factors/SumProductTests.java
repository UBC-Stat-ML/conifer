package bayonet.factors;

import java.util.Random;

import org.ejml.ops.RandomMatrices;
import org.ejml.simple.SimpleMatrix;
import org.jgrapht.UndirectedGraph;

import bayonet.factors.algo.SumProduct;
import bayonet.factors.discrete.DiscreteFactorGraph;
import bayonet.graphs.GraphUtils;



public class SumProductTests
{
  
  public static DiscreteFactorGraph<Integer> buildRandomMarkov(Random rand, int nStates, int length)
  {
    // build topology
    UndirectedGraph<Integer, ?> topology = GraphUtils.newUndirectedGraph();
    for (int i = 0; i < length; i++)
    {
      topology.addVertex(i);
      if (i > 0)
        topology.addEdge(i-1, i);
    }
    
    // build potentials
    DiscreteFactorGraph<Integer> result = new DiscreteFactorGraph<Integer>(topology);
    for (int i = 0; i < length; i++)
    {
      // create unary
      double [] data = //new double[nStates];//
        RandomMatrices.createRandom(1, nStates, 0, 1.0/nStates, rand).data;
//      for (int s = 0; s < nStates; s++)
//        data[s] = 1.0/nStates/2.0;
      result.setUnary(i, data);
      
      if (i > 0) 
      {
        // create binary
        SimpleMatrix matrix = //new SimpleMatrix(nStates, nStates);
          new SimpleMatrix(RandomMatrices.createRandom(nStates, nStates, 0, 1.0/nStates, rand));
        
//        for (int s1 = 0; s1 < nStates; s1++)
//          for (int s2 = 0; s2 < nStates; s2++)
//            matrix.set(s1, s2, 1.0);
        
        result.setBinary(i-1, i, matrix);
      }
        
    }
    
    return result;
  }
  
  public static void main(String [] args)
  {
    Random rand = new Random(1);
    
    int len = 8000;
    DiscreteFactorGraph<Integer> markov = buildRandomMarkov(rand, 2, len);
    long start = System.currentTimeMillis();
    SumProduct<Integer> sp = new SumProduct<Integer>(markov);
    System.out.println(sp.logNormalization());
//    System.out.println(-len * Math.log(2.0));
    System.out.println("Time: " + (System.currentTimeMillis() - start)/1000.0);
  }
}
