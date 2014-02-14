package bayonet.factors.discrete;

import fig.basic.UnorderedPair;
import goblin.Taxon;

import java.util.ArrayList;
import java.util.List;

import ma.RateMatrixLoader;
import nuts.math.GMFct;
import nuts.math.Graphs;

import org.apache.commons.lang3.tuple.Pair;
import org.ejml.simple.SimpleMatrix;
import org.jgrapht.UndirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;

import pty.learn.DiscreteBP;
import pty.smc.models.CTMC.SimpleCTMC;
import bayonet.factors.BaseFactorGraph;
import bayonet.factors.BinaryFactor;
import bayonet.factors.FactorOperation;
import bayonet.factors.UnaryFactor;
import bayonet.factors.algo.SumProduct;

import com.google.common.collect.Lists;

import conifer.ml.data.PhylogeneticHeldoutDataset;
import conifer.ml.data.PhylogeneticHeldoutDataset.PhylogeneticHeldoutDatasetOptions;



public class DiscreteFactorGraph<V> extends BaseFactorGraph<V>
{
  public DiscreteFactorGraph(UndirectedGraph<V, ?> topology)
  {
    super(topology);
  }

  @Override
  public FactorOperation<V> marginalizationOperation()
  {
    return discreteFactorGraphOperations;
  }
  
  public void setBinary(V mNode, V oNode, SimpleMatrix m2oPotentials)
  { 
    binaries.put(Pair.of(mNode, oNode), new DiscreteBinaryFactor<V>(m2oPotentials.transpose(), mNode, oNode));
    binaries.put(Pair.of(oNode, mNode), new DiscreteBinaryFactor<V>(m2oPotentials, oNode, mNode));
  }
  
  public void setUnary(V node, double [] values)
  {
    setUnaries(node, new SimpleMatrix(1, values.length, true, values));
  }
  
  private int nSites = -1;
  public void setUnaries(V node, SimpleMatrix site2ValuePotentials)
  {
    checkNSites(site2ValuePotentials.numRows());
    unaries.put(node, new DiscreteUnaryFactor<V>(node, site2ValuePotentials, new int[site2ValuePotentials.numRows()]));
  }
  private void checkNSites(int tested)
  {
    if (nSites == -1)
      nSites = tested;
    else if (nSites != tested)
      throw new RuntimeException();
  }



  private final FactorOperation<V> discreteFactorGraphOperations = new FactorOperation<V>() 
  {
    @SuppressWarnings({ "rawtypes", "unchecked" })
    @Override
    public UnaryFactor<V> pointwiseProduct(final List<UnaryFactor<V>> unaries)
    {
      final int nFactors = unaries.size();
      final DiscreteUnaryFactor [] cast = new DiscreteUnaryFactor[nFactors];
      for (int factorIndex = 0; factorIndex < nFactors; factorIndex++)
        cast[factorIndex] = (DiscreteUnaryFactor) unaries.get(factorIndex);
      
      final int nSites = cast[0].nSites();
      final int nVariableValues = cast[0].nVariableValues();
      
      final int [] newScales = new int[nSites];
      final SimpleMatrix newMatrix = new SimpleMatrix(nSites, nVariableValues);
      
      for (int site = 0; site < nSites; site++)
      {
        int sumScales = 0;
        for (int factor = 0; factor < nFactors; factor++)
          sumScales += cast[factor].scales[site];
        newScales[site] = sumScales;
      }
      
      for (int site = 0; site < nSites; site++)
        for (int varValue = 0; varValue < nVariableValues; varValue++)
        {
          double prodUnnorm = 1.0;
          for (int factor = 0; factor < nFactors; factor++)
            prodUnnorm *= cast[factor].site2valuePotentials.get(site, varValue);
          newMatrix.set(site, varValue, prodUnnorm);
        }
      
      return new DiscreteUnaryFactor(cast[0].node, newMatrix, newScales);
    }

    @Override
    @SuppressWarnings({ "unchecked", "rawtypes" })
    public UnaryFactor<V> marginalize(
        final BinaryFactor<V> _binary,
        final List<UnaryFactor<V>> unariesOnMarginalized)
    {
      final int maxDegree = 2;
      if (unariesOnMarginalized.size() <= maxDegree)
      {
        final int degree = unariesOnMarginalized.size();
        final DiscreteUnaryFactor<V> 
          dbf0 = degree >= 1 ? (DiscreteUnaryFactor) unariesOnMarginalized.get(0) : null,
          dbf1 = degree == 2 ? (DiscreteUnaryFactor) unariesOnMarginalized.get(1) : null;
        
        final SimpleMatrix 
          site2valuePotentials0 = degree >= 1 ? dbf0.site2valuePotentials : null,
          site2valuePotentials1 = degree == 2 ? dbf1.site2valuePotentials : null;
          
        final int [] 
          scales0 = degree >= 1 ? dbf0.scales : null,
          scales1 = degree == 2 ? dbf1.scales : null;
        
        final DiscreteBinaryFactor<V> binary = (DiscreteBinaryFactor) _binary;
        final SimpleMatrix o2mPotentials = binary.o2mPotentials;
        
        final SimpleMatrix newMatrix = new SimpleMatrix(nSites, binary.nOtherVariableValues());
        final int [] newScales = new int[nSites];
        
        final int nOtherValues = binary.nOtherVariableValues();
        final int nMarginalizedValues = binary.nMarginalizedVariableValues();
        
        // Warning: this part of the code is less readable and easy to maintain
        // because it is in the inner loop of phylogenetic computations
             if (degree == 0) ;
        else if (degree == 1) for (int site = 0; site < nSites; site++) newScales[site] = scales0[site];
        else                  for (int site = 0; site < nSites; site++) newScales[site] = scales0[site] + scales1[site];
             
         if (degree == 0) 
           for (int site = 0; site < nSites; site++)
             for (int otherIndex = 0; otherIndex < nOtherValues; otherIndex++)
             {
               double sum = 0.0;
               for (int margIndex = 0; margIndex < nMarginalizedValues; margIndex++)
                 sum += o2mPotentials.get(otherIndex, margIndex);
               newMatrix.set(site, otherIndex, sum);
             }
         else if (degree == 1) 
           for (int site = 0; site < nSites; site++)
             for (int otherIndex = 0; otherIndex < nOtherValues; otherIndex++)
             {
               double sum = 0.0;
               for (int margIndex = 0; margIndex < nMarginalizedValues; margIndex++)
                 sum += o2mPotentials.get(otherIndex, margIndex) 
                       * site2valuePotentials0.get(site, margIndex);
               newMatrix.set(site, otherIndex, sum);
             }
         else 
           for (int site = 0; site < nSites; site++)
             for (int otherIndex = 0; otherIndex < nOtherValues; otherIndex++)
             {
               double sum = 0.0;
               for (int margIndex = 0; margIndex < nMarginalizedValues; margIndex++)
                 sum += o2mPotentials.get(otherIndex, margIndex) 
                       * site2valuePotentials0.get(site, margIndex) 
                       * site2valuePotentials1.get(site, margIndex);
               newMatrix.set(site, otherIndex, sum);
             }
        
        return new DiscreteUnaryFactor<V>(binary.otherNode(), newMatrix, newScales);
      }
      else
        return marginalizeOnReducedUnariesDegree(this, maxDegree, _binary, unariesOnMarginalized);
    }
  };
  
  private static <V> UnaryFactor<V> marginalizeOnReducedUnariesDegree(
      FactorOperation<V> operation, 
      int maxDegree,
      final BinaryFactor<V> binary,
      final List<UnaryFactor<V>> unariesOnMarginalized)
  {
    System.out.println("hell");
    if (unariesOnMarginalized.size() <= maxDegree)
      throw new RuntimeException();
    ArrayList<UnaryFactor<V>> reducedList = Lists.newArrayList();
    // first items stay as is
    for (int i = 0; i < maxDegree - 1; i++)
      reducedList.add(unariesOnMarginalized.get(i));
    // last one is obtained by reducing the rest
    reducedList.add(operation.pointwiseProduct(unariesOnMarginalized.subList(maxDegree - 1,  unariesOnMarginalized.size())));
    if (reducedList.size() != maxDegree)
      throw new RuntimeException();
    return operation.marginalize(binary, reducedList);
  }
  
  private static class DiscreteBinaryFactor<V> implements BinaryFactor<V>
  {
    private final SimpleMatrix o2mPotentials;
    private final V m, o;
    
    private DiscreteBinaryFactor(SimpleMatrix o2mPotentials, V m, V o)
    {
      this.o2mPotentials = o2mPotentials;
      this.m = m;
      this.o = o;
    }

    public int nOtherVariableValues()
    {
      return o2mPotentials.numRows();
    }
    
    public int nMarginalizedVariableValues()
    {
      return o2mPotentials.numCols();
    }

    @Override public V marginalizedNode() { return m; }
    @Override public V otherNode() { return o; }
  }
  
  private static final class DiscreteUnaryFactor<V> implements UnaryFactor<V>
  {
    private static final int MIN_SCALE = -50;
    private static final double UNDERFLOW_THRESHOLD = Math.exp(MIN_SCALE);
    private static final double UNDERFLOW_THRESHOLD_INVERSE = Math.exp(-MIN_SCALE);
    
    private final SimpleMatrix site2valuePotentials;
    private final int [] scales; // base e
    private final double logNormalization;
    private final V node;
    
    public DiscreteUnaryFactor(V node, SimpleMatrix site2valuePotentials, int [] scales)
    {
      this.node = node;
      this.site2valuePotentials = site2valuePotentials;
      this.scales = scales;
      
      double logNorm = 0.0;
      double tempProd = 1.0;
      for (int site = 0; site < nSites(); site++)
      {
        double currentNorm = norm(site, site2valuePotentials);
        final int currentScale = scales[site];
        
        // update normalization
        logNorm = logNorm - currentScale;
        tempProd *= currentNorm;
        
        while (tempProd < UNDERFLOW_THRESHOLD)
        {
          tempProd *= UNDERFLOW_THRESHOLD_INVERSE;
          logNorm += MIN_SCALE;
        }
        
        // rescale if needed
        while (currentNorm < UNDERFLOW_THRESHOLD)
        {
          scales[site] = scales[site] - MIN_SCALE;
          for (int valueIndex = 0; valueIndex < nVariableValues(); valueIndex++)
            site2valuePotentials.set(site, valueIndex, site2valuePotentials.get(site,valueIndex) * UNDERFLOW_THRESHOLD_INVERSE);
          currentNorm = norm(site, site2valuePotentials); //scales[site];
        }
      }
      logNorm += Math.log(tempProd);
      
      this.logNormalization = logNorm;
    }

    public int nVariableValues()
    {
      return site2valuePotentials.numCols();
    }

    @Override public V connectedVariable() { return node; }
    
    public double logNormalization(int site)
    {
      return Math.log(norm(site, site2valuePotentials)) - scales[site];
    }
    
    private static double norm(int site, SimpleMatrix m)
    {
      double sum = 0.0;
      for (int valueIndex = 0; valueIndex < m.numCols(); valueIndex++)
        sum += m.get(site, valueIndex);
      return sum;
    }

    @Override
    public double logNormalization()
    {
      return logNormalization;
    }

    public int nSites()
    {
      return scales.length;
    }
    
  }
  
  public static void main(String [] args)
  {
    PhylogeneticHeldoutDatasetOptions phyloOptions = new PhylogeneticHeldoutDatasetOptions();
//    phyloOptions.alignmentFile = "/Users/bouchard/Documents/data/utcs/23S.E/R0/cleaned.alignment.fasta";
//    phyloOptions.treeFile = "/Users/bouchard/Documents/data/utcs/23S.E.raxml.nwk"; 
    phyloOptions.maxNSites = 100;
    phyloOptions.minFractionObserved = 0.9;
    PhylogeneticHeldoutDataset phyloData = PhylogeneticHeldoutDataset.loadData(phyloOptions);
    
    SimpleCTMC ctmc = new SimpleCTMC(RateMatrixLoader.k2p(), 1);
    GMFct<Taxon> pots = DiscreteBP.toGraphicalModel(phyloData.rootedTree, ctmc, phyloData.obs, 0);
    
    DiscreteFactorGraph<Taxon> converted = fromGM(pots);
    
    for (int i = 0; i < 100; i++)
    {
      long start = System.currentTimeMillis();
//      TreeSumProd<Taxon> tsp = new TreeSumProd<Taxon>(
//          pots);
//      
//      
//      System.out.println("method1 = " + tsp.logZ());
//      System.out.println("time = " + (System.currentTimeMillis()-start));
//      System.out.println();
//      
      
      
      start = System.currentTimeMillis();
      SumProduct<Taxon> sp = new SumProduct<Taxon>(converted);
      System.out.println("method2 = " + sp.computeMarginal(new Taxon("internal66")).logNormalization());
      System.out.println("time = " + (System.currentTimeMillis()-start));
    }
  }
  
  public static <V> DiscreteFactorGraph<V> fromGM(GMFct<V> model)
  {
    // create graph
    UndirectedGraph<V, DefaultEdge> ug = new SimpleGraph<V, DefaultEdge>(DefaultEdge.class);
    DiscreteFactorGraph<V> newFG = new DiscreteFactorGraph<V>(ug);
    
    // add vertex
    for (V vertex : model.graph().vertexSet())
    {
      ug.addVertex(vertex);
     
      int nValues = model.nStates(vertex);
      SimpleMatrix newMatrix = new SimpleMatrix(1, nValues);
      
      boolean shouldAdd = false;
      for (int valueIndex = 0; valueIndex < nValues; valueIndex++)
      {
        double value = model.get(vertex, valueIndex);
        newMatrix.set(0, valueIndex, value);
        if (value != 1.0)
          shouldAdd = true;
      }
      if (shouldAdd)
        newFG.setUnaries(vertex, newMatrix);
    }
    
    // add edges
    for (UnorderedPair<V, V> e : Graphs.edgeSet(model.graph()))
    {
      V n0 = e.getFirst(), 
        n1 = e.getSecond();
      ug.addEdge(n0, n1);
      
      int nValues0 = model.nStates(n0),
          nValues1 = model.nStates(n1);
      
      SimpleMatrix trans = new SimpleMatrix(nValues0, nValues1);
      
      for (int s0 = 0; s0 < nValues0; s0++)
        for (int s1 = 0; s1 < nValues1; s1++)
          trans.set(s0, s1, model.get(n0, n1, s0, s1));
      
      newFG.setBinary(n0, n1, trans);
    }
    
    return newFG;
  }
  
}
