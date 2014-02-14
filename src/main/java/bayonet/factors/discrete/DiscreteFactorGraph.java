package bayonet.factors.discrete;

import java.util.ArrayList;
import java.util.List;

import ma.RateMatrixLoader;
import nuts.math.GMFct;
import nuts.math.Graphs;
import nuts.math.TreeSumProd;

import org.apache.commons.lang3.tuple.Pair;
import org.ejml.simple.SimpleMatrix;
import org.jgrapht.UndirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;

import pty.SumT;
import pty.learn.DiscreteBP;
import pty.smc.models.CTMC.SimpleCTMC;

import com.google.common.collect.Lists;

import conifer.ml.data.PhylogeneticHeldoutDataset;
import conifer.ml.data.PhylogeneticHeldoutDataset.PhylogeneticHeldoutDatasetOptions;
import conifer.ml.tests.TestRealData;

import fig.basic.UnorderedPair;
import goblin.Taxon;

import bayonet.factors.BaseFactorGraph;
import bayonet.factors.BinaryFactor;
import bayonet.factors.FactorOperation;
import bayonet.factors.UnaryFactor;
import bayonet.factors.algo.EdgeSorter;
import bayonet.factors.algo.SumProduct;



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
  
  public void setUnary(V node, SimpleMatrix site2ValuePotentials)
  {
    unaries.put(node, new DiscreteUnaryFactor<V>(node, site2ValuePotentials, new int[site2ValuePotentials.numRows()]));
  }

  private final FactorOperation<V> discreteFactorGraphOperations = new FactorOperation<V>() {

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
      if (unariesOnMarginalized.size() == 0)
        throw new RuntimeException(); // TODO: might be able to return identity in some cases
      
      final int maxDegree = 2;
      if (unariesOnMarginalized.size() <= maxDegree)
      {
        final boolean hasTwo = unariesOnMarginalized.size() == 2;
        final DiscreteUnaryFactor<V> 
          dbf0 =           (DiscreteUnaryFactor) unariesOnMarginalized.get(0),
          dbf1 = hasTwo ? ((DiscreteUnaryFactor) unariesOnMarginalized.get(1)) : null;
        
        final SimpleMatrix 
          site2valuePotentials0 =          dbf0.site2valuePotentials,
          site2valuePotentials1 = hasTwo ? dbf1.site2valuePotentials : null;
          
        final int [] 
          scales0 =          dbf0.scales,
          scales1 = hasTwo ? dbf1.scales : null;
        
        final DiscreteBinaryFactor<V> binary = (DiscreteBinaryFactor) _binary;
        final SimpleMatrix o2mPotentials = binary.o2mPotentials;
        
        final int nSites = dbf0.nSites();
        final SimpleMatrix newMatrix = new SimpleMatrix(dbf0.nSites(), binary.nOtherVariableValues());
        final int [] newScales = new int[nSites];
        
        final int nOtherValues = binary.nOtherVariableValues();
        final int nMarginalizedValues = binary.nMarginalizedVariableValues();
        
        for (int site = 0; site < nSites; site++)
        {
          // scales are in log scale, so simply add them
          newScales[site] = hasTwo ? scales0[site] + scales1[site] : scales0[site];
        }
        
        for (int site = 0; site < nSites; site++)
        {
          // do the pointwise and marginalization
          for (int otherIndex = 0; otherIndex < nOtherValues; otherIndex++)
          {
            double sum = 0.0;
            if (hasTwo)
              for (int margIndex = 0; margIndex < nMarginalizedValues; margIndex++)
                sum += o2mPotentials.get(otherIndex, margIndex) 
                      * site2valuePotentials0.get(site, margIndex) 
                      * site2valuePotentials1.get(site, margIndex);
            else
              for (int margIndex = 0; margIndex < nMarginalizedValues; margIndex++)
                sum += o2mPotentials.get(otherIndex, margIndex) 
                      * site2valuePotentials0.get(site, margIndex);
            newMatrix.set(site, otherIndex, sum);
          }
        }
        
        return new DiscreteUnaryFactor<V>(binary.otherNode(), newMatrix, newScales);
      }
      else
        return marginalizeOnReducedUnariesDegree(this, maxDegree, _binary, unariesOnMarginalized);
    }
  };
  
  public static <V> UnaryFactor<V> marginalizeOnReducedUnariesDegree(
      FactorOperation<V> operation, 
      int maxDegree,
      final BinaryFactor<V> binary,
      final List<UnaryFactor<V>> unariesOnMarginalized)
  {
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
  
  public static class DiscreteBinaryFactor<V> implements BinaryFactor<V>
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
  
  public static final class DiscreteUnaryFactor<V> implements UnaryFactor<V>
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
          scales[site] += scales[site] + MIN_SCALE;
          for (int valueIndex = 0; valueIndex < nVariableValues(); valueIndex++)
            site2valuePotentials.set(site, valueIndex, site2valuePotentials.get(site,valueIndex) * UNDERFLOW_THRESHOLD_INVERSE);
          currentNorm = scales[site];
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
    phyloOptions.alignmentFile = "/Users/bouchard/Documents/data/utcs/23S.E/R0/cleaned.alignment.fasta";
    phyloOptions.treeFile = "/Users/bouchard/Documents/data/utcs/23S.E.raxml.nwk"; 
    PhylogeneticHeldoutDataset phyloData = PhylogeneticHeldoutDataset.loadData(phyloOptions);
    
    SimpleCTMC ctmc = new SimpleCTMC(RateMatrixLoader.k2p(), 1);
    GMFct<Taxon> pots = DiscreteBP.toGraphicalModel(phyloData.rootedTree, ctmc, phyloData.obs, 0);
    TreeSumProd<Taxon> tsp = new TreeSumProd<Taxon>(
        pots);
    
    System.out.println("method1 = " + tsp.logZ());
    
    DiscreteFactorGraph<Taxon> converted = fromGM(pots);
    
    System.out.println("Messages Fwd:" + EdgeSorter.newEdgeSorter(converted.getTopology(), new Taxon("internal66")).forwardMessages());
//    System.out.println("Messages Bwd:" + EdgeSorter.newEdgeSorter(converted.getTopology(), new Taxon("internal66")).backwardMessages());
    
    
    SumProduct<Taxon> sp = new SumProduct<Taxon>(converted);
    System.out.println(converted.getTopology().toString());
    System.out.println("method2 = " + sp.computeMarginal(new Taxon("internal66")).logNormalization());
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
      
      for (int valueIndex = 0; valueIndex < nValues; valueIndex++)
        newMatrix.set(0, valueIndex, model.get(vertex, valueIndex));
      
      newFG.setUnary(vertex, newMatrix);
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
