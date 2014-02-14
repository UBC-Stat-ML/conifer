package bayonet.factors;

import java.util.Map;

import org.apache.commons.lang3.tuple.Pair;
import org.jgrapht.UndirectedGraph;

import com.google.common.collect.Maps;



public abstract class BaseFactorGraph<V> implements FactorGraph<V>
{
  protected final Map<V, UnaryFactor<V>> unaries = Maps.newHashMap();
  protected final Map<Pair<V,V>, BinaryFactor<V>> binaries = Maps.newHashMap();
  protected final UndirectedGraph<V, ?> topology;
  
  public BaseFactorGraph(UndirectedGraph<V, ?> topology)
  {
    this.topology = topology;
  }

  @Override
  public UndirectedGraph<V, ?> getTopology()
  {
    return topology;
  }
  
  @Override
  public UnaryFactor<V> getUnary(V node)
  {
    if (!topology.containsVertex(node))
      throw new RuntimeException();
    return unaries.get(node);
  }
  
  @Override
  public BinaryFactor<V> getBinary(V marginalizedNode, V otherNode)
  {
    if (!topology.containsEdge(marginalizedNode, otherNode))
      throw new RuntimeException();
    return binaries.get(Pair.of(marginalizedNode, otherNode));
  }
}
