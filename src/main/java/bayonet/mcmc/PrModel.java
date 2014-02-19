package bayonet.mcmc;

import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import org.jgrapht.Graphs;
import org.jgrapht.UndirectedGraph;

import com.google.common.collect.Sets;

import bayonet.graphs.GraphUtils;
import briefj.ReflexionUtils;



/**
 * A semi-automatic MCMC framework focused on extensibility to 
 * rich data structures. 
 * 
 * Difference from previous art:
 * - vs. STAN: STAN is currently limited to continuous r.v., bayonet supports many combinatorial ones
 * - vs. pymc: pymc does not support in-place modification, an important requirement in combinatorial spaces
 * - vs. JAGS: does not support custom data types for the random variables
 * 
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class PrModel
{
  @SuppressWarnings("rawtypes")
  private final UndirectedGraph graph = GraphUtils.newUndirectedGraph();
  
  @SuppressWarnings("rawtypes")
  private final Set 
    stochasticVariables = Sets.newLinkedHashSet(),
    observedVariables = Sets.newLinkedHashSet();
  
  public PrModel(Object specification)
  {
    parse(specification);
  }
  
  public void parse(Object o) 
  { 
    try {  _parse(o); } 
    catch (Exception e) { throw new RuntimeException(e); }
  }
  
  @SuppressWarnings("unchecked")
  public void addRelation(Factor f)
  {
    graph.addVertex(f);
    try { 
      if (!addLinks(f, f))
        throw new RuntimeException("A factor should contain at least one FactorArgument " +
        		"defining a random variable, or at least one FactorComponent recursively satisfying that property. " +
        		"See FactorArgument.defines()");
    } 
    catch (Exception e) { throw new RuntimeException(); }
  }
  
  @SuppressWarnings({ "rawtypes", "unchecked" })
  public Set getLatentStochasticVariables()
  {
    Set result = Sets.newLinkedHashSet(stochasticVariables);
    result.removeAll(observedVariables);
    return result;
  }
  
  @SuppressWarnings({ "rawtypes", "unchecked" })
  public Set getObservedVariables()
  {
    return Collections.unmodifiableSet(observedVariables);
  }
  
  @SuppressWarnings("unchecked")
  public void setObserved(Object variable)
  {
    observedVariables.add(variable);
  }
  
  @SuppressWarnings({ "rawtypes", "unchecked" })
  public <T> Set<T> getStochasticVariables(Class<T> ofType)
  {
    Set result = new LinkedHashSet();
    for (Object variable : getLatentStochasticVariables())
      if (ofType.isAssignableFrom(variable.getClass()))
        result.add(variable);
    return result;
  }
  
  @SuppressWarnings("unchecked")
  public List<Factor> neighbors(Object variable)
  {
    return Graphs.neighborListOf(graph, variable);
  }
  
  @SuppressWarnings("unchecked")
  private boolean addLinks(Factor f, Object o) throws IllegalArgumentException, IllegalAccessException
  {
    boolean variableDefined = false;
    for (Field argumentField : ReflexionUtils.getAnnotatedDeclaredFields(f.getClass(), FactorArgument.class, true))
    {
      if (!Modifier.isFinal(argumentField.getModifiers()))
        throw new RuntimeException("Fields annotated with FactorArgument should be final.");
      Object value = argumentField.get(o);
      if (value != null)
      {
        if (argumentField.getAnnotation(FactorArgument.class).defines())
        {
          stochasticVariables.add(value);
          variableDefined = true;
        }
        graph.addVertex(value);
        graph.addEdge(value, f);
      }
    }
    for (Field componentField : ReflexionUtils.getAnnotatedDeclaredFields(f.getClass(), FactorComponent.class, true))
    {
      if (!Modifier.isFinal(componentField.getModifiers()))
        throw new RuntimeException("Fields annotated with FactorComponent should be final.");
      Object value = componentField.get(o);
      if (value != null)
        if (addLinks(f, value))
          variableDefined = true;
    }
    return variableDefined;
  }
  
  private boolean _parse(Object o) throws IllegalArgumentException, IllegalAccessException  
  {
    List<Field> fields = ReflexionUtils.getAnnotatedDeclaredFields(o.getClass(), AddToModel.class, true);
    for (Field field : fields)
    {
      Object toAdd = field.get(o); 
      boolean isFactor = toAdd instanceof Factor;
      if (isFactor)
        addRelation((Factor) toAdd);
      boolean recDidSomething = _parse(toAdd);
      if (!isFactor && !recDidSomething)
        throw new RuntimeException("An annotation @AddToModel should be on a field of type Factor, or " +
            "should recursively contain fields with the same annotation and that property.");
    }
    return !fields.isEmpty();
  }
}
