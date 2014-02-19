package bayonet.mcmc.moves;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import com.google.common.collect.Lists;


import bayonet.mcmc.Factor;
import bayonet.mcmc.Move;
import bayonet.mcmc.MoveFactory;
import bayonet.mcmc.PrModel;
import bayonet.mcmc.RealVariable;


/**
 * A very simple MH move for real random variables, using a standard 
 * normal to propose.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class SimpleMHMove implements Move, MoveFactory
{
//  final SummaryStatistics acceptanceProbabilities = new SummaryStatistics();
  private final MoveableVariable variable;
  private final Collection<Factor> connectedFactors;
  
  public SimpleMHMove(MoveableVariable variable,
      Collection<Factor> connectedFactors)
  {
    this.variable = variable;
    this.connectedFactors = connectedFactors;
  }

  @Override
  public void execute(Random rand)
  {
    final double logDensityBefore = computeLogUnnormalizedPotentials();
    final double propRatio = variable.proposeInPlace(rand);
    final double logDensityAfter = computeLogUnnormalizedPotentials();
    final double ratio = Math.exp(propRatio + logDensityAfter - logDensityBefore);
    final boolean accept = rand.nextDouble() < ratio;
    variable.acceptRejectInPlace(accept);
//    acceptanceProbabilities.addValue(accept ? 1.0 : 0.0);
  }
  
  @Override
  public List<Move> build(PrModel model)
  {
    List<Move> result = Lists.newArrayList();
    
    for (MoveableVariable mv : model.getStochasticVariables(MoveableVariable.class))
      result.add(new SimpleMHMove(mv, model.neighbors(mv)));
    
    for (RealVariable realv : model.getStochasticVariables(RealVariable.class))
      result.add(new SimpleMHMove(new MoveableFromReal(realv), model.neighbors(realv)));
    
    return result;
  }
  
  private static class MoveableFromReal implements MoveableVariable
  {
    private final RealVariable real;
    private double old = Double.NaN;

    private MoveableFromReal(RealVariable real)
    {
      super();
      this.real = real;
    }

    @Override
    public double proposeInPlace(Random rand)
    {
      if (!Double.isNaN(old))
        throw new RuntimeException();
      old = real.getValue();
      final double newValue = old + rand.nextGaussian();
      real.setValue(newValue);
      return 0.0;
    }

    @Override
    public void acceptRejectInPlace(boolean accept)
    {
      if (!accept)
        real.setValue(old);
      old = Double.NaN;
    }
  }
  
  /**
   * Compute the part of the density that will be affected by 
   * chaning the variable held in this object.
   * 
   * @return
   */
  private double computeLogUnnormalizedPotentials()
  {
    double result = 0.0;
    for (Factor f : connectedFactors)
      result += f.logDensity();
    return result;
  }

  @SuppressWarnings("rawtypes")
  @Override
  public List variablesCovered()
  {
    return Collections.singletonList(variable);
  }
}