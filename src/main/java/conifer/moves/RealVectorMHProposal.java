package conifer.moves;

import java.util.List;
import java.util.Random;

import blang.factors.Factor;
import blang.mcmc.ConnectedFactor;
import blang.mcmc.MHProposalDistribution;
import blang.mcmc.SampledVariable;
import blang.variables.RealVectorInterface;
import briefj.opt.Option;

// TODO: this should go in blang!

public class RealVectorMHProposal implements MHProposalDistribution
{
  public static double bandWidth = 0.01;
  
  @SampledVariable RealVectorInterface variable;
  
  @ConnectedFactor List<Factor> connectedFactors;
  
  private double [] savedValue = null;

  @Override
  public Proposal propose(Random rand)
  {
    //System.out.println("Computing RealVectorHMProposal");
    if (savedValue != null)
      throw new RuntimeException();
    double [] variableArray = variable.getVector();
    savedValue = variableArray.clone();
    
    for (int i = 0; i < variableArray.length; i++)
    {
      variableArray[i] += bandWidth * rand.nextGaussian();
    }
    variable.setVector(variableArray);
    
    return new ProposalRealization();
  }
  
  private class ProposalRealization implements Proposal
  {

    @Override
    public double logProposalRatio()
    {
      return 0;
    }

    @Override
    public void acceptReject(boolean accept)
    {
      if (!accept)
      {
        double [] variableArray = variable.getVector();
        for (int i = 0; i < variableArray.length; i++)
          variableArray[i] = savedValue[i];
        variable.setVector(variableArray);
      }
      savedValue = null;
    }
    
  }


}
