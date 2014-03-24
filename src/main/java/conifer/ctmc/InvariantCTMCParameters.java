package conifer.ctmc;



public class InvariantCTMCParameters implements CTMCParameters
{
  private final RateMatrixToEmissionModel emissionModel;
  private final double[] initialDistribution;
  
  public InvariantCTMCParameters(double[] initialDistribution,
      RateMatrixToEmissionModel emissionModel)
  {
    this.emissionModel = emissionModel;
    this.initialDistribution = initialDistribution;
  }

  @Override
  public CTMC getProcess()
  {
    return new InvariantCTMC(initialDistribution);
  }

  @Override
  public RateMatrixToEmissionModel getEmissionModel()
  {
    return emissionModel;
  }

  @Override
  public double[][] getRateMatrix()
  {
    return new double[initialDistribution.length][initialDistribution.length];
  }
}
