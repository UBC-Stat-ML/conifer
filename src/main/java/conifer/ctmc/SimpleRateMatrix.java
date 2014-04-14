package conifer.ctmc;

import java.io.File;
import java.util.Arrays;

import com.google.gson.Gson;


import briefj.BriefIO;



public class SimpleRateMatrix implements CTMCParameters
{
  private final double [][] rateMatrix;
  private final RateMatrixToEmissionModel emissionModel;
  
  public SimpleRateMatrix(double[][] rateMatrix,
      RateMatrixToEmissionModel emissionModel)
  {
    this.rateMatrix = rateMatrix;
    this.emissionModel = emissionModel;
  }
  
  @Override
  public double[][] getRateMatrix()
  {
    return rateMatrix;
  }

  @Override
  public CTMC getProcess()
  {
    return new EigenCTMC(rateMatrix);
  }

  @Override
  public RateMatrixToEmissionModel getEmissionModel()
  {
    return emissionModel;
  }
  
  public int nStates() { return rateMatrix.length; }
  
  public static SimpleRateMatrix fromJSONFile(File jsonFile)
  {
    String jsonString = BriefIO.fileToString(jsonFile);
    return fromJSONString(jsonString);
  }
  
  private static SimpleRateMatrix fromJSONString(String jsonString)
  {
    return new Gson().fromJson(jsonString, SimpleRateMatrix.class);
  }
  public static SimpleRateMatrix fromResource(String resourceURL)
  {
    String jsonString = BriefIO.resourceToString(resourceURL); 
    return fromJSONString(jsonString);
  }
  
  @Override
  public String toString()
  {
    return "SimpleRateMatrix [rateMatrix=" + Arrays.deepToString(rateMatrix)
        + ", emissionModel=" + emissionModel + "]";
  }


}
