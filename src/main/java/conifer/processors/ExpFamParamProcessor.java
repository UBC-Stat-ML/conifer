package conifer.processors;

import java.io.File;
import java.util.List;

import conifer.ctmc.expfam.CTMCExpFam;
import conifer.ctmc.expfam.CTMCState;
import conifer.ctmc.expfam.ExpFamParameters;

import bayonet.coda.CodaParser;
import bayonet.coda.SimpleCodaPlots;
import blang.processing.NodeProcessor;
import blang.processing.ProcessorContext;
import briefj.OutputManager;
import briefj.collections.Counter;
import briefj.tomove.Results;


/**
 * Store and plot samples from the exponential family-based representation of
 * rate matrices.
 * 
 * Two types of statistics are stored: the weights (natural parameterization) 
 * used in the exponential family, and the values of the rates in the rate matrix.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class ExpFamParamProcessor implements NodeProcessor<ExpFamParameters>
{
  private ExpFamParameters parameters;
  private OutputManager output = null;
  private int interval = 2;
  private int current = 0;

  @Override
  public void process(ProcessorContext context)
  {
    ensureInitialized(context);
    
    CTMCExpFam<CTMCState>.LearnedReversibleModel model = parameters.getModel();
    
    // TODO: split into different folders
    // TODO: better viz (remove dep on coda at same time)
    // TODO: add stationary distribution estimates
    
    // print weights
    Counter<String> weights = model.getWeights();
    for (String key : weights.keySet())
      output.write(key, "mcmcIter", context.getMcmcIteration(), key, weights.getCount(key));
       

    Counter<CTMCState> stationary = model.getStationaryDistribution();
    for (CTMCState state0 : stationary.keySet())
      output.write("stationary("+state0.toString()+")", "mcmcIter", context.getMcmcIteration(), state0, stationary.getCount(state0));
       
    
    List<CTMCState> states = parameters.globalExponentialFamily.stateIndexer.objectsList();
    for (CTMCState state0 : states)
    {
      Counter<CTMCState> rates = model.getRates(state0);
      for (CTMCState state1 : rates.keySet())
      {
        String key = "q(" + state0 + "," + state1 + ")";
        output.write(key, "mcmcIter", context.getMcmcIteration(), key, rates.getCount(state1));
      }
    }
    output.flush();
    current++;
    if (current == interval)
    {
      interval = interval * 2;
      current = 0;
      File 
        indexFile = new File(output.getOutputFolder(), "CODAindex.txt"),
        chainFile = new File(output.getOutputFolder(), "CODAchain1.txt");
      CodaParser.CSVToCoda(indexFile, chainFile, output.getOutputFolder());
      SimpleCodaPlots codaPlots = new SimpleCodaPlots(chainFile, indexFile);
      codaPlots.toPDF(new File(output.getOutputFolder(), "codaPlots.pdf"));
    }
  }

  private void ensureInitialized(ProcessorContext context)
  {
    if (output != null)
      return;
    output = new OutputManager();
    String varName = context.getModel().getName(parameters);
    File csvSamples = new File(Results.getResultFolder(), varName + "-csv");
    output.setOutputFolder(csvSamples);
  }

  @Override
  public void setReference(ExpFamParameters variable)
  {
    this.parameters = variable;
  }

}
