package conifer.processors;

import java.io.File;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.time.StopWatch;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import com.beust.jcommander.internal.Maps;

import conifer.ctmc.expfam.CTMCExpFam;
import conifer.ctmc.expfam.CTMCState;
import conifer.ctmc.expfam.ExpFamParameters;
import bayonet.coda.CodaParser;
import bayonet.coda.EffectiveSize;
import bayonet.coda.SimpleCodaPlots;
import blang.processing.NodeProcessor;
import blang.processing.ProcessorContext;
import briefj.OutputManager;
import briefj.collections.Counter;
import briefj.run.Results;


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
  private OutputManager samplesOutput = null, output = null;
  private int interval = 10;
  private int current = 0;
  private boolean progressInfo = false;

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
      samplesOutput.write(key, "mcmcIter", context.getMcmcIteration(), key, weights.getCount(key));
       

//    Counter<CTMCState> stationary = model.getStationaryDistribution();
//    for (CTMCState state0 : stationary.keySet())
//      samplesOutput.write("stationary("+state0.toString()+")", "mcmcIter", context.getMcmcIteration(), state0, stationary.getCount(state0));
//
//
//    List<CTMCState> states = parameters.globalExponentialFamily.stateIndexer.objectsList();
//    for (CTMCState state0 : states)
//    {
//      Counter<CTMCState> rates = model.getRates(state0);
//      for (CTMCState state1 : rates.keySet())
//      {
//        String key = "q(" + state0 + "," + state1 + ")";
//        samplesOutput.write(key, "mcmcIter", context.getMcmcIteration(), key, rates.getCount(state1));
//      }
//    }
    
    current++;
    if (context.getOptions().CODA && ((progressInfo && current == interval) || context.isLastProcessCall()))
      try
      {
        samplesOutput.flush();
        interval = interval * 2;
        current = 0;
       // File
       //   indexFile = new File(samplesOutput.getOutputFolder(), "CODAindex.txt"),
       //   chainFile = new File(samplesOutput.getOutputFolder(), "CODAchain1.txt");
        //CodaParser.CSVToCoda(indexFile, chainFile, samplesOutput.getOutputFolder());
        //SimpleCodaPlots codaPlots = new SimpleCodaPlots(chainFile, indexFile);
        //codaPlots.toPDF(new File(samplesOutput.getOutputFolder(), "codaPlots.pdf"));
        
        //EffectiveSize essCalculator = new EffectiveSize(chainFile, indexFile);
        //List<Double> essValues = essCalculator.getESSValues();
        //DescriptiveStatistics stats = new DescriptiveStatistics();
        //for (double val : essValues) stats.addValue(val);
        //double time = watch.getTime() / 1000.0;
        //String variableName = context.getModel().getName(parameters);
        //Map<String,Double> summaries = Maps.newHashMap();
        //summaries.put("max", stats.getMax());
        //summaries.put("median", stats.getPercentile(50));
        //summaries.put("min", stats.getMin());
        // These following two lines are added by Tingting
        //summaries.put("25th Quantile", stats.getPercentile(25));
        //summaries.put("75th Quantile", stats.getPercentile(75));
//        for (String statKey : summaries.keySet())
//        {
//          double ess = summaries.get(statKey);
//          double essPerSec = ess/time;
//          output.printWrite(variableName  + "-ess", "iteration", context.getMcmcIteration(), "statistic", statKey, "ess", ess, "time", time, "essPerSec", essPerSec);
//
//        }
      }
      catch (Exception e)
      {
        System.err.println("Warning: plot not working because R not found.");
      }
  }
  
  public StopWatch watch = new StopWatch();

  private void ensureInitialized(ProcessorContext context)
  {
    if (samplesOutput != null)
      return;
    progressInfo = context.getOptions().progressCODA && context.getOptions().CODA;
    watch.start();
    samplesOutput = new OutputManager();
    output = new OutputManager();
    String varName = context.getModel().getName(parameters);
    File csvSamples = new File(Results.getResultFolder(), varName + "-csv");
    samplesOutput.setOutputFolder(csvSamples);
    output.setOutputFolder(Results.getResultFolder());
  }

    @Override
    public void setReference(ExpFamParameters variable)
    {
        this.parameters = variable;
    }

}



