package conifer.processors;

import java.io.File;
import java.util.List;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.time.StopWatch;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import bayonet.coda.CodaParser;
import bayonet.coda.EffectiveSize;
import bayonet.coda.SimpleCodaPlots;
import blang.processing.NodeProcessor;
import blang.processing.ProcessorContext;
import blang.variables.RealValued;
import briefj.OutputManager;
import briefj.run.Results;
import conifer.io.CopyNumberTreeObservation;

public class CopyNumberTreeProcessor implements NodeProcessor<CopyNumberTreeObservation>
{

    private CopyNumberTreeObservation variable;
    private OutputManager output = null;
    private OutputManager samplesOutput = null;
    private int interval = 4;
    private int current = 0;
    private String variableName = null;
    private boolean progress;
    private boolean CODA;
        
   
    @Override
    public void process(ProcessorContext context)
    {
        ensureInitialized(context);
        int iteration = context.getMcmcIteration();
        
        Parsimony parsimony = variable.parsimony;
        
        samplesOutput.write(variableName, "mcmcIter", iteration, variableName, variable.getValue());
        
        if (CODA)
        {
          current++;
          if (current == interval && progress)
          {
            interval = interval * 2;
            current = 0;
            generateCODA(iteration);       
          }

          if (context.isLastProcessCall())
            generateCODA(iteration);
        }
    }

    private void ensureInitialized(ProcessorContext context)
    {
        if (output != null)
            return;
          
          output = new OutputManager();
          samplesOutput = new OutputManager();
          if (variableName == null)
            variableName = context.getModel().getName(variable);
          csvSamples = new File(Results.getResultFolder(), variableName + "-csv");
          output.setOutputFolder(Results.getResultFolder());
          samplesOutput.setOutputFolder(csvSamples);
          progress = context.getOptions().progressCODA;
          CODA = context.getOptions().CODA;
    }

    @Override
    public void setReference(CopyNumberTreeObservation variable)
    {
        this.variable = variable;
    }

    public StopWatch watch = new StopWatch();
    private File csvSamples;
    
    private void generateCODA(int iteration)
    {
      samplesOutput.flush();
      File 
      indexFile = new File(csvSamples, "CODAindex.txt"),
      chainFile = new File(csvSamples, "CODAchain1.txt");
      CodaParser.CSVToCoda(indexFile, chainFile, csvSamples);
      SimpleCodaPlots codaPlots = new SimpleCodaPlots(chainFile, indexFile);
      codaPlots.toPDF(new File(csvSamples, "codaPlots.pdf"));
      
      EffectiveSize essCalculator = new EffectiveSize(chainFile, indexFile);
      List<Double> essValues = null;
      try {
      essValues =  essCalculator.getESSValues();
      } catch (Exception e)
      {
        return;
      }
      if (essValues.size() != 1)
        throw new RuntimeException();
      double ess = essValues.get(0);
      double time = watch.getTime() / 1000.0;
      double essPerSec = ess/time;
      output.printWrite(variableName + "-ess", "iteration", iteration, "ess", ess, "time", time, "essPerSec", essPerSec);
    }
    
}
