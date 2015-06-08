package conifer.cn;

import static blang.variables.IntegerVariable.intVar;

import java.io.File;
import java.util.Random;

import bayonet.distributions.BetaBinomial;
import bayonet.distributions.BetaBinomial.MeanPrecisionParameterization;
import bayonet.rplot.PlotHistogram;

public class TestBetaBinomial
{
    public static void main(String args[])
    {
        double mean = 0.5; 
        double precision = 100; 
        int trials = 10000; 
        int nSamples = 10000; 
        BetaBinomial<MeanPrecisionParameterization> beta = BetaBinomial.on(intVar(1)).withMeanPrecision(mean, precision, trials);
        
        Random random = new Random(1);

        double[] samples = new double[nSamples];
        for (int i = 0; i < nSamples; i++)
        {
            beta.generate(random);    
            samples[i] = beta.getRealization().getValue();
        }
        
        PlotHistogram plt = PlotHistogram.from(samples);
        plt.toPDF(new File("/Users/jewellsean/Desktop/test/hist.pdf"));
    }
}
