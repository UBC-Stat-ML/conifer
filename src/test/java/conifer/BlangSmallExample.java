package conifer;

import java.io.File;
import java.util.List;

import tutorialj.Tutorial;

import com.google.common.collect.Lists;

import bayonet.distributions.Exponential;
import bayonet.distributions.Exponential.MeanParameterization;
import bayonet.distributions.Uniform;
import bayonet.distributions.Uniform.MinMaxParameterization;
import bayonet.rplot.PlotHistogram;
import blang.MCMCRunner;
import blang.annotations.DefineFactor;
import blang.processing.ProcessorContext;
import  blang.variables.RealVariable;
import briefj.tomove.Results;

import static  blang.variables.RealVariable.*;

/**
 * A small example to demonstrate the basic features of Blang.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class BlangSmallExample extends MCMCRunner
{
  /* ******************************************************************* *
     Specification of the probability model used for this small example: 
   * ******************************************************************* */
  
  RealVariable observation = real(0.06);
  
  @DefineFactor(onObservations = true) 
  Uniform<MinMaxParameterization> likelihood = Uniform.on(observation);
  
  @DefineFactor 
  Exponential<MeanParameterization> prior = Exponential.on(likelihood.parameters.max).withMean(10.0);

  
  /* ******************************************************************* */ 

  @Tutorial(showSource = false, showLink = true, linkPrefix = "src/test/java/")
  public static void main(String [] args)
  {
    BlangSmallExample example = new BlangSmallExample();
    example.factory.mcmcOptions.nMCMCSweeps = 100000;
    example.run();
    File histogramFile = Results.getFileInResultFolder("posteriorOnNumberOfFutureHumans.pdf");
    PlotHistogram.from(example.samples).toPDF(histogramFile);
  }
  
  /**
   * Blang provides a way of building a factor graph declaratively. 
   * 
   * Recall that factor graphs are bipartite undirected graphs with two types 
   * of nodes: **factors** and **variables**. Each factor computes one factor
   * in a target unnormalized density broken into a product of factors. 
   * Each variable holds one coordinate of a state. 
   * 
   * For now, think of a factor as a prior or conditional probability density 
   * (although Blang can be used in other contexts such as for Markov random fields, 
   * we will focus on factor graphs
   * induced from directed models in this tutorial).
   * 
   * **Factors** are declared 
   * using the ``@DefineFactor`` annotation (an annotation is a mechanism designed 
   * to extend the functionality of the Java language). Blang will 
   * look for fields with a DefineFactor annotations and build model composed of 
   * these factors. 
   * 
   * We will look at the inner working of factors later on,
   * for now we focus on how to use them. For example, in our example, there are two 
   * factors, one for the prior (an exponential density), and one for the likelihood
   * (a uniform density).
   *  
   * **Variables** are fields within factors classes that have a 
   * ``@FactorArgument`` and ``@FactorComponent`` annotations (the differences between 
   * the two will be explained shortly). These fields can be of any type (more on that
   * later), and they typically correspond to either the parameters or the realization 
   * of a probability density.
   * 
   * In order to use existing factors, the only thing to do with variables is
   * to *connect* them across factors. In our example, this is achieved via
   * ``Exponential.on(likelihood.parameters.max)``. Here ``Exponential.on(x)`` is just a 
   * a static function short-hand for creating a unit mean exponential density with 
   * realization x, 
   * ``new Exponential(x, ...)``. The consequence of doing this is that the realization field
   * of the exponential prior is ``==`` (object equality, i.e. two references to the same
   * location in the heap) to the realization field of 
   * the maximum parameter of the uniform density factor. Blang automatically discovers 
   * these relationships and use them to build edges in the factor graph.
   * 
   * #### Building an MCMC sampler
   * 
   * Obtaining an MCMC sampler from a model as specified above is very simple:
   * 
   * 1. Make the class in which you declared the factor extend ``MCMCRunner``.
   * 2. In the main function, create an instance of that class, and call ``.run()``.
   * 
   * If you run the sampler, you will see an histogram created in ``posteriorOnMax.pdf``,
   * showing the approximate posterior over the number of future members of the human 
   * species under our simple model.
   * 
   * #### Differentiating parameters from realisations
   * 
   * As mentioned earlier, there are two types of factor arguments in our setup: 
   * parameters and realisations. As we will see shortly, it 
   * is useful to differentiate the two. This is done via the syntax 
   * ``@FactorArgument(makeStochastic = true)``, which indicates that the argument
   * is a realisation (therefore making the variable stochastic/random).
   * 
   * See ``bayonet.distributions.Exponential.java`` and ``bayonet.distributions.Uniform.java``
   * for examples of this construction.
   * 
   * #### Specifying what should be sampled
   * 
   * All the latent (non-observed) random variables are sampled by default. Blang uses the 
   * following strategy to determine the random variables: the variable should appear as annotated by  
   * ``@FactorArgument(makeStochastic = true)`` in one of the factors (see above).
   * 
   * A random variable can be either latent or observed. To mark a variable as observed,
   * simply use ``@DefineFactor(onObservations = true)`` on the factor in which the 
   * observed variable is declared as stochastic.
   * 
   * #### Processing the samples
   * 
   * The simplest way of process the samples is to override the function ``process()``
   * as follows:
   * 
   * 
   */
  @Tutorial(showSource = true, showSignature = true, showLink = true, linkPrefix = "src/test/java/")
  @Override
  protected void process(ProcessorContext context)
  {
    samples.add(likelihood.parameters.max.getValue() - observation.getValue());
  }
  
  List<Double> samples = Lists.newArrayList();
}
