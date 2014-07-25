package conifer;

import java.util.Random;

import conifer.ctmc.CTMC;
import conifer.moves.SingleBranchScaling;
import conifer.moves.SingleNNI;

import bayonet.distributions.Exponential;
import bayonet.distributions.Uniform;
import bayonet.distributions.Exponential.MeanParameterization;
import bayonet.distributions.Uniform.MinMaxParameterization;
import blang.annotations.DefineFactor;
import blang.annotations.Samplers;
import tutorialj.Tutorial;



public class ConiferTutorial
{
  /**
   * # conifer
   * 
   * Installing from source
   * ----------------------
   * 
   * Requires: gradle, git, eclipse
   * 
   * - Clone the repository
   * - Type ``gradle eclipse`` from the root of the repository
   * - From eclipse:
   *   - ``Import`` in ``File`` menu
   *   - ``Import existing projects into workspace``
   *   - Select the root of the newly created repo
   *   - Deselect ``Copy projects into workspace`` to avoid having duplicates
   * 
   */
  @Tutorial(startTutorial = "README.md", showSource = false)
  public void usage() {}
  
  

}
