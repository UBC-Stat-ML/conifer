package conifer;

import java.util.Collection;
import java.util.List;
import java.util.Random;

import bayonet.math.SpecialFunctions;
import blang.core.RealConstant;
import conifer.ctmc.RateMatrices;
import conifer.ctmc.SimpleRateMatrix;
import conifer.factors.NonClockTreePriorUtils;
import conifer.models.DiscreteGammaMixture;
import conifer.models.MultiCategorySubstitutionModel;

public class Utils 
{
  public static UnrootedTree unrootedTree(
      Collection<TreeNode> leaves)
  {
    return NonClockTreePriorUtils.unrootedTree(leaves);
  }
  
  @SuppressWarnings({ "rawtypes", "unchecked" })
  public static MultiCategorySubstitutionModel kimura(int nSites)
  { 
    return new MultiCategorySubstitutionModel(
      new DiscreteGammaMixture(
        new RealConstant(0.0), 
        new RealConstant(1.0), 
        SimpleRateMatrix.fromResource("/conifer/ctmc/kimura1980.txt"), 
        1
      ), 
      nSites
    );
  }
  
  
  @SuppressWarnings({ "rawtypes", "unchecked" })
  public static MultiCategorySubstitutionModel codon(int nSites)
  { 
	Random rand = new Random();  
    return new MultiCategorySubstitutionModel(
      new DiscreteGammaMixture(
        new RealConstant(0.0), 
        new RealConstant(1.0), 
        RateMatrices.randomGTR(rand, 61),
        1
      ), 
      nSites
    );
  }
  
  public static double sampleUnitRateExponential(Random random)
  {
    return - Math.log(random.nextDouble());
  }
  
  public static double sampleExponential(Random random, double rate)
  {
    return sampleUnitRateExponential(random) / rate;
  }
  
  public static <S> S sample(List<S> items, Random rand)
  {
      return items.get(rand.nextInt(items.size()));
  }
  
  public static double gammaQuantile(double y, double shape, double scale) 
  {
    return 0.5 * scale * pointChi2(y, 2.0 * shape);
  }

  private static double pointChi2(double prob, double v) 
  {
      // Returns z so that Prob{x<z}=prob where x is Chi2 distributed with df
      // = v
      // RATNEST FORTRAN by
      // Best DJ & Roberts DE (1975) The percentage points of the
      // Chi2 distribution. Applied Statistics 24: 385-388. (AS91)
  
      final double e = 0.5e-6, aa = 0.6931471805, p = prob;
      double ch, a, q, p1, p2, t, x, b, s1, s2, s3, s4, s5, s6;
      double epsi = .01;
      if (p < 0.000002 || p > 1 - 0.000002) {
          epsi = .000001;
      }
      // if (p < 0.000002 || p > 0.999998 || v <= 0) {
      //      throw new IllegalArgumentException("Arguments out of range p" + p + " v " + v);
      //  }
      double g = SpecialFunctions.lnGamma(v / 2);
      double xx = v / 2;
      double c = xx - 1;
      if (v < -1.24 * Math.log(p)) {
          ch = Math.pow((p * xx * Math.exp(g + xx * aa)), 1 / xx);
          if (ch - e < 0) {
              return ch;
          }
      } else {
          if (v > 0.32) {
              x = normalQuantile(p, 0, 1);
              p1 = 0.222222 / v;
              ch = v * Math.pow((x * Math.sqrt(p1) + 1 - p1), 3.0);
              if (ch > 2.2 * v + 6) {
                  ch = -2 * (Math.log(1 - p) - c * Math.log(.5 * ch) + g);
              }
          } else {
              ch = 0.4;
              a = Math.log(1 - p);
  
  
              do {
                  q = ch;
                  p1 = 1 + ch * (4.67 + ch);
                  p2 = ch * (6.73 + ch * (6.66 + ch));
                  t = -0.5 + (4.67 + 2 * ch) / p1
                          - (6.73 + ch * (13.32 + 3 * ch)) / p2;
                  ch -= (1 - Math.exp(a + g + .5 * ch + c * aa) * p2 / p1)
                          / t;
              } while (Math.abs(q / ch - 1) - epsi > 0);
          }
      }
      do {
          q = ch;
          p1 = 0.5 * ch;
          if ((t = SpecialFunctions.incompleteGammaP(xx, p1, g)) < 0) {
              throw new IllegalArgumentException(
                      "Arguments out of range: t < 0");
          }
          p2 = p - t;
          t = p2 * Math.exp(xx * aa + g + p1 - c * Math.log(ch));
          b = t / ch;
          a = 0.5 * t - b * c;
  
          s1 = (210 + a * (140 + a * (105 + a * (84 + a * (70 + 60 * a))))) / 420;
          s2 = (420 + a * (735 + a * (966 + a * (1141 + 1278 * a)))) / 2520;
          s3 = (210 + a * (462 + a * (707 + 932 * a))) / 2520;
          s4 = (252 + a * (672 + 1182 * a) + c * (294 + a * (889 + 1740 * a))) / 5040;
          s5 = (84 + 264 * a + c * (175 + 606 * a)) / 2520;
          s6 = (120 + c * (346 + 127 * c)) / 5040;
          ch += t
                  * (1 + 0.5 * t * s1 - b
                  * c
                  * (s1 - b
                  * (s2 - b
                  * (s3 - b
                  * (s4 - b * (s5 - b * s6))))));
      } while (Math.abs(q / ch - 1) > e);
  
      return (ch);
  }
  
  public static double normalQuantile(double z, double m, double sd)
  {
    return m + Math.sqrt(2.0) * sd * SpecialFunctions.inverseErf(2.0 * z - 1.0);
  }

}
