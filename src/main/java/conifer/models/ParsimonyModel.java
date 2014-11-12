package conifer.models;

import static blang.variables.RealVariable.real;
import blang.annotations.FactorArgument;
import blang.factors.Factor;
import blang.variables.RealVariable;
import conifer.Parsimony;
import conifer.ParsimonyVector;

/**
 * 
 * @author jewellsean
 *
 */
public class ParsimonyModel implements Factor
{
    @FactorArgument 
    public final RealVariable betaBinomialPrecision; 
    
    @FactorArgument(makeStochastic = true)
    public final Parsimony parsimony;  
    
    public ParsimonyModel(RealVariable betaBinomialPrecision, Parsimony parsimony)
    {
        this.betaBinomialPrecision = betaBinomialPrecision;
        this.parsimony = parsimony; 
    }
    
    public static ParsimonyModel on(Parsimony parsimony)
    {
        RealVariable betaBinomial = real(1); 
        return new ParsimonyModel(betaBinomial, parsimony);
    }

    
    
    @Override
    public double logDensity()
    {
        return 0;
    }
    
    
}
