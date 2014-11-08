package conifer.models;

import blang.annotations.FactorArgument;
import blang.variables.RealVariable;
import static blang.variables.RealVariable.real;
import conifer.Parsimony;
import conifer.ParsimonyVector;

/**
 * 
 * @author jewellsean
 *
 */
public class ParsimonyModel
{
    @FactorArgument 
    public final RealVariable betaBinomialPrecision; 
    
    @FactorArgument(makeStochastic=true)
    public final Parsimony parsimony;  
    
    public ParsimonyModel(RealVariable betaBinomialPrecision, Parsimony parsimony)
    {
        this.betaBinomialPrecision = betaBinomialPrecision;
        this.parsimony = parsimony; 
    }
    
    public static ParsimonyModel initalization(int nSites)
    {
        RealVariable betaBinomial = real(1); 
        Parsimony parsimony = new Parsimony(ParsimonyVector.oneInit(nSites));
        return new ParsimonyModel(betaBinomial, parsimony);
    }
    
    
}
