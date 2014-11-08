package conifer.models;

import conifer.Parsimony;
import blang.annotations.FactorArgument;
import blang.annotations.FactorComponent;

public class ParsimonyModel
{

    @FactorArgument public final double betaBinomialPrecision; 
    @FactorComponent public final Parsimony parsimony;  
    
    public ParsimonyModel(double betaBinomialPrecision, Parsimony parsimony)
    {
        this.betaBinomialPrecision = betaBinomialPrecision;
        this.parsimony = parsimony; 
    }
    
    
}
