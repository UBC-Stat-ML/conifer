package conifer.models;

import conifer.Parsimony;
import blang.annotations.FactorArgument;
import blang.annotations.FactorComponent;

public class ParsimonyModel
{

    @FactorArgument public final double precision; 
    @FactorComponent public final Parsimony parsimony;  
    
    public ParsimonyModel(double precision, Parsimony parsimony)
    {
        this.precision = precision;
        this.parsimony = parsimony; 
    }
    
    
}
