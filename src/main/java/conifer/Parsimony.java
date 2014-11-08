package conifer;

import blang.annotations.Processors;
import blang.annotations.Samplers;
import conifer.moves.ParsimonySampler;
import conifer.processors.ParsimonyProcessor;

@Samplers({ParsimonySampler.class})
@Processors({ParsimonyProcessor.class})

/**
 * 
 * @author jewellsean
 *
 */
public class Parsimony
{

/*
 * Store the min-value for each site in a vector
 * this is the current state of affairs 
 * 
 */
    private ParsimonyVector M; 
        
    public Parsimony(ParsimonyVector M)
    {
        this.M = M; 
    }
  
   
    public ParsimonyVector getM()
    {
        return M; 
    }

}
