package conifer;

import blang.annotations.Samplers;
import conifer.moves.ParsimonySampler;

@Samplers({ParsimonySampler.class})
//@Processors({})
public class Parsimony
{

/*
 * Store the min-value for each site in a vector
 * this is the current state of affairs 
 * 
 */
    private int[] M; 

    public Parsimony(int[] M)
    {
        this.M = M; 
    }
    

}
