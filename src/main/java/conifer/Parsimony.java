package conifer;

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
