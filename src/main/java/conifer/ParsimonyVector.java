package conifer;

/**
 * 
 * @author jewellsean
 *
 */
public class ParsimonyVector
{
    private int[] M; 
    
    
    public ParsimonyVector(int[] M)
    {
        this.M = M; 
    }
    
    public int getIndex(int i)
    {
        return M[i];
    }
    
    public void setIndex(int i, int val)
    {
        M[i] = val;
    }
    
    public int[] getParsimony()
    {
        return M; 
    }
    
    public void setParsimony(int[] M)
    {
        this.M = M; 
    }
}
