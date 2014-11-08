package conifer;

/**
 * 
 * @author jewellsean
 *
 */
public class ParsimonyVector
{
    private int[] M; 
    
    public static ParsimonyVector oneInit(int nSites)
    {
        int[] M = new int[nSites];
        for (int i = 0; i < nSites; i++)
            M[i] = 1;
        return new ParsimonyVector(M);
    }
            
    
    
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
