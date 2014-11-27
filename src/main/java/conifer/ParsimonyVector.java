package conifer;

import java.util.Map;

/**
 * 
 * @author jewellsean
 *
 */
public class ParsimonyVector
{
    private int[] M; 
    public Map<Integer, String> leafString;
    
    public static ParsimonyVector oneInit(int nSites, Map<Integer, String> leafString)
    {
        int[] M = new int[nSites];
        for (int i = 0; i < nSites; i++)
            M[i] = 1;
        return new ParsimonyVector(M, leafString);
    }
            
    
    
    public ParsimonyVector(int[] M, Map<Integer, String> leafString)
    {
        this.M = M; 
        this.leafString = leafString;
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
    
    public String getLeafString(Integer leaf)
    {
        return leafString.get(leaf);
    }
    
}
