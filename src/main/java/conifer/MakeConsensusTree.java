package conifer;

import java.io.File;

import briefj.opt.Option;
import briefj.run.Mains;

public class MakeConsensusTree implements Runnable
{
    @Option(required = true) 
    public String treeFilePath;
    @Option(required = true)
    public String outputFile;
    
    public static void main(String args[])
    {
        Mains.instrumentedRun(args, new MakeConsensusTree());
    }

    @Override
    public void run()
    {
        MajorityRuleTree.buildAndWriteConsensusTree(new File(treeFilePath),
                new File(outputFile));
        
    }
    
}
