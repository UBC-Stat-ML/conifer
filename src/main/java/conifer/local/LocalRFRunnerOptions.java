package conifer.local;

import java.util.Random;

import conifer.rejfreeutil.RFSamplerOptions;
import briefj.opt.Option;
import briefj.opt.OptionSet;


public class LocalRFRunnerOptions {
    @Option
    public long maxRunningTimeMilli = Long.MAX_VALUE;

    @Option
    public int maxSteps = 1000;

    @Option
    public double maxTrajectoryLength = Double.POSITIVE_INFINITY;

    @Option
    public Random samplingRandom = new Random(1);

    @OptionSet(name = "rfOptions")
    public RFSamplerOptions rfOptions = new RFSamplerOptions();

    @Option
    public boolean silent = false;
}