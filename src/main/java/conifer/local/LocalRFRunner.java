package conifer.local;

import java.util.Collection;
import java.util.List;
import java.util.concurrent.TimeUnit;

import com.google.common.base.Stopwatch;

import conifer.rejfreeprocessors.MomentRayProcessor;
import conifer.rejfreeprocessors.SaveRaysProcessor;
import conifer.rejfreeprocessors.SaveSamplesProcessor;
import blang.ProbabilityModel;
import blang.variables.RealVariable;
import briefj.OutputManager;
import briefj.opt.OptionSet;
import briefj.run.Results;


public class LocalRFRunner {
    @OptionSet(name = "runnerOptions")
    public final LocalRFRunnerOptions options;

    public LocalRFRunner(LocalRFRunnerOptions options) {
        this.options = options;
    }

    public LocalRFRunner() {
        this.options = new LocalRFRunnerOptions();
    }

    public ProbabilityModel model;
    public LocalRFSampler sampler;

    public OutputManager output = Results.getGlobalOutputManager();

    public MomentRayProcessor momentRayProcessor = null;
    public SaveRaysProcessor saveRaysProcessor = null;
    public SaveSamplesProcessor saveSamplesProcessor = null;

    public void init(Object modelSpec) {
        if (isInit())
            throw new RuntimeException("Model already initialized.");
        model = new ProbabilityModel(modelSpec);
        sampler = new LocalRFSampler(model, options.rfOptions);
    }

    private boolean isInit() {
        return model != null;
    }

    private void checkInit() {
        if (!isInit())
            throw new RuntimeException("Model should first be initialized.");
    }

    public void addMomentRayProcessor() {
        checkInit();
        if (momentRayProcessor != null)
            throw new RuntimeException("Already added");
        momentRayProcessor = new MomentRayProcessor();
        sampler.addRayProcessor(momentRayProcessor);
    }

    public void addSaveRaysProcessor(Collection<RealVariable> variables) {
        checkInit();
        if (saveRaysProcessor != null)
            throw new RuntimeException("Already added");
        saveRaysProcessor = new SaveRaysProcessor(variables);
        sampler.addRayProcessor(saveRaysProcessor);
    }

    public void addSaveAllRaysProcessor() {
        checkInit();
        List<RealVariable> all = model.getLatentVariables(RealVariable.class);
        addSaveRaysProcessor(all);
    }

    public void addSaveSamplesProcessor(List<RealVariable> variables) {
        checkInit();
        if (saveSamplesProcessor != null)
            throw new RuntimeException("Already added");
        saveSamplesProcessor = new SaveSamplesProcessor(variables);
        sampler.addPointProcessor(saveSamplesProcessor);
    }

    public Stopwatch watch = null;

    public void run() {
        checkInit();

        watch = Stopwatch.createStarted();
        sampler.iterate(options.samplingRandom, options.maxSteps, options.maxTrajectoryLength, options.maxRunningTimeMilli);
        watch.stop();
        long elapsed = watch.elapsed(TimeUnit.MILLISECONDS);

        if (!options.silent)
            output.printWrite("general-sampler-diagnostic",
                    "wallClockTimeMilli", elapsed,
                    "trajectoryLength", sampler.getTrajectoryLength(),
                    "nCollisions", sampler.getNCollisions(),
                    "nCollidedVariables", sampler.getNCollidedVariables(),
                    "nRefreshments", sampler.getNRefreshments(),
                    "nRefreshedVariables", sampler.getNRefreshedVariables());

        output.flush();
    }


}
