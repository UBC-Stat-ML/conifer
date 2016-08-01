package conifer.rejfreeprocessors;

import conifer.local.LocalRFSampler;
import conifer.local.TrajectoryRay;
import blang.variables.RealVariable;


public interface RayProcessor {
    public void init(LocalRFSampler sampler);

    /**
     * Process a ray that bounced (closed) at the given time.
     *
     * @param var
     * @param ray
     * @param time
     * @param sampler
     */
    public void processRay(RealVariable var, TrajectoryRay ray, double time, LocalRFSampler sampler);
}