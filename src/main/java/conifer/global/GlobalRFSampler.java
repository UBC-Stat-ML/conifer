package conifer.global;

import java.util.*;

import blang.ProbabilityModel;
import blang.factors.Factor;
import conifer.local.CollisionContext;
import conifer.local.CollisionFactor;
import conifer.local.EventQueue;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.jblas.DoubleMatrix;

import conifer.rejfreeutil.RFSamplerOptions;
import conifer.rejfreeutil.StaticUtils;
import conifer.rejfreeutil.RFSamplerOptions.RefreshmentMethod;

import com.google.common.collect.Lists;

import bayonet.distributions.Exponential;
import bayonet.opt.DifferentiableFunction;
import bayonet.opt.LBFGSMinimizer;

import static conifer.rejfreeutil.StaticUtils.*;

public class GlobalRFSampler {
    public static interface CollisionSolver {
        public double collisionTime(final DoubleMatrix initialPoint, final DoubleMatrix velocity, DifferentiableFunction energy, final double exponential);
    }

    /**
     * *Negative* log density of the target distribution.
     */

    private final EventQueue<CollisionFactor> _collisionQueue = new EventQueue<>();
    private final Map<CollisionFactor, Boolean> isCollisionMap = new LinkedHashMap<>();

    private final DifferentiableFunction energy;
    private final CollisionSolver solver;
    private final RFSamplerOptions options;

    private List<DoubleMatrix> trajectory = Lists.newArrayList();
    private List<DoubleMatrix> samples = Lists.newArrayList();

    private DoubleMatrix currentPosition, currentVelocity;

    private SummaryStatistics collisionToRefreshmentRatio = new SummaryStatistics();
    private SummaryStatistics collectedPerEvent = new SummaryStatistics();

    public final List<CollisionFactor> allFactors = new ArrayList<>();
    public final ProbabilityModel model;


    /**
     * @param energy The negative log density of the target, assumed to be convex
     */
    public GlobalRFSampler(DifferentiableFunction energy, DoubleMatrix initialPosition, RFSamplerOptions options, CollisionSolver solver, ProbabilityModel model) {
        if (options.refreshmentMethod != RefreshmentMethod.GLOBAL && options.refreshmentMethod != RefreshmentMethod.LOCAL)
            throw new RuntimeException();
        this.solver = solver;
        this.energy = energy;
        this.options = options;
        this.currentPosition = initialPosition;
        this.currentVelocity = null;
        this.model = model;
        for(Factor f: model.linearizedFactors())
            allFactors.add((CollisionFactor) f);

    }


    public  GlobalRFSampler(DifferentiableFunction energy, DoubleMatrix initialPosition, RFSamplerOptions options, ProbabilityModel model) {
        this(energy, initialPosition, options, new PegasusConvexCollisionSolver(), model);
    }

    public GlobalRFSampler(DifferentiableFunction energy, DoubleMatrix initialPosition, ProbabilityModel model) {
        this(energy, initialPosition, new RFSamplerOptions(), model);
    }

    public static GlobalRFSampler initializeRFWithLBFGS(DifferentiableFunction energy, RFSamplerOptions options, ProbabilityModel model) {
        return new GlobalRFSampler(energy, optimizePosition(energy), options, model);
    }

    public static GlobalRFSampler initializeRFWithLBFGS(DifferentiableFunction energy, ProbabilityModel model) {
        return initializeRFWithLBFGS(energy, new RFSamplerOptions(), model);
    }

    private static DoubleMatrix optimizePosition(DifferentiableFunction energy) {
        double[] min = new LBFGSMinimizer().minimize(energy, new DoubleMatrix(energy.dimension()).data, 1e-2);
        return new DoubleMatrix(min);
    }

    private DoubleMatrix
            mean,
            variance;

    public DoubleMatrix getMean() {
        return mean;
    }

    public DoubleMatrix getVariance() {
        return variance;
    }

    //Todo: remove the redundancy of this method compared with iterate() method
    public void iterateWithSuperPosition(Random rand, int numberOfIterations, boolean useSuperPosition) {
        if (currentVelocity == null)
            currentVelocity = uniformOnUnitBall(energy.dimension(), rand);

        trajectory.add(currentPosition);
        double totalTime = 0.0;
        mean = new DoubleMatrix(dimensionality());
        variance = new DoubleMatrix(dimensionality(), dimensionality());

        //initialize collision queue first
        initCollisionQueue(rand, 0.0);

        for (int iter = 0; iter < numberOfIterations; iter++) {
            // obtain the collision factor, smallest collision time in the queue and indicator of whether
            // it is a true collision or a virtual collision
            final Map.Entry<Double, CollisionFactor> collision = _collisionQueue.pollEvent();

            //this collisionTime starts from zero and it is not the collision delta time.
            final double collisionTime = collision.getKey();
            final CollisionFactor collisionFactor = collision.getValue();
            final boolean isActualCollision = isCollisionMap.get(collisionFactor);
            double collisionDeltaTime = collisionTime - totalTime;

            // compare collision delta time with refreshTime
            double refreshTime = options.refreshRate == 0 ? Double.POSITIVE_INFINITY : Exponential.generate(rand, options.refreshRate);
            double eventTime = Math.min(collisionDeltaTime, refreshTime);

            totalTime += eventTime;
            collisionToRefreshmentRatio.addValue(eventTime / refreshTime);

            // collect state
            collectSamples(currentPosition, currentVelocity, eventTime, rand);

            // update state
            boolean collisionOccurs = isActualCollision && (eventTime < refreshTime);
            currentPosition = position(currentPosition, currentVelocity, eventTime);
            allFactors.get(0).setPosition(currentPosition);
            trajectory.add(currentPosition);

            //update velocity
            if (collisionOccurs)
                currentVelocity = StaticUtils.bounce(currentVelocity, gradient(currentPosition));
            if(collisionDeltaTime >= refreshTime)
                currentVelocity = refreshVelocity(currentPosition, currentVelocity, rand);

            if(isActualCollision){
                // clear _collisionQueue and collision map
                _collisionQueue.clear();
                isCollisionMap.clear();

                // initialize collision queue again using the updated velocity and updated position
                // initialize collision map
                initCollisionQueue(rand, totalTime);
            }else{
                updateCandidateCollision(rand, collisionFactor, totalTime);
            }

        }
        mean.divi(totalTime);
        variance.divi(totalTime);
    }

    public void iterate(Random rand, int numberOfIterations, boolean useSuperPosition, ProbabilityModel model){
        if(!useSuperPosition){

            iterate(rand, numberOfIterations);

        }else{

            iterateWithSuperPosition(rand, numberOfIterations, true);
        }
    }

    public void iterate(Random rand, int numberOfIterations) {
        if (currentVelocity == null)
            currentVelocity = uniformOnUnitBall(energy.dimension(), rand);

        trajectory.add(currentPosition);
        double totalTime = 0.0;
        mean = new DoubleMatrix(dimensionality());
        variance = new DoubleMatrix(dimensionality(), dimensionality());
        for (int iter = 0; iter < numberOfIterations; iter++) {
            // simulate event
            final double exponential = StaticUtils.generateUnitRateExponential(rand);
            double collisionTime = solver.collisionTime(currentPosition, currentVelocity, energy, exponential);
            double refreshTime = options.refreshRate == 0 ? Double.POSITIVE_INFINITY : Exponential.generate(rand, options.refreshRate);
            double eventTime = Math.min(collisionTime, refreshTime);
            totalTime += eventTime;
            collisionToRefreshmentRatio.addValue(collisionTime / refreshTime);

            // collect state
            collectSamples(currentPosition, currentVelocity, eventTime, rand);

            // update state
            boolean collisionOccurs = collisionTime < refreshTime;
            currentPosition = position(currentPosition, currentVelocity, eventTime);
            trajectory.add(currentPosition);
            if (collisionOccurs)
                currentVelocity = StaticUtils.bounce(currentVelocity, gradient(currentPosition));
            else
                currentVelocity = refreshVelocity(currentPosition, currentVelocity, rand);
        }
        mean.divi(totalTime);
        variance.divi(totalTime);
    }

    public void setVelocity(DoubleMatrix velocity) {
        this.currentVelocity = velocity.dup();
    }

    private DoubleMatrix refreshVelocity(DoubleMatrix currentPosition,
                                         DoubleMatrix currentVelocity, Random rand) {
        return uniformOnUnitBall(currentVelocity.length, rand);
    }

    public void initCollisionQueue(Random rand, double currentTime){
        CollisionContext context = new CollisionContext(rand, currentVelocity);
        double candidateCollisionTime = 0;
        for(int i=0; i< allFactors.size();i++){
            Pair<Double, Boolean> collisionInfo = allFactors.get(i).getLowerBoundForCollisionDeltaTime(context);
            candidateCollisionTime = currentTime + collisionInfo.getLeft();
            isCollisionMap.put(allFactors.get(i), collisionInfo.getRight());
             _collisionQueue.add(allFactors.get(i), candidateCollisionTime);
        }
    }




    private static final double epsilon = 1;

    private void updateCandidateCollision(Random rand, CollisionFactor factor, double currentTime) {
        _collisionQueue.remove(factor);

        CollisionContext context = new CollisionContext(rand, currentVelocity);
        Pair<Double, Boolean> collisionInfo = factor.getLowerBoundForCollisionDeltaTime(context);

        double candidateCollisionTime = currentTime + collisionInfo.getLeft();
        isCollisionMap.put(factor, collisionInfo.getRight());

        if (_collisionQueue.containsTime(candidateCollisionTime)) {
            System.err.println("The sampler has hit an event of probability zero: two collisions scheduled exactly at the same time.");
            System.err.println("Because of numerical precision, this could possibly happen, but very rarely.");

            System.err.println("For internal implementation reasons, one of the collisions at time " + candidateCollisionTime + " was moved to " + (candidateCollisionTime + epsilon));
            candidateCollisionTime += epsilon;
        }

        _collisionQueue.add(factor, candidateCollisionTime);
    }

    private void collectSamples(DoubleMatrix initialPosition,
                                DoubleMatrix velocity, double eventTime, Random rand) {
        mean.addi(initialPosition.mul(eventTime)
                .addi(velocity.mul(eventTime * eventTime / 2.0)));
        variance
                .addi(initialPosition.mmul(initialPosition.transpose()).mul(eventTime))
                .addi(initialPosition.mmul(velocity.transpose()).mul(eventTime * eventTime))
                .addi(velocity.mmul(velocity.transpose()).mul(eventTime * eventTime * eventTime / 3.0));

        if (options.collectRate == 0.0)
            return;
        double timeConsumed = Exponential.generate(rand, options.collectRate);
        int nCollected = 0;
        while (timeConsumed < eventTime) {
            nCollected++;
            samples.add(position(initialPosition, velocity, timeConsumed));
            timeConsumed += Exponential.generate(rand, options.collectRate);
        }
        collectedPerEvent.addValue(nCollected);
    }

    private DoubleMatrix gradient(DoubleMatrix position) {
        return new DoubleMatrix(energy.derivativeAt(position.data));
    }

    public List<DoubleMatrix> getTrajectory() {
        return trajectory;
    }

    public SummaryStatistics getCollisionToRefreshmentRatio() {
        return collisionToRefreshmentRatio;
    }

    public List<DoubleMatrix> getSamples() {
        return samples;
    }

    public SummaryStatistics getCollectedPerEvent() {
        return collectedPerEvent;
    }

    public int dimensionality() {
        return energy.dimension();
    }

    public DoubleMatrix getCurrentPosition() {
        return currentPosition;
    }

    public void setCurrentPosition(DoubleMatrix currentPosition) {
        this.currentPosition = currentPosition;
    }
}
