package conifer.global;

import java.util.*;
import java.util.concurrent.TimeUnit;

import blang.ProbabilityModel;
import blang.factors.Factor;
import com.google.common.base.Stopwatch;
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

    private int nBounces = 0;
    private int nRefreshments = 0;
    private double totalTime = 0;
    private double averageTime = 0;
    private int nVirtualCollisions =0;

    private long collisionCalculationTime = 0;
    private long candidateCollisionCalculationTime = 0;
    private long gradientCalculationTime =0;
    private long refreshmentTime = 0;

    public Stopwatch watch = null;

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
        this(energy, initialPosition, options, new MultipleConvexCollisionSolver(SolverNames.Pegasus), model);
    }

    public  GlobalRFSampler(DifferentiableFunction energy, DoubleMatrix initialPosition, RFSamplerOptions options, ProbabilityModel model, SolverNames solverNames) {
        this(energy, initialPosition, options, new MultipleConvexCollisionSolver(solverNames), model);
    }

    public GlobalRFSampler(DifferentiableFunction energy, DoubleMatrix initialPosition, ProbabilityModel model) {
        this(energy, initialPosition, new RFSamplerOptions(), model);
    }

    public static GlobalRFSampler initializeRFWithLBFGS(DifferentiableFunction energy, RFSamplerOptions options, ProbabilityModel model) {
        return new GlobalRFSampler(energy, optimizePosition(energy), options, model);
    }

    public static GlobalRFSampler initializeRFWithLBFGS(DifferentiableFunction energy, RFSamplerOptions options, ProbabilityModel model, SolverNames solverNames) {
        return new GlobalRFSampler(energy, optimizePosition(energy), options, new MultipleConvexCollisionSolver(solverNames), model);
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
        totalTime = 0.0;
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
            if (collisionOccurs){
                nBounces ++;
                //watch = Stopwatch.createStarted();
                currentVelocity = StaticUtils.bounce(currentVelocity, gradient(currentPosition));
                //watch.stop();
                //gradientCalculationTime = gradientCalculationTime + watch.elapsed(TimeUnit.NANOSECONDS);

            }

            if(collisionDeltaTime >= refreshTime){
                currentVelocity = refreshVelocity(currentPosition, currentVelocity, rand);
                nRefreshments++;

            }

            if(isActualCollision){
                // clear _collisionQueue and collision map
                _collisionQueue.clear();
                isCollisionMap.clear();

                //watch = Stopwatch.createStarted();
                // initialize collision queue again using the updated velocity and updated position
                // initialize collision map
                initCollisionQueue(rand, totalTime);
                //watch.stop();
                collisionCalculationTime = collisionCalculationTime + watch.elapsed(TimeUnit.NANOSECONDS);
            }else{
                //watch = Stopwatch.createStarted();
                updateCandidateCollision(rand, collisionFactor, totalTime);
                //watch.stop();
                nVirtualCollisions++;
                //candidateCollisionCalculationTime = candidateCollisionCalculationTime + watch.elapsed(TimeUnit.NANOSECONDS);
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
            //watch = Stopwatch.createStarted();
            double collisionTime = solver.collisionTime(currentPosition, currentVelocity, energy, exponential);
            //watch.stop();
            //collisionCalculationTime = collisionCalculationTime + watch.elapsed(TimeUnit.NANOSECONDS);
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
            if (collisionOccurs){
                //watch = Stopwatch.createStarted();
                currentVelocity = StaticUtils.bounce(currentVelocity, gradient(currentPosition));
                //watch.stop();
                //gradientCalculationTime = gradientCalculationTime + watch.elapsed(TimeUnit.NANOSECONDS);
                nBounces++;
            }

            else{
                //watch = Stopwatch.createStarted();
                currentVelocity = refreshVelocity(currentPosition, currentVelocity, rand);
                //watch.stop();
                //refreshmentTime = refreshmentTime + watch.elapsed(TimeUnit.NANOSECONDS);
                nRefreshments++;
            }

        }
        mean.divi(totalTime);
        variance.divi(totalTime);
        averageTime = totalTime/numberOfIterations;
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

    public int getNBounces(){

        return nBounces;
    }

    public int getNRefreshments(){

        return nRefreshments;
    }

    public double getTotalTimeWithSuperposition(){

        return totalTime;

    }

    public double getAverageTimeWithoutSuperposition(){
        return averageTime;
    }

    public long getCollisionCalculationTime(){
        return collisionCalculationTime;
    }

    public long getCandidateCollisionCalculationTime(){
        return candidateCollisionCalculationTime;
    }

    public long getGradientCalculationTime(){

        return gradientCalculationTime;
    }

    public int getNVirtualCollisions(){
        return nVirtualCollisions;
    }

    public long getRefreshmentTime(){
        return refreshmentTime;
    }

    public void setCurrentPosition(DoubleMatrix currentPosition) {
        this.currentPosition = currentPosition;
    }
}
