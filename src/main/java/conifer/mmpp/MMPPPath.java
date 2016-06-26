package conifer.mmpp;

import com.google.common.collect.Lists;

import java.util.List;

/**
 * Created by crystal on 2016-04-26.
 */
public class MMPPPath
{
    private final List<Double> times = Lists.newArrayList();
    private final List<Integer> states = Lists.newArrayList();
    private final List<Integer> poissonStates = Lists.newArrayList();
    private final List<Double> poissonTimes = Lists.newArrayList();


    public void addSegment(int currentPoint, double time)
    {
        if (times.size() != states.size())
            throw new RuntimeException();
        if (!isEmpty() && lastState() == currentPoint)
            times.set(times.size()-1, times.get(states.size()-1) + time);
        else
        {
            times.add(time);
            states.add(currentPoint);
        }
    }

    public List<Integer> getPoissonEventState(){ return poissonStates;}
    public List<Double> getPoissonTimes(){ return poissonTimes;}

    public int nSojourns()
    {
        return times.size();
    }

    public boolean isEmpty()
    {
        return nSojourns() == 0;
    }

    public int lastState()
    {
        return states.get(states.size() - 1);
    }

    public int firstState()
    {
        return states.get(0);
    }

    public double totalLength()
    {
        double sum = 0.0;
        for (double t : times)
            sum += t;
        return sum;
    }

    public void addPoissonEventTimesOfMMPP(double eventTime)
    {
        poissonTimes.add(eventTime);
    }

    public void addPoissonEventStatesOfMMPP(int eventState){

        poissonStates.add(eventState);
    }
}