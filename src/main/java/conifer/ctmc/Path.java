package conifer.ctmc;

import java.util.List;

import com.google.common.collect.Lists;



public class Path
{
  private final List<Double> times = Lists.newArrayList();
  private final List<Integer> states = Lists.newArrayList();

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
  
}