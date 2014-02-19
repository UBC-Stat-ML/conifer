package bayonet.mcmc;

import java.util.List;
import java.util.Random;




public interface Move
{
  public void execute(Random rand);
  @SuppressWarnings("rawtypes")
  public List variablesCovered();
}
