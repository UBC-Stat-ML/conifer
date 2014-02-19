package bayonet.mcmc;

import java.util.List;



public interface MoveFactory
{
  public List<Move> build(PrModel model);
}
