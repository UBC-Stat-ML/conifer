package bayonet.mcmc;

import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;




public class Sampler
{
  private final List<Move> moves;
  private final PrModel model;
  
  @SuppressWarnings({ "rawtypes", "unchecked" })
  private List<Move> init(List<MoveFactory> factories)
  {
    if (model.getObservedVariables().isEmpty())
      throw new RuntimeException("At least one node needs to be set to observed. Use PrModel.setObserved().");
    
    List<Move> result = Lists.newArrayList();
    Set coveredNode = Sets.newHashSet();
    
    for (MoveFactory factory : factories)
    {
      List<Move> moves = factory.build(model);
      for (Move move : moves)
        result.addAll(move.variablesCovered());
      coveredNode.addAll(moves);
    }
    
    if (!coveredNode.containsAll(model.getLatentStochasticVariables()))
      throw new RuntimeException();
    
    return result;
  }
  
  private void sweep(Random rand)
  {
    Collections.shuffle(moves, rand);
    for (Move move : moves)
      move.execute(rand);
    
    // TODO: some sort of logging/checkpointing
  }
}
