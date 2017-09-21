package conifer.RandomUtils;

import java.util.List;
import java.util.Random;

/**
 * Created by crystal on 2017-09-21.
 */



public class DiscreteUniform
{
    public static <S> S sample(List<S> items, Random rand)
    {
        return items.get(rand.nextInt(items.size()));
    }
}
