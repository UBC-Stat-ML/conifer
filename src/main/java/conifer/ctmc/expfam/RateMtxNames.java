package conifer.ctmc.expfam;

import static conifer.ctmc.expfam.SerializedExpFamMixture.*;
import static conifer.io.Indexers.*;
import static conifer.io.PhylogeneticObservationFactory.*;
import static conifer.ctmc.RateMatrices.*;

import java.util.Random;

import conifer.ctmc.RateMatrices;
import conifer.ctmc.SimpleRateMatrix;
import conifer.io.PhylogeneticObservationFactory;
import briefj.Indexer;


public enum RateMtxNames
{

    KIMURA1980 {
      @Override
      public SerializedExpFamMixture getSerialized()
      {
        return SerializedExpFamMixture.kimura1980();
      }
      
      public Indexer<String> getIndexer()
      {
        return dnaIndexer();
      }
      
      public PhylogeneticObservationFactory getFactory()
      {
        return nucleotidesFactory();
      }
      
      public SimpleRateMatrix getRateMtx()
      {
        return RateMatrices.kimura1980();
      }
      
      
    },
    ACCORDANCE {
      @Override
      public SerializedExpFamMixture getSerialized()
      {
        return SerializedExpFamMixture.accordance();
      }
      
      public Indexer<String> getIndexer()
      {
        return proteinIndexer();
      }
      
      public PhylogeneticObservationFactory getFactory()
      {
        return proteinFactory();
      }
      
      public SimpleRateMatrix getRateMtx()
      {
        return RateMatrices.accordance();
      }
    },
    PAIR {
      @Override
      public SerializedExpFamMixture getSerialized()
      {
        return pair();
      }
      
      public Indexer<String> getIndexer()
      {
        return proteinPairIndexer();
      }
      
      public PhylogeneticObservationFactory getFactory()
      {
        return proteinPairFactory();
      }
      public SimpleRateMatrix getRateMtx()
      {
        Random rand = new Random();
        return RateMatrices.randomGTR(rand, 400);
      }
    },
    POLARITY {
      @Override
      public SerializedExpFamMixture getSerialized()
      {
        return SerializedExpFamMixture.polarity();
      }
      
      public Indexer<String> getIndexer()
      {
        return proteinIndexer();
      }
      
      public PhylogeneticObservationFactory getFactory()
      {
        return proteinFactory();
      }
      
      public SimpleRateMatrix getRateMtx()
      {
        return RateMatrices.polarity();
      }
      
      
    },
    POLARITYSIZE {
      @Override
      public SerializedExpFamMixture getSerialized()
      {
        return SerializedExpFamMixture.polaritySize();
      }
      
      public Indexer<String> getIndexer()
      {
        return proteinIndexer();
      }
      
      public PhylogeneticObservationFactory getFactory()
      {
        return proteinFactory();
      }
      
      public SimpleRateMatrix getRateMtx()
      {
        return RateMatrices.polaritySize();
      }
    };
    
    public abstract  SerializedExpFamMixture getSerialized();


    public abstract Indexer<String> getIndexer();


    public abstract PhylogeneticObservationFactory getFactory();


    public abstract SimpleRateMatrix getRateMtx();
    
}
