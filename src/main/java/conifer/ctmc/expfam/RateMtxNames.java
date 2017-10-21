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

    DNAGTR{

      @Override
      public SerializedExpFamMixture getSerialized()
      {
        return SerializedExpFamMixture.dnaGTR();
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
        Random rand = new Random(2);
        return RateMatrices.randomGTR(rand, 4);
      }

    },


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
    },
    POLARITYSIZEGTR {
      @Override
      public SerializedExpFamMixture getSerialized()
      {
        return SerializedExpFamMixture.polaritySizeGTR();
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
        return RateMatrices.polaritySizeGTR();
      }
    },
    PROTEINSIMPLEGTR{
      @Override
      public SerializedExpFamMixture getSerialized()
      {
        return SerializedExpFamMixture.proteinSimpleGTR();
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
        Random rand = new Random(1);
        return RateMatrices.randomGTR(rand, 20);
      }   
      
    };
    
    public abstract  SerializedExpFamMixture getSerialized();


    public abstract Indexer<String> getIndexer();


    public abstract PhylogeneticObservationFactory getFactory();


    public abstract SimpleRateMatrix getRateMtx();
    
}
