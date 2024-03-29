package conifer.ctmc.expfam;

import java.io.File;
import java.util.List;
import java.util.Map;

import org.jgrapht.UndirectedGraph;

import bayonet.graphs.GraphUtils;
import blang.inits.DesignatedConstructor;
import blang.inits.Input;
import briefj.BriefIO;
import briefj.Indexer;
import briefj.collections.Counter;
import briefj.collections.UnorderedPair;

import com.google.common.collect.Maps;
import com.google.gson.Gson;

import conifer.ctmc.RateMatrixToEmissionModel;
import conifer.io.FeatureFactory;




public class SerializedExpFamMixture
{
    public int nCategories;
    public List<String> orderedLatents;
    public List<List<String>> supportEdges;
    public List<UnaryFeature> unaryFeatures;
    public List<BinaryFeature> binaryFeatures;
    public boolean fullSupport;

    public static final Object partition = null;

    public static void main(String [] args)
    {
        SerializedExpFamMixture s = fromResource("/conifer/ctmc/expfam/dnaGTR-expfam.txt");
        System.out.println(s);
        System.out.println(s.getCTMCStateIndexer());
        System.out.println(s.getSupport());
    }
    
    
    public SerializedExpFamMixture(int nCategories, List<String> orderedLatents, List<List<String>> supportEdges, List<UnaryFeature> unaryFeatures, List<BinaryFeature> binaryFeatures,  boolean fullSupport){
    	this.nCategories = nCategories;
    	this.orderedLatents = orderedLatents;
    	this.supportEdges = supportEdges;
    	this.unaryFeatures = unaryFeatures;
    	this.binaryFeatures = binaryFeatures;
    	this.fullSupport = fullSupport;
    	}
  
 
    public BivariateFeatureExtractor<CTMCState> getBivariateFeatureExtractor()
    {
        final Map<UnorderedPair<CTMCState, CTMCState>, Map<String, Double>> binaryFeaturesMap = binaryFeaturesMap(binaryFeatures);
        return new BivariateFeatureExtractor<CTMCState>() {

            @SuppressWarnings("unchecked")
            @Override
            public void extract(@SuppressWarnings("rawtypes") Counter counts, CTMCState state1, CTMCState state2)
            {
                Map<String,Double> feats = binaryFeaturesMap.get(UnorderedPair.of(state1, state2));
                if (feats != null)
                    for (String key : feats.keySet())
                    {
                        if (counts.containsKey(key) && counts.getCount(key) != feats.get(key))
                            throw new RuntimeException("Duplicate key with different values:" + key);
                        counts.setCount(key, feats.get(key));
                    }
            }
        };
    }

    public UnivariateFeatureExtractor<CTMCState> getUnivariateFeatureExtractor()
    {
        final Map<CTMCState, Map<String, Double>> binaryFeaturesMap = unaryFeaturesMap(unaryFeatures);
        return new UnivariateFeatureExtractor<CTMCState>() {

            @SuppressWarnings("unchecked")
            @Override
            public void extract(@SuppressWarnings("rawtypes") Counter counts, CTMCState state)
            {
                Map<String,Double> feats = binaryFeaturesMap.get(state);
                if (feats != null)
                    for (String key : feats.keySet())
                    {
                        if (counts.containsKey(key) && counts.getCount(key) != feats.get(key))
                            throw new RuntimeException("Duplicate key with different values:" + key);
                        counts.setCount(key, feats.get(key));
                    }
            }
        };
    }

    public Indexer<CTMCState> getCTMCStateIndexer()
    {
        Indexer<CTMCState> result = new Indexer<CTMCState>();

        for (int category = 0; category < nCategories; category++)
            for (String latent : orderedLatents)
                result.addToIndex(new CTMCState(category, latent, partition));

        return result;
    }

    public UndirectedGraph<CTMCState, ?> getSupport()
    {
        UndirectedGraph<CTMCState, ?> result = GraphUtils.newUndirectedGraph();

        for (CTMCState state : getCTMCStateIndexer().objects())
            result.addVertex(state);

        if (fullSupport && (supportEdges != null) || (!fullSupport && (supportEdges == null)))
            throw new RuntimeException("Only one of fullSupport or supportEdges should be specified");


        for (int category = 0; category < nCategories; category++)
            if (fullSupport)
            {
                for (String latent0 : orderedLatents)
                    for (String latent1 : orderedLatents)
                        if (latent0 != latent1)
                            result.addEdge(new CTMCState(category, latent0, partition), new CTMCState(category, latent1, partition));
            }
            else
            {
                for (List<String> edge : supportEdges)
                {
                    if (edge.size() != 2)
                        throw new RuntimeException("Edges should contain exactly two latent states. Problem:" + edge);

                    CTMCState
                            state0 = new CTMCState(category, edge.get(0), partition),
                            state1 = new CTMCState(category, edge.get(1), partition);

                    result.addEdge(state0, state1);
                }
            }

        return result;
    }

    @DesignatedConstructor
    public static SerializedExpFamMixture fromJSON(@Input String file)
    {
        return fromJSONString(BriefIO.fileToString(new File(file)));
    }
    
    public static SerializedExpFamMixture fromJSON(File file)
    {
        return fromJSONString(BriefIO.fileToString(file));
    }
    
    public static SerializedExpFamMixture fromJSONString(String string)
    {
        return new Gson().fromJson(string, SerializedExpFamMixture.class);
    }

    public static SerializedExpFamMixture fromResource(String resourceURL)
    {
        String jsonString = BriefIO.resourceToString(resourceURL, SerializedExpFamMixture.class);
        return fromJSONString(jsonString);
    }

    public static Map<CTMCState,Map<String,Double>> unaryFeaturesMap(List<UnaryFeature> features)
    {
        Map<CTMCState,Map<String,Double>> result = Maps.newHashMap();
        for (UnaryFeature feature : features)
            result.put(feature.state, feature.features);
        return result;
    }

    public static Map<UnorderedPair<CTMCState, CTMCState>, Map<String,Double>> binaryFeaturesMap(List<BinaryFeature> features)
    {
        Map<UnorderedPair<CTMCState, CTMCState>, Map<String,Double>> result = Maps.newHashMap();
        for (BinaryFeature feature : features)
            result.put(UnorderedPair.of(feature.state0, feature.state1), feature.features);
        return result;
    }

    public static class UnaryFeature
    {
        public CTMCState state;
        public Map<String,Double> features;
        @Override
        public String toString()
        {
            return "UnaryFeature [state=" + state + ", features=" + features + "]";
        }
        // create a constructor for UnaryFeature
        public UnaryFeature(CTMCState state, Map<String,Double> features){
        	this.state = state;
        	this.features = features;
        }
        
        
    }

    public static class BinaryFeature
    {
        public CTMCState state0, state1;
        public Map<String,Double> features;
        // create a constructor for BinaryFeature
        public BinaryFeature(CTMCState state0, CTMCState state1, Map<String, Double> features){
        	this.state0 = state0;
        	this.state1 = state1;
        	this.features = features;
        }
        @Override
        public String toString()
        {
            return "BinaryFeature [state0=" + state0 + ", state1=" + state1
                    + ", features=" + features + "]";
        }
    }

    public RateMatrixToEmissionModel getEmissionModel()
    {
        // TODO: functionality not yet implemented
        return null;
    }

    public Indexer<String> getLatentIndexer()
    {
        Indexer<String> result = new Indexer<String>();
        for (String latent : orderedLatents)
            result.addToIndex(latent);
        return result;
    }



    @Override
    public String toString()
    {
        return "SerializedExpFamMixture [nCategories=" + nCategories
                + ", orderedLatents=" + orderedLatents + ", supportEdges="
                + supportEdges + ", unaryFeatures=" + unaryFeatures
                + ", binaryFeatures=" + binaryFeatures + ", fullSupport=" + fullSupport
                + "]";
    }


    public static SerializedExpFamMixture rateMtxModel(final RateMtxNames selectedRateMtx) {
        if (selectedRateMtx == null) {
            throw new IllegalArgumentException("model is null!");
        }

        return selectedRateMtx.getSerialized();
    }





    public static SerializedExpFamMixture kimura1980()
    {
        return FeatureFactory.kimura1980();
    }

    public static SerializedExpFamMixture dnaGTR()
    {
        return FeatureFactory.dnaGTR();
    }


    public static SerializedExpFamMixture polarity()
    {
    	return FeatureFactory.polarity();
    }



	public static SerializedExpFamMixture polaritySize()
    {
        return  FeatureFactory.polaritySize();
        	
    }

    public static SerializedExpFamMixture polaritySizeGTR()
    {
        return  FeatureFactory.polaritySizeGTR();
    }

    public static SerializedExpFamMixture proteinSimpleGTR(){

        return FeatureFactory.proteinSimpleGTR();
        		
    }


}