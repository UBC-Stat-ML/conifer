package conifer.rejfreemodels.phylo;

import blang.annotations.DefineFactor;
import blang.annotations.FactorArgument;
import blang.annotations.FactorComponent;
import blang.factors.FactorList;
import blang.variables.RealVariable;
import conifer.ctmc.expfam.CTMCExpFam;
import conifer.ctmc.expfam.CTMCState;
import conifer.ctmc.expfam.ExpFamParameters;
import conifer.local.CollisionFactor;
import conifer.rejfreemodels.normal.NormalFactor;
import org.jblas.DoubleMatrix;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by crystal on 2016-07-28.
 */
public class ExpectedCompleteReversibleModel {


    @DefineFactor
    public final List<CollisionFactor> localFactors;

    @FactorComponent
    public final ExpFamParameters parameters;

    @FactorArgument(makeStochastic = true)
    public final List<RealVariable> weights;

    public CTMCExpFam<CTMCState> ctmcExpFam;

    public CTMCExpFam<CTMCState>.ExpectedCompleteReversibleObjective auxObjective;

    public ExpectedCompleteReversibleModel( ExpFamParameters parameters, CTMCExpFam<CTMCState>.ExpectedCompleteReversibleObjective objective,
                                            CTMCExpFam<CTMCState> ctmcExpFam) {
        this.parameters = parameters;
        this.ctmcExpFam = ctmcExpFam;
        this.auxObjective = objective;
        this.weights = getWeights();
        this.localFactors = localFactors();
    }

    public List<RealVariable> getWeights(){
        return transformWeightIntoReal(parameters);
    }

    public List<RealVariable> transformWeightIntoReal(ExpFamParameters parameters){

        List<RealVariable> result = new ArrayList<>();
        for(int i=0; i<parameters.getVector().length;i++){
            result.add(RealVariable.real(parameters.getVector()[i]));
        }
        return result;
    }

    private List<CollisionFactor> localFactors() {

        List<CollisionFactor> result = new ArrayList<>();


        // check to see if I need to add the Normal factor for weights
        DoubleMatrix precision = new DoubleMatrix().eye(parameters.getVector().length);

        CollisionFactor g = new NormalFactor(precision, weights);
        result.add(g);


        //add initial count factors first
        for(int state=0; state< parameters.globalExponentialFamily.nStates; state++){
            CollisionFactor f = new InitialCountFactor(parameters,auxObjective,ctmcExpFam, state, weights);
            result.add(f);
        }

        //add sojourn time factors
        for(int state =0; state< parameters.globalExponentialFamily.nStates;state++){
            final int [] curSupports = ctmcExpFam.supports[state];
            for(int endIdx =0; endIdx < curSupports.length; endIdx++){
                CollisionFactor f = new SojournTimeFactor(parameters, auxObjective, ctmcExpFam, state, curSupports[endIdx],weights);
                result.add(f);
            }
        }
        //add transition count factors
        for(int state=0; state < parameters.globalExponentialFamily.nStates;state++){
            final int [] curSupports = ctmcExpFam.supports[state];
            for(int endIdx =0; endIdx < curSupports.length; endIdx++){
                CollisionFactor f = new TransitionCountFactor(parameters, auxObjective, ctmcExpFam, state, curSupports[endIdx],weights);
                result.add(f);
            }

        }

        return result;
    }
}


