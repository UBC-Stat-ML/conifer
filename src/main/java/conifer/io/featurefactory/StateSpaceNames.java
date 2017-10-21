package conifer.io.featurefactory;

import java.util.List;

import conifer.io.PhylogeneticObservationFactory;

public enum StateSpaceNames {
	DNA{

		@Override
		public List<String> getStateSpace() {
			// TODO Auto-generated method stub
			return PhylogeneticObservationFactory.nucleotidesFactory().orderedSymbols;
		}
		
		
		
		
	},
	
	PROTEIN{

		@Override
		public List<String> getStateSpace() {
			// TODO Auto-generated method stub
			return PhylogeneticObservationFactory.proteinFactory().orderedSymbols;
		}
		
		
		
		
		
	}, 
	CODON{

		@Override
		public List<String> getStateSpace() {
			// TODO Auto-generated method stub
			return PhylogeneticObservationFactory.codonFactory().orderedSymbols;
		}
		
		
		
		
		
		
	};
	
	
	public abstract List<String> getStateSpace();
	

}
