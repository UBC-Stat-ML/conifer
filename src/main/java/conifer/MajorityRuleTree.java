package conifer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import org.apache.commons.io.FilenameUtils;

import briefj.BriefIO;
import briefj.run.Results;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.io.FastaUtils;
import jebl.evolution.coalescent.ExponentialGrowth;
import jebl.evolution.io.ImportException;
import jebl.evolution.io.NewickExporter;
import jebl.evolution.io.NewickImporter;
import jebl.evolution.io.NexusExporter;
import jebl.evolution.io.TreeExporter;
import jebl.evolution.taxa.Taxon;
import jebl.evolution.trees.ConsensusTreeBuilder;
import jebl.evolution.trees.RootedTree;
import jebl.evolution.trees.Tree;
import jebl.evolution.trees.TreeBuilderFactory;
import jebl.evolution.treesimulation.CoalescentIntervalGenerator;
import jebl.evolution.treesimulation.IntervalGenerator;
import jebl.evolution.treesimulation.TreeSimulator;

public class MajorityRuleTree {

	public static void buildAndWriteConsensusTree(File inputFile, File outFile) {
		Tree[] trees = null;
		try {
			trees = loadTrees(inputFile);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ImportException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

//		System.out.println("Working Directory = "
//				+ System.getProperty("user.dir"));

		ConsensusTreeBuilder<?> treeBuilder = TreeBuilderFactory.buildRooted(
				trees, .5, TreeBuilderFactory.ConsensusMethod.GREEDY);
		Tree consensusTree = treeBuilder.build();

		// write the tree
//		System.out.println("Working Directory = "
//				+ System.getProperty("user.dir"));

		try {
			writeConsensusTree(consensusTree, outFile, false);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private static Tree[] loadTrees(File inputFile) throws IOException,
			ImportException {

		List<RootedTree> treeList = new ArrayList<RootedTree>(100);

		NewickImporter importer = new NewickImporter(new BufferedReader(
				new FileReader(inputFile)), true);
		int maxTrees = 100;
		int burnIn = 0;
		int cnt = 0;
		int cntBurnin = 0;
		while (importer.hasTree() && cnt < maxTrees) {
			RootedTree tree = (RootedTree) importer.importNextTree();
			if (cntBurnin >= burnIn) {
				treeList.add(tree);
				cnt++;
			}
			cntBurnin++;
		}

		return treeList.toArray(new Tree[treeList.size()]);
	}

	public static void writeConsensusTree(Tree consensusTree, File outputFile,
			boolean isFileTypeNewick) throws IOException {
		Writer writer = new FileWriter(outputFile);
		TreeExporter treeWriter = (isFileTypeNewick ? new NewickExporter(writer)
				: new NexusExporter(writer));
		treeWriter.exportTree(consensusTree);
		writer.close();
	}
	
	
	/**
	 * Constructs a nwk topology with branch lengths from the given fastaFile, basically, 
	 * gets the fasta file, parses it and extracts species names and creates a tree 
	 * and augment with branch lengths and writes to file.
	 */ 
	public static File randomUnrootedPhylogenyFromMultipleAlignmentFastaFile(File fastaFile) {
		
		Map<TreeNode,CharSequence> data = FastaUtils.readFasta(fastaFile);
		UnrootedTree tree = UnrootedTreeLikelihood.defaultTree(data.keySet());
		String outputFileName = FilenameUtils.removeExtension(fastaFile.getName()) + "_" + System.currentTimeMillis() + ".nwk";
		File outputFile = Results.getFileInResultFolder(outputFileName);
		final PrintWriter treeWriter = BriefIO.output(outputFile);
		treeWriter.println(tree.toNewick());
		treeWriter.flush();
		return outputFile;
	}
	
	/**
	 * Use JEBL instead
	 * @param fastaFile
	 * @return
	 */
	//TODO: this is an inelegant implementation, doesn't have to be so coupled. Refactor in due time.
	public static File treeFromFastaFile(File fastaFile) { 
		Map<TreeNode,CharSequence> data = FastaUtils.readFasta(fastaFile);
		Collection<Taxon> taxa = new Vector<Taxon>();
		
		for (TreeNode item : data.keySet()) {
			Taxon taxon = Taxon.getTaxon(item.toString());
			taxon.setAttribute("heightAttribute", 0.0);
			taxa.add(taxon);
		}
	
		// 2. Now given the taxa, simulate a tree
		TreeSimulator t = new TreeSimulator(taxa, "heightAttribute");
		ExponentialGrowth exponentialGrowth = new ExponentialGrowth();
		exponentialGrowth.setN0(10);
		exponentialGrowth.setGrowthRate(0.1);

		IntervalGenerator intervals = new CoalescentIntervalGenerator(exponentialGrowth);
		RootedTree tree = t.simulate(intervals);
		
		// 3. Write the tree down to the file
		String outputFileName = FilenameUtils.removeExtension(fastaFile.getName()) + "_" + System.currentTimeMillis() + ".nwk";
		File outputFile = Results.getFileInResultFolder(outputFileName);
		Writer writer = null;
		try {
			writer = new FileWriter(outputFile);
			TreeExporter treeWriter = new NewickExporter(writer);
			treeWriter.exportTree(tree);
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			throw new RuntimeException("Where is the file?");
		}
		
		return outputFile;		
	}
	
	
	public static File FESape4Tree() {
		return new File("src/main/resources/conifer/sampleInput/FES.ape.4.nwk");
	}
	
	public static File UTYape4Tree() {
		return new File("src/main/resources/conifer/sampleInput/UTY.ape.4.nwk");
	}
}
