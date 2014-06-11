package conifer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.List;

import jebl.evolution.io.ImportException;
import jebl.evolution.io.NewickExporter;
import jebl.evolution.io.NewickImporter;
import jebl.evolution.io.NexusExporter;
import jebl.evolution.io.TreeExporter;
import jebl.evolution.trees.ConsensusTreeBuilder;
import jebl.evolution.trees.RootedTree;
import jebl.evolution.trees.Tree;
import jebl.evolution.trees.TreeBuilderFactory;

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
			if (cntBurnin > burnIn) {
				treeList.add(tree);
				cnt++;
			}
			cntBurnin++;
		}

		return treeList.toArray(new Tree[treeList.size()]);
	}

	private static void writeConsensusTree(Tree consensusTree, File outputFile,
			boolean isFileTypeNewick) throws IOException {
		Writer writer = new FileWriter(outputFile);
		TreeExporter treeWriter = (isFileTypeNewick ? new NewickExporter(writer)
				: new NexusExporter(writer));
		treeWriter.exportTree(consensusTree);
		writer.close();
	}
}
