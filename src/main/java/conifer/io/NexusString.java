package conifer.io;

import jebl.evolution.graphs.Node;
import jebl.evolution.io.TreeExporter;
import jebl.evolution.taxa.Taxon;
import jebl.evolution.trees.RootedTree;
import jebl.evolution.trees.Tree;
import jebl.evolution.trees.Utils;
import jebl.util.Attributable;

import java.awt.*;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Writer;
import java.util.*;
import java.util.List;

/**
 * Export trees to Nexus format.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Joseph Heled
 * @author Sohrab Salehi (sohrab.salehi@gmail.com)
 * @version $Id: NexusExporter.java 1060 2010-06-01 08:52:55Z rambaut $
 */

public class NexusString implements TreeExporter {

	public NexusString(Writer writer) {
		this(writer, false);
	}

	/**
     *
     * @param writer where export text goes
     */
    public NexusString(Writer writer, boolean writeMetaComments) {
		this(writer, writeMetaComments, false);
    }

    /**
     *
     * @param writer where export text goes
     */
    public NexusString(Writer writer, boolean writeMetaComments, boolean interleave) {
		this.writeMetaComments = writeMetaComments;
        this.writer = new PrintWriter(writer);
    }

    
    /**
     * Export a single tree
     *
     * @param tree
     * @throws java.io.IOException
     */
    public void exportTree(Tree tree) throws IOException {
        List<Tree> trees = new ArrayList<Tree>();
        trees.add(tree);
        exportTrees(trees);
    }

    private void writeTrees(Collection<? extends Tree> trees, boolean checkTaxa) throws IOException {
        for( Tree t : trees ) {
            if( checkTaxa && !establishTreeTaxa(t) ) {
                throw new IllegalArgumentException();
            }
            final boolean isRooted = t instanceof RootedTree;
            final RootedTree rtree = isRooted ? (RootedTree)t : Utils.rootTheTree(t);


            StringBuilder builder = new StringBuilder("");

            // TREE & UTREE are depreciated in the NEXUS format in favour of a metacomment
            // [&U] or [&R] after the TREE command. Andrew.
            // TT: The [&U], [&R] should actually come *after* the " = " and be uppercase, see
            // e.g. tree_rest in http://www.cs.nmsu.edu/~epontell/nexus/nexus_grammar .
            // Before 2008-05-05 we incorrectly inserted it before the treeName.

            appendAttributes(rtree, exportExcludeKeys, builder);

            appendTree(rtree, rtree.getRootNode(), builder);
            builder.append(";");

            writer.println(builder);
        }
    }

    public void exportTrees(Collection<? extends Tree> trees) throws IOException {
        exportTrees(trees, false);
    }

    public void exportTrees(Collection<? extends Tree> trees, boolean writeTaxa) throws IOException {
        if (writeTaxa) {
            TreeSet<Taxon> taxa = new TreeSet<Taxon>();
            for (Tree tree : trees) {
                taxa.addAll(tree.getTaxa());
            }

            establishTaxa(taxa);
        }

        //writer.println("begin trees;");
        writeTrees(trees, false);
        //writer.println("end;");
    }


    /**
     * Write a new taxa block and record them for later reference.
     * @param taxonArray
     */
    private void setTaxa(Taxon[] taxonArray) {
        this.taxa = new HashSet<Taxon>();

        writer.println("begin taxa;");
        writer.println("\tdimensions ntax=" + taxonArray.length + ";");
        writer.println("\ttaxlabels");

        for (Taxon taxon : taxonArray) {
            taxa.add(taxon);

            StringBuilder builder = new StringBuilder("\t");
            appendTaxonName(taxon, builder);
            appendAttributes(taxon, null, builder);
            writer.println(builder);
        }
        writer.println(";\nend;\n");
    }

    final private String nameRegex = "^(\\w|-)+$";

    /**
     * name suitable for printing - quotes if necessary
     * @param taxon
     * @param builder
     * @return
     */
    private StringBuilder appendTaxonName(Taxon taxon, StringBuilder builder) {
        String name = taxon.getName();
        if (!name.matches(nameRegex)) {
            // JEBL way of quoting the quote character
            name = name.replace("\'", "\'\'");
            builder.append("\'").append(name).append("\'");
            return builder;
        }
        return builder.append(name);
    }



    private boolean establishTreeTaxa(Tree tree) {
        return establishTaxa(tree.getTaxa());
    }

    private boolean establishTaxa(Collection<? extends Taxon> ntaxa) {
        if( taxa != null && taxa.size() == ntaxa.size()  && taxa.containsAll(ntaxa)) {
            return false;
        }

        setTaxa(ntaxa.toArray(new Taxon[]{}));
        return true;
    }

    /**
     * Prepare for writing a tree. If a taxa block exists and is suitable for tree,
     * do nothing. If not, write a new taxa block.
     * @param tree
     * @param node
     * @param builder
     */
    private void appendTree(RootedTree tree, Node node, StringBuilder builder) {
        if (tree.isExternal(node)) {
            appendTaxonName(tree.getTaxon(node), builder);

            appendAttributes(node, null, builder);

            if( tree.hasLengths() ) {
                builder.append(':');
                builder.append(roundDouble(tree.getLength(node), 6));
            }
        } else {
            builder.append('(');
            List<Node> children = tree.getChildren(node);
            final int last = children.size() - 1;
            for (int i = 0; i < children.size(); i++) {
                appendTree(tree, children.get(i), builder);
                builder.append(i == last ? ')' : ',');
            }

            appendAttributes(node, null, builder);

            Node parent = tree.getParent(node);
            // Don't write root length. This is ignored elsewhere and the nexus importer fails
            // whet it is present.
            if (parent != null) {
                if (tree.hasLengths()) {
                    builder.append(":").append(roundDouble(tree.getLength(node), 6));
                }
            }
        }
    }

    public static double roundDouble(double value, int decimalPlace) {
        double power_of_ten = 1;
        while (decimalPlace-- > 0)
            power_of_ten *= 10.0;
        return Math.round(value * power_of_ten) / power_of_ten;
    }

    private StringBuilder appendAttributes(Attributable item, String[] excludeKeys, StringBuilder builder) {
	    if (!writeMetaComments) {
		    return builder;
	    }

        boolean first = true;
        for( String key : item.getAttributeNames() ) {
            // we should replace the explicit check for name by something more general.
            // Like a reserved character at the start (here &). however we have to worry about backward
            // compatibility so no change yet with name.
            boolean exclude = false;
            if(excludeKeys != null) {
                for(String eKey : excludeKeys) {
                    if(eKey.equals(key)) {
                        exclude = true;
                    }
                }
            }
            Object value = item.getAttribute(key);

            if( !exclude && !key.startsWith("&") && value != null ) {
                if (first) {
                    builder.append("[&");
                    first = false;
                } else {
                    builder.append(",");
                }

                if( key.indexOf(' ') < 0 ) {
                    builder.append(key);
                } else {
                    builder.append("\"").append(key).append("\"");
                }

                builder.append('=');

                appendAttributeValue(value, builder);
            }
        }
        if (!first) {
            builder.append("]");
        }

        return builder;
    }

    private StringBuilder appendAttributeValue(Object value, StringBuilder builder) {
        if (value instanceof Object[]) {
            builder.append("{");
            Object[] elements = ((Object[])value);

            if (elements.length > 0) {
                appendAttributeValue(elements[0], builder);
                for (int i = 1; i < elements.length; i++) {
                    builder.append(",");
                    appendAttributeValue(elements[i], builder);
                }
            }
            return builder.append("}");
        }

        if (value instanceof Color) {
            return builder.append("#").append(Integer.toHexString(((Color)value).getRGB()).substring(2));
        }

        if (value instanceof String) {
            return builder.append("\"").append(value).append("\"");
        }

        return builder.append(value);
    }

    public static final String treeNameAttributeKey = "name";
    public static final String[] exportExcludeKeys = new String[] {treeNameAttributeKey, "R", "U"};

    static public boolean isGeneratedTreeName(String name) {
        return name != null && name.matches("tree_[0-9]+");
    }

    private Set<Taxon> taxa = null;
    protected final PrintWriter writer;
	private boolean writeMetaComments;
    public static final int MAX_ROW_LENGTH = 60;
}
