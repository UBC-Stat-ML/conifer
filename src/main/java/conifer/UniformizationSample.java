package conifer;

import java.io.*;
import java.util.Arrays;
import java.util.Random;

import briefj.BriefArrays;
import briefj.opt.Option;
import com.fasterxml.jackson.databind.JsonMappingException;
import conifer.ctmc.EndPointSampler;
import conifer.ctmc.expfam.ExpFamMixture;
import conifer.ctmc.expfam.RateMtxNames;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.MultiCategorySubstitutionModel;
import briefj.run.Mains;
import briefj.run.Results;
import briefj.BriefIO;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.JsonMappingException;
import conifer.processors.FileNameString;
import org.apache.commons.io.FileUtils;
import org.jblas.DoubleMatrix;


public class UniformizationSample implements Runnable
{
    @Option(gloss="File of provided alignment")
    public File inputFile;


    @Option(gloss="File of the tree topology")
    public File treeFile;

    @Option(gloss="ESS Experiment Number")
    public int rep = 1;

    @Option(gloss="Rate Matrix Method")
    public RateMtxNames selectedRateMtx = RateMtxNames.POLARITYSIZEGTR;

    @Option(gloss = "True rate matrix generating data")
    public File rateMtxFile;

    @Option(gloss="Use cache or not")
    public boolean cached;

    private final PrintWriter detailWriter = BriefIO.output(Results.getFileInResultFolder("experiment.details.txt"));


    public void run()  {
        ObjectMapper mapper = new ObjectMapper();
        double[][] array;
        EndPointSampler.cached=cached;

        try (FileInputStream in = new FileInputStream(rateMtxFile)) {
            array = mapper.readValue(in, double[][].class);
            long startTime = System.currentTimeMillis();
            UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood1 =
                    UnrootedTreeLikelihood
                            .fromFastaFile(inputFile, selectedRateMtx)
                            .withSingleRateMatrix(array)
                            .withExpFamMixture(ExpFamMixture.rateMtxModel(selectedRateMtx))
                            .withTree(treeFile);
            Random rand = new Random(1);
            likelihood1.evolutionaryModel.samplePosteriorPaths(rand, likelihood1.observations, likelihood1.tree);

            String fileName = inputFile.getName();
            FileNameString fileNameString = new FileNameString(fileName);
            String curBr = fileNameString.subStringBetween(fileName, "br", "seed");
            String curSeed = fileNameString.subStringBetween(fileName, "seed", "align");

            logToFile("Total time in seconds: " + ((System.currentTimeMillis() - startTime) / 1000.0));
            logToFile("branch length:"+ curBr);
            logToFile("used seed:"+ curSeed);
            logToFile("cached" + cached);

            File newDirectory = new File(Results.getResultFolder().getParent() + "br"+ curBr+ "seed" + curSeed+
            "cached"+ cached);
            newDirectory.mkdir();
            try
            {
                FileUtils.copyDirectory(Results.getResultFolder(), newDirectory);
            } catch (IOException e)
            {
                e.printStackTrace();
            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (JsonMappingException e) {

        } catch (IOException e) {
            e.printStackTrace();
        }


    }

    public static void main(String [] args)
    {
        Mains.instrumentedRun(args, new UniformizationSample());
    }


    public void logToFile(String someline) {
        this.detailWriter.println(someline);
        this.detailWriter.flush();
    }


}
