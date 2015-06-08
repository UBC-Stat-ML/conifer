package conifer.rgen;

import java.io.File;

import bayonet.rplot.PlotHistogram;
import bayonet.rplot.RJavaBridge;
import bayonet.rplot.RUtils;
import briefj.BriefIO;

public class GenerateIupacEncoding extends RJavaBridge {
	public final double maxCopyNumber;

	private GenerateIupacEncoding(int maxCopyNumber)
	{
		//cn-ctmc-iupac-encoding
		this.maxCopyNumber = maxCopyNumber;
	}

	public static GenerateIupacEncoding withCopyNumber(int maxCopyNumber)
	{
		return new GenerateIupacEncoding(maxCopyNumber);
	}

	public String getResourceDir() { return System.getProperty("user.dir") + "/src/main/resources/"; }

	@Override public String rTemplateResourceURL() { return "/conifer/rgen/makeUPACCNCTMC.txt"; }
	
	 public String getOutput()
	 {
	    return RUtils.escapeQuote(output);
	 }
	 
	 private transient String output;
	 
	 public void toUPAC(String output)
	 {
	    this.output = output;
	    RUtils.callRBridge(this);
	 }
	
	public static void main(String[] args) {
		int maxCopyNumber = 6;
		GenerateIupacEncoding.withCopyNumber(maxCopyNumber).toUPAC("cn-ctmc-iupac-encoding-" + maxCopyNumber + ".txt");
	}
}
