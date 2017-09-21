package conifer.local;

import org.jblas.DoubleMatrix;

/**
 * Created by crystal on 2016-10-18.
 */
public class ARFSDataStruct {
    public DoubleMatrix next_q = null;
    public DoubleMatrix q = null;
    public double trajectoryLength = -1;

    public ARFSDataStruct(DoubleMatrix next_q,
                          DoubleMatrix q, double trajectoryLength) {
        this.next_q = next_q;
        this.q = q;
        this.trajectoryLength = trajectoryLength;

    }
}
