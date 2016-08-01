package conifer.local;


/**
 * A ray is parameterized by a starting point (time, position at that time),
 * and the velocity just after that time (i.e. that starting point is
 * assumed to be a collision, so the we store the velocity just after the
 * bounce).
 *
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 */
public class TrajectoryRay {
    /*
     * last collision time
     */
    public final double t;

    /*
     * position and velocity right after that last collision
     */
    public final double position_t;
    public final double velocity_t;

    public TrajectoryRay(double t, double position_t, double velocity_t) {
        this.t = t;
        this.position_t = position_t;
        this.velocity_t = velocity_t;
    }

    /**
     * @param time
     * @return Position at some time larger or equal to the last collision time, assuming
     * no collision affecting that coordinate until that time.
     */
    public double position(double time) {
        if (time < t)
            throw new RuntimeException("Current time cannot be smaller than start t (time=" + time + ",t=" + t + ")");

        return position_t + (time - t) * velocity_t;
    }

    @Override
    public String toString() {
        return "TrajectoryRay [t=" + t + ", position_t=" + position_t
                + ", velocity_t=" + velocity_t + "]";
    }
}