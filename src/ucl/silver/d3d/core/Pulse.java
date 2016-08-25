package ucl.silver.d3d.core;

/**
 * <p>Title: D3D</p>
 * <p>Description: 3D Diffusion-Reaction Simulator</p>
 * <p>Copyright: Copyright (c) 2016</p>
 * <p>Company: The Silver Lab at University College London</p>
 * @author Jason S. Rothman
 * @version 1.6
 */

public class Pulse extends ParamVector {

    public double time = 0;
    public double duration = 0;
    public double amplitude = 1;

    public Pulse(Project p, double TIME, double DURATION, double AMPLITUDE) {
        super(p);
        time = TIME;
        duration = DURATION;
        amplitude = AMPLITUDE;
    }

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("time")) {
            return project.timeUnits;
        }
        if (name.equalsIgnoreCase("duration")) {
            return project.timeUnits;
        }
        return super.units(name);
    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        if (super.setMyParams(o, v)) {
            return true;
        }

        String n = o.getName();

        if (n.equalsIgnoreCase("time")) {
            time = v;
            return true;
        }
        if (n.equalsIgnoreCase("duration")) {
            duration = v;
            return true;
        }
        if (n.equalsIgnoreCase("amplitude")) {
            amplitude = v;
            return true;
        }
        return false;
    }
}
