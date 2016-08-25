package ucl.silver.d3d.core;

/**
 * <p>Title: D3D</p>
 * <p>Description: 3D Diffusion-Reaction Simulator</p>
 * <p>Copyright: Copyright (c) 2016</p>
 * <p>Company: The Silver Lab at University College London</p>
 * @author Jason S. Rothman
 * @version 1.6
 */

public class PulseTimer extends ParamVector {

    public int numPulses = 0;

    public boolean impulse = false;

    public Pulse[] pulses = null;

    public double[] timer = null; // pulse array used during simulations

    @Override
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("numPulses")) {
            return false;
        }
        return true;
    }

    public PulseTimer(Project p, double TIME) {
        super(p);
        add(TIME);
    }

    public PulseTimer(Project p, double TIME, double DURATION) {
        super(p);
        add(TIME, DURATION);
    }

    public PulseTimer(Project p, double TIME, double DURATION, double AMPLITUDE) {
        super(p);
        add(TIME, DURATION, AMPLITUDE);
    }

    public PulseTimer(Project p, double TIME, double DURATION, double AMPLITUDE, int numPulses, double pulseInterval) {
        super(p);
        addTrain(TIME, DURATION, AMPLITUDE, numPulses, pulseInterval);
    }

    public final int add(double TIME) {
        return addPulse(new Pulse(project, TIME, 0, 1));
    }

    public final int add(double TIME, double DURATION) {
        return addPulse(new Pulse(project, TIME, DURATION, 1));
    }

    public final int add(double TIME, double DURATION, double AMPLITUDE) {
        return addPulse(new Pulse(project, TIME, DURATION, AMPLITUDE));
    }

    public final int addTrain(double TIME, double DURATION, double AMPLITUDE, int numPulses, double pulseInterval) {

        int j = -1;

        for (int i = 0; i < numPulses; i++) {
            j = addPulse(new Pulse(project, TIME + i * pulseInterval, DURATION, AMPLITUDE));
        }

        return j;

    }

    private int addPulse(Pulse newPulse) {

        int i = 0;

        if (pulses != null) {
            i = pulses.length;
        }

        Pulse[] newArray = new Pulse[i+1];

        if (i > 0) {
            System.arraycopy(pulses, 0, newArray, 0, i);
        }

        newArray[i] = newPulse;

        pulses = newArray; // replace old array with new one

        numPulses = pulses.length;

        //Master.log("Added pulse #" + i);

        return i;

    }

    public final boolean initTimer() {

        int it;
        double time, duration, sumPulses = 0;

        int itmax = (int) (project.simTime / project.dt) + 1;

        if ((pulses == null) || (itmax <= 0)) {
            timer = null;
            return false;
        }

        timer = new double[itmax];

        if (impulse) {
            
            for (int i = 0; i < pulses.length; i += 1) {

                it = (int) (pulses[i].time / project.dt);

                if ((it < 0) || (it >= itmax)) {
                    continue;
                }

                timer[it] = 1;

            }
            
        } else {

            for (it = 0; it < itmax; it++) {

                sumPulses = 0;

                for (int i = 0; i < pulses.length; i += 1) {

                    time = it * project.dt;
                    duration = pulses[i].duration;

                    if (duration <= 0) {
                        duration = project.dt;
                    }

                    if ((time >= pulses[i].time) && (time <= pulses[i].time + duration)) {
                        sumPulses += pulses[i].amplitude;
                    }

                }
                
                timer[it] = sumPulses;
                
            }

        }

        return true;

    }

    @Override
    public boolean addUser(ParamVector pv) {

        super.addUser(pv);

        if (pulses != null) {
            for (int i = 0; i < pulses.length; i++) {
                if (pulses[i] != null) {
                    pulses[i].addUser(pv);
                }
            }
        }

        return true;

    }

    @Override
    public boolean createVector(boolean close) {

        if (!super.createVector(false)) {
            return false;
        }

        if (pulses != null) {
            for (int i = 0; i < pulses.length; i++) {
                if (pulses[i] != null) {
                    addBlankParam();
                    pulses[i].createVector(true);
                    addVector(pulses[i].getVector());
                    pulses[i].addUser(this);
                }
            }
        }

        if (close) {
            closeVector();
        }

        return true;

    }

    @Override
    public void updateVector(ParamObject[] v) {

        super.updateVector(v);

        if (pulses != null) {
            for (int i = 0; i < pulses.length; i++) {
                if (pulses[i] != null) {
                    pulses[i].updateVector(v);
                }
            }
        }

    }
}
