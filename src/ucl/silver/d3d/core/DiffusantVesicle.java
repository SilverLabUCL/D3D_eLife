package ucl.silver.d3d.core;

import ucl.silver.d3d.utils.*;

/**
 * <p>Title: D3D</p>
 * <p>Description: 3D Diffusion-Reaction Simulator</p>
 * <p>Copyright: Copyright (c) 2016</p>
 * <p>Company: The Silver Lab at University College London</p>
 * @author Jason S. Rothman
 * @version 1.6
 */

public class DiffusantVesicle extends ParamVector {

    public double radius = 0; // um
    public double fluorescence = 1.0;

    public double D; // diffusion coefficient (um^2/ms)
    public double step3; // um/sqrt(3)

    public double localD; // local diffusion coefficient (um^2/ms)
    public double localStep3; // um/sqrt(3)
    public double DsDff;

    public double x0, y0, z0; // starting center location
    public double x, y, z; // current center location
    public double d20; // distance to center of geometry

    public double sqrDisplacement = -1;
    public double firstCollision = 0;

    public boolean mobile = true;
    public boolean insideGeometry = true;

    public Voxel voxel = null; // voxel where vesicle resides

    public DiffusantVesicle nextReady = null;

    public DiffusantVesicle[] connectTo = null;
    public double[] connectorOffTime = null; // ms
    public double connectorLifeTime = 0;

    private final double invsqrt3 = 1.0 / Math.sqrt(3);

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("radius")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("D")) {
            return project.spaceUnits + "^2/" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("step3")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("extraX")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("extraY")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("extraZ")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("x0")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("y0")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("z0")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("x")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("y")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("z")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("d20")) {
            return project.spaceUnits;
        }
        return super.units(name);
    }

    public DiffusantVesicle(Project p, String Name, double RADIUS, double d, double X, double Y, double Z) {

        super(p);

        name = Name;

        radius = RADIUS;
        D = d;
        x = X;
        y = Y;
        z = Z;
        x0 = X;
        y0 = Y;
        z0 = Z;

        createVector(true);

        setD(D); // computes mean step

    }

    public void setLocation(double X, double Y, double Z, boolean initStartLocation) {

        x = X;
        y = Y;
        z = Z;

        if (initStartLocation) {
            x0 = X;
            y0 = Y;
            z0 = Z;
        }

    }

    public void initStartLocation() {
        x0 = x;
        y0 = y;
        z0 = z;
    }

    public double squareDisplacement() {

        double dx = x - x0;
        double dy = y - y0;
        double dz = z - z0;

        return dx * dx + dy * dy + dz * dz;

    }

    public double squareDistance(DiffusantVesicle dv) {

        double dx = x - dv.x;
        double dy = y - dv.y;
        double dz = z - dv.z;

        return dx * dx + dy * dy + dz * dz;

    }

    public boolean overlap(DiffusantVesicle dv, double extraDistance) {

        double dx = x - dv.x;
        double dy = y - dv.y;
        double dz = z - dv.z;

        double sqrDistance = dx * dx + dy * dy + dz * dz;
        double minDistance = radius + dv.radius + extraDistance;

        return (sqrDistance < minDistance * minDistance);

    }

    public boolean copy(DiffusantVesicle dv) {

        if (dv == null) {
            return false;
        }

        name = dv.name;

        radius = dv.radius;
        fluorescence = dv.fluorescence;

        D = dv.D;
        step3 = dv.step3;
        localStep3 = dv.localStep3;
        DsDff = dv.DsDff;

        x0 = dv.x0;
        y0 = dv.y0;
        z0 = dv.z0;

        x = dv.x;
        y = dv.y;
        z = dv.z;

        d20 = dv.d20;

        mobile = dv.mobile;

        insideGeometry = dv.insideGeometry;

        voxel = dv.voxel;

        nextReady = dv.nextReady;

        return true;

    }

    public void init() {
        setD(D); // computes mean step
    }

    public double getVolume() {
        return volume(radius);
    }

    public static double volume(double RADIUS) {
        if (RADIUS > 0) {
            return 4.0 * Math.PI * Math.pow(RADIUS, 3.0) / 3.0;
        } else {
            return Double.NaN;
        }
    }

    public static double radius(double VOLUME) {
        if (VOLUME > 0) {
            return Math.pow(3.0 * VOLUME / (4.0 * Math.PI), (1.0 / 3.0));
        } else {
            return Double.NaN;
        }
    }

    public static double step(double diffusionConstant, double dt) {
        if ((diffusionConstant > 0) && (dt > 0)) {
            return Math.sqrt(6.0 * diffusionConstant * dt);
        } else {
            return Double.NaN;
        }
    }

    public static double step3(double diffusionConstant, double dt) { // step/sqrt(3)
        if ((diffusionConstant > 0) && (dt > 0)) {
            return Math.sqrt(2.0 * diffusionConstant * dt);
        } else {
            return Double.NaN;
        }
    }

    public static double D(double STEP, double dt) {
        if ((STEP > 0) && (dt > 0)) {
            return STEP * STEP / (6.0 * dt);
        } else {
            return Double.NaN;
        }
    }

    public static double D3(double STEP3, double dt) {
        if ((STEP3 > 0) && (dt > 0)) {
            return STEP3 * STEP3 / (2.0 * dt);
        } else {
            return Double.NaN;
        }
    }

    public static double dt(double diffusionConstant, double STEP) {
        if ((diffusionConstant > 0) && (STEP > 0)) {
            return STEP * STEP / (6.0 * diffusionConstant);
        } else {
            return Double.NaN;
        }
    }

    public static double Dratio_short(double mobileVolumeFraction, double immobileVolumeFraction) { // D/D0 ratio from Tokuyama and Oppenheim 1994

        double Dratio;
        double Dratio_mobile = Dratio_short_Tokuyama(mobileVolumeFraction);
        double Dratio_immobile = Dratio_short_Freed(immobileVolumeFraction);
        double dff = Dff_short_Banchio(mobileVolumeFraction);

        if (immobileVolumeFraction > 0) {
            Dratio = Dratio_mobile / (1 + (Dratio_mobile / dff) * ((1 / Dratio_immobile) - 1)); // Michailidou et al. 2009 Eq. 3
            //Dratio = Dratio_mobile * Dratio_immobile; // simple method
        } else {
            Dratio = Dratio_mobile;
        }

        return Dratio;

    }

    public static double Dff_short_Banchio(double mobileVolumeFraction) {
        return 1 - 1.5 * mobileVolumeFraction + 0.75 * mobileVolumeFraction * mobileVolumeFraction; // Banchio & Brady 2003
    }

    public static double Dratio_short_Tokuyama(double mobileVolumeFraction) { // Tokuyama and Oppenheim 1994

        double b = Math.sqrt(9.0 * mobileVolumeFraction / 8.0);
        double c = 11.0 * mobileVolumeFraction / 16.0;
        double HI_mobile = (2 * b * b / (1 - b)) - (c / (1 + 2 * c)) - (b * c * (2 + c) / ((1 + c) * (1 - b + c)));

        HI_mobile = 1 / (1 + HI_mobile);

        return HI_mobile;
    }

    public static double Dratio_short_Freed(double immobileVolumeFraction) { // Ds/D0 ratio from Freed and Muthukumar 1978

        double viscosity = 0.0007;
	double radius = 0.025; // um
	double volume = 4 * Math.PI * radius * radius * radius / 3;
        double concentration = immobileVolumeFraction / volume;
        double drag0 = 6 * Math.PI * viscosity * radius;

        double dratio, dratioFound = Double.NaN;
        double dmax = 1;
        double dstep = 0.001;
        double isteps = 1 + dmax / dstep;

        double eqLeft, eqRight, ak;

        for (int i = 0; i < isteps; i++) {

            dratio = dmax - i * dstep;
            eqLeft = drag0 / dratio;
            ak = radius * Math.sqrt(concentration * eqLeft / viscosity);
            eqRight = drag0 * (1 + ak + (ak * ak / 3));

            if (eqRight < eqLeft) {
                dratioFound = dratio;
                break;
            }

        }

        //Master.log("Freed Ds/D0 = " + dratioFound);

        return dratioFound;

    }

    public final void setD(double d) {

        D = d;
        step3 = step3(D, project.dt);

        setParamObject("D", D);
        setParamObject("step3", step3);

    }

    public void setStep(double s) {

        step3 = s * invsqrt3;
        D = D(s, project.dt);

        setParamObject("D", D);
        setParamObject("step3", step3);

    }

    public void setStep3(double s3) {

        step3 = s3;
        D = D(s3 * Math.sqrt(3), project.dt);

        setParamObject("D", D);
        setParamObject("step3", step3);

    }

    public void unconnectAll(){

        if (connectTo == null) {
            return;
        }

        for (int i = 0; i < connectTo.length; i++) {

            if (connectTo[i] == null) {
                continue;
            }

            //Master.log("unconnected #" + i + ": " + this + " from " + connectToNew[i]);

            connectTo[i].unconnectFrom(this);

        }

        connectTo = null;
        connectorOffTime = null;

    }

    public int numberOfConnections(boolean madeByThisVesicle){

        int count = 0;

        if (connectTo == null) {
            return 0;
        }

        if (madeByThisVesicle) {

            for (int i = 0; i < connectTo.length; i++) {

                if ((connectTo[i] != null) && (connectorOffTime[i] > 0)) {
                    count++;
                }

            }

            return count;

        }

        for (int i = 0; i < connectTo.length; i++) {

            if (connectTo[i] != null) {
                count++;
            }

        }

        return count;

    }

    public boolean unconnectFrom(DiffusantVesicle dv){

        boolean removedConnection = false;

        if (connectTo == null) {
            return false;
        }

        for (int i = 0; i < connectTo.length; i++) {

            if (connectTo[i] == dv) {
                connectTo[i] = null;
                connectorOffTime[i] = 0;
                removedConnection = true;
            }

        }

        return removedConnection;

    }

    public boolean connectToNew(DiffusantVesicle dv, MersenneTwisterFast mt, double unconnectRate, double time) {

        if (connectTo == null) {
            return false;
        }

        for (int i = 0; i < connectTo.length; i++) {

            if (connectTo[i] == dv) {
                return false; // already exists, not a new connection
            }

            if (connectTo[i] == null) {

                connectTo[i] = dv;

                if (Double.isInfinite(unconnectRate)) {
                    connectorOffTime[i] = Double.POSITIVE_INFINITY;
                } else {
                    connectorLifeTime = mt.nextPoisson(unconnectRate);
                    //lifetime = (1.0 / unconnectRate);
                    //Master.log("" + lifetime);
                    connectorOffTime[i] = time + connectorLifeTime;
                }

                //Master.log("connect #" + i + ": " + this + " to " + dv);

                return true;

            }

        }

        return false;

    }

    public boolean connectTo(DiffusantVesicle dv) {

        if (connectTo == null) {
            return false;
        }

        for (int i = 0; i < connectTo.length; i++) {

            if (connectTo[i] == dv) {
                return true; // already exists
            }

            if (connectTo[i] == null) {

                connectTo[i] = dv;
                connectorOffTime[i] = 0;

                //Master.log("connect #" + i + ": " + this + " to " + dv);

                return true;

            }

        }

        return false;

    }

    public boolean isConnected() {

        if (connectTo == null) {
            return false;
        }

        for (int i = 0; i < connectTo.length; i++) {
            if (connectTo[i] != null) {
                return true;
            }
        }

        return false;

    }

    public boolean isConnectedTo(DiffusantVesicle dv) {

        if (connectTo == null) {
            return false;
        }

        for (int i = 0; i < connectTo.length; i++) {
            if (connectTo[i] == dv) {
                return true;
            }
        }

        return false;

    }

    public boolean noMoreConnectors() {

        if (connectTo == null) {
            return true;
        }

        for (int i = 0; i < connectTo.length; i++) {

            if (connectTo[i] == null) {
                return false;
            }

        }

        return true;

    }

    public boolean testConnectorSeperation(double connectorLength) {

        if (connectTo == null) {
            return false;
        }

        for (int i = 0; i < connectTo.length; i++) {

            if (connectTo[i] == null) {
                continue;
            }

            if (!overlap(connectTo[i], connectorLength)) {
                return true;
            }

        }

        return false;

    }

    public String xyz() {
        return "x = " + x + ", y = " + y + ", z = " + z;
    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        String n = o.getName();

        if (!(o.paramVector instanceof DiffusantVesicle)) {
            return false;
        }

        if (n.equalsIgnoreCase("radius")) {
            if (v <= 0) {
                return false;
            }
            radius = v;
        }
        if (n.equalsIgnoreCase("fluorescence")) {
            if (v < 0) {
                return false;
            }
            fluorescence = v;
        }
        if (n.equalsIgnoreCase("D")) {
            if (v < 0) {
                return false;
            }
            setD(v);
        }
        if (n.equalsIgnoreCase("step3")) {
            if (v < 0) {
                return false;
            }
            setStep3(v);
        }
        if (n.equalsIgnoreCase("x")) {
            if (v <= 0) {
                return false;
            }
            x = v;
        }
        if (n.equalsIgnoreCase("y")) {
            if (v <= 0) {
                return false;
            }
            y = v;
        }
        if (n.equalsIgnoreCase("z")) {
            if (v <= 0) {
                return false;
            }
            z = v;
        }
        if (n.equalsIgnoreCase("mobile")) {
            if (v == 1) {
                mobile = true;
            } else {
                mobile = false;
            }
            return true;
        }
        if (n.equalsIgnoreCase("insideGeometry")) {
            if (v == 1) {
                insideGeometry = true;
            } else {
                insideGeometry = false;
            }
            return true;
        }
        return false;
    }

}
