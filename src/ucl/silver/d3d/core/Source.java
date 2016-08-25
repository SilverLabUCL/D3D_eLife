package ucl.silver.d3d.core;

import java.awt.Color;
import java.io.*;

/**
 * <p>Title: D3D</p>
 * <p>Description: 3D Diffusion-Reaction Simulator</p>
 * <p>Copyright: Copyright (c) 2016</p>
 * <p>Company: The Silver Lab at University College London</p>
 * @author Jason S. Rothman
 * @version 1.6
 */

public class Source extends ParamVector {

    public String diffusantName = null;
    public int diffusantNum = -1; // diffusant array number in Project

    private int initFlag = -1; // compute init values with (0) C0 (1) molecules (2) Qvoxel (3) Qtotal
    public double C0 = 0; // source concentration per voxel (mM)
    public double molecules = 0; // number of molecules per voxel
    public double I0 = 0; // current per voxel (pA)
    public double Qvoxel = 0; // flux per voxel (mol/ms)
    public double Qtotal = 0; // flux (mol/ms)
    public int valence = 0; // valence of ion, for conversion of flux to pA

    public double krate = 0; // rate of release/uptake (1/ms)

    public boolean clamp = false; // source "clamps" voxels to C0
    public boolean writeCurrent = false;
    public boolean useGeometryCoordinates = true;
    public boolean useProjectArray = false;
    private boolean useArray = false;

    public double C0_conversionFactor = 1, I0_conversionFactor = 1;

    private Color defaultColor = new Color(153, 0, 51);

    public ColorD3D color = null;
    CoordinatesVoxels coordinates = null; // coordinates of source
    public Save save = null;
    public PulseTimer pulseTimer = null; // timer for source

    private RandomAccessFile raFile = null;
    public double[] array = null;

    private final double avogadro = 6.0221415E23;
    private final double farady = 96500;

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("c0")) {
            return project.concUnits;
        }
        if (name.equalsIgnoreCase("I0")) {
            return "pA";
        }
        if (name.equalsIgnoreCase("Qvoxel")) {
            return "mol/" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("Qtotal")) {
            return "mol/" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("krate")) {
            return "1/" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("fileName")) {
            return "DIR";
        }
        return super.units(name);
    }

    @Override
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("c0_conversionFactor")) {
            return false;
        }
        if (name.equalsIgnoreCase("i0_conversionFactor")) {
            return false;
        }
        return true;
    }

    public Source(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c, double concentration, PulseTimer pt) {

        super(p);

        if (NAME.length() > 0) {
            name = NAME;
        } else {
            name = "Source";
        }

        color = new ColorD3D(name + "_color", defaultColor);

        diffusantNum = DiffusantNum;
        initFlag = 0;
        C0 = concentration;
        pulseTimer = pt;
        save = new Save(p);

        if (c == null) {
            useGeometryCoordinates = true;
            coordinates = null;
        } else {
            useGeometryCoordinates = false;
            coordinates = new CoordinatesVoxels(p);
            coordinates.matchVoxels(c);
        }

        createVector(true); // sets ParamVector

    }

    public Source(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c, PulseTimer pt, int Molecules) {

        super(p);

        if (NAME.length() > 0) {
            name = NAME;
        } else {
            name = "Source";
        }

        color = new ColorD3D(name + "_color", defaultColor);

        diffusantNum = DiffusantNum;
        initFlag = 1;
        molecules = Molecules;
        pulseTimer = pt;
        save = new Save(p);

        if (c == null) {
            useGeometryCoordinates = true;
            coordinates = null;
        } else {
            useGeometryCoordinates = false;
            coordinates = new CoordinatesVoxels(p);
            coordinates.matchVoxels(c);
        }

        createVector(true); // sets ParamVector

    }

    public Source(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c, PulseTimer pt, double qtotal) {

        super(p);

        if (NAME.length() > 0) {
            name = NAME;
        } else {
            name = "Source";
        }

        color = new ColorD3D(name + "_color", defaultColor);

        diffusantNum = DiffusantNum;
        initFlag = 3;
        Qtotal = qtotal;
        pulseTimer = pt;
        save = new Save(p);

        if (c == null) {
            useGeometryCoordinates = true;
            coordinates = null;
        } else {
            useGeometryCoordinates = false;
            coordinates = new CoordinatesVoxels(p);
            coordinates.matchVoxels(c);
        }

        createVector(true); // sets ParamVector

    }

    public CoordinatesVoxels coordinates() {
        if ((useGeometryCoordinates) || (coordinates == null)) {
            return project.geometry;
        } else  {
            return coordinates;
        }
    }

    public void updateCoordinates() {
        // nothing to do
    }

    public boolean SetDiffusantNum(int arrayNum) {

        Diffusant d = project.getDiffusant(arrayNum);

        if (d != null) {
            diffusantNum = arrayNum;
            diffusantName = d.name;
            updateVectors();
            return true;
        } else {
            Master.log("Source SetDiffusantNum error: failed to find diffusant #" + arrayNum);
            updateVectors();
            return false;
        }

    }

    public boolean SetDiffusantName(String dName) {

        int i = project.getDiffusantNum(dName);

        if (i >= 0) {
            diffusantName = dName;
            diffusantNum = i;
            updateVectors();
            return true;
        } else {
            Master.log("Source SetDiffusantName error: failed to find diffusant " + dName);
            updateVectors();
            return false;
        }

    }

    public void init() {

        double molesPerVoxel, litersPerVoxel;

        if (save != null) {
            save.init();
            saveFileName();
            saveDimensions();
        }

        setParamError("diffusantNum", null);

        diffusantName = project.diffusantName(diffusantNum);

        if (!project.checkDiffusantNum(diffusantNum)) {
            error("init", "diffusantNum", "out of range");
        }

        coordinates().update();

        litersPerVoxel = Math.pow(project.dx, 3) * 1e-18 * 1000;

        switch (initFlag) {

            case 0: // start with C0 (mM)
                molesPerVoxel = C0 * litersPerVoxel * 1e-3;
                molecules = molesPerVoxel * avogadro;
                Qvoxel = molesPerVoxel / project.dt; // mols/msec
                Qtotal = Qvoxel * coordinates().voxels; // mols/msec
                break;

            case 1: // start with molecules
                molesPerVoxel = molecules / avogadro;
                C0 = 1e3 * molesPerVoxel / litersPerVoxel; // mM
                Qvoxel = molesPerVoxel / project.dt; // mols/msec
                Qtotal = Qvoxel * coordinates().voxels; // mols/msec
                break;

            case 2: // start with Qvoxel
                Qtotal = Qvoxel * coordinates().voxels;
                molesPerVoxel = Qvoxel * project.dt;
                molecules = molesPerVoxel * avogadro;
                C0 = 1e3 * molesPerVoxel / litersPerVoxel; // mM
                break;

            case 3: // start with Qtotal
                Qvoxel = Qtotal / coordinates().voxels; // spread flux evenly amongst voxels
                molesPerVoxel = Qvoxel * project.dt;
                molecules = molesPerVoxel * avogadro;
                C0 = 1e3 * molesPerVoxel / litersPerVoxel; // mM
                //Master.log("Qvoxel=" + Qvoxel);
                break;

            case 4: // import from file which is in units of pA

                // compute factor for converting pA/totalVoxels to mM/voxel

                Qtotal = 1 / (valence * farady * 1e3 * 1e12); // total flux
                Qvoxel = Qtotal / coordinates().voxels; // flux per voxel
                molesPerVoxel = Qvoxel * project.dt;
                C0_conversionFactor = 1.0e3 * molesPerVoxel / litersPerVoxel; // mM/voxel

                // compute factor for converting mM/voxel to pA/voxel

                molesPerVoxel = 1 * litersPerVoxel * 1e-3;
                Qvoxel = molesPerVoxel / project.dt; // mols/msec
                I0_conversionFactor = Qvoxel * valence * farady * 1e3 * 1e12; // convert molar flux to pA

                break;

            default:

                error("init: bad initFlag: " + initFlag);

        }

        if ((valence > 0) && (initFlag < 4)) {
            I0 = Qvoxel * valence * farady * 1.0e3; // convert molar flux to Amps
            I0 *= 1.0e12; // amps to pA
        }

        if (C0 < 0) {
            error("init: negative C0: " + C0);
        }

        updateVectors();

    }

    public void initPulseTimer() {
        if (pulseTimer != null) {
            pulseTimer.initTimer();
        }
    }

    public double release(RunFiniteDifference fd, Geometry geometry) {

        int index, numVoxels;
        double c0, i0, c1, c2, newconc, avgC = 0, avgI = 0;
        boolean timerHigh = false;

        if (coordinates().spaceVoxels == 0) {
            return 0;
        }

        if ((fd.diffus == null) || (fd.it >= fd.itmax)) {
            return 0;
        }

        if ((pulseTimer == null) || (pulseTimer.timer == null)) {
            timerHigh = true; // on for all time
        } else if (pulseTimer.timer[fd.it] > 0) {
            timerHigh = true;
        }

        //c2 = C0 * timer;
        //current = I0 * timer;

        if (useArray) {

            if (useProjectArray) {

                if (fd.it < project.sourceArray.length) {

                    c0 = project.sourceArray[(int) fd.it];

                    if (c0 < 0) {
                        c0 = 0;
                    } else {
                        c0 *= C0_conversionFactor;
                    }

                    i0 = c0 * I0_conversionFactor;

                } else {
                    c0 = 0;
                    i0 = 0;
                }

            } else {

                if (fd.it < array.length) {

                    c0 = array[(int) fd.it];

                    if (c0 < 0) {
                        c0 = 0;
                    } else {
                        c0 *= C0_conversionFactor;
                    }

                    i0 = c0 * I0_conversionFactor;

                } else {
                    c0 = 0;
                    i0 = 0;
                }

            }

        } else {
            c0 = C0;
            i0 = I0;
        }

        if (timerHigh) {

            if (fd.diffus[0].length == 1) { // single compartment
                
                if (clamp) {

                    c2 = c0 - fd.diffus[diffusantNum][0];
                    fd.diffus[diffusantNum][0] = c0;

                    avgC = c2;

                } else if (krate > 0) {

                    c1 = fd.diffus[diffusantNum][0];
                    c2 = c1 * krate * project.dt;
                    newconc = c1 + c2;

                    if (newconc >= 0) {
                        fd.diffus[diffusantNum][0] = newconc;
                        avgC = c2;
                    }

                } else if (krate < 0) {

                    c1 = fd.diffus[diffusantNum][0];
                    c2 = (c1 - c0) * krate * project.dt;
                    newconc = c1 + c2;

                    if (newconc >= c0) {
                        fd.diffus[diffusantNum][0] = newconc;
                        avgC = c2;
                    }

                } else {

                    fd.diffus[diffusantNum][0] += c0;
                    avgC = c0;
                    avgI = i0;

                }

            } else {

                if (coordinates().index == null) {
                    //Master.exit("Source error: coordinates index array has not be initialized.");
                    return 0;
                }

                numVoxels = coordinates().index.length;

                for (int i = 0; i < numVoxels; i++) {

                    index = coordinates().index[i];

                    if (clamp) {

                        c2 = c0 - fd.diffus[diffusantNum][index];
                        fd.diffus[diffusantNum][index] = c0;

                        avgC += c2;

                    } else if (krate > 0) {

                        c1 = fd.diffus[diffusantNum][index];
                        c2 = c1 * krate * project.dt;
                        newconc = c1 + c2;

                        if (newconc >= 0) {
                            fd.diffus[diffusantNum][index] = newconc;
                            avgC += c2;
                        }

                    } else if (krate < 0) {

                        c1 = fd.diffus[diffusantNum][index];
                        c2 = (c1 - c0) * krate * project.dt;
                        newconc = c1 + c2;

                        if (newconc >= c0) {
                            fd.diffus[diffusantNum][index] = newconc;
                            avgC += c2;
                        }

                    } else {

                        fd.diffus[diffusantNum][index] += c0;
                        avgC += c0;
                        avgI += i0;

                    }

                }

                if (numVoxels > 1) {
                    avgC /= 1.0 * numVoxels;
                    avgI /= 1.0 * numVoxels;
                }

            }

        }

        if (writeCurrent) {
            save.saveData(avgI);
        } else {
            save.saveData(avgC);
        }

        return avgC;

    }

    public double displayValue(int xVoxel, int yVoxel, int zVoxel) {
        if (coordinates().isInside(xVoxel, yVoxel, zVoxel)) {
            return C0;
        }
        return -1;
    }

    public Color displayColor(int xVoxel, int yVoxel, int zVoxel, double min, double max) {

        if (color == null) {
            return null;
        }

        if (!coordinates().isInside(xVoxel, yVoxel, zVoxel)) {
            return null;
        }

        return color.getColor(C0, min, max);

    }

    public boolean check() {

        boolean ok = true;

        if (coordinates().voxels != coordinates().spaceVoxels) {
            Master.log("Warning: " + name + " is located at non-space voxel(s)");
        }
        if (useProjectArray && project.sourceArray != null) {
            useArray = true;
        } else if (array != null) {
            useArray = true;
        } else {
            useArray = false;
        }
        return ok;
    }

    

    public boolean saveFileName(){

        String dname = "";

        Diffusant[] diffusants = project.diffusants;

        if (save == null) {
            return false;
        }

        if ((diffusantNum >= 0) && (diffusantNum < diffusants.length)) {
            dname = diffusants[diffusantNum].name;
        }

        save.fileName(name, dname);

        return true;

    }

    public void saveDimensions() {

        if ((save == null) || (!save.autoDimensions)) {
            return;
        }

        save.xdim = project.timeUnits;

        if ((diffusantName == null) || (diffusantName.length() == 0)) {
            save.ydim = project.concUnits;
        } else {
            save.ydim = diffusantName + " (" + project.concUnits + ")";
        }

    }

    public boolean saveInit() {

        int dataPoints = 1;

        if (save == null) {
            return false;
        }

        return save.init(name, coordinates, -1, dataPoints);

    }

    public boolean saveFinish() {

        if (save == null) {
            return false;
        }

        return save.finish(name, coordinates, -1);

    }

    @Override
    public boolean addUser(ParamVector pv) {

        super.addUser(pv);

        if (color != null) {
            color.addUser(pv);
        }

        if (coordinates != null) {
            coordinates().addUser(pv);
        }

        if (save != null) {
            save.addUser(pv);
        }

        if (pulseTimer != null) {
            pulseTimer.addUser(pv);
        }

        return true;

    }

    @Override
    public boolean createVector(boolean close) {

        if (!super.createVector(false)) {
            return false;
        }

        if (color != null) {
            addBlankParam();
            color.createVector(true);
            addVector(color.getVector());
            color.addUser(this);
        }

        if (coordinates != null) {
            addBlankParam();
            coordinates.createVector(true);
            addVector(coordinates.getVector());
            coordinates.addUser(this);
        }

        if (save != null) {
            addBlankParam();
            save.createVector(true);
            addVector(save.getVector());
            save.addUser(this);
        }

        if (pulseTimer != null) {
            addBlankParam();
            pulseTimer.createVector(true);
            addVector(pulseTimer.getVector());
            pulseTimer.addUser(this);
        }

        if (close){
            this.closeVector();
        }

        return true;

    }

    @Override
    public void updateVector(ParamObject[] v) {

        super.updateVector(v);

        if (color != null) {
            color.updateVector(v);
        }

        if (coordinates != null) {
            coordinates.updateVector(v);
        }

        if (save != null) {
            save.updateVector(v);
        }

        if (pulseTimer != null) {
            pulseTimer.updateVector(v);
        }

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        String n = o.getName();

        if ((o.paramVector instanceof CoordinatesVoxels) && (coordinates != null)) {
            if (coordinates.setMyParams(o, v)) {
                updateCoordinates();
                return true;
            }
            return false;
        }

        if (!(o.paramVector instanceof Source)) {
            return false;
        }

        if (n.equalsIgnoreCase("diffusantNum")) {
            if (v < 0) {
                return false;
            }
            return SetDiffusantNum((int) v);
        }
        if (n.equalsIgnoreCase("krate")) {
            krate = v;
            return true;
        }
        if (n.equalsIgnoreCase("C0")) {
            if (v < 0) {
                return false;
            }
            C0 = v;
            initFlag = 0;
            init();
            return true;
        }
        if (n.equalsIgnoreCase("molecules")) {
            if (v < 0) {
                return false;
            }
            molecules = (int) v;
            initFlag = 1;
            init();
            return true;
        }
        if (n.equalsIgnoreCase("Qvoxel")) {
            if (v < 0) {
                return false;
            }
            Qvoxel = v;
            initFlag = 2;
            init();
            return true;
        }
        if (n.equalsIgnoreCase("Qtotal")) {
            if (v < 0) {
                return false;
            }
            Qtotal = v;
            initFlag = 3;
            init();
            return true;
        }
        if (n.equalsIgnoreCase("valence")) {
            if (v < 0) {
                return false;
            }
            valence = (int) v;
            return true;
        }
        if (n.equalsIgnoreCase("writeCurrent")) {
            if (v == 1) {
                writeCurrent = true;
            } else {
                writeCurrent = false;
            }
            return true;
        }
        if (n.equalsIgnoreCase("clamp")) {
            if (v == 1) {
                clamp = true;
            } else {
                clamp = false;
            }
            return true;
        }
        return false;
    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, String s) {

        if ((o == null) || (s == null)) {
            return false;
        }

        String n = o.getName();

        if ((o.paramVector instanceof CoordinatesVoxels) && (coordinates != null)) {
            if (coordinates.setMyParams(o, s)) {
                updateCoordinates();
                return true;
            }
            return false;
        }

        if (!(o.paramVector instanceof Source)) {
            return false;
        }

        if (n.equalsIgnoreCase("name")) {
            name = s;
            return true;
        }
        if (n.equalsIgnoreCase("diffusantName")) {
            return SetDiffusantName(s);
        }
        
        return false;
    }
}
