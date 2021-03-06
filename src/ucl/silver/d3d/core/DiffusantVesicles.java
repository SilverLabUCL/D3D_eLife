package ucl.silver.d3d.core;

import ucl.silver.d3d.utils.MersenneTwisterFast;

import java.awt.Color;

/**
 * <p>Title: D3D</p>
 * <p>Description: 3D Diffusion-Reaction Simulator</p>
 * <p>Copyright: Copyright (c) 2016</p>
 * <p>Company: The Silver Lab at University College London</p>
 * @author Jason S. Rothman
 * @version 1.6
 */

public class DiffusantVesicles extends Diffusant {

    public double minRadius; // um
    public double maxRadius; // um
    public double meanRadius; // um
    public double setMeanRadius; // um

    public double meanVolume; // um^3
    public double totalVolume; // um^3

    public double meanFluorescence;

    public double meanStep3; // um

    public double maxX0, maxY0, maxZ0; // um

    public transient DiffusantVesicle[] vesicle; // the array of vesicles
    public transient double xyz[][] = null;

    public double setDensity; // vesicles / um^3
    public double setVolumeFraction;

    public double density; // vesicles / um^3
    public double volumeFraction;

    public int numVesicles; // total number of vesicles in array

    public int mobileVesicles, immobileVesicles;
    public double mobileVolumeFraction, immobileVolumeFraction;
    public double mobilePercent, immobilePercent; // e.g. 0.25
    public double setImmobilePercent = 0; // set immobile fraction

    public double kPhoto; // k of photolysis reaction (1/ms)
    public double kPhotoSave; // variable for saving kPhoto to output file

    public ColorD3D colorReady = new ColorD3D( "colorReady", Color.white);
    public ColorD3D colorImmobile = new ColorD3D( "colorImmobile", Color.gray);
    public ColorD3D colorConnected = new ColorD3D( "colorConnected", new Color(102,102,0));

    public boolean PV_includeVesicleArray = false;

    public boolean saveXYZ = false; // save vesicles xyz position
    public boolean saveMSD = false; // save mean square distance

    public Save save_XYZ = null;
    public Save save_MSD = null;

    public MersenneTwisterFast mt = new MersenneTwisterFast(); // create/init random number generator

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("minRadius")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("maxRadius")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("meanRadius")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("setMeanRadius")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("meanVolume")) {
            return project.spaceUnits + "^3";
        }
        if (name.equalsIgnoreCase("totalVolume")) {
            return project.spaceUnits + "^3";
        }
        if (name.equalsIgnoreCase("meanStep3")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("meanExtraX")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("meanExtraY")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("meanExtraZ")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("maxX0")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("maxY0")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("maxZ0")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("density")) {
            return "vesicles/" + project.spaceUnits + "^3";
        }
        if (name.equalsIgnoreCase("setDensity")) {
            return "vesicles/" + project.spaceUnits + "^3";
        }
        if (name.equalsIgnoreCase("kPhoto")) {
            return "1/" + project.timeUnits;
        }
        return super.units(name);
    }

    public DiffusantVesicles(Project p, String NAME, double InitialConcentration, double DiffusionConstant, CoordinatesVoxels c,
            double radius, PulseTimer pt, PSF PSF) {

        super(p, NAME, InitialConcentration, DiffusionConstant, c);

        setMeanRadius = radius;
        pulseTimer = pt;
        reaction = true;
        psf = PSF;

        if (save != null) {
            save.xdim = project.timeUnits;
            save.ydim = "kPhoto (1/" + project.timeUnits + ")";
        }

        createVector(true);

    }

    @Override
    public void init() {

        super.init();

        if (vesicle != null) {
            for (int i = 0; i < vesicle.length; i++) {
                if (vesicle[i] != null) {
                    vesicle[i].init();
                }
            }
        }

        vesicleStats();

        if (save != null) {
            save.init();
            saveFileName();
            saveDimensions();
            save.updateVectors();
        }

        if (saveXYZ) {

            if (save_XYZ == null) {
                save_XYZ = new Save(project);
                save_XYZ.samples2save = 0;
                save_XYZ.skipSamples = 0;
                save_XYZ.sampleRate = -1;
                save_XYZ.sampleInterval = 0;
                save_XYZ.saveWhileComputing = true;
                save_XYZ.save2BinaryFile = false;
                save_XYZ.save2TextFile = true;
            }

            save_XYZ.init();
            saveXYZFileName(-1);
            saveXYZDimensionsInit();
            save_XYZ.updateVectors();

        }

        if (saveMSD) {

            if (save_MSD == null) {
                save_MSD = new Save(project);
                save_MSD.saveWhileComputing = false;
                save_MSD.save2BinaryFile = true;
                save_MSD.save2TextFile = false;
            }

            save_MSD.init();
            saveMSDFileName();
            saveMSDDimensions();
            save_MSD.updateVectors();

        }

    }

    public void vesicleStats() {

        int count = 0;
        double volume;
        double spaceVolume;

        if (project.monteCarlo == null) {
            return;
        }

        if (vesicle == null) {
            return;
        }

        spaceVolume = project.monteCarlo.spaceVolume();

        minRadius = Double.POSITIVE_INFINITY;
        maxRadius = 0;
        meanRadius = 0;

        meanVolume = 0;
        totalVolume = 0;

        meanFluorescence = 0;

        D = 0;
        meanStep3 = 0;

        maxX0 = 0;
        maxY0 = 0;
        maxZ0 = 0;

        mobileVesicles = 0;
        immobileVesicles = 0;
        mobilePercent = 0;
        immobilePercent = 0;

        density = 0;
        volumeFraction = 0;
        mobileVolumeFraction = 0;
        immobileVolumeFraction = 0;

        for (int i = 0; i < vesicle.length; i++) {

            if (vesicle[i] == null) {
                continue;
            }

            if (!vesicle[i].insideGeometry) {
                continue;
            }

            count++;

            volume = vesicle[i].getVolume();

            minRadius = Math.min(vesicle[i].radius, minRadius);
            maxRadius = Math.max(vesicle[i].radius, maxRadius);

            meanRadius += vesicle[i].radius;
            meanVolume += volume;
            totalVolume += volume;
            meanFluorescence += vesicle[i].fluorescence;

            D += vesicle[i].D;
            meanStep3 += vesicle[i].step3;

            maxX0 = Math.max(maxX0, Math.abs(vesicle[i].x0));
            maxY0 = Math.max(maxY0, Math.abs(vesicle[i].y0));
            maxZ0 = Math.max(maxZ0, Math.abs(vesicle[i].z0));

            if (vesicle[i].mobile) {
                mobileVesicles++;
            } else {
                immobileVesicles++;
            }

        }

        numVesicles = immobileVesicles + mobileVesicles;

        meanRadius /= numVesicles;
        meanVolume /= numVesicles;
        meanFluorescence /= numVesicles;
        D /= numVesicles;
        meanStep3 /= numVesicles;

        immobilePercent = immobileVesicles / ( 1.0 * numVesicles);
        mobilePercent = mobileVesicles / ( 1.0 * numVesicles);

        density = (1.0 * numVesicles) / spaceVolume;
        volumeFraction = totalVolume / spaceVolume;
        mobileVolumeFraction = mobilePercent * volumeFraction;
        immobileVolumeFraction = immobilePercent * volumeFraction;

        setParamObject("minRadius", minRadius);
        setParamObject("maxRadius", maxRadius);
        setParamObject("meanRadius", meanRadius);

        setParamObject("meanVolume", meanVolume);
        setParamObject("totalVolume", totalVolume);

        setParamObject("meanFluorescence", meanFluorescence);

        setParamObject("D", D);
        setParamObject("meanStep3", meanStep3);

        setParamObject("maxX0", maxX0);
        setParamObject("maxY0", maxY0);
        setParamObject("maxZ0", maxZ0);

        setParamObject("numVesicles", numVesicles);
        setParamObject("mobileVesicles", mobileVesicles);
        setParamObject("immobileVesicles", immobileVesicles);
        setParamObject("mobileVesicleFraction", mobilePercent);
        setParamObject("immobilePercent", immobilePercent);

        setParamObject("density", density);
        setParamObject("volumeFraction", volumeFraction);
        setParamObject("mobileVolumeFraction", mobileVolumeFraction);
        setParamObject("immobileVolumeFraction", immobileVolumeFraction);

    }

    public void computeDensity() {

        double sumVolume = 0.0, count = 0.0;

        double spaceVolume = project.monteCarlo.spaceVolume();

        for (int i = 0; i < vesicle.length; i++) {

            if ((vesicle[i] != null) && vesicle[i].insideGeometry) {
                count++;
                sumVolume += vesicle[i].getVolume();
            }

        }

        density = count / spaceVolume;
        volumeFraction = sumVolume / spaceVolume;

        setParamObject("density", density);
        setParamObject("volumeFraction", volumeFraction);

    }

    public static int numVesiclesPossible(double volumeGeometry, double vesicleRadius, double vFraction) {

        double vesicleVolume = 4.0 * Math.PI * Math.pow(vesicleRadius, 3) / 3.0;

        int nVesicles = (int) (vFraction * volumeGeometry / vesicleVolume);

        return nVesicles;

    }

    @Override
    public double displayValue(int i, int j, int k) {

        if (concentration != null) {
            return concentration[i][j][k];
        }

        if (psf == null) {
            return C0;
        }

        psf.checkExists();

        return psf.getArrayValue(i, j, k);

    }

    @Override
    public boolean saveInit() {

        int dataPoints = 1;

        if (save != null) {
            if (!save.init(name, coordinates(), -1, dataPoints)) {
                return false;
            }
        }
        
        if (save_MSD != null) {
            if (!save_MSD.init(name, coordinates(), -1, dataPoints)) {
                return false;
            }
        }

        return true;

    }

    @Override
    public boolean saveFinish() {

        if (save != null) {
            save.finish(name, coordinates(), -1);
        }

        if (save_MSD != null) {
            save_MSD.finish(name, coordinates(), -1);
        }

        return true;

    }

    @Override
    public boolean save() {

        double svalue = -1;

        if (save != null) {
            save.saveData(kPhotoSave);
        }

        if (save_MSD != null) {
            if (save_MSD.skipCounter == 0) {
                svalue = meanSquareDisplacement();
            }
            save_MSD.saveData(svalue);
        }

        return true;

    }

    public boolean saveXYZFileName(int msec) {

        if (save_XYZ == null) {
            return false;
        }

        //save_XYZ.fileName("XYZF", name);
        save_XYZ.fileName("XYZ", name);

        if (msec >= 0) {
            save_XYZ.outputFile += "_" + Integer.toString(msec) + project.timeUnits;
        }

        return true;

    }

    public void saveXYZDimensionsInit() {

        if ((save_XYZ == null) || (!save_XYZ.autoDimensions)) {
            return;
        }

        save_XYZ.xdim = "XYZ positions";

        if ((name == null) || (name.length() == 0)) {
            save_XYZ.ydim = project.spaceUnits;
        } else {
            save_XYZ.ydim = name + " (" + project.spaceUnits + ")";
        }

    }

    public boolean saveXYZ(int msec) {

        double f;
        String ostr;

        int dataPoints = 1;

        if (!saveXYZ || (vesicle == null)) {
            return false;
        }

        saveXYZFileName(msec);

        if (!save_XYZ.init(name, coordinates(), msec, dataPoints)) {
            return false;
        }

        for (int i = 0; i < vesicle.length; i++) {

            if (vesicle[i] == null) {
                continue;
            }

            if (!vesicle[i].insideGeometry) {
                continue;
            }

            f = vesicle[i].fluorescence;
            //f = vesicles[i].fluorescence * vesicles[i].voxel.PSFd;
            //f = vesicles[i].localD;

            //ostr = Double.toString(vesicle[i].x) + '\t' + Double.toString(vesicle[i].y) + '\t' + Double.toString(vesicle[i].z) + '\t' + Double.toString(f);
            ostr = Double.toString(vesicle[i].x) + '\t' + Double.toString(vesicle[i].y) + '\t' + Double.toString(vesicle[i].z);

            save_XYZ.writeString(ostr);

        }

        save_XYZ.finish(name, coordinates(), msec);

        //Master.log("Saved vesicles positions at " + msec + " " + project.timeUnits);

        return false;

    }

    public boolean saveMSDFileName() {

        if (save_MSD == null) {
            return false;
        }

        save_MSD.fileName("MSD", name);

        return true;

    }

    public void saveMSDDimensions() {

        if ((save_MSD == null) || (!save_MSD.autoDimensions)) {
            return;
        }

        save_MSD.xdim = "ms";

        if ((name == null) || (name.length() == 0)) {
            save_MSD.ydim = project.spaceUnits;
        } else {
            save_MSD.ydim = name + " (" + project.spaceUnits + "^2)";
        }

    }

    public boolean checkVesicleNum(int i) {
        if ((vesicle != null) && (i >= 0) && (i < vesicle.length)) {
            return true;
        }
        return false;
    }

    public void initVesicles() {

        double spaceVolume = project.monteCarlo.spaceVolume();

        if ((setVolumeFraction > 0) && (spaceVolume > 0)) {
            numVesicles = numVesiclesPossible(spaceVolume, setMeanRadius, setVolumeFraction);
        } else if ((setDensity > 0) && (spaceVolume > 0)) {
            numVesicles = (int) (setDensity * spaceVolume);
        }

        Master.log("numVesicles = " + numVesicles);

        if (numVesicles <= 0) {
            Master.exit("DiffusantVesicles: initVesicles: bad number of vesicles: " + numVesicles);
        }

        vesicle = new DiffusantVesicle[numVesicles];

        for (int i = 0; i < vesicle.length; i++) {
            vesicle[i] = new DiffusantVesicle(project, "ready", setMeanRadius, D, Double.NaN, Double.NaN, Double.NaN);
        }

    }

    public void shuffleVesicles() {

        int r, size = vesicle.length;

        DiffusantVesicle temp;

        for (int i = 0; i < size; i++) {

            r = (int) (mt.nextDouble() * size);

            if ((i != r) && (r < size)) { // swap
                temp = vesicle[i];
                vesicle[i] = vesicle[r];
                vesicle[r] = temp;
            }

        }

    }

    public void initStartLocation() {

        for (int i = 0; i < vesicle.length; i++) {
            vesicle[i].initStartLocation();
        }

    }

    public double meanSquareDisplacement() {

        double sumSD = 0.0, count = 0.0;

        for (int i = 0; i < vesicle.length; i++) {
            if (vesicle[i].mobile && vesicle[i].insideGeometry) {
                sumSD += vesicle[i].squareDisplacement();
                count += 1.0;
            }
        }

        return sumSD / count;

    }

    public boolean initImmobileVesicles(Coordinates keepClear) {

        long immobile = Math.round(setImmobilePercent * numVesicles);

        mobileVesicles = 0;
        immobileVesicles = 0;

        for (int i = 0; i < vesicle.length; i++) {

            if (vesicle[i] == null) {
                continue;
            }

            vesicle[i].mobile = true;

            if (vesicle[i].insideGeometry) {
                mobileVesicles++;
            }

        }

        if (immobile == 0) {
            return false; // nothing to do
        }

        shuffleVesicles();

        for (int i = 0; i < vesicle.length; i++) {

            if (vesicle[i] == null) {
                continue;
            }

            if (!vesicle[i].insideGeometry) {
                continue;
            }

            if ((keepClear != null) && keepClear.isInside(vesicle[i].x, vesicle[i].y, vesicle[i].z)) {
                continue;
            }

            vesicle[i].mobile = false;
            immobileVesicles++;
            mobileVesicles--;

            if (immobileVesicles >= immobile) {
                break; // finished
            }

        }

        if (immobileVesicles != immobile) {
            Master.exit("DiffusantVesicles: initImmobileVesicles: failed to initialize immobile vesicles: " + (immobile - immobileVesicles));
        }

        immobilePercent = 1.0 * immobileVesicles / (immobileVesicles + mobileVesicles) ;

        return false;

    }

    int countVesicles(String vesicleType) {

        int count = 0;

        if (vesicleType.equalsIgnoreCase("mobile")) {

            for (int i = 0; i < vesicle.length; ++i) {
                if (vesicle[i].mobile) {
                    count++;
                }
            }

        } else if (vesicleType.equalsIgnoreCase("immobile")) {

            for (int i = 0; i < vesicle.length; ++i) {
                if (!vesicle[i].mobile) {
                    count++;
                }
            }

        } else { // "ready" or "docked" or "reserve"

            for (int i = 0; i < vesicle.length; ++i) {
                if (vesicle[i].name.equalsIgnoreCase(vesicleType)) {
                    count++;
                }
            }

        }

        return count;

    }

    public void react(int it) {

        double pt, photolysis;

        if ((pulseTimer == null) || (it >= pulseTimer.timer.length)) {
            return;
        }

        pt = pulseTimer.timer[it];
        kPhotoSave = kPhoto * pt;

        if (pt == 0) {
            return;
        }

        for (int i = 0; i < vesicle.length; ++i) {

            if (vesicle[i] == null) {
                continue;
            }

            if (vesicle[i].voxel == null) {
                continue;
            }

            if (!vesicle[i].insideGeometry) {
                continue;
            }

            if (vesicle[i].voxel.PSFi >= 0) {
                photolysis = vesicle[i].fluorescence * kPhotoSave * vesicle[i].voxel.PSFi * project.dt;
                vesicle[i].fluorescence -= photolysis;
            }

        }

        return;

    }

    @Override
    public boolean addUser(ParamVector pv) {

        super.addUser(pv);

        if (save_XYZ != null) {
            save_XYZ.addUser(pv);
        }

        if (save_MSD != null) {
            save_MSD.addUser(pv);
        }

        if (colorReady != null) {
            colorReady.addUser(pv);
        }

        if (colorImmobile != null) {
            colorImmobile.addUser(pv);
        }

        if (colorConnected != null) {
            colorConnected.addUser(pv);
        }

        if (PV_includeVesicleArray && (vesicle != null)) {
            for (int i = 0; i < vesicle.length; i++) {
                if (vesicle[i] != null) {
                    vesicle[i].addUser(pv);
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

        if (save_XYZ != null) {
            addBlankParam();
            save_XYZ.createVector(true);
            addVector(save_XYZ.getVector());
            save_XYZ.addUser(this);
        }

        if (save_MSD != null) {
            addBlankParam();
            save_MSD.createVector(true);
            addVector(save_MSD.getVector());
            save_MSD.addUser(this);
        }

        if (colorReady != null) {
            addBlankParam();
            colorReady.createVector(true);
            addVector(colorReady.getVector());
            colorReady.addUser(this);
        }

        if (colorImmobile != null) {
            addBlankParam();
            colorImmobile.createVector(true);
            addVector(colorImmobile.getVector());
            colorImmobile.addUser(this);
        }

        if (colorConnected != null) {
            addBlankParam();
            colorConnected.createVector(true);
            addVector(colorConnected.getVector());
            colorConnected.addUser(this);
        }

        if (PV_includeVesicleArray && (vesicle != null)) {
            for (int i = 0; i < vesicle.length; i++) {
                if (vesicle[i] != null) {
                    addBlankParam();
                    vesicle[i].createVector(true);
                    addVector(vesicle[i].getVector());
                    vesicle[i].addUser(this);
                }
            }
        }

        if (close) {
            closeVector();
        }

        return true;

    }

    @Override
    public void clearVector() {

        super.clearVector();

        if (PV_includeVesicleArray && (vesicle != null)) {
            for (int i = 0; i < vesicle.length; i++) {
                if (vesicle[i] != null) {
                    vesicle[i].clearVector();
                }
            }
        }

    }

    @Override
    public void updateVector(ParamObject[] v) {
        
        super.updateVector(v);

        if (save_XYZ != null) {
            save_XYZ.updateVector(v);
        }

        if (save_MSD != null) {
            save_MSD.updateVector(v);
        }

        if (colorReady != null) {
            colorReady.updateVector(v);
        }

        if (colorImmobile != null) {
            colorImmobile.updateVector(v);
        }

        if (colorConnected != null) {
            colorConnected.updateVector(v);
        }

        if (PV_includeVesicleArray && (vesicle != null)) {
            for (int i = 0; i < vesicle.length; i++) {
                if (vesicle[i] != null) {
                    vesicle[i].updateVector(v);
                }
            }
        }

    }

    public boolean setVesicles(String varName, double value) {

        boolean OK = false;
        boolean atLeastOne = false;

        if (vesicle == null) {
            return false;
        }

        for (int i = 0; i < vesicle.length; i++) {
            if (vesicle[i] != null) {
                OK = vesicle[i].set(varName, value);
                if (OK) {
                    atLeastOne = true;
                }
            }
        }

        vesicleStats();

        return atLeastOne;

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {
        
        String n = o.getName();

        if (o == null) {
            return false;
        }

        if (!(o.paramVector instanceof DiffusantVesicles)) {
            return false;
        }

        if (super.setMyParams(o, v)) {
            return true;
        }

        if (n.equalsIgnoreCase("setMeanRadius")) {
            if (v <= 0) {
                return false;
            }
            return setVesicles("radius",v);
        }
        if (n.equalsIgnoreCase("meanFluorescence")) {
            if (v < 0) {
                return false;
            }
            return setVesicles("fluorescence",v);
        }
        if (n.equalsIgnoreCase("D")) {
            if (v < 0) {
                return false;
            }
            D = v;
            return setVesicles("D",v);
        }
        if (n.equalsIgnoreCase("meanStep3")) {
            if (v < 0) {
                return false;
            }
            return setVesicles("step3",v);
        }
        if (n.equalsIgnoreCase("meanExtraX")) {
            return setVesicles("extraX",v);
        }
        if (n.equalsIgnoreCase("meanExtraY")) {
            return setVesicles("extraY",v);
        }
        if (n.equalsIgnoreCase("meanExtraZ")) {
            return setVesicles("extraZ",v);
        }
        if (n.equalsIgnoreCase("meanDockRefractoryPeriod")) {
            if (v < 0) {
                return false;
            }
            return setVesicles("dockRefractoryPeriod",v);
        }
        if (n.equalsIgnoreCase("setDensity")) {
            if (v < 0) {
                return false;
            }
            setDensity = v;
            return true;
        }
        if (n.equalsIgnoreCase("setVolumeFraction")) {
            if (v < 0) {
                return false;
            }
            setVolumeFraction = v;
            return true;
        }
        if (n.equalsIgnoreCase("kPhoto")) {
            if (v < 0) {
                return false;
            }
            kPhoto = v;
            return true;
        }
        if (n.equalsIgnoreCase("setImmobilePercent")) {
            if (v < 0) {
                return false;
            }
            setImmobilePercent = v;
            return true;
        }
        if (n.equalsIgnoreCase("savePositions")) {
            if (v == 1) {
                saveXYZ = true;
            } else {
                saveXYZ = false;
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

        if (!(o.paramVector instanceof DiffusantVesicles)) {
            return false;
        }

        if (super.setMyParams(o, s)) {
            return true;
        }

        return false;
    }

}
