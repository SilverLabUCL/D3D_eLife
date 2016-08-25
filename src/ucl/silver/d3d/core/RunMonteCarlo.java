package ucl.silver.d3d.core;

import ucl.silver.d3d.gui.*;
import ucl.silver.d3d.utils.*;

/**
 * <p>Title: D3D</p>
 * <p>Description: 3D Diffusion-Reaction Simulator</p>
 * <p>Copyright: Copyright (c) 2016</p>
 * <p>Company: The Silver Lab at University College London</p>
 * @author Jason S. Rothman
 * @version 1.6
 */

public class RunMonteCarlo
        extends ParamVector implements Runnable {

    public double minVesicleRadius; // um
    public double maxVesicleRadius; // um
    public double minAllowedDistanceBetweenVesicles; // um
    public double minDistanceBetweenVesicles; // um

    public int numVesicles; // total number of vesicles
    public double vesicleDensity; // vesicles / um^3
    public double vesicleVolumeFraction;
    public double totalVesicleVolume; // um^3

    public boolean initVesicleRandom = true;

    public double Dcyto = 0; // um^2/ms

    public double voxelVolume = 0; // um^3

    public int mobileVesicles = 0;
    public int immobileVesicles = 0;
    public double mobilePercent = 0;
    public double immobilePercent = 0;
    public double mobileVolumeFraction = 0;
    public double immobileVolumeFraction = 0;

    public boolean connectVesicles = false;
    public boolean connectorsBinomial = false;
    public int maxNumConnectors = 5;
    public double meanNumConnectors = 1.5;
    public double connectorLength = 0.01; // um
    public double connectRate = 0.7/0.1; // P/ms
    public double unconnectRate = 0.175/0.1; // P/ms
    public double avgConnectorLifeTime = 0;
    public int connectorLifeTimeCounter = 0;

    public double minVesicleStep; // um
    public double removeVesicleOverlapStep3 = 0.001/Math.sqrt(3); // um

    public boolean PBC = false; // periodic boundary conditions
    public boolean freeDiffusion = false; // vesicles are allowed to overlap
    public boolean removeVesicleOverlap = true;

    public boolean frapOn = false;
    public boolean saveFluorescence = false;

    public boolean driftOn = false;
    public double driftRateX = 0.0046/1000.0; // um/ms
    public double driftRateY = 0.0046/1000.0; // um/ms
    public double driftRateZ = 0.0116/1000.0; // um/ms
    public double driftOnset = 0; // ms
    private double driftDX, driftDY, driftDZ; // um

    public Coordinates[] mito = null; // mito coordinates
    public boolean mitoCoordinatesOn = false;
    public double mitoRadius = (0.25 / 0.89) / 2.0; // um Palay
    public double mitoAxialRatio = 2.0 / 0.25; // Palay
    public double setMitoVolumeFraction = 0.28; // Zoltan average
    public double mitoVolumeFraction = 0;
    public double mitoVolume = 0; // um^3
    public String mitoShape = "cylinder";
    public int mito_ijkselect = -2;

    // simulation variables

    public transient Grid grid = null;
    public transient DiffusantVesicles[] diffusant = null;
    public transient DetectorPSF[] detector = null;
    public transient Thread thread = null;
    public transient Geometry geometry;

    public transient PSF PSFi = null; // illumination PSF
    public transient PSF PSFd = null; // detection PSF

    private int iPSF = -1;
    private int dPSF = -1;

    public int itime;
    public double time;
    public boolean timer = true;
    public double stepx, stepy, stepz;

    public boolean batchesExist = false;
    public boolean autoInit = true;
    public boolean runSimulation = false;
    public boolean cancel = false;
    public boolean initialized = false;
    public boolean preview = false;
    public boolean removingVesicleOverlap = false;

    public boolean hydrodynamicsLocalD = false; // hydrodynamic interactions local density computation
    public boolean hydroWallZ = false;
    public boolean hydrodynamicsLocalDVoxels = false;

    public int overlapTrialLimit = (int) 1e4;
    public int absLimit = (int) 1e6;

    public transient StopWatch timer1 = new StopWatch();
    public transient StopWatch timer2 = new StopWatch();

    public transient MersenneTwisterFast mt = new MersenneTwisterFast(); // create/init random number generator

    private double saveRanGauss = 999999; // for Gaussian random number generator

    public boolean saveMSDspatial = false;
    private Save[] MSDspatial = null;
    public double MSDspatialTbgn = -1;
    private double MSDspatial_t1 = -1, MSDspatial_t2 = -1;
    public double MSDspatial_win = 100; // ms
    public int MSDspatial_numBins = 9;
    public double MSDspatial_binWidth = 0.05; // um

    private final DiffusantVesicle testVesicle = new DiffusantVesicle(project, "ready", 0, 0, 0, 0, 0);

    public double vesicleVolume = Double.NaN;

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("minVesicleRadius")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("maxVesicleRadius")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("minAllowedDistanceBetweenVesicles")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("minDistanceBetweenVesicles")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("minVesicleStep")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("totalVesicleVolume")) {
            return project.volumeUnits;
        }
        if (name.equalsIgnoreCase("vesicleDensity")) {
            return "vesicles/" + project.volumeUnits;
        }
        if (name.equalsIgnoreCase("Dcyto")) {
            return project.diffusionUnits;
        }
        if (name.equalsIgnoreCase("connectorLength")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("connectRate")) {
            return "1/" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("unconnectRate")) {
            return "1/" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("avgConnectorLifeTime")) {
            return project.timeUnits;
        }
        if (name.equalsIgnoreCase("mitoRadius")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("mitoVolume")) {
            return project.volumeUnits;
        }
        if (name.equalsIgnoreCase("driftRateX")) {
            return project.spaceUnits + "/" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("driftRateY")) {
            return project.spaceUnits + "/" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("driftRateZ")) {
            return project.spaceUnits + "/" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("driftOnset")) {
            return project.timeUnits;
        }
        if (name.equalsIgnoreCase("time")) {
            return project.timeUnits;
        }
        if (name.equalsIgnoreCase("stepx")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("stepy")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("stepz")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("MSDspatialTbgn")) {
            return project.timeUnits;
        }
        if (name.equalsIgnoreCase("MSDspatial_win")) {
            return project.timeUnits;
        }
        if (name.equalsIgnoreCase("MSDspatial_binWidth")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("vesicleVolume")) {
            return project.volumeUnits;
        }
        return super.units(name);
    }

    @Override
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("minVesicleRadius")) {
            return false;
        }
        if (name.equalsIgnoreCase("maxVesicleRadius")) {
            return false;
        }
        if (name.equalsIgnoreCase("minAllowedDistanceBetweenVesicles")) {
            return false;
        }
        if (name.equalsIgnoreCase("minDistanceBetweenVesicles")) {
            return false;
        }
        if (name.equalsIgnoreCase("totalVesicleVolume")) {
            return false;
        }
        if (name.equalsIgnoreCase("VesicleDensity")) {
            return false;
        }
        if (name.equalsIgnoreCase("vesicleVolumeFraction")) {
            return false;
        }
        if (name.equalsIgnoreCase("numVesicles")) {
            return false;
        }
        if (name.equalsIgnoreCase("mobileVesicles")) {
            return false;
        }
        if (name.equalsIgnoreCase("immobileVesicles")) {
            return false;
        }
        if (name.equalsIgnoreCase("mobilePercent")) {
            return false;
        }
        if (name.equalsIgnoreCase("immobilePercent")) {
            return false;
        }
        if (name.equalsIgnoreCase("mobileVolumeFraction")) {
            return false;
        }
        if (name.equalsIgnoreCase("immobileVolumeFraction")) {
            return false;
        }
        if (name.equalsIgnoreCase("mitoVolumeFraction")) {
            return false;
        }
        if (name.equalsIgnoreCase("mitoVolume")) {
            return false;
        }
        if (name.equalsIgnoreCase("time")) {
            return false;
        }
        if (name.equalsIgnoreCase("batchesExist")) {
            return false;
        }
        if (name.equalsIgnoreCase("runSimulation")) {
            return false;
        }
        if (name.equalsIgnoreCase("cancel")) {
            return false;
        }
        if (name.equalsIgnoreCase("initialized")) {
            return false;
        }
        if (name.equalsIgnoreCase("preview")) {
            return false;
        }
        return super.canEdit(name);
    }

    public RunMonteCarlo(Project p) {

        super(p);

        createVector(true);

    }

    public void init() {

        if (project.numBatches() > 0) {
            batchesExist = true;
        }

    }

    public void setOutputRate(double newRate) {

    }

    public void startSimulation(boolean PREVIEW) {

        preview = PREVIEW;

        grid = Master.grid();

        if ((grid != null) && preview) {
            grid.preview = true;
            //voxelGrid.diffusant = diffus;
            Master.mainframe.panel2D.displayModeSet("Diffusant.0"); // set display to "Diffusant"
            Master.mainframe.panel2D.displayModeUpdate();
        }

        runSimulation = true;
        cancel = false;
        thread = new Thread(this);
        thread.start(); // this calls run() below

    }

    public void pauseSimulation() {

        runSimulation = false;

        timer1.stop();
        timer2.stop();

        Master.log("paused Monte Carlo simulation at time " + time + " " + project.timeUnits);

    }

    public void cancelSimulation() {
        runSimulation = false;
        cancel = true;
        Master.log("cancelled Monte Carlo simulation at time " + time + " " + project.timeUnits);
    }

    public boolean initSimulation() {

        if (checkVariables()) {
            return true;
        }

        if (autoInit) {
            if (initAll()) {
                return true;
            }
        }

        itime = 0;
        time = 0;

        if (driftOn) {
            driftDX = driftRateX * project.dt;
            driftDY = driftRateY * project.dt;
            driftDZ = driftRateZ * project.dt;
        }

        timer1.start();

        initialized = true;

        printParameters();

        initSave();

        Master.log("initialized Monte Carlo simulation");

        return false;

    }

    public void printParameters() {

        double voxels = geometry.voxels;
        double spaceVoxels = geometry.spaceVoxels;

        Master.log("Geometry xWidth = " + geometry.xWidth);
        Master.log("Geometry yWidth = " + geometry.yWidth);
        Master.log("Geometry zWidth = " + geometry.zWidth);

        Master.log("Geometry % nonspace = " + ( 100 * (voxels - spaceVoxels) / voxels ) );

        for (int i = 0; i < diffusant.length; i++) {

            Master.log("DiffusantVesicle #" + i + ", avg D = " + diffusant[i].D);
            Master.log("DiffusantVesicle #" + i + ", estimated meanStep = " + Math.sqrt(6.0 * diffusant[i].D * project.dt));
            Master.log("DiffusantVesicle #" + i + ", avg meanStep = " + (diffusant[i].meanStep3 * Math.sqrt(3)));
            Master.log("DiffusantVesicle #" + i + ", volume fraction = " + diffusant[i].volumeFraction);
            Master.log("DiffusantVesicle #" + i + ", immobilePercent = " + diffusant[i].immobilePercent);

            //for (int jmito = 0; jmito < diffusants[imito].vesicles.length; jmito++) {
            //    Master.log("DiffusantVesicle #" + jmito + ", D = " + diffusants[imito].vesicles[jmito].step);
            //}

        }

    }

    public void finishSimulation(){
        
        finishSave();
        vesicleStats();

        if (connectVesicles) {
            avgConnectorLifeTime /= connectorLifeTimeCounter;
            Master.log("average connector life time (ms): " + avgConnectorLifeTime + " (n=" + Integer.toString(connectorLifeTimeCounter) + ")");
        }

    }

    public void run() {

        double simTime = project.simTime;

        while (runSimulation) { // runSimulation thru batches

            if (!initialized) {

                if (project.simulationInit(preview)) {
                    cancelSimulation(); // ERROR
                }

                saveVesiclePositions(0);

            }

            if (time < simTime) {
                Master.log("starting Monte Carlo simulation at time " + time + " " + project.timeUnits);
            }

            while (runSimulation && !cancel && (time < simTime)) {

                for (int s = 0; s < diffusant.length; s++) {
                    diffusant[s].save(); // save diffusant variables
                }

                runMSDspatial();

                if (frapOn) {

                    detectFluorescence();

                    for (int k = 0; k < diffusant.length; k++) {
                        diffusant[k].react(itime);
                    }

                }

                moveVesiclesCichocki(false);

                if (hydrodynamicsLocalD) {
                    localDensityAll(false);
                }

                if (connectVesicles) {
                    connectVesicles();
                    unconnectVesicles();
                }

                if (driftOn && (time > driftOnset)) {
                    drift();
                }

                if ((grid != null) && preview) {
                    grid.repaint();
                }

                itime += 1;
                time += project.dt;

                if (timer) {
                    timer2.timer(time);
                }

            }

            if (cancel || (time >= simTime)) {
                finish();
            }

        }

    }

    public void finish() {

        if (!initialized) {
            return;
        }

        project.simulationFinish();

        saveVesiclePositions((int) project.simTime);

        //avgFirstCollision();

        if (grid != null) {
            grid.repaint();
        }

        if (batchesExist) {
            initialized = false;
        } else {
            runSimulation = false;
        }

        if ((grid != null) && (!runSimulation)) {
            grid.preview = false;
        }

        timer1.stop();

        Master.log("finished Monte Carlo simulation. time = " + timer1.toString());

        if (timer) {
            timer2.stop();
        }

    }

    public boolean initAll(){

        if (checkVariables()) {
            return true;
        }

        if (initDiffusantVesiclesArray()) {
            return true;
        }

        if (initDetectors()){
            return true;
        }

        if (initPSFs()) {
            return true;
        }

        if (initVoxels()) {
            return true;
        }

        if (initMitochondria()) {
            return true;
        }

        if (initVesicles()) {
            return true;
        }

        if (checkDX()) {
            return true;
        }

        if (initDT()) {
            return true;
        }

        if (initVesicleRandom) {
            if (initVesiclesRandom()) {
                return true;
            }
        }

        if (initVesiclesImmobile()) {
            return true;
        }

        if (removeVesicleOverlap && removeVesicleOverlap()) {
            return true;
        }

        initVesicleStartLocations();

        initConnectors();

        vesicleStats();

        countVesiclesWithinPSF();

        if (PSFi != null) {
            PSFi.sum(geometry);
        }

        if (PSFd != null) {
            if (geometry.simpleCuboid) {
                PSFd.sum(geometry);
            } else {
                sumPSF(PSFd, geometry.voxelSpace);
            }
        }

        meanFluorescence(true);

        if (hydrodynamicsLocalD) {
            localDensityAll(true);
        }

        if (maxNumConnectors == 0) {
            connectVesicles = false;
        }

        Master.updatePanel2D();

        return false;

    }

    public boolean checkVariables(){

        geometry = project.geometry;

        return false;

    }

    public boolean initDiffusantVesiclesArray(){

        int count = 0;

        if (project.diffusants == null) {
            error("No Diffusant Vesicles!");
            return true;
        }

        for (int i = 0; i < project.diffusants.length; i++) {
            if (project.diffusants[i] instanceof DiffusantVesicles) {
                count++;
            }
        }

        if (count == 0) {
            error("No Diffusant Vesicles!");
            return true;
        }

        diffusant = new DiffusantVesicles[project.diffusants.length];

        for (int i = 0; i < project.diffusants.length; i++) {

            if (project.diffusants[i] instanceof DiffusantVesicles) {
                diffusant[i] = (DiffusantVesicles) project.diffusants[i];
            } else {
                diffusant[i] = null;
            }
        }

        return false;

    }

    public boolean initDetectors() {

        int count = 0;

        Detector[] detect = project.detectors;

        if (detect == null) {
            return false; // nothing to do
        }

        for (int i = 0; i < detect.length; i++) {
            if (detect[i] instanceof DetectorPSF) {
                count++;
            }
        }

        if (count == 0) {
            return false; // nothing to do
        }

        detector = new DetectorPSF[detect.length];

        for (int i = 0; i < detect.length; i++) {

            if (detect[i] instanceof DetectorPSF) {
                detector[i] = (DetectorPSF) detect[i];
            } else {
                detector[i] = null;
            }

        }

        return false;

    }

    public void sumPSF(PSF psf, Voxel[][][] voxels){

        double avg = 0, sum = 0, count = 0;

        if (voxels == null){
            return;
        }

        if (psf.array == null) {
            return;
        }

        if (psf.array.length != voxels.length) {
            return;
        }

        if (psf.array[0].length != voxels[0].length) {
            return;
        }

        if (psf.array[0][0].length != voxels[0][0].length) {
            return;
        }

        for (int i = 0; i < psf.array.length; i++) {
            for (int j = 0; j < psf.array[0].length; j++) {
                for (int k = 0; k < psf.array[0][0].length; k++) {
                    if (voxels[i][j][k].isSpace) {
                        sum += psf.array[i][j][k] * voxels[i][j][k].PSFweight;
                        count++;
                    }
                }
            }
        }

        if (count > 0) {
            avg = sum / count;
        }

        Master.log("weighted dPSF sum: " + sum);

        psf.sum = sum;
        psf.avg = avg;

        psf.setParamObject("sum", sum);
        psf.setParamObject("avg", avg);

    }

    public void setPSFSelect(int IPSF, int DPSF) {
        iPSF = IPSF;
        dPSF = DPSF;
    }

    public boolean initPSFs() {

        int count = 0;
        double e2 = 0.135335;

        if (!frapOn) {
            return false; // nothing to do
        }

        if (geometry.PSFs == null) {
            return true;
        }

        if ((iPSF < 0) || (iPSF >= geometry.PSFs.length)) {
            Master.exit("MonteCarlo : initPSFs : bad value for iPSF : " + iPSF);
        }

        if ((dPSF < 0) || (dPSF >= geometry.PSFs.length)) {
            Master.exit("MonteCarlo : initPSFs : bad value for dPSF : " + dPSF);
        }

        PSFi = geometry.PSFs[iPSF];

        PSFi.useGeometryCoordinates = true;
        PSFi.array = null;
        PSFi.checkExists();

        PSFd = geometry.PSFs[dPSF];

        PSFd.useGeometryCoordinates = true;
        PSFd.array = null;
        PSFd.checkExists();

        for (int k = 0; k < PSFi.array[0][0].length; k++) {
            for (int j = 0; j < PSFi.array[0].length; j++) {
                for (int i = 0; i < PSFi.array.length; i++) {
                    if (PSFi.array[i][j][k] > e2) {
                        count++;
                    }
                }
            }
        }

        Master.log("PSFi(voxel > e^-2) n = " + count);

        count = 0;

        for (int k = 0; k < PSFd.array[0][0].length; k++) {
            for (int j = 0; j < PSFd.array[0].length; j++) {
                for (int i = 0; i < PSFd.array.length; i++) {
                    if (PSFd.array[i][j][k] > e2) {
                        count++;
                    }
                }
            }
        }

        Master.log("PSFc(voxel > e^-2) n = " + count);

        return false;

    }

    public boolean initVoxels() {

        if (PBC) {
            geometry.initVoxelsPBC(iPSF, dPSF);
        } else {
            geometry.initVoxels(iPSF, dPSF);
        }

        voxelVolume = project.dx * project.dx * project.dx;

        return false;
    }

    public boolean initVesicles() {

        double hydroDsDff;

        if (diffusant == null) {
            return true;
        }

        for (int i = 0; i < diffusant.length; i++) {

            if (diffusant[i] != null) {

                diffusant[i].initVesicles();

                if ((diffusant[i].xyz != null) && (diffusant[i].vesicle.length == diffusant[i].xyz.length)) {
                    for (int j = 0; j < diffusant[i].vesicle.length; j++) {
                        setVesicleLocation(diffusant[i].vesicle[j], diffusant[i].xyz[j][0], diffusant[i].xyz[j][1], diffusant[i].xyz[j][2], true);
                        addToVoxelList(diffusant[i].vesicle[j]);
                    }
                }

            }

        }

        if (hydroWallZ) {

            for (int i = 0; i < diffusant.length; i++) {

                if ((diffusant[i] == null) || (diffusant[i].vesicle == null)) {
                    continue;
                }

                hydroDsDff = DiffusantVesicle.Dratio_short(diffusant[i].mobileVolumeFraction, diffusant[i].immobileVolumeFraction);
                hydroDsDff /= DiffusantVesicle.Dff_short_Banchio(diffusant[i].setVolumeFraction * (1 - diffusant[i].setImmobilePercent));

                for (int j = 0; j < diffusant[i].vesicle.length; j++) {
                    diffusant[i].vesicle[j].DsDff = hydroDsDff;
                }

            }

        }

        return false;

    }

    public boolean initVesiclesImmobile() {

        if (diffusant == null) {
            return true;
        }

        for (int i = 0; i < diffusant.length; i++) {

            if (diffusant[i] == null) {
                continue;
            }

            if (diffusant[i].initImmobileVesicles(null)){
                return true;
            }

        }

        return false;

    }

    public boolean saveVesicleDensity(String saveTag) {
        return false;
    }

    public boolean saveVesiclePositions(int msec) {

        if (preview || (diffusant == null)) {
            return true;
        }

        for (int i = 0; i < diffusant.length; i++) {

            if (diffusant[i] == null) {
                continue;
            }

            if (diffusant[i].saveXYZ){
                diffusant[i].saveXYZ(msec);
            }

        }

        return false;

    }

    public double spaceVolume() {
        return project.geometry.spaceVolume - mitoVolume;
    }

    public void vesicleStats() {

        double distance;
        double spaceVolume = spaceVolume();

        minVesicleRadius = Double.POSITIVE_INFINITY;
        maxVesicleRadius = 0;
        minAllowedDistanceBetweenVesicles = Double.POSITIVE_INFINITY;

        totalVesicleVolume = 0;

        numVesicles = 0;
        mobileVesicles = 0;
        immobileVesicles = 0;

        vesicleDensity = 0;
        vesicleVolumeFraction = 0;

        if (diffusant == null) {
            return;
        }

        for (int i = 0; i < diffusant.length; i++) {

            if (diffusant[i] == null) {
                continue;
            }

            //diffusants[k].init();
            diffusant[i].vesicleStats();

            minVesicleRadius = Math.min(minVesicleRadius, diffusant[i].minRadius);
            maxVesicleRadius = Math.max(maxVesicleRadius, diffusant[i].maxRadius);

            totalVesicleVolume += diffusant[i].totalVolume;

            numVesicles += diffusant[i].numVesicles;
            mobileVesicles += diffusant[i].mobileVesicles;
            immobileVesicles += diffusant[i].immobileVesicles;

        }

        mobilePercent = 1.0 * mobileVesicles / (immobileVesicles + mobileVesicles);
        immobilePercent = 1.0 * immobileVesicles / (immobileVesicles + mobileVesicles);

        for (int i = 0; i < diffusant.length; i++) {

            if (diffusant[i] == null) {
                continue;
            }

            for (int j = 0; j < diffusant.length; j++) {

                if (diffusant[j] == null) {
                    continue;
                }

                distance = diffusant[i].minRadius + diffusant[j].minRadius;
                minAllowedDistanceBetweenVesicles = Math.min(distance, minAllowedDistanceBetweenVesicles);

            }
        }

        vesicleVolume = 4 * Math.PI * minVesicleRadius * minVesicleRadius * minVesicleRadius / 3.0;

        vesicleDensity = (numVesicles * 1.0) / spaceVolume;
        vesicleVolumeFraction = totalVesicleVolume / spaceVolume;

        mobileVolumeFraction = mobilePercent * vesicleVolumeFraction;
        immobileVolumeFraction = immobilePercent * vesicleVolumeFraction;

        setParamObject("minVesicleRadius", minVesicleRadius);
        setParamObject("maxVesicleRadius", maxVesicleRadius);
        setParamObject("minAllowedDistanceBetweenVesicles", minAllowedDistanceBetweenVesicles);

        setParamObject("vesicleVolume", totalVesicleVolume);

        setParamObject("numVesicles", numVesicles);
        setParamObject("mobileVesicles", mobileVesicles);
        setParamObject("immobileVesicles", immobileVesicles);
        setParamObject("mobilePercent", mobilePercent);
        setParamObject("immobilePercent", immobilePercent);
        setParamObject("mobileVolumeFraction", mobileVolumeFraction);
        setParamObject("immobileVolumeFraction", immobileVolumeFraction);

        setParamObject("vesicleDensity", vesicleDensity);
        setParamObject("vesicleVolumeFraction", vesicleVolumeFraction);

        minDistance();

    }

    public boolean checkDX() {

        double maxRadius = 0;

        if (diffusant == null) {
            return true;
        }

        for (int i = 0; i < diffusant.length; i++) {

            if ( (diffusant[i] == null) || (diffusant[i].vesicle == null)) {
                continue;
            }

            for (int j = 0; j < diffusant[i].vesicle.length; j++) {

                if (diffusant[i].vesicle[j] == null) {
                    continue;
                }

                maxRadius = Math.max(maxRadius, diffusant[i].vesicle[j].radius);

            }

        }

        if (maxRadius <= 0) {
            Master.exit("MonteCarlo: checkDX: bad maxRadius: " + maxRadius);
        }

        if ((int) (100 * project.dx) < (int) (100 * 2 * maxRadius)) {
            //Master.exit("MonteCarlo: checkDX: voxel size is smaller than maximum vesicle diameter: " + (2 * maxRadius));
        }

        return false;

    }

    public boolean initDT() {

        double dt, maxD = 0;

        if (diffusant == null) {
            return true;
        }

        for (int i = 0; i < diffusant.length; i++) {

            if ( (diffusant[i] == null) || (diffusant[i].vesicle == null)) {
                continue;
            }

            for (int j = 0; j < diffusant[i].vesicle.length; j++) {
                maxD = Math.max(maxD, diffusant[i].vesicle[j].D);
            }

        }

        if (maxD <= 0) {
            Master.exit("initDT: negative maxD: " + maxD);
            return true;
        }

        dt = minVesicleStep * minVesicleStep / ( 6.0 * maxD );
        project.dt = dt;
        project.stability = dt * 3 * maxD / (project.dx * project.dx);

        // can only compute stability since maxD, dt and dx are set elsewhere

        return false;

    }

    public boolean initMitochondria() {

        int count = 0, iMitoTrials = 20, itrialMax = 100, radiusTrials = 2000;
        double d, x = 0, y = 0, z = 0;
        double radius, rstep, vf, volume;
        boolean overlap, finished = false;

        mitoVolume = 0;

        if (!mitoCoordinatesOn || (setMitoVolumeFraction <= 0)) {
            return false;
        }

        if (mito != null) {

            for (int i = 0 ; i < mito.length; i++) {
                mitoVolume += mito[i].geometryVolume();
            }

            mitoVolumeFraction = mitoVolume / geometry.volume;

            Master.log("final mito volume fraction = " + mitoVolumeFraction);

            return false;

        }

        Coordinates test = new Coordinates(project);
        Coordinates[] temp = new Coordinates[iMitoTrials];

        radius = mitoRadius;
        rstep = radius / radiusTrials;

        if (mito_ijkselect == -2) {

            d = mt.nextDouble();

            if (d < 0.3333333333) {
                mito_ijkselect = 0; // xy plane
            } else if (d < 0.6666666666) {
                mito_ijkselect = 1; // yz plane
            } else {
                mito_ijkselect = 2; // zx plane
            }

        }

        Master.log("initializing mitochondria...");

        for (int imito = 0; imito < iMitoTrials; imito++) {

            overlap = true;

            for (int itrial = 0; itrial < itrialMax; itrial++) {

                x = geometry.x1 + mt.nextDouble() * (geometry.x2 - geometry.x1);
                y = geometry.y1 + mt.nextDouble() * (geometry.y2 - geometry.y1);
                z = geometry.z1 + mt.nextDouble() * (geometry.z2 - geometry.z1);

                if (mito_ijkselect == -1) {

                    d = mt.nextDouble(); // random orientation

                    if (d < 0.3333333333) {
                        mito_ijkselect = 0;
                    } else if (d < 0.6666666666) {
                        mito_ijkselect = 1;
                    } else {
                        mito_ijkselect = 2;
                    }

                }

                setMitoPosition(test, x, y, z, radius, mitoAxialRatio, mito_ijkselect);

                overlap = false;

                for (int jmito = 0; jmito < imito; jmito++) {
                    if (temp[jmito].intersectsCuboid(test)) {
                        overlap = true;
                        break;
                    }
                }

                if (!overlap) {
                    break;
                }

            }

            if (overlap) {
                Master.exit("MonteCarlo: initMitochondria : failed to place mitochondria #" + imito);
            }

            volume = test.geometryVolume();

            vf = (mitoVolume + volume) / geometry.volume;

            if (vf == setMitoVolumeFraction) {

                finished = true;

            } else if (vf > setMitoVolumeFraction) {

                finished = false;

                for (int i = 1; i < radiusTrials; i++) {

                    radius -= i * rstep;

                    if (radius <= 0) {
                        break;
                    }

                    setMitoPosition(test, x, y, z, radius, mitoAxialRatio, mito_ijkselect);

                    volume = test.geometryVolume();

                    vf = (mitoVolume + volume) / geometry.volume;

                    if (vf <= setMitoVolumeFraction) {
                        finished = true;
                        break;
                    }

                }

                if (!finished) {
                    Master.exit("MonteCarlo: initMitochondria : failed to place mitochondria #" + imito);
                }

            }

            temp[imito] = new Coordinates(project, test);
            temp[imito].updateDimensions();
            temp[imito].name = "mito #" + imito;
            temp[imito].color.setColor("[r=155,g=155,b=155]");

            mitoVolume += volume;
            count++;

            Master.log("mito #" + imito + " radius = " + radius);
            //Master.log("mito shape " + temp[imito].shape);
            //Master.log("mito shapeNum " + temp[imito].shapeNum);
            //Master.log("mitoCoordinatesOn % volume = " + vf * 100.0);

            if (finished) {
                break;
            }

        }

        if (count == 0) {
            return true;
        }

        mito = new Coordinates[count];

        System.arraycopy(temp, 0, mito, 0, count);

        mitoVolumeFraction = mitoVolume / geometry.volume;

        Master.log("final mito volume fraction = " + mitoVolumeFraction);

        return false;

    }

    public void setMitoPosition(Coordinates c, double x, double y, double z, double radius, double mitoAxialRatio, int ijkSelect) {

        String xyz = "";

        double halfLength = 0.5 * 2.0 * radius * mitoAxialRatio;

        switch (ijkSelect) {
            case 0:
                c.x1 = x - halfLength; // geometry.x1;
                c.y1 = y - radius;
                c.z1 = z - radius;
                c.x2 = x + halfLength; // geometry.x2;
                c.y2 = y + radius;
                c.z2 = z + radius;
                xyz = "x";
                break;
            case 1:
                c.x1 = x - radius;
                c.y1 = y - halfLength; // geometry.y1;
                c.z1 = z - radius;
                c.x2 = x + radius;
                c.y2 = y + halfLength; // geometry.y2;
                c.z2 = z + radius;
                xyz = "y";
                break;
            case 2: // xy
                c.x1 = x - radius;
                c.y1 = y - radius;
                c.z1 = z - halfLength; // geometry.z1;
                c.x2 = x + radius;
                c.y2 = y + radius;
                c.z2 = z + halfLength; // geometry.z2;
                xyz = "z";
                break;
        }

        if (mitoShape.equalsIgnoreCase("cuboid")) {
            c.setShape("cuboid");
        } else if (mitoShape.equalsIgnoreCase("ellipsoid")) {
            c.setShape("ellipsoid");
        } else if (mitoShape.equalsIgnoreCase("cylinder")) {
            c.setShape("cylinder" + xyz);
        }

        c.updateDimensions();

    }

    public boolean testMitochondriaOverlap(DiffusantVesicle dv) {

        if (mito == null) {
            return false;
        }

        for (int i = 0; i < mito.length; i++) {
            if (mito[i] == null) {
                break;
            }
            if (mito[i].isInside(dv.x, dv.y, dv.z, dv.radius)) {
                return true;
            }
        }

        return false;

    }

    public static double hydroWall_ll(double radius, double z) {

        double x = radius / z;

        // z = shortest distance from sphere center to wall

        x = Math.min(x, 1);

        double b1 = 1 - (9 * x / 16) + (x * x * x / 8) - (x * x * x * x * 45 / 256) - (x * x * x * x * x / 16);

        return b1;

    }

    public static double hydroWall_T(double radius, double h) {

        double b2;

        //double h = z - radius; // distance between wall and edge of sphere

        h = Math.max(h, 0);
        b2 = 6 * h * h + 2 * radius * h;
        b2 /= 6 * h * h + 9 * radius * h + 2 * radius * radius;

        return b2;

    }

    public boolean setVesicleLocation(DiffusantVesicle dv, double x, double y, double z, boolean initStartLocation) {

        int xVoxel, yVoxel, zVoxel;

        if (PBC) {

            dv.x = x;
            dv.y = y;
            dv.z = z;

            if (initStartLocation) {
                dv.x0 = x;
                dv.y0 = y;
                dv.z0 = z;
            }

            xVoxel = (int) geometry.computeVoxelX(x);
            yVoxel = (int) geometry.computeVoxelY(y);
            zVoxel = (int) geometry.computeVoxelZ(z);

            if ((xVoxel < 0) || (xVoxel >= geometry.voxelSpacePBC.length)) {
                //return false;
                return true;
            }

            if ((yVoxel < 0) || (yVoxel >= geometry.voxelSpacePBC[0].length)) {
                //return false;
                return true;
            }

            if ((zVoxel < 0) || (zVoxel >= geometry.voxelSpacePBC[0][0].length)) {
                //return false;
                return true;
            }

            dv.voxel = geometry.voxelSpacePBC[xVoxel][yVoxel][zVoxel];

            return true;

        } else {

            xVoxel = (int) geometry.computeVoxelX(x);
            yVoxel = (int) geometry.computeVoxelY(y);
            zVoxel = (int) geometry.computeVoxelZ(z);

            if ((xVoxel < 0) || (xVoxel >= geometry.voxelSpace.length)) {
                return false;
            }

            if ((yVoxel < 0) || (yVoxel >= geometry.voxelSpace[0].length)) {
                return false;
            }

            if ((zVoxel < 0) || (zVoxel >= geometry.voxelSpace[0][0].length)) {
                return false;
            }

            dv.x = x;
            dv.y = y;
            dv.z = z;

            if (initStartLocation) {
                dv.x0 = x;
                dv.y0 = y;
                dv.z0 = z;
            }

            dv.voxel = geometry.voxelSpace[xVoxel][yVoxel][zVoxel];

            return true;

        }

    }

    public boolean outOfBounds(DiffusantVesicle dv) {

        if (geometry.simpleCuboid) {
            return beyondLimits(dv); // simple computation
        } else {
            if (beyondLimits(dv)) {
                return true;
            }
            return insideVoxelNonSpace(dv); // check nonspace voxels
        }

    }

    public boolean beyondLimits(DiffusantVesicle dv) {

        double extra;

        if (PBC || freeDiffusion) {
            extra = 0;
        } else {
            extra = dv.radius;
        }

        if ((dv.x < geometry.x1 + extra) || (dv.x > geometry.x2 - extra)) {
            return true;
        }

        if ((dv.y < geometry.y1 + extra) || (dv.y > geometry.y2 - extra)) {
            return true;
        }

        if ((dv.z < geometry.z1 + extra) || (dv.z > geometry.z2 - extra)) {
            return true;
        }

        return false;

    }

    public boolean insideVoxelNonSpace(DiffusantVesicle dv) {

        double dlimit;
        Voxel ivoxel;

        if (dv.voxel == null) {
            return false;
        }

        if (!dv.voxel.isSpace) {
            return true;
        }

        if (freeDiffusion) {
            return false;
        }

        dlimit = project.dx * 0.5 + dv.radius; // um

        for (int i = 0; i < dv.voxel.numNonSpaceNeighbors; i++) {

            ivoxel = dv.voxel.nonSpaceNeighbors[i];

            if ((dv.x > ivoxel.x - dlimit) && (dv.x < ivoxel.x + dlimit)) {
                if ((dv.y > ivoxel.y - dlimit) && (dv.y < ivoxel.y + dlimit)) {
                    if ((dv.z > ivoxel.z - dlimit) && (dv.z < ivoxel.z + dlimit)) {
                        return true;
                    }
                }
            }

        }

        return false;

    }

    public void initVesicleStartLocations() {

        for (int i = 0; i < diffusant.length; i++) {

            if ((diffusant[i] == null) || (diffusant[i].vesicle == null)) {
                continue;
            }

            diffusant[i].initStartLocation();

        }

    }

    public boolean initVesicleRandom(DiffusantVesicle dv, Coordinates c, boolean allowOverlaps, int overlapTrialLimit, int absLimit) {

        int trial = 0;
        double x, y, z;
        boolean overlap;

        if ((dv == null) || (c == null)) {
            return true;
        }

        while (true) {

            trial++;

            if (trial >= absLimit) {
                return true;
            }

            x = (mt.nextDouble() * (c.x2 - c.x1)) + c.x1;
            y = (mt.nextDouble() * (c.y2 - c.y1)) + c.y1;
            z = (mt.nextDouble() * (c.z2 - c.z1)) + c.z1;

            if (!setVesicleLocation(dv, x, y, z, true)){
                continue;
            }

            if (outOfBounds(dv)) {
                continue;
            }

            if (testMitochondriaOverlap(dv)) {
                continue;
            }

            overlap = false;

            if (!freeDiffusion) {
                overlap = (testVesicleOverlap(dv, null) != null);
            }

            if (!overlap || (allowOverlaps && (trial > overlapTrialLimit))) {

                if (!addToVoxelList(dv)) {
                    return true;
                }

                return false;

            }

        }

    }

    public boolean initVesiclesRandom() {

        int nVesicles;

        boolean allowOverlaps = true;

        Coordinates c = null;

        if (diffusant == null) {
            return true;
        }

        for (int i = 0; i < diffusant.length; i++) {

            if ((diffusant[i] == null) || (diffusant[i].vesicle == null)) {
                continue;
            }

            if (diffusant[i].coordinates == null) {
                c = new Coordinates(project, geometry);
            } else {
                c = new Coordinates(project, diffusant[i].coordinates);
            }

            nVesicles = diffusant[i].vesicle.length;

            Master.log("randomly placing ready vesicles (" + nVesicles + ")");

            for (int j = 0; j < nVesicles; j++) {

                if (!diffusant[i].vesicle[j].insideGeometry) {
                    continue;
                }

                if (initVesicleRandom(diffusant[i].vesicle[j], c, allowOverlaps, overlapTrialLimit, absLimit)) {
                    return true;
                }

                if ((j > 0) && (Math.IEEEremainder(j, 10000) == 0)) {
                    Master.log("placed vesicle " + j + " / " + nVesicles);
                }

            }

        }

        return false;

    }
    
    public boolean removeVesicleOverlap() {

        int radiusNM, lastRadiusNM = 0, count = 0;
        double halfDistance = 0;

        if (freeDiffusion) {
            return false; // OK
        }

        if (!testVesicleOverlap()){
            return false; // OK
        }

        if (diffusant == null) {
            return true;
        }

        vesicleStats(); // will update vesicles step size

        removingVesicleOverlap = true;

        Master.log("elimination of vesicle overlap...");

        while (halfDistance < minVesicleRadius) {

            halfDistance = 0.5 * minDistance();

            radiusNM = (int) (1000 * halfDistance);

            if (radiusNM > lastRadiusNM) {
                Master.log("current radius: " + radiusNM + " nm");
                lastRadiusNM = radiusNM;
                count = 0;
            }

            // set current diameter to minimal distance between vesicles

            for (int i = 0; i < diffusant.length; i++) {

                //if ((diffusant[k] == null) || (diffusant[k].vesicle == null)) {
                //    continue;
                //}

                for (int j = 0; j < diffusant[i].vesicle.length; j++) {
                    if (diffusant[i].vesicle[j].insideGeometry) {
                        //diffusant[k].vesicle[j].radius = Math.min(halfDistance, diffusant[k].meanRadius);
                        diffusant[i].vesicle[j].radius = halfDistance;
                    }
                }

            }

            moveVesiclesCichocki(true);

            count++;

            if (count > absLimit) {
                break;
            }

        }

        for (int i = 0; i < diffusant.length; i++) {

            if ((diffusant[i] == null) || (diffusant[i].vesicle == null)) {
                continue;
            }

            //radius = diffusant[k].meanRadius;

            for (int j = 0; j < diffusant[i].vesicle.length; j++) {
                if (diffusant[i].vesicle[j].insideGeometry) {
                    diffusant[i].vesicle[j].radius = diffusant[i].meanRadius;
                    diffusant[i].vesicle[j].x0 = diffusant[i].vesicle[j].x;
                    diffusant[i].vesicle[j].y0 = diffusant[i].vesicle[j].y;
                    diffusant[i].vesicle[j].z0 = diffusant[i].vesicle[j].z;
                }
            }

        }

        removingVesicleOverlap = false;

        if (testVesicleOverlap()){
            error("removeVesicleOverlap: failed to remove vesicle overlap.");
            return true;
        }

        if (initVesiclesImmobile()) {
            return true;
        }

        return false;

    }

    public double minDistanceSlow() {

        int countOverlaps = 0;
        double sqrDistance, minSqrDist, minSqrDist2;
        double DX, DY, DZ;

        minSqrDist = Double.POSITIVE_INFINITY;
        minDistanceBetweenVesicles = Double.NaN;

        boolean found = false;

        if (diffusant == null) {
            return Double.NaN;
        }

        for (int i = 0; i < diffusant.length; i++) {

            if ((diffusant[i] == null) || (diffusant[i].vesicle == null)) {
                continue;
            }

            for (int j = 0; j < diffusant[i].vesicle.length - 1; j++) {

                if (!diffusant[i].vesicle[j].insideGeometry) {
                    continue;
                }

                minSqrDist2 = Double.POSITIVE_INFINITY;

                for (int k = j+1; k < diffusant[i].vesicle.length; k++) {

                    //if (j == k) {
                    //    continue;
                    //}

                    if (!diffusant[i].vesicle[k].insideGeometry) {
                        continue;
                    }

                    DX = diffusant[i].vesicle[j].x - diffusant[i].vesicle[k].x;
                    DY = diffusant[i].vesicle[j].y - diffusant[i].vesicle[k].y;
                    DZ = diffusant[i].vesicle[j].z - diffusant[i].vesicle[k].z;

                    sqrDistance = DX * DX + DY * DY + DZ * DZ;

                    if (sqrDistance == 0) {
                        Master.log("warning: zero sqrDistance " + j + ", " + k);
                    }

                    if (Math.sqrt(sqrDistance) < 0.04) {
                        if ((Math.abs(DZ) <= 0.03) && (Math.abs(DZ) > 0)) {
                            if ((Math.abs(DX) < 0.01) && (Math.abs(DY) < 0.01)) {
                                //Master.log("" + DX + ", " + DY + ", " + DZ);
                                //Master.log("" + (Math.sqrt(DX * DX + DY * DY)));
                                diffusant[i].vesicle[k].insideGeometry = false;
                                diffusant[i].vesicle[k].x = Double.NaN;
                                diffusant[i].vesicle[k].y = Double.NaN;
                                diffusant[i].vesicle[k].z = Double.NaN;
                                countOverlaps++;
                            }
                        }
                    }

                    if (sqrDistance < minSqrDist) {
                        minSqrDist = sqrDistance;
                        found = true;
                    }

                    if (sqrDistance < minSqrDist2) {
                        minSqrDist2 = sqrDistance;
                    }

                }

                //Master.log("" + Math.sqrt(minSqrDist2));

            }

        }

        if (countOverlaps > 0) {
            Master.log("removed overlaps = " + countOverlaps);
        }

        if (found){
            minDistanceBetweenVesicles = Math.sqrt(minSqrDist);
        } else {
            minDistanceBetweenVesicles = Double.NaN;
        }

        return minDistanceBetweenVesicles;

    }

    public double minDistance() {

        double min, minSqrDist = Double.MAX_VALUE;

        boolean found = false;

        minDistanceBetweenVesicles = Double.NaN;

        if (diffusant == null) {
            return Double.NaN;
        }

        for (int i = 0; i < diffusant.length; i++) {

           if ((diffusant[i] == null) || (diffusant[i].vesicle == null)) {
                continue;
            }

            for (int j = 0; j < diffusant[i].vesicle.length; ++j) {

                if (diffusant[i].vesicle[j] == null) {
                    continue;
                }

                if (diffusant[i].vesicle[j].insideGeometry) {

                    min = minSqrDistance(diffusant[i].vesicle[j], diffusant[i].vesicle[j]);

                    if (min >= 0) {

                        if (min < minSqrDist) {
                            minSqrDist = min;
                            found = true;
                        }

                    } else {
                        //Master.log("failed to find min distance for vesicles " + j);
                    }

                }

            }

        }

        if (found) {
            minDistanceBetweenVesicles = Math.sqrt(minSqrDist);
        } else {
            minDistanceBetweenVesicles = minDistanceSlow();
        }

        return minDistanceBetweenVesicles;

    }

    public double minSqrDistance(DiffusantVesicle dv, DiffusantVesicle ignoreVesicle) {

        double dx, dy, dz, sqrDistance = 0, min = Double.MAX_VALUE;
        int ii, jj, kk;
        boolean found = false;

        DiffusantVesicle ivesicle;

        Voxel voxel = dv.voxel;
        VoxelPBC voxelPBC = null;
        Voxel ivoxel;

        if (voxel == null) {
            return -1;
        }

        if (PBC && (voxel instanceof VoxelPBC)) {
            voxelPBC = (VoxelPBC) voxel;
        }

        for (int i = 0; i < voxel.numNeighbors; i++) {

            ivoxel = voxel.neighbors[i];
            ivesicle = ivoxel.firstReady;

            while (ivesicle != null) {

                if ((ivesicle != dv) && (ivesicle != ignoreVesicle) && (!Double.isNaN(ivesicle.x)) && ivesicle.insideGeometry ) {

                    dx = dv.x - ivesicle.x;
                    dy = dv.y - ivesicle.y;
                    dz = dv.z - ivesicle.z;

                    sqrDistance = dx * dx + dy * dy + dz * dz;

                    if (sqrDistance < min) {
                        min = sqrDistance;
                        found = true;
                    }

                    if (sqrDistance <= 0) {
                        Master.log("warning: zero sqrDistance " + dv + " and " + ivesicle);
                    }

                }

                ivesicle = ivesicle.nextReady;

            }

        }

        if (PBC && (voxelPBC != null)) {

            for (int i = 0; i < voxelPBC.numPBCneighbors; i++) {

                ivoxel = voxelPBC.PBCneighbors[i];
                ivesicle = ivoxel.firstReady;

                while (ivesicle != null) {

                    if ((ivesicle != dv) && (ivesicle != ignoreVesicle) && (!Double.isNaN(ivesicle.x)) && ivesicle.insideGeometry ) {

                        ii = voxelPBC.PBCi[i];
                        jj = voxelPBC.PBCj[i];
                        kk = voxelPBC.PBCk[i];

                        dx = dv.x - (ii * 2 * geometry.x2 + ivesicle.x);
                        dy = dv.y - (jj * 2 * geometry.y2 + ivesicle.y);
                        dz = dv.z - (kk * 2 * geometry.z2 + ivesicle.z);

                        sqrDistance = dx * dx + dy * dy + dz * dz;

                        if (sqrDistance < min) {
                            min = sqrDistance;
                            found = true;
                        }

                    }

                    ivesicle = ivesicle.nextReady;

                }

            }

        }

        if (found) {
            return min;
        } else {
            return -1;
        }

    }

    public boolean testVesicleOverlap(){

        if (diffusant == null) {
            return false;
        }

        for (int i = 0; i < diffusant.length; i++) {

            if ((diffusant[i] == null) || (diffusant[i].vesicle == null)) {
                continue;
            }

            for (int j = 0; j < diffusant[i].vesicle.length; ++j) {

                if (diffusant[i].vesicle[j] == null) {
                    continue;
                }

                if (!diffusant[i].vesicle[j].insideGeometry) {
                    continue;
                }

                if (testVesicleOverlap(diffusant[i].vesicle[j], diffusant[i].vesicle[j]) != null) {
                    return true;
                }

            }
        }

        return false;

    }
    
    public DiffusantVesicle testVesicleOverlap(DiffusantVesicle dv, DiffusantVesicle ignoreVesicle) {

        double dx, dy, dz, sqrDistance = 0, minDBV;
        int ii, jj, kk;

        DiffusantVesicle ivesicle;

        Voxel voxel = dv.voxel;
        VoxelPBC voxelPBC = null;
        Voxel ivoxel;

        if (voxel == null) {
            return null;
        }

        if (PBC && (voxel instanceof VoxelPBC)) {
            voxelPBC = (VoxelPBC) voxel;
        }

        for (int i = 0; i < voxel.numNeighbors; i++) {

            ivoxel = voxel.neighbors[i];
            ivesicle = ivoxel.firstReady;

            while (ivesicle != null) {

                if ((ivesicle != dv) && (ivesicle != ignoreVesicle) && (!Double.isNaN(ivesicle.x)) && ivesicle.insideGeometry && (ivesicle.radius > 0 )) {

                    dx = dv.x - ivesicle.x;
                    dy = dv.y - ivesicle.y;
                    dz = dv.z - ivesicle.z;

                    sqrDistance = dx * dx + dy * dy + dz * dz;

                    minDBV = dv.radius + ivesicle.radius;

                    if (sqrDistance < minDBV * minDBV) {
                        return ivesicle;
                    }

                }

                ivesicle = ivesicle.nextReady;

            }

        }

        if (PBC && (voxelPBC != null)) {

            for (int i = 0; i < voxelPBC.numPBCneighbors; i++) {

                ivoxel = voxelPBC.PBCneighbors[i];
                ivesicle = ivoxel.firstReady;

                while (ivesicle != null) {

                    if ((ivesicle != dv) && (ivesicle != ignoreVesicle) && (!Double.isNaN(ivesicle.x)) && ivesicle.insideGeometry ) {

                        ii = voxelPBC.PBCi[i];
                        jj = voxelPBC.PBCj[i];
                        kk = voxelPBC.PBCk[i];

                        dx = dv.x - (ii * 2 * geometry.x2 + ivesicle.x);
                        dy = dv.y - (jj * 2 * geometry.y2 + ivesicle.y);
                        dz = dv.z - (kk * 2 * geometry.z2 + ivesicle.z);

                        sqrDistance = dx * dx + dy * dy + dz * dz;

                        minDBV = dv.radius + ivesicle.radius;

                        if (sqrDistance < minDBV * minDBV) {
                            return ivesicle;
                        }

                    }

                    ivesicle = ivesicle.nextReady;

                }

            }

        }

        return null;

    }

    public double avgMinSqrDistance() {

        double minSD;
        double avgMinSD = 0, count = 0;

        for (int i = 0; i < diffusant.length; i++) {

            if ((diffusant[i] == null) || (diffusant[i].vesicle == null)) {
                continue;
            }

            for (int j = 0; j < diffusant[i].vesicle.length; ++j) {

                if (diffusant[i].vesicle[j] == null) {
                    continue;
                }

                minSD = minSqrDistance(diffusant[i].vesicle[j], diffusant[i].vesicle[j]);

                if (minSD > 0) {
                    avgMinSD += minSD;
                    count += 1;
                }

            }

        }

        return (avgMinSD/count);

    }

    public void runMSDspatial() {

        double svalue;

        if (!saveMSDspatial || (MSDspatialTbgn < 0)) {
            return;
        }

        if ((MSDspatial == null) && (time >= MSDspatialTbgn)) {
            initMSDspatial(MSDspatialTbgn);
            initD20();
            MSDspatial_t1 = MSDspatialTbgn;
            MSDspatial_t2 = MSDspatial_t1 + MSDspatial_win;
        }

        if ((MSDspatial_t1 < 0) || (MSDspatial_t2 < 0)) {
            return;
        }
        
        if ((MSDspatial != null) && (time > MSDspatial_t2)) {
            
            for (int i = 0; i < MSDspatial.length; i++) {
                MSDspatial[i].finish("Vesicles", project.geometry, -1);
            }
            
            MSDspatial = null;
            MSDspatial_t1 = -1;
            MSDspatial_t2 = -1;
            MSDspatialTbgn = -1;
            
            return;
            
        }

        if ((time >= MSDspatial_t1) && (time <= MSDspatial_t2)) {
            for (int i = 0; i < MSDspatial.length; i++) {
                if (MSDspatial[i].skipCounter == 0) {
                    svalue = MSDspatial(i);
                } else {
                    svalue = -1; // does not matter
                }
                MSDspatial[i].saveData(svalue);
            }
        }

    }

    public void initMSDspatial(double time) {

        int bin, it;
        int dataPoints = 1;
        String fname;

        MSDspatial = new Save[MSDspatial_numBins];

        for (int i = 0; i < MSDspatial.length; i++) {
            MSDspatial[i] = new Save(project);
            MSDspatial[i].saveWhileComputing = false;
            MSDspatial[i].save2BinaryFile = true;
            MSDspatial[i].save2TextFile = false;
            MSDspatial[i].init();
            bin = (int) ((i + 1) * MSDspatial_binWidth * 1000);
            it = (int) time;
            fname = "MSD_t" + Integer.toString(it) + "_d" + Integer.toString(bin);
            Master.log("init save " + fname);
            MSDspatial[i].fileName(fname, "");
            MSDspatial[i].xdim = "ms";
            MSDspatial[i].ydim = "Vesicles" + " (" + project.spaceUnits + "^2)";
            MSDspatial[i].updateVectors();
            MSDspatial[i].init("Vesicles", project.geometry, -1, dataPoints);
        }

    }

    public void initD20() {

        DiffusantVesicle dv;
        double dx, dy, dz;

        for (int i = 0; i < diffusant.length; i++) {

            if ((diffusant[i] == null) || (diffusant[i].vesicle == null)) {
                continue;
            }

            for (int j = 0; j < diffusant[i].vesicle.length; ++j) {

                dv = diffusant[i].vesicle[j];

                if ((dv == null) || !dv.mobile || !dv.insideGeometry) {
                    continue;
                }

                dv.x0 = dv.x;
                dv.y0 = dv.y;
                dv.z0 = dv.z;

                dx = dv.x - 0;
                dy = dv.y - 0;
                dz = dv.z - 0;

                dv.d20 = Math.sqrt(dx * dx + dy * dy + dz * dz);

            }

        }

    }

    public double MSDspatial(int ibin) {

        double sumSD = 0.0, count = 0.0;
        double x1 = ibin * MSDspatial_binWidth;
        double x2 = x1 + MSDspatial_binWidth;

        DiffusantVesicle dv;

        for (int i = 0; i < diffusant.length; i++) {

            if ((diffusant[i] == null) || (diffusant[i].vesicle == null)) {
                continue;
            }

            for (int j = 0; j < diffusant[i].vesicle.length; ++j) {

                dv = diffusant[i].vesicle[j];

                if ((dv == null) || !dv.mobile || !dv.insideGeometry) {
                    continue;
                }

                if ((dv.d20 >= x1) && (dv.d20 < x2)) {
                    sumSD += dv.squareDisplacement();
                    count += 1.0;
                    //Master.log("" + dv + " " + dv.d2AZ0);
                }

            }

        }

        return (sumSD / count);

    }

    public boolean initConnectors() {

        if (!connectVesicles) {
            return false;
        }

        if (connectorsBinomial) {
            return initConnectorBinomialDistribution();
        }

        Master.log("init " + maxNumConnectors + " connectors / vesicle");

        for (int i = 0; i < diffusant.length; i++) {

            if ((diffusant[i] == null) || (diffusant[i].vesicle == null)) {
                continue;
            }

            for (int j = 0; j < diffusant[i].vesicle.length; j++) {

                if (!(diffusant[i].vesicle[j] instanceof DiffusantVesicle)) {
                    continue;
                }

                diffusant[i].vesicle[j].connectTo = new DiffusantVesicle[maxNumConnectors];
                diffusant[i].vesicle[j].connectorOffTime = new double[maxNumConnectors];

            }

        }

        return false;

    }

    private boolean initConnectorBinomialDistribution() {

        double avg = 0, count = 0;
        int nConnectors, numConnected = 0;
        //double meanNumConnectors = 0.9;//1.5; // 1.0
        double probability = meanNumConnectors / (1.0 * maxNumConnectors); // Fernandez-Busnadiego 2013

        DiffusantVesicle dv;

        int[] numConnectorsHisto = new int[maxNumConnectors + 1];

        Master.log("max number of connectors = " + maxNumConnectors);
        Master.log("connector probability = " + probability);
        
        for (int i = 0; i < diffusant.length; i++) {

            if ((diffusant[i] == null) || (diffusant[i].vesicle == null)) {
                continue;
            }

            for (int j = 0; j < diffusant[i].vesicle.length; j++) {

                if (!(diffusant[i].vesicle[j] instanceof DiffusantVesicle)) {
                    continue;
                }

                dv = diffusant[i].vesicle[j];

                nConnectors = 0;

                for (int k = 0; k < maxNumConnectors; k++) {
                    if (mt.nextDouble() < probability) {
                        nConnectors++;
                    }
                }

                if (nConnectors > 0) {
                    dv.connectTo = new DiffusantVesicle[nConnectors]; // need twice as many for tracking connectors
                    dv.connectorOffTime = new double[nConnectors]; // need twice as many for tracking connectors
                }

                numConnectorsHisto[nConnectors]++;

            }

        }

        for (int i = 0; i < numConnectorsHisto.length; i++) {
            Master.log("" + i + " connectors: " + numConnectorsHisto[i]);
            avg += i * numConnectorsHisto[i];
            count += numConnectorsHisto[i];
            if (i > 0) {
                numConnected += numConnectorsHisto[i];
            }
        }

        avg /= count;

        Master.log("average = " + avg + " connectors / vesicle");
        Master.log("fraction connected = " + (numConnected/count));

        return false;

    }

    public void connectVesicles() {

        int ii, jj, kk;
        double dx, dy, dz, sqrDistance, minDBV;

        DiffusantVesicle dv, kvesicle;

        VoxelPBC voxelPBC = null;

        for (int i = 0; i < diffusant.length; i++) {

            if ((diffusant[i] == null) || (diffusant[i].vesicle == null)) {
                continue;
            }

            for (int j = 0; j < diffusant[i].vesicle.length; j++) {

                dv = diffusant[i].vesicle[j];

                for (int k = 0; k < dv.voxel.numNeighbors; k++) {

                    kvesicle = (DiffusantVesicle) dv.voxel.neighbors[k].firstReady;

                    while (kvesicle != null) {

                        if (dv.noMoreConnectors()) {
                            break;
                        }

                        if ((kvesicle != dv) && !dv.isConnectedTo(kvesicle) && !kvesicle.noMoreConnectors()) {

                            if (dv.overlap(kvesicle, connectorLength)) {

                                if (mt.nextDouble() < connectRate * project.dt) {

                                    if (kvesicle.connectToNew(dv, mt, unconnectRate, time)) {

                                        avgConnectorLifeTime += dv.connectorLifeTime;
                                        connectorLifeTimeCounter += 1;

                                        //Master.log("connected vesicles " + dv.connectorLifeTime);

                                        if (!dv.connectTo(kvesicle)) {
                                            Master.exit("connection failure");
                                        }

                                    } else {
                                        Master.exit("connection failure");
                                    }

                                }

                            }

                        }

                        kvesicle = kvesicle.nextReady;

                    }

                }

                if (dv.noMoreConnectors()) {
                    continue;
                }

                if (!PBC) {
                    continue;
                }

                if (dv.voxel instanceof VoxelPBC) {
                    voxelPBC = (VoxelPBC) dv.voxel;
                } else {
                    continue;
                }

                if (voxelPBC == null) {
                    continue;
                }

                if (true) {
                    continue;
                }

                for (int k = 0; k < voxelPBC.numPBCneighbors; k++) {

                    kvesicle = (DiffusantVesicle) voxelPBC.PBCneighbors[k].firstReady;

                    while (kvesicle != null) {

                        if (dv.noMoreConnectors()) {
                            break;
                        }

                        if ((kvesicle != dv) && !dv.isConnectedTo(kvesicle) && !kvesicle.noMoreConnectors()) {

                            ii = voxelPBC.PBCi[i];
                            jj = voxelPBC.PBCj[i];
                            kk = voxelPBC.PBCk[i];

                            dx = dv.x - (ii * 2 * geometry.x2 + kvesicle.x);
                            dy = dv.y - (jj * 2 * geometry.y2 + kvesicle.y);
                            dz = dv.z - (kk * 2 * geometry.z2 + kvesicle.z);

                            sqrDistance = dx * dx + dy * dy + dz * dz;

                            minDBV = dv.radius + kvesicle.radius + connectorLength;

                            if (sqrDistance < minDBV * minDBV) {

                                if (mt.nextDouble() < connectRate * project.dt) {

                                    if (kvesicle.connectToNew(dv, mt, unconnectRate, time)) {

                                        avgConnectorLifeTime += dv.connectorLifeTime;
                                        connectorLifeTimeCounter += 1;

                                        //Master.log("connected vesicles PBC " + dv.connectorLifeTime);

                                        if (!dv.connectTo(kvesicle)) {
                                            Master.exit("connection failure");
                                        }

                                    } else {
                                        Master.exit("connection failure");
                                    }

                                }

                            }

                        }

                        kvesicle = kvesicle.nextReady;

                    }

                }

            }

        }

    }

    public void unconnectVesicles() {

        DiffusantVesicle dv, dv2;

        if (unconnectRate <= 0) {
            return;
        }

        for (int i = 0; i < diffusant.length; i++) {

            if ((diffusant[i] == null) || (diffusant[i].vesicle == null)) {
                continue;
            }

            for (int j = 0; j < diffusant[i].vesicle.length; j++) {

                dv = diffusant[i].vesicle[j];

                if (dv.connectTo == null) {
                    continue;
                }

                for (int k = 0; k < dv.connectTo.length; k++) {
                    
                    if (dv.connectTo[k] == null) {
                        continue;
                    }

                    if ((dv.connectorOffTime[k] > 0) && (time >= dv.connectorOffTime[k])) {
                        dv2 = dv.connectTo[k];
                        dv2.unconnectFrom(dv);
                        dv.connectTo[k] = null;
                        dv.connectorOffTime[k] = 0;
                    }

                }

            }

        }

    }

    public double[] computeNumberOfConnectors() {

        double avg = 0, connected = 0;
        int n, ntotal = 0;

        double[] results = new double[2];

        for (int i = 0; i < diffusant.length; i++) {

            if ((diffusant[i] == null) || (diffusant[i].vesicle == null)) {
                continue;
            }

            for (int j = 0; j < diffusant[i].vesicle.length; j++) {

                n = diffusant[i].vesicle[j].numberOfConnections(true);
                avg += n;
                ntotal++;

                if (n > 0) {
                    connected++;
                }

            }

        }

        if (ntotal > 0) {
            avg /= ntotal;
            connected /= ntotal;
        } else {
            avg = 0;
            connected = 0;
        }

        //Master.log("" + avg + " connectors / vesicle");
        //Master.log("" + connected + " fraction connected");

        results[0] = avg;
        results[1] = connected;

        return results;

    }

    public double avgFirstCollision() {

        double avg = 0, count = 0, sqrd = 0;

        for (int i = 0; i < diffusant.length; i++) {

            if ((diffusant[i] == null) || (diffusant[i].vesicle == null)) {
                continue;
            }

            for (int j = 0; j < diffusant[i].vesicle.length; ++j) {

                if ((diffusant[i].vesicle[j] == null) || !diffusant[i].vesicle[j].mobile ) {
                    continue;
                }

                if (diffusant[i].vesicle[j].firstCollision > 0) {
                    avg += diffusant[i].vesicle[j].firstCollision;
                    sqrd += diffusant[i].vesicle[j].sqrDisplacement;
                    count++;
                    //Master.log("" + diffusant[i].vesicle[j].firstCollision);
                }

            }

        }

        if (count == 0) {
            Master.log("no collisions");
        } else {

            avg /= count;
            sqrd /= count;

            Master.log("avg first collision = " + avg + " ms (n=" + count + ")");
            Master.log("square displacement = " + sqrd + " um^2");
            Master.log("D = " + (sqrd / (6 * avg)) + " um^2/ms");

        }

        return avg;

    }

    public double localDensityAll(boolean print) {

        double d, avg = 0, avgD = 0, avgDD0 = 0, avgStep = 0, count = 0;
        double min = 99999, max = 0;

        DiffusantVesicle dv;

        //if (immobileVesicleFraction > 0) {
        //    Master.exit("aborted MC simulation: immobile fraction not allowed with local density");
        //}

        if (diffusant == null) {
            return Double.NaN;
        }

        for (int i = 0; i < diffusant.length; i++) {

            if ((diffusant[i] == null) || (diffusant[i].vesicle == null)) {
                continue;
            }

            for (int j = 0; j < diffusant[i].vesicle.length; ++j) {

                if (diffusant[i].vesicle[j] == null) {
                    continue;
                }

                if (diffusant[i].vesicle[j].insideGeometry) {

                    dv = diffusant[i].vesicle[j];

                    if (PBC) {
                        d = localDensityPBC(dv);
                    } else if (hydrodynamicsLocalDVoxels) {
                        d = localDensityVoxels(dv);
                    } else {
                        d = localDensity(dv);
                    }

                    if (print && (d > 0)) {
                        avg += d;
                        avgStep += dv.localStep3;
                        count++;
                        min = Math.min(min, d);
                        max = Math.max(max, d);
                    }

                }

            }

        }

        avg /= count;
        avgStep /= count;

        if (print) {
            Master.log("local density avg = " + avg);
            Master.log("local density min = " + min);
            Master.log("local density max = " + max);
            //Master.log("local D/D0 avg = " + avgDD0 / count);
            Master.log("local Dshort avg = " + DiffusantVesicle.D3(avgStep, project.dt));
            Master.log("local r-step avg = " + avgStep);
        }

        return avg;

    }

    public double localDensity(DiffusantVesicle dv) {

        int xVoxel, yVoxel, zVoxel, nVoxels, nRadiiExt;
        int x0, y0, z0, x1, y1, z1;
        double dx, dy, dz;
        double d, h, Rext, r, v, totalV, sumV, volumeCap = 0;
        double localMobileVolumeFraction, Ds;

        Voxel voxel;
        DiffusantVesicle ivesicle;

        if ((dv == null) || (geometry.voxelSpace == null)) {
            return Double.NaN;
        }

        dz = geometry.z2 - dv.z;

        nRadiiExt = 4; // Rext
        Rext = nRadiiExt * dv.radius;
        h = Rext - dz;

        if (h > 0) {
            volumeCap = (Math.PI * h * h / 3) * (3 * Rext - h);
        }

        totalV = 4.0 * Math.PI * Rext * Rext * Rext / 3.0;
        sumV = vesicleVolume;

        xVoxel = (int) geometry.computeVoxelX(dv.x);
        yVoxel = (int) geometry.computeVoxelY(dv.y);
        zVoxel = (int) geometry.computeVoxelZ(dv.z);

        nVoxels = (int) Math.ceil(Rext / project.dx);

        x0 = xVoxel - nVoxels;
        y0 = yVoxel - nVoxels;
        z0 = zVoxel - nVoxels;

        x1 = xVoxel + nVoxels;
        y1 = yVoxel + nVoxels;
        z1 = zVoxel + nVoxels;

        x0 = Math.max(x0, 0);
        y0 = Math.max(y0, 0);
        z0 = Math.max(z0, 0);

        x1 = Math.min(x1, geometry.voxelSpace.length - 1);
        y1 = Math.min(y1, geometry.voxelSpace[0].length - 1);
        z1 = Math.min(z1, geometry.voxelSpace[0][0].length - 1);

        for (int i = x0; i <= x1; i++) {
            for (int j = y0; j <= y1; j++) {
                for (int k = z0; k <= z1; k++) {

                    voxel = geometry.voxelSpace[i][j][k];
                    ivesicle = voxel.firstReady;

                    while (ivesicle != null) {

                        if ((ivesicle != dv) && ivesicle.insideGeometry && ivesicle.mobile) {

                            dx = dv.x - ivesicle.x;
                            dy = dv.y - ivesicle.y;
                            dz = dv.z - ivesicle.z;

                            d = Math.sqrt(dx * dx + dy * dy + dz * dz);
                            r = ivesicle.radius;

                            if (d < Rext - r) {
                                sumV += 4.0 * Math.PI * r * r * r / 3.0;
                            } else if (d < Rext + r) {
                                v = Math.PI * (Rext + r - d) * (Rext + r - d) * (d * d + 2 * d * r - 3 * r * r + 2 * d * Rext + 6 * r * Rext - 3 * Rext * Rext) / (12 * d);
                                sumV += v;
                            }

                        }

                        ivesicle = ivesicle.nextReady;

                    }

                }
            }
        }

        localMobileVolumeFraction =  sumV / (totalV - volumeCap);
        
        Ds = DiffusantVesicle.Dratio_short(localMobileVolumeFraction, immobileVolumeFraction);
        dv.localD = Dcyto * Ds;
        dv.localStep3 = DiffusantVesicle.step3(dv.localD, project.dt);
        dv.DsDff = Ds / DiffusantVesicle.Dff_short_Banchio(localMobileVolumeFraction);

        return localMobileVolumeFraction;

    }

    public double localDensityVoxels(DiffusantVesicle dv) {

        int xVoxel, yVoxel, zVoxel, nVoxels, nRadiiExt;
        int x0, y0, z0, x1, y1, z1;
        double dx, dy, dz;
        double d, h, Rext, r, totalV = 0, sumV_mobile, sumV_immobile = 0;
        double localMobileVolumeFraction, localImmobileVolumeFraction, Ds;

        Voxel voxel;
        DiffusantVesicle ivesicle;

        if ((dv == null) || (geometry.voxelSpace == null)) {
            return Double.NaN;
        }

        nRadiiExt = 4; // Rext
        Rext = nRadiiExt * dv.radius;
        
        sumV_mobile = vesicleVolume;

        xVoxel = (int) geometry.computeVoxelX(dv.x);
        yVoxel = (int) geometry.computeVoxelY(dv.y);
        zVoxel = (int) geometry.computeVoxelZ(dv.z);

        nVoxels = (int) Math.ceil(Rext / project.dx);

        x0 = xVoxel - nVoxels;
        y0 = yVoxel - nVoxels;
        z0 = zVoxel - nVoxels;

        x1 = xVoxel + nVoxels;
        y1 = yVoxel + nVoxels;
        z1 = zVoxel + nVoxels;

        x0 = Math.max(x0, 0);
        y0 = Math.max(y0, 0);
        z0 = Math.max(z0, 0);

        x1 = Math.min(x1, geometry.voxelSpace.length - 1);
        y1 = Math.min(y1, geometry.voxelSpace[0].length - 1);
        z1 = Math.min(z1, geometry.voxelSpace[0][0].length - 1);

        for (int i = x0; i <= x1; i++) {
            for (int j = y0; j <= y1; j++) {
                for (int k = z0; k <= z1; k++) {

                    voxel = geometry.voxelSpace[i][j][k];

                    dx = dv.x - voxel.x;
                    dy = dv.y - voxel.y;
                    dz = dv.z - voxel.z;

                    d = Math.sqrt(dx * dx + dy * dy + dz * dz);

                    if (d > Rext) {
                        continue;
                    }

                    if (voxel.isSpace) {
                        totalV += voxelVolume;
                    }

                    ivesicle = voxel.firstReady;

                    while (ivesicle != null) {

                        if ((ivesicle != dv) && ivesicle.insideGeometry) {

                            if (ivesicle.mobile) {
                                sumV_mobile += vesicleVolume;
                            } else {
                                sumV_immobile += vesicleVolume;
                            }

                        }

                        ivesicle = ivesicle.nextReady;

                    }

                }
            }
        }

        //localDensity =  sumV / (totalV - volumeCap);
        localMobileVolumeFraction = sumV_mobile / totalV;
        localImmobileVolumeFraction = sumV_immobile / totalV;

        Ds = DiffusantVesicle.Dratio_short(localMobileVolumeFraction, immobileVolumeFraction);
        
        dv.localD = Dcyto * Ds;
        dv.localStep3 = DiffusantVesicle.step3(dv.localD, project.dt);
        dv.DsDff = Ds / DiffusantVesicle.Dff_short_Banchio(localMobileVolumeFraction);

        if (Double.isNaN(dv.localStep3)) {
            dv.localStep3 = dv.step3;
            //Master.log("NaN localSteps: " + localMobileVolumeFraction + "," + dv.localStep3);
        }

        return localMobileVolumeFraction;

    }

    public double localDensityPBC(DiffusantVesicle dv) {

        int xVoxel, yVoxel, zVoxel, nVoxels, nRadiiExt;
        int x0, y0, z0, x1, y1, z1, ii, jj, kk;
        double dx, dy, dz;
        double d, Rext, r, v, totalV, sumV;
        double xdir, ydir, zdir;
        double localMobileVolumeFraction, Ds;

        Voxel voxel;
        DiffusantVesicle ivesicle;

        if ((dv == null) || (geometry.voxelSpacePBC == null)) {
            return Double.NaN;
        }

        nRadiiExt = 4; // Rext
        Rext = nRadiiExt * dv.radius;
        totalV = 4.0 * Math.PI * Rext * Rext * Rext / 3.0;
        sumV = vesicleVolume; // all vesicles with same radius

        xVoxel = (int) geometry.computeVoxelX(dv.x);
        yVoxel = (int) geometry.computeVoxelY(dv.y);
        zVoxel = (int) geometry.computeVoxelZ(dv.z);

        nVoxels = (int) Math.ceil(Rext / project.dx);

        x0 = xVoxel - nVoxels;
        y0 = yVoxel - nVoxels;
        z0 = zVoxel - nVoxels;

        x1 = xVoxel + nVoxels;
        y1 = yVoxel + nVoxels;
        z1 = zVoxel + nVoxels;

        for (int i = x0; i <= x1; i++) {
            for (int j = y0; j <= y1; j++) {
                for (int k = z0; k <= z1; k++) {

                    if (i < 0) {
                        ii = geometry.voxelSpacePBC.length + i;
                        xdir = -1;
                    } else if (i >= geometry.voxelSpacePBC.length) {
                        ii = i - geometry.voxelSpacePBC.length;
                        xdir = 1;
                    } else {
                        ii = i;
                        xdir = 0;
                    }

                    if (j < 0) {
                        jj = geometry.voxelSpacePBC[0].length + j;
                        ydir = -1;
                    } else if (j >= geometry.voxelSpacePBC[0].length) {
                        jj = j - geometry.voxelSpacePBC[0].length;
                        ydir = 1;
                    } else {
                        jj = j;
                        ydir = 0;
                    }

                    if (k < 0) {
                        kk = geometry.voxelSpacePBC[0][0].length + k;
                        zdir = -1;
                    } else if (k >= geometry.voxelSpacePBC[0][0].length) {
                        kk = k - geometry.voxelSpacePBC[0][0].length;
                        zdir = 1;
                    } else {
                        kk = k;
                        zdir = 0;
                    }

                    voxel = geometry.voxelSpacePBC[ii][jj][kk];
                    ivesicle = voxel.firstReady;

                    while (ivesicle != null) {

                        if ((ivesicle != dv) && ivesicle.insideGeometry && ivesicle.mobile) {

                            dx = dv.x - (xdir * 2 * geometry.x2 + ivesicle.x);
                            dy = dv.y - (ydir * 2 * geometry.y2 + ivesicle.y);
                            dz = dv.z - (zdir * 2 * geometry.z2 + ivesicle.z);

                            d = Math.sqrt(dx * dx + dy * dy + dz * dz);
                            r = ivesicle.radius;

                            if (d < Rext - r) {
                                sumV += 4.0 * Math.PI * r * r * r / 3.0;
                            } else if (d < Rext + r) {
                                v = Math.PI * (Rext + r - d) * (Rext + r - d) * (d * d + 2 * d * r - 3 * r * r + 2 * d * Rext + 6 * r * Rext - 3 * Rext * Rext) / (12 * d);
                                sumV += v;
                            }

                        }

                        ivesicle = ivesicle.nextReady;

                    }

                }
            }
        }

        localMobileVolumeFraction = sumV / totalV;
        Ds = DiffusantVesicle.Dratio_short(localMobileVolumeFraction, immobileVolumeFraction);
        
        dv.localD = Dcyto * Ds;
        dv.localStep3 = DiffusantVesicle.step3(dv.localD, project.dt);
        dv.DsDff = Ds / DiffusantVesicle.Dff_short_Banchio(localMobileVolumeFraction);

        return localMobileVolumeFraction;

    }

    public int countVesiclesWithinPSF(){

        int dnum;
        int count = 0;
        double e2 = 0.135335;

        if ( !frapOn || (detector == null) || (dPSF == -1)) {
            return 0; // nothing to do
        }

        for (int i = 0; i < detector.length; i++) {

            if (detector[i] == null) {
                continue;
            }

            dnum = detector[i].diffusantNum;

            if ((dnum < 0) || (dnum >= diffusant.length)) {
                continue;
            }

            if (diffusant[dnum].vesicle == null) {
                continue;
            }

            for (int j = 0; j < diffusant[dnum].vesicle.length; j++) {

                if (!diffusant[dnum].vesicle[j].insideGeometry) {
                    continue;
                }

                if (diffusant[dnum].vesicle[j].voxel.PSFd > e2) {
                    count++;
                }

            }

        }

        Master.log("vesicles within cPSF n = " + count);

        return count;

    }

    public void detectFluorescence(){

        int dnum;
        double w, avg = 0, count = 0;

        if ( !frapOn || !saveFluorescence || (detector == null) || (dPSF == -1)) {
            return; // nothing to do
        }

        for (int i = 0; i < detector.length; i++) {

            if (detector[i] == null) {
                continue;
            }

            dnum = detector[i].diffusantNum;

            if ((dnum < 0) || (dnum >= diffusant.length)) {
                continue;
            }

            if (diffusant[dnum].vesicle == null) {
                continue;
            }

            for (int j = 0; j < diffusant[dnum].vesicle.length; j++) {

                if (!diffusant[dnum].vesicle[j].insideGeometry) {
                    continue;
                }

                if (freeDiffusion || PBC) {
                    w = 1.0;
                } else {
                    w = diffusant[i].vesicle[j].voxel.PSFweight; // extra weighting due to non-space voxels
                }

                avg += diffusant[dnum].vesicle[j].fluorescence * diffusant[dnum].vesicle[j].voxel.PSFd * w;
                count++;

            }

            if ((PSFd.avg > 0) && (count > 0)) {
                //avg /= (PSFd.avg * count);
                avg /= PSFd.avg;
            } else {
                avg = 0;
            }

            detector[i].save.saveData(avg);

        }

    }

    public double meanFluorescence(boolean printResults) {

        double w, avg = 0, count = 0;

        if ((!frapOn) || (diffusant == null)) {
            return 0;
        }

        for (int i = 0; i < diffusant.length; i++) {

            if ((diffusant[i] == null) || (diffusant[i].vesicle == null)) {
                continue;
            }

            for (int j = 0; j < diffusant[i].vesicle.length; j++) {

                if (!diffusant[i].vesicle[j].insideGeometry) {
                    continue;
                }

                if (freeDiffusion || PBC) {
                    w = 1.0;
                } else {
                    w = diffusant[i].vesicle[j].voxel.PSFweight;
                }

                avg += diffusant[i].vesicle[j].fluorescence * diffusant[i].vesicle[j].voxel.PSFd * w;
                count++;

            }

        }

        if ((PSFd.avg > 0) && (count > 0)) {
            //avg /= (PSFd.avg * count);
            avg /= PSFd.avg;
        } else {
            avg = 0;
        }

        if (printResults) {
            Master.log("mean fluorescence: " + avg);
        }

        return avg;

    } 

    public void saveArray(double[] data, String fxn, String saveTag) {

        Master.writeToFiles = true;

        Save save = new Save(project);

        save.saveWhileComputing = true;
        save.save2BinaryFile = false;
        save.save2TextFile = true;

        save.fileName(fxn, saveTag);
        save.xdim = project.spaceUnits;
        save.ydim = project.timeUnits;

        save.samples2save = data.length;

        if (!save.init(fxn, geometry, -1, 1)) {
            return;
        }

        for (int i = 0; i < data.length; i++) {
            save.writeData(data[i]);
        }

        if (save != null) {
            save.finish(fxn, geometry, -1);
        }

    }

    public boolean addToVoxelList(DiffusantVesicle dv) {

        if (dv.voxel == null){
            return false;
        }

        dv.nextReady = dv.voxel.firstReady; // save existing vesicles
        dv.voxel.firstReady = dv; // replace with new vesicles, creating chain

        return true;

    }

    public boolean removeFromVoxelList(DiffusantVesicle dv, Voxel voxel) {

        if (voxel == null){
            return false;
        }

        boolean found = false;

        DiffusantVesicle vtest = voxel.firstReady;
        DiffusantVesicle vnext;
        DiffusantVesicle vhold = null;

        while (vtest != null) {

            vnext = vtest.nextReady;

            if (vtest == dv) {
                if (vhold == null) {
                    voxel.firstReady = vnext; // dv was first in list
                } else {
                    vhold.nextReady = vnext;
                }
                found = true;
            }

            vhold = vtest;
            vtest = vnext;

        }

        return found;

    }

    double ranGauss() { // 0.0 mean, 1.0 stdv

        double v1, v2, w, step;

        if (saveRanGauss != 999999) {
            step = saveRanGauss;
            saveRanGauss = 999999;
            return step;
        }

        do {
            v1 = (2.0 * mt.nextDouble() - 1.0);
            v2 = (2.0 * mt.nextDouble() - 1.0);
            w = v1 * v1 + v2 * v2;
        } while (w >= 1.0);

        w = Math.sqrt((-2.0 * Math.log(w)) / w);

        step = v1 * w; // first value
        saveRanGauss = v2 * w; // save second value

        return step;

    }

    public boolean moveVesicleGauss(DiffusantVesicle dv) {

        double step3;

        if (removingVesicleOverlap) {
            step3 = removeVesicleOverlapStep3;
        } else if (hydrodynamicsLocalD){
            step3 = dv.localStep3;
        } else {
            step3 = dv.step3;
        }

        stepx = ranGauss() * step3;
        stepy = ranGauss() * step3;
        stepz = ranGauss() * step3;

        return setVesicleLocation(dv, dv.x + stepx, dv.y + stepy, dv.z + stepz, false);

    }

    public boolean moveVesicleGaussHydroWallz(DiffusantVesicle dv) {

        double step3;
        double z, b_ll = 1, b_T = 1;

        if (hydrodynamicsLocalD) {
            step3 = dv.localStep3;
        } else {
            step3 = dv.step3;
        }

        if (removingVesicleOverlap) {

            step3 = removeVesicleOverlapStep3;

        } else {

            z = geometry.z2 - dv.z;
            
            b_ll = hydroWall_ll(dv.radius, z);
            b_T = hydroWall_T(dv.radius, Math.abs(z - dv.radius));

            b_ll = 1 / (1 + dv.DsDff * ((1 / b_ll) - 1)); // Michailidou et al. 2009
            b_T = 1 / (1 + dv.DsDff * ((1 / b_T) - 1)); // Michailidou et al. 2009

            b_ll = Math.sqrt(b_ll);
            b_T = Math.sqrt(b_T);

        }

        stepx = ranGauss() * step3 * b_ll;
        stepy = ranGauss() * step3 * b_ll;
        stepz = ranGauss() * step3 * b_T;
        
        return setVesicleLocation(dv, dv.x + stepx, dv.y + stepy, dv.z + stepz, false);

    }

    public void moveVesiclesCichocki(boolean moveAll) {

        boolean outOfBounds, overlap, sameVoxel, moveConnected, ok;

        DiffusantVesicle dv, odv;

        if (diffusant == null) {
            return;
        }

        for (int i = 0; i < diffusant.length; i++) {

            if ((diffusant[i] == null) || (diffusant[i].vesicle == null)) {
                continue;
            }

            diffusant[i].shuffleVesicles();

            for (int j = 0; j < diffusant[i].vesicle.length; j++) {

                if (!diffusant[i].vesicle[j].insideGeometry) {
                    continue;
                }

                if (!moveAll && !diffusant[i].vesicle[j].mobile) {
                    continue;
                }

                dv = diffusant[i].vesicle[j];
                testVesicle.copy(dv);

                if (hydroWallZ) {
                    if (!moveVesicleGaussHydroWallz(testVesicle)) {
                        continue; // dont move
                    }
                } else {
                    if (!moveVesicleGauss(testVesicle)) {
                        continue; // dont move
                    }
                }

                outOfBounds = outOfBounds(testVesicle);

                if (PBC && outOfBounds) {
                    outOfBounds = wrapAtBorder(testVesicle);
                }

                if (outOfBounds) {
                    continue; // dont move
                }

                if (testMitochondriaOverlap(testVesicle)) {
                    if ((time > 0) && (dv.firstCollision == 0)) {
                        dv.firstCollision = time;
                        dv.sqrDisplacement = dv.squareDisplacement();
                    }
                    continue; // dont move
                }

                if (!freeDiffusion) {

                    odv = testVesicleOverlap(testVesicle, dv);

                    overlap = odv != null;

                    if (overlap) {
                        if ((time > 0) && (dv.firstCollision == 0)) {
                            dv.firstCollision = time;
                            dv.sqrDisplacement = dv.squareDisplacement();
                        }
                        continue;
                    }

                }

                if (connectVesicles && (dv.connectTo != null)) {

                    moveConnected = true;

                    for (int k = 0; k < dv.connectTo.length; k++) {

                        if (dv.connectTo[k] == null) {
                            continue;
                        }

                        if (!testVesicle.overlap(dv.connectTo[k], connectorLength)) {
                            moveConnected = false;
                            break;
                        }

                    }

                    if (!moveConnected) {
                        continue;
                    }

                }

                sameVoxel = false;

                if (testVesicle.voxel == diffusant[i].vesicle[j].voxel) {
                    sameVoxel = true;
                }

                if (!sameVoxel) {
                    removeFromVoxelList(diffusant[i].vesicle[j], diffusant[i].vesicle[j].voxel);
                }

                diffusant[i].vesicle[j].copy(testVesicle);

                if (!sameVoxel) {
                    addToVoxelList(diffusant[i].vesicle[j]);
                }

            }

        }

    }

    public void drift() {

        double x, y, z;
        boolean outOfBounds, sameVoxel;

        if (diffusant == null) {
            return;
        }

        if (!PBC) {
            Master.exit("drift error: PBC is not on");
        }

        if (mito != null) {
            for (int i = 0; i < mito.length; i++) {
                x = mito[i].xCenter + driftDX;
                y = mito[i].yCenter + driftDY;
                z = mito[i].zCenter + driftDZ;
                mito[i].setCenter(x, y, z);
            }
        }

        for (int i = 0; i < diffusant.length; i++) {

            if ((diffusant[i] == null) || (diffusant[i].vesicle == null)) {
                continue;
            }

            for (int j = 0; j < diffusant[i].vesicle.length; j++) {

                if (!diffusant[i].vesicle[j].insideGeometry) {
                    continue;
                }

                testVesicle.copy(diffusant[i].vesicle[j]);

                setVesicleLocation(testVesicle, testVesicle.x + driftDX, testVesicle.y  + driftDY, testVesicle.z + driftDZ, false);

                outOfBounds = outOfBounds(testVesicle);

                if (outOfBounds) {
                    outOfBounds = wrapAtBorder(testVesicle);
                    testVesicle.fluorescence = 1.0; // reset F to simulate non-frapped vesicles moving into PSF
                }

                if (outOfBounds) {
                    Master.exit("drift error: wrap at border error");
                }

                sameVoxel = false;

                if (testVesicle.voxel == diffusant[i].vesicle[j].voxel) {
                    sameVoxel = true;
                }

                if (!sameVoxel) {
                    removeFromVoxelList(diffusant[i].vesicle[j], diffusant[i].vesicle[j].voxel);
                }

                diffusant[i].vesicle[j].copy(testVesicle);

                if (!sameVoxel) {
                    addToVoxelList(diffusant[i].vesicle[j]);
                }

            }

        }

        return;

    }

    public boolean wrapAtBorder(DiffusantVesicle dv) {

        double x = dv.x;
        double y = dv.y;
        double z = dv.z;

        if (x > geometry.x2) {
            x -= 2 * geometry.x2;
            dv.x0 -= 2 * geometry.x2;
        } else if (x < geometry.x1) {
            x -= 2 * geometry.x1;
            dv.x0 -= 2 * geometry.x1;
        }

        if (y > geometry.y2) {
            y -= 2 * geometry.y2;
            dv.y0 -= 2 * geometry.y2;
        } else if (y < geometry.y1) {
            y -= 2 * geometry.y1;
            dv.y0 -= 2 * geometry.y1;
        }

        if (z > geometry.z2) {
            z -= 2 * geometry.z2;
            dv.z0 -= 2 * geometry.z2;
        } else if (z < geometry.z1) {
            z -= 2 * geometry.z1;
            dv.z0 -= 2 * geometry.z1;
        }

        setVesicleLocation(dv, x, y, z, false);

        return outOfBounds(dv);

    }

    public boolean initSave() {

        project.checkDirectory();

        return true;

    }

    public boolean finishSave() {

        return true;
    }

    @Override
    public boolean addUser(ParamVector pv) {

        if (pv == null) {
            return false;
        }

        super.addUser(pv);

        return true;

    }

    @Override
    public boolean createVector(boolean close) {

        if (!super.createVector(false)) {
            return false;
        }

        if (mito != null) {

            addBlankParam();

            for (int i = 0; i < mito.length; i++) {
                mito[i].createVector(true);
                addVector(mito[i].getVector());
                mito[i].addUser(this);
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
    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        String n = o.getName();

        if (!(o.paramVector instanceof RunMonteCarlo)) {
            return false;
        }

        if (n.equalsIgnoreCase("minVesicleStep")) {
            if (v <= 0) {
                return false;
            }
            minVesicleStep = v;
            return true;
        }
        if (n.equalsIgnoreCase("freeDiffusion")) {
            if (v == 1) {
                freeDiffusion = true;
            } else {
                freeDiffusion = false;
            }
            return true;
        }
        if (n.equalsIgnoreCase("PBC")) {
            if (v == 1) {
                PBC = true;
            } else {
                PBC = false;
            }
            return true;
        }
        if (n.equalsIgnoreCase("connectVesicles")) {
            if (v == 1) {
                connectVesicles = true;
            } else {
                connectVesicles = false;
            }
            return true;
        }
        if (n.equalsIgnoreCase("maxNumConnectors")) {
            if (v < 0) {
                return false;
            }
            maxNumConnectors = (int) v;
            return true;
        }
        if (n.equalsIgnoreCase("connectRate")) {
            if (v < 0) {
                return false;
            }
            connectRate = v;
            return true;
        }
        if (n.equalsIgnoreCase("unconnectRate")) {
            if (v < 0) {
                return false;
            }
            unconnectRate = v;
            return true;
        }
        if (n.equalsIgnoreCase("hydrodynamicsLocalD")) {
            if (v == 1) {
                hydrodynamicsLocalD = true;
            } else {
                hydrodynamicsLocalD = false;
            }
            return true;
        }
        if (n.equalsIgnoreCase("hydroWallZ")) {
            if (v == 1) {
                hydroWallZ = true;
            } else {
                hydroWallZ = false;
            }
            return true;
        }
        if (n.equalsIgnoreCase("MSDspatialTbgn")) {
            MSDspatialTbgn = v;
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

        if (!(o.paramVector instanceof RunMonteCarlo)) {
            return false;
        }

        if (n.equalsIgnoreCase("name")) {
            name = s;
            return true;
        }
        return false;
    }
}
