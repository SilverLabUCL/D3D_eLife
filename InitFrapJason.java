package ucl.silver.d3d.init;

import ucl.silver.d3d.core.*;
import java.util.Random;

/**
 * <p>Title: D3D</p>
 *
 * <p>Description: 3D Diffusion Simulator</p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: University College London</p>
 *
 * @author Jason Rothman
 * @version 1.0
 */
public class InitFrapJason extends InitProject {

    public String directory = "/Jason/D3D/Simulations/";
    public String folder = "Testing";

    public double simTime = 8.0; // msec
    public double saveRate = 20; // kHz
    public double printRate = 0.002; // kHz

    public boolean saveWhileComputing = false;
    
    public double temperature = 37.0; // C
    public double stabilityFD = 0.4;
    
    public double vesicleRadius = 0.025; // um
    
    public String d1_name = "Vesicles";
    public double d1_C = 1.0;
    public double d1_D = 0.060e-3;

    public String d2_name = "";
    public double d2_C = -1;
    public double d2_D = -1;

    public boolean FRAP_on = false;
    public double FRAP_kPhoto = 2.0; // photolysis k factor
    public double FRAP_onset = 100; // ms
    public double FRAP_bleach_length = 0.5; // ms
    public boolean FRAP_bleach_on = true; // bleach pulse
    public double FRAP_small_length = 2.0; // ms
    //public double FRAP_small_amp_ratio = 0.005; // WRONG
    public double FRAP_small_amp_ratio = 0.0014; // ratio measured laser power after objective
    public boolean FRAP_small_on = true;

    //public boolean mitoNonSpaceVoxels = false;
    //public double mitoRadius = (0.25 / 0.89) / 2.0; // um Palay
    //public double mitoAxialRatio = 2.0 / 0.25; // Palay
    //public double mitoVolumeFraction = 0.28; // Zoltan average
    //public double mitoVolumeFractionTolerance = 0.02;
    //public double mitoVolumeFractionActual = 0;
    //public boolean mitoVolumeFractionExact = false;
    //public boolean mitoPreset = false;
    //public boolean mitoCylinder = true;
    //public int mito_ijkselect = -1;

    public int iPSFselect = 1; // (0) Gauss (1) Wilson (2) Torok (3) Impulse
    public int iPSF, dPSF;

    public String psf_file = "D3D_Wilson_PSF.d3d";

    public boolean writePSF = false;
    public boolean readPSF = false;

    public double refrac_water = 1.338; // at 37C
    public double lambda_illum = 0.488; // um
    public double NA = 0.90; // fit to Zoltan's iPSF
    public double rscale = 1;

    public double excite_fwhm_xy = 0.343; // Gauss PSF
    public double excite_fwhm_z = 1.285;

    public double detect_fwhm_xy = 0.257; // measured Gauss PSF
    public double detect_fwhm_z = 0.886; // measured Gauss PSF

    //public double detect_fwhm_xy = 0.213; // theoretical Gauss PSF
    //public double detect_fwhm_z = 0.710; // theoretical Gauss PSF

    int d1_num = -1;
    int d2_num = -1;
    int xVoxels, yVoxels, zVoxels;

    public double xvcenter, yvcenter, zvcenter;

    double zGradDensity[] = {6356, 2864, 1409, 673, 276};
    //double zGradDensity[] = {6652, 3243, 1254, 1254, 236};
    double zGradVolumeFraction[] = {0, 0, 0, 0, 0};
    int zGradNumVesicles[] = {0, 0, 0, 0, 0};
    CoordinatesVoxels[] zGradCoordinates;

    PSFgauss gpsf = null;
    PSFwilson wpsf = null;
    PSFtorok tpsf = null;
    PSF epsf = null;

    DiffusantPhoto frap = null;
    DiffusantPhoto frap2 = null;

    PulseTimer timer = null;

    public long seed = 8682522807148012L + System.nanoTime();
    
    Random ran = new Random(seed);

    InitGeometryMito igm;

    public InitFrapJason(Project p) {
        super(p);
        igm = new InitGeometryMito(p);
        String[] flist = {"initMonteCarloFrap", "initMonteCarloMSD", "initFiniteDifferenceFrap"};
        initFuncList = flist;
        createVector(true);
    }

    @Override
    public boolean initFunction(String initSelect) {

        double rf, tau;

        int i = initFunctionNum(initSelect);

        if (!Master.foundMainStartUpArguments) {
            i = 2;
        }

        switch (i) {

            case 0:
                return initMonteCarloFrap();
                //return initMonteCarloFrapMovie();

            case 1:
                return initMonteCarloMSD();

            case 2:
                return initFiniteDifferenceFrap();
                
            case 3:
                return initCichocki();

            case 4:
                return InitFrapJasonMisc.writeBatchesFRAP();
                //return InitFrapJasonMisc.writeBatchesMSD();
                //return InitFrapJasonMisc.writeBatchesFD();
                //return InitFrapJasonMisc.writeQdel();

            case 5:

                //InitFrapJasonMisc.callWriteBatchFiles2();

                
                tau = 0.06;

                for (int j = 0; j < 9; j++) {

                    rf = 0.01;

                    for (int k = 0; k < 9; k++) {

                        //Master.log("" + tau + ", " + rf);
                        InitFrapJasonMisc.callWriteBatchFiles2(tau, rf);

                        rf += 0.01;

                    }

                    tau += 0.01;

                }

                return true;

            case 6:
                //InitFrapJasonMisc.cleanUpFrap();
                //InitFrapJasonMisc.cleanUpFrap3("23p8", "23p75");
                //InitFrapJasonMisc.cleanUpFrap3("26p3", "26p25");
                //InitFrapJasonMisc.cleanUpFrapRenumber();
                //return InitFrapJasonMisc.cleanUpFrap2("d060");
                //InitFrapJasonMisc.cleanUpMSD("MSD", "MSD");
                //InitFrapJasonMisc.cleanUpMSD("ResidenceTime", "RT");
                return true;

            case 10:
                InitFrapJasonMisc.frapAllConfigs2();
                return true;

            default:
                Master.log("InitFrapJason.initFunction error: failed to find init function " + initSelect);
                return true; // error

        }

    }

    private boolean initCichocki() {

        double dx, lambda, step, t0;

        DiffusantVesicles dv;

        directory = "/Jason/FrapVGlut1/Simulations/Sims16_Cichocki/";

        switch(0) {
            case 0:
                lambda = 0.035;
                folder = "dx035";
                break;
            case 1:
                lambda = 0.025;
                folder = "dx025";
                break;
            case 2:
                lambda = 0.018;
                folder = "dx018";
                break;
            default:
                return true;
        }

        project.newMonteCarlo();
        project.name = "Monte Carlo FRAP";

        if (!Master.foundMainStartUpArguments) {
            project.directory = directory;
        }

        saveRate = 500;
        printRate = 0.2;

        RunMonteCarlo mc = project.monteCarlo;

        if (mc == null) {
            return true; // error
        }
        
        double cubeWidth = 1.0;

        double vesicleVolumeFraction = 0.20; // 40

        project.folder = "den" + Integer.toString((int)(vesicleVolumeFraction * 100)) + "_" + folder;

        dx = 2 * vesicleRadius;

        d1_name = "Vesicles";

        d1_D = 0.001;

        igm.xdim = cubeWidth;
        igm.ydim = cubeWidth;
        igm.zdim = cubeWidth;

        stabilityFD = 0.005;

        step = lambda * vesicleRadius;

        System.out.println("dx = " + step);

        t0 = vesicleRadius * vesicleRadius / ( 6 * d1_D);

        System.out.println("t0 = " + t0);

        simTime = 16 * t0 * 20;

        if (initProject(dx)) {
            return true;
        }

        if (igm.initGeometryMC(-1)) {
            return true;
        }

        dv = initVesicles(d1_D);

        dv.setVolumeFraction = vesicleVolumeFraction;
        dv.setImmobilePercent = 0; // 0.275
        dv.saveXYZ = false;
        dv.saveMSD = true;

        mc.minVesicleStep = step;
        mc.freeDiffusion = false;
        mc.PBC = true;

        if (!Master.foundMainStartUpArguments) {
            mc.initAll(); // // do not run when running batches
        }

        Master.addBatchList(26, "geometry.cubeWidth", "1.0,1.0,1.0,1.0");

        return false;

    }

    private boolean initMonteCarloFrap() {

        double dx;
        DiffusantVesicles dv;

        //directory = "/Jason/FrapVGlut1/Simulations/Sims12_FRAP/";
        directory = "/Jason/FrapVGlut1/Simulations/Testing/";
        //directory = "/Jason/FrapVGlut1/Simulations/Sims12_FRAP_MSD/Den44/";

        //folder = "testing";
        //folder = "hindered_noimmobile";
        //folder = "MinStep";
        //folder = "cubeWidth";
        //folder = "DT";
        //folder = "kPhoto";
        //folder = "D";
        //folder = "ImmobilePercent_1000nm";
        //folder = "ImmobilePercent_1200nm";
        //folder = "ImmobilePercent_1500nm";
        //folder = "ImmobilePercent_1800nm";
        //folder = "SmallPlusBleach";
        //folder = "SmallPulsesOnly/p005";
        //folder = "zLong";
        //folder = "zPSF_Long";

        project.newMonteCarlo();
        project.name = "Monte Carlo FRAP";

        if (!Master.foundMainStartUpArguments) {
            project.directory = directory;
            project.folder = folder;
        }

        RunMonteCarlo mc = project.monteCarlo;

        if (mc == null) {
            return true; // error
        }

        printRate = 0.004;

        //NA = 0.80;

        FRAP_on = true;
        FRAP_kPhoto = 2.0; // USE THIS ONE

        double baselinetime = 500;
        //double baselinetime = 1;

        FRAP_onset = baselinetime + 70; // first pulse is 70 ms before bleaching pulse
        simTime = 10000 + FRAP_onset;

        //FRAP_onset = 15;
        //simTimeStr = 1000 + FRAP_onset;

        //FRAP_onset = 1;
        //simTime = 10 + FRAP_onset;

        //FRAP_small_on = false;
        //FRAP_bleach_on = false;

        //double vesicleVolumeFraction = 0.26; // Zoltan (old)
        double vesicleVolumeFraction = 0.17; // Zoltan (new 02/09/2015)
        //double vesicleVolumeFraction = 0.10; // estimate for dispersed + clustered
        
        //double cubeWidth = 2.5;
        double cubeWidth = 3;

        //ellipsoid = true;
        
        igm.mitoVolumeFraction = 0.28;

        igm.mito_ijkselect = -2;
        //mito_ijkselect = 1;
        
        igm.mitoNonSpaceVoxels = true; // mito are non-space voxels
        //mitoPreset = true;

        //mc.mitoCoordinatesOn = true; // mito are coordinate objects like vesicles
        mc.setMitoVolumeFraction = igm.mitoVolumeFraction;
        //mc.mito = getMitochondriaPreset0();
        //mc.mito = getMitochondriaPreset1();
        //mc.mito = getMitochondriaPreset2();
        mc.mitoRadius = igm.mitoRadius;
        mc.mitoAxialRatio = igm.mitoAxialRatio;

        vesicleRadius = 0.049 / 2.0;
        dx = 0.05;

        d1_name = "Vesicles";

        //d1_D = 0.090e-3;
        d1_D = 0.060e-3;

        igm.xdim = cubeWidth;
        igm.ydim = cubeWidth;
        igm.zdim = cubeWidth;

        stabilityFD = 0.005;

        if (Master.foundMainStartUpArguments) {
            iPSFselect = 2;
        } else {
            iPSFselect = 0;//2; // (0) Gauss (1) Wilson (2) Torok (3) Impulse
        }

        //writePSF = true;
        //readPSF = true;
        //psf_file = "PSF_Torok.d3d";

        detect_fwhm_xy = 0.2551; // measured, bead deconvolved
        detect_fwhm_z = 0.9157; // measured, bead deconvolved

        //detect_fwhm_xy = 0.230021; // theoretical computed from Torok fit
        //detect_fwhm_z = 1.023993; // theoretical computed from Torok fit

        if (initProject(dx)) {
            return true;
        }
        if (igm.initGeometryMC(-1)) {
            return true;
        }
        if (initExcitationPSF()) {
            return true;
        }

        project.simTime = simTime;

        dv = initVesicles(d1_D);

        dv.colorReady.setColor("[r=147,g=210,b=21]");
        //dv.colorImmobile.setColor("[r=102,g=153,b=0]");
        //dv.colorReady.setColor("[r=255,g=255,b=255]");
        dv.colorImmobile.setColor("[r=150,g=150,b=150]");

        dv.color.setColor("[r=0,g=0,b=153]");
        dv.color.setGradientColor("[r=255,g=255,b=255]");

        dv.setVolumeFraction = vesicleVolumeFraction;
        dv.setImmobilePercent = 0.25; // 0.175; // 0.275;
        //dv.saveXYZ = true;
        //dv.saveMSD = true;

        //dv = initVesicles(d1_D * 0.5);
        //dv.numVesiclesReady = -1;

        if (initDetectors()) {
            return true;
        }

        //mc.minVesicleStep = 0.0024494897427831783; // can use this for all D
        //mc.minVesicleStep = 0.001;
        mc.minVesicleStep = 0.002;
        //mc.minVesicleStep = 0.005;

        //mc.freeDiffusion = true;
        //mc.PBC = true;

        mc.frapOn = true;
        mc.saveFluorescence = true;

        mc.setPSFSelect(iPSF, dPSF);

        boolean driftOn = false;

        if (driftOn) {
            mc.driftOn = true;
            mc.driftOnset = baselinetime;
            mc.PBC = true;
        }
        
        boolean cluster = false;

        if (cluster) {
            mc.clusterImmobileVesicles = true;
            mc.clusterImmobilePercent = dv.setImmobilePercent;
            dv.setImmobilePercent = 0.01; // start sparse
        }

        //sims2do();

        //fileDat2Bin();

        if (!Master.foundMainStartUpArguments) {
            mc.initAll(); // // do not run when running batches
        }
        
        //mc.testMoveVesicles(false, false);
        //mc.testMinDistance();

        //project.writeParamFile();

        //Master.addBatchList(0, "shape.cubeWidth", "0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6");
        //Master.addBatchList(0, "montecarlo.minVesicleStep", "0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001");

        //Master.addBatchList(0, "Vesicles.kPhoto", "3.5,3.5,3.5,3.5,3.5,3.5,3.5,3.5,3.5,3.5");

        //Master.addBatchList(0, "Vesicles.D", "0.008e-3,0.008e-3,0.008e-3,0.008e-3,0.008e-3,0.008e-3,0.008e-3,0.008e-3,0.008e-3,0.008e-3");
        //Master.addBatchList(10, "Vesicles.D", "0.008e-3,0.008e-3,0.008e-3,0.008e-3,0.008e-3,0.008e-3,0.008e-3,0.008e-3,0.008e-3,0.008e-3");

        //Batch[] batches = Master.project.batches;
        
        //if (Master.project.batches != null && Master.project.batches.length > 0) {
        //    for (int i = 0; i < batches.length; i++) {
        //        batches[i].folder = "c20_d010_im17p5_" + batches[i].batchNum;
        //    }
        //}

        return false;

    }

    private boolean initMonteCarloMSD() {

        double dx;
        DiffusantVesicles dv;

        //directory = "/Jason/FrapVGlut1/Simulations/Sims23_MSD_D090/NoMitoHI/Den26/";
        directory = "/Jason/FrapVGlut1/Simulations/Testing/";

        //folder = "c10_den26_im27p5_dx0020_0";
        folder = "msd";

        project.newMonteCarlo();
        project.name = "Monte Carlo MSD";

        if (!Master.foundMainStartUpArguments) {
            project.directory = directory;
            project.folder = folder;
        }

        RunMonteCarlo mc = project.monteCarlo;

        if (mc == null) {
            return true; // error
        }

        printRate = 0.02;
        saveRate = 20; // kHz

        //simTime = 12000;//300;
        simTime = 10;

        double cubeWidth = 2.0;
        //double cubeWidth = 1.0;

        igm.mitoNonSpaceVoxels = true; // mito are non-space voxels
        //mitoPreset = true;
        
        igm.mitoVolumeFraction = 0.28;
        //mitoVolumeFraction = 0;

        igm.mito_ijkselect = -2; // 1

        //mitoNonSpaceVoxels = true; // mito are non-space voxels
        //mitoPreset = true;

        //mc.mitoCoordinatesOn = true; // mito are coordinate objects like vesicles (SLOWER)
        mc.setMitoVolumeFraction = igm.mitoVolumeFraction;
        mc.mitoRadius = igm.mitoRadius;
        mc.mitoAxialRatio = igm.mitoAxialRatio;

        //cubeWidth = cubeWidth/2.0;

        vesicleRadius = 0.049 / 2.0;
        dx = 0.05;

        d1_name = "Vesicles";

        d1_D = 0.060e-3;

        igm.xdim = cubeWidth;
        igm.ydim = cubeWidth;
        igm.zdim = cubeWidth;

        if (initProject(dx)) {
            return true;
        }

        if (igm.initGeometryMC(-1)) {
            return true;
        }

        dv = initVesicles(d1_D);

        //dv.colorReady.setColor("[r=204,g=255,b=0]");
        //dv.colorImmobile.setColor("[r=102,g=153,b=0]");

        //dv.setVolumeFraction = 0.19;//0.26; // Zoltan
        dv.setVolumeFraction = 0.17; // Zoltan (new 02/09/2015)
        dv.setImmobilePercent = 0.25;
        //dv.saveXYZ = true;
        dv.saveMSD = true;
        //dv.saveDisplacementAfterCollisions = true;
        dv.init();

        //mc.minVesicleStep = 0.0024494897427831783; // can use this for all D
        mc.minVesicleStep = 0.0005;
        //mc.minVesicleStep = 0.0010;
        //mc.minVesicleStep = 0.0015;

        //mc.freeDiffusion = true;
        mc.PBC = true;
        //mc.wrapz = false;
        //mc.saveResidenceTime = true;
        //mc.hydroWallZ = true;
        //mc.hydrodynamicsDscale = true;
        //mc.hydrodynamicsLocalD = true;
        //mc.vesicleLattice = true;

        mc.Dcyto = d1_D;

        //mc.saveMSDtimes = new double[10];
        //mc.saveMSDtimes[0] = 50;
        //mc.saveMSDtimes[1] = 100;
        //mc.saveMSDtimes[2] = 500;
        //mc.saveMSDtimes[3] = 1000;
        //mc.saveMSDtimes[4] = 2000;
        //mc.saveMSDtimes[5] = 3000;
        //mc.saveMSDtimes[6] = 4000;
        //mc.saveMSDtimes[7] = 6000;
        //mc.saveMSDtimes[8] = 8000;
        //mc.saveMSDtimes[9] = 10000;

        boolean cluster = false;
        boolean connectVesicles = false;
        //boolean connectVesicles = true;

        if (mc.hydrodynamicsLocalD) {
            //Master.log("Dh = " + (d1_D * DiffusantVesicle.DD0short(dv.setVolumeFraction)));
        }

        if (cluster) {
            mc.clusterImmobileVesicles = true;
            mc.clusterImmobilePercent = dv.setImmobilePercent;
            dv.setImmobilePercent = 0.05; // start sparse
        }

        if (connectVesicles) {
            mc.connectVesicles = true;
            mc.maxNumConnectors = 3;
            mc.connectRate = 10; // P/ms
            mc.unconnectRate = 0.01; // 1.75 kHz
        }

        if (dv.setVolumeFraction > 0.3) {
            mc.overlapTrialLimit = 10;
        }

        if (!Master.foundMainStartUpArguments) {
            mc.initAll(); // // do not run when running batches
        }
        
        //double msd = mc.avgMinSqrDistance();
        //double tauShort = (msd - vesicleDiameter * vesicleDiameter)/d1_D;
        //msd = Math.sqrt(msd);

        //Master.log("average min distance between vesicles = " + msd);
        //Master.log("tau short = " + tauShort + " ms");

        return false;

    }

    private boolean initMonteCarloFrapMovie() {

        double dx;
        DiffusantVesicles dv;

        //directory = "/Jason/FrapVGlut1/Simulations/Sims12_FRAP/";
        directory = "/Jason/FrapVGlut1/Simulations/";

        folder = "Testing";

        project.newMonteCarlo();
        project.name = "Monte Carlo FRAP";

        if (!Master.foundMainStartUpArguments) {
            project.directory = directory;
            project.folder = folder;
        }

        RunMonteCarlo mc = project.monteCarlo;

        if (mc == null) {
            return true; // error
        }

        printRate = 0.004;

        FRAP_on = true;
        FRAP_kPhoto = 2.0;

        double baselinetime = 500;

        //FRAP_onset = baselinetime + 70; // first pulse is 70 ms before bleaching pulse
        //simTime = 10000 + FRAP_onset;

        //FRAP_onset = 15;
        //simTimeStr = 1000 + FRAP_onset;

        FRAP_onset = 2;
        simTime = 60 + FRAP_onset;
        FRAP_small_on = false;

        //double cubeWidth = 2.5;
        //double cubeWidth = 2.0;
        double cubeWidth = 1.1;

        //ellipsoid = true;

        //mitoNonSpaceVoxels = true; // mito are non-space voxels
        //mitoPreset = true;

        //mc.mitoCoordinatesOn = true; // mito are coordinate objects like vesicles
        //mc.mito = getMitochondriaPreset();

        //mc.setMitoDensity = 5;//28;
        //mc.mitoShape = "ellipsoid";

        dx = 2 * vesicleRadius;

        d1_name = "Vesicles";

        //d1_D = 0.090e-3;
        d1_D = 0.060e-3;

        igm.xdim = cubeWidth;
        igm.ydim = cubeWidth;
        igm.zdim = cubeWidth * 1.5;

        if (initProject(dx)) {
            return true;
        }
        if (igm.initGeometryMC(-1)) {
            return true;
        }
        if (initExcitationPSF()) {
            return true;
        }

        dv = initVesicles(d1_D);

        //dv.colorReady.setColor("Heat");
        //dv.colorReady.setColor("BTC");
        //dv.colorReady.setColor("[r=153,g=204,b=0]");
        dv.colorReady.setColor("[r=102,g=255,b=0]");
        dv.colorReady.setGradientColor("[r=255,g=0,b=0]");
        //dv.colorReady.setColor("Magenta");
        //dv.colorReady.invert = true;
        //dv.colorImmobile.setColor("[r=102,g=153,b=0]");
        //dv.color.setColor("[r=0,g=0,b=0]");
        //dv.color.setColor("[r=153,g=0,b=0]");
        //dv.color.setGradientColor("[r=255,g=255,b=255]");
        //dv.color.invert = true;

        dv.setVolumeFraction = 0.26; // Zoltan
        dv.setImmobilePercent = 0;//0.275

        if (initDetectors()) {
            return true;
        }

        mc.minVesicleStep = 0.002; // can use this for all D

        mc.frapOn = true;
        mc.saveFluorescence = true;
        mc.PBC = true;

        mc.setPSFSelect(iPSF, dPSF);

        mc.initAll(); // do not run when running batches

        return false;

    }

    public boolean initLaciLongPulse() {

        directory = "/Jason/FrapVGlut1/D3D/Laci/";
        folder = "Cube8x8x10um";

        double dx = 0.05;

        if (initProject(dx)) {
            return true;
        }
        if (igm.initGeometryMC(-1)) {
            return true;
        }
        if (initExcitationPSF()) {
            return true;
        }
        if (initDiffusantsFD()) {
            return true;
        }
        if (initDetectors()) {
            return true;
        }
        if (initSources()) {
            return true;
        }

        return false;

    }

    public boolean initFiniteDifferenceFrap() {

        String cstr, dstr, ipstr;
        double immobileVesicleFraction;

        double dx = 0.05;

        printRate = 0.001;

        double dim = 3;

        igm.xdim = dim;
        igm.ydim = dim;
        igm.zdim = dim;

        int figSelect = 2; // 2 or 5

        if (figSelect == 2) {
            directory = "/Jason/FrapVGlut1/Simulations/Sims25_FRAP_FiniteDifference/Figure2C/";
            stabilityFD = 0.00185;
            //immobileVesicleFraction = 0;
            immobileVesicleFraction = 0.30; // from fits
            d1_D = 0.027e-3;
        } else if (figSelect == 5) {
            directory = "/Jason/FrapVGlut1/Simulations/Sims25_FRAP_FiniteDifference/Figure5B/";
            stabilityFD = 0.0024;
            immobileVesicleFraction = 0.25;
            //d1_D = 0.090e-3;
            d1_D = 0.060e-3 * 0.295503;
        } else {
            return false;
        }

        //ellipsoid = false;
        igm.eighth_geometry = true;

        d1_name = "Mobile";
        d1_C = 1.0 - immobileVesicleFraction;
        //d1_C = 1.0;

        d2_name = "Immobile";
        d2_C = 1.0 - d1_C; // Laci
        d2_D = 0;
        //d2_C = -1;
        //d2_D = -1;

        FRAP_on = true;
        FRAP_onset = 100;
        //simTime = 100000 + FRAP_onset;
        simTime = 10000 + FRAP_onset;

        //FRAP_small_on = false;
        //FRAP_bleach_on = false;

        iPSFselect = 2; // (0) Gauss (1) Wilson (2) Torok (3) Impulse

        cstr = "c" + String.format("%02d", (int) (dim * 10));
        dstr = InitFrapJasonMisc.diffusionString(d1_D * 1e6);
        ipstr = InitFrapJasonMisc.IPstring(immobileVesicleFraction * 100);

        if (iPSFselect == 0) { // Axelrod Gauss

            excite_fwhm_z = -1;
            excite_fwhm_xy = 0.2708;

            detect_fwhm_z = -1;
            detect_fwhm_xy = excite_fwhm_xy;

            folder = "FD_Gauss_";

            FRAP_kPhoto = 1.0;

        } else if (iPSFselect == 2) {

            detect_fwhm_xy = 0.2551; // measured, bead deconvolved
            detect_fwhm_z = 0.9157; // measured, bead deconvolved

            //folder = "FD_Torok_";
            folder = "FD_FRAP_";

            FRAP_kPhoto = 2.0;

        }

        folder += dstr + "_" + ipstr;

        if (FRAP_small_on) {
            folder += "_probes";
        } else {
            folder += "_noprobes";
        }

        folder += "_" + Integer.toString((int) (simTime/1000)) + "s";
        
        project.name = "Finite Difference FRAP";

        if (!Master.foundMainStartUpArguments) {
            project.directory = directory;
            project.folder = folder;
        }

        if (initProject(dx)) {
            return true;
        }
        if (igm.initGeometryMC(-1)) {
            return true;
        }
        if (initExcitationPSF()) {
            return true;
        }
        if (initDiffusantsFD()) {
            return true;
        }
        if (initDetectors()) {
            return true;
        }

        return false;

    }

    public boolean initMFsmall() {

        directory = "/Jason/FrapVGlut1/D3D/";
        folder = "MFsmall";

        double dim = 2;
        double dx = 0.05;

        igm.xdim = dim;
        igm.ydim = dim;
        igm.zdim = dim;

        simTime = 3000.0;

        stabilityFD = 0.00001;

        igm.ellipsoid = false;
        igm.eighth_geometry = true;

        d1_name = "Mobile";
        d1_C = 0.7;
        d1_D = 0.00001;

        d2_name = "Immobile";
        d2_C = 0.3;
        d2_D = 0;

        iPSFselect = 0; // Gaussian

        FRAP_on = true;
        FRAP_kPhoto = 3.85;
        FRAP_onset = 0.1; // ms
        FRAP_bleach_length = 0.5; // ms

        //writePSF = true;
        //readPSF = true;

        if (initProject(dx)) {
            return true;
        }
        if (igm.initGeometryMC(-1)) {
            return true;
        }
        if (initExcitationPSF()) {
            return true;
        }
        if (initDiffusantsFD()) {
            return true;
        }
        if (initDetectors()) {
            return true;
        }

        return false;

    }

    public boolean initProject(double dx) {

        if (!Master.foundMainStartUpArguments) {
           project.simTime = simTime;
        }
        
        project.dx = dx;
        project.set("saveRate", saveRate);
        project.set("printRate", printRate);

        project.stability = stabilityFD;

        return false;

    }

    public boolean initExcitationPSF() {

        double sdxy, sdz;
        
        if (readPSF) {

            switch (iPSFselect) {
                case 0:
                    gpsf = (PSFgauss) Master.readObject(directory + psf_file);
                    break;
                case 1:
                    wpsf = (PSFwilson) Master.readObject(directory + psf_file);
                    break;
                case 2:
                    tpsf = (PSFtorok) Master.readObject(directory + psf_file);
                    break;
                case 3:
                    epsf = (PSF) Master.readObject(directory + psf_file);
                    break;
                default:
                    return true;
            }

        } else {

            switch (iPSFselect) {

                case 0:

                    sdxy = PSFgauss.computeSTDV(excite_fwhm_xy);
                    sdz = PSFgauss.computeSTDV(excite_fwhm_z);
                    
                    gpsf = new PSFgauss(project, null, sdxy, sdxy, sdz);

                    if (igm.eighth_geometry) {
                        gpsf.arbitraryVoxelCenter = true;
                        gpsf.xySymmetric = false;
                    } else {
                        gpsf.xySymmetric = true;
                    }

                    gpsf.set("xVoxelCenter", xvcenter);
                    gpsf.set("yVoxelCenter", yvcenter);
                    gpsf.set("zVoxelCenter", zvcenter);
                    gpsf.name = "Illumination PSF";
                    
                    iPSF = geometry.addPSF(gpsf);

                    break;

                case 1:

                    wpsf = new PSFwilson(project, null);

                    wpsf.numericalAperture = NA;
                    wpsf.waveLength = lambda_illum;
                    wpsf.refractiveIndex = refrac_water;

                    if (igm.eighth_geometry) {
                        wpsf.arbitraryVoxelCenter = true;
                        wpsf.xySymmetric = false;
                    } else {
                        wpsf.xySymmetric = true;
                    }

                    wpsf.set("xVoxelCenter", xvcenter);
                    wpsf.set("yVoxelCenter", yvcenter);
                    wpsf.set("zVoxelCenter", zvcenter);
                    wpsf.name = "Illumination PSF";
                    wpsf.rscale = rscale;

                    iPSF = geometry.addPSF(wpsf);

                    break;

                case 2:

                    tpsf = new PSFtorok(project, null);
                    tpsf.setFitZoltan3();
                    //tpsf.setFitJason();
                    //tpsf.setFitDavid();

                    if (igm.eighth_geometry) {
                        tpsf.arbitraryVoxelCenter = true;
                        tpsf.xySymmetric = false;
                    } else {
                        tpsf.xySymmetric = true;
                    }

                    tpsf.set("xVoxelCenter", xvcenter);
                    tpsf.set("yVoxelCenter", yvcenter);
                    tpsf.set("zVoxelCenter", zvcenter);
                    tpsf.name = "Illumination PSF";
                    tpsf.rscale = rscale;

                    iPSF = geometry.addPSF(tpsf);

                    break;

                case 3:
                    
                    epsf = new PSF(project, null);

                    if (igm.eighth_geometry) {
                        epsf.arbitraryVoxelCenter = true;
                        epsf.xySymmetric = false;
                    } else {
                        epsf.xySymmetric = true;
                    }

                    epsf.set("xVoxelCenter", xvcenter);
                    epsf.set("yVoxelCenter", yvcenter);
                    epsf.set("zVoxelCenter", zvcenter);
                    epsf.name = "Illumination PSF";

                    iPSF = geometry.addPSF(epsf);

                    break;

                default:
                    return true;
            }

            if (writePSF) {
                switch (iPSFselect) {
                    case 0:
                        gpsf.compute();
                        gpsf.sum(geometry);
                        Master.writeObject(directory + psf_file, gpsf);
                        break;
                    case 1:
                        wpsf.compute();
                        wpsf.sum(geometry);
                        Master.writeObject(directory + psf_file, wpsf);
                        break;
                    case 2:
                        tpsf.compute();
                        tpsf.sum(geometry);
                        Master.writeObject(directory + psf_file, tpsf);
                        break;
                    case 3:
                        epsf.compute();
                        epsf.sum(geometry);
                        Master.writeObject(directory + psf_file, epsf);
                        break;
                    default:
                        return true;
                }
            }

        }

        //c.set(0, tpsf.d.y0, 0, tpsf.d.xVox-1, tpsf.d.y0, tpsf.d.zVox-1);
        //Master.export3DArray(directory + "PSFtorok_xz_densityE.txt", c, tpsf.p);

        //wpsf.checkExists();

        //c.set(0, 0, 0, geom.xVoxels()-1, geom.yVoxels()-1, 0);
        //Master.export3DArrayLaci(directory + "PSF_Wilson.txt", c, wpsf.p);

        //c.set(0, 61, 81, geom.xVoxels()-1, 61, 81);
        //Master.export3DArrayLaci(directory + "PSF_Gauss.txt", c, dpsf.p);

        return false;

    }

    public boolean initDiffusantsFD() {

        int dnum = -1;

        //int dnum = Master.addDiffusant("DyeBleached", 0, dye_D); // bleached dye

        initPulseTimer();

        switch (iPSFselect) {
            default:
                frap = new DiffusantPhoto(project, d1_name, d1_C, d1_D, null, dnum, timer, gpsf, FRAP_kPhoto);
                //frap2 = new DiffusantPhoto("Dye2", C2, D2, dnum, timer, gpsf);
                break;
            case 1:
                frap = new DiffusantPhoto(project, d1_name, d1_C, d1_D, null, dnum, timer, wpsf, FRAP_kPhoto);
                //frap2 = new DiffusantPhoto("Dye2", C2, D2, dnum, timer, wpsf);
                break;
            case 2:
                frap = new DiffusantPhoto(project, d1_name, d1_C, d1_D, null, dnum, timer, tpsf, FRAP_kPhoto);
                //frap2 = new DiffusantPhoto("Dye2", C2, D2, dnum, timer, tpsf);
                break;
            case 3:
                frap = new DiffusantPhoto(project, d1_name, d1_C, d1_D, null, dnum, timer, epsf, FRAP_kPhoto);
                //frap2 = new DiffusantPhoto("Dye2", C2, D2, dnum, timer, epsf);
                break;
        }

        //frap.save.save2BinaryFile = true;

        d1_num = project.addDiffusant(frap);

        if ((d2_C >= 0) && (d2_D >= 0)) {

            switch (iPSFselect) {
                default:
                    frap2 = new DiffusantPhoto(project, d2_name, d2_C, d2_D, null, dnum, timer, gpsf, FRAP_kPhoto);
                    break;
                case 1:
                    frap2 = new DiffusantPhoto(project, d2_name, d2_C, d2_D, null, dnum, timer, wpsf, FRAP_kPhoto);
                    break;
                case 2:
                    frap2 = new DiffusantPhoto(project, d2_name, d2_C, d2_D, null, dnum, timer, tpsf, FRAP_kPhoto);
                    break;
                case 3:
                    frap2 = new DiffusantPhoto(project, d2_name, d2_C, d2_D, null, dnum, timer, epsf, FRAP_kPhoto);
                    break;
            }

            frap2.save.save2TextFile = false;

            d2_num = project.addDiffusant(frap2);

        }

        return false;

    }

    public void initPulseTimer() {

        if (!FRAP_on) {
            return;
        }

        if (FRAP_bleach_on) {
            timer = new PulseTimer(project, FRAP_onset, FRAP_bleach_length, 1.0);
        } else {
            timer = new PulseTimer(project, FRAP_onset, FRAP_bleach_length, 0);
        }

        if (FRAP_small_on) {
            timer.add(FRAP_onset - 70, FRAP_small_length, FRAP_small_amp_ratio);
            timer.add(FRAP_onset - 40, FRAP_small_length, FRAP_small_amp_ratio);
            timer.add(FRAP_onset - 10, FRAP_small_length, FRAP_small_amp_ratio);
            timer.add(FRAP_onset + 20, FRAP_small_length, FRAP_small_amp_ratio);
            timer.add(FRAP_onset + 50, FRAP_small_length, FRAP_small_amp_ratio);
            timer.add(FRAP_onset + 80, FRAP_small_length, FRAP_small_amp_ratio);
            timer.add(FRAP_onset + 110, FRAP_small_length, FRAP_small_amp_ratio);
            timer.add(FRAP_onset + 150, FRAP_small_length, FRAP_small_amp_ratio);
            timer.add(FRAP_onset + 200, FRAP_small_length, FRAP_small_amp_ratio);
            timer.add(FRAP_onset + 260, FRAP_small_length, FRAP_small_amp_ratio);
            timer.add(FRAP_onset + 310, FRAP_small_length, FRAP_small_amp_ratio);
            timer.add(FRAP_onset + 390, FRAP_small_length, FRAP_small_amp_ratio);
            timer.add(FRAP_onset + 490, FRAP_small_length, FRAP_small_amp_ratio);
            timer.add(FRAP_onset + 810, FRAP_small_length, FRAP_small_amp_ratio);
            timer.add(FRAP_onset + 1110, FRAP_small_length, FRAP_small_amp_ratio);
            timer.add(FRAP_onset + 1710, FRAP_small_length, FRAP_small_amp_ratio);
            timer.add(FRAP_onset + 2710, FRAP_small_length, FRAP_small_amp_ratio);
            timer.add(FRAP_onset + 3910, FRAP_small_length, FRAP_small_amp_ratio);
            timer.add(FRAP_onset + 4880, FRAP_small_length, FRAP_small_amp_ratio);
            timer.add(FRAP_onset + 5910, FRAP_small_length, FRAP_small_amp_ratio);
            timer.add(FRAP_onset + 6910, FRAP_small_length, FRAP_small_amp_ratio);
            timer.add(FRAP_onset + 7910, FRAP_small_length, FRAP_small_amp_ratio);
            timer.add(FRAP_onset + 8910, FRAP_small_length, FRAP_small_amp_ratio);
            timer.add(FRAP_onset + 9880, FRAP_small_length, FRAP_small_amp_ratio);
        }

    }

    public DiffusantVesicles initVesicles(double D) {

        DiffusantVesicles dv;

        initPulseTimer();

        switch (iPSFselect) {
            default:
                dv = new DiffusantVesicles(project, d1_name, 0, D, null, vesicleRadius, timer, gpsf);
                break;
            case 1:
                dv = new DiffusantVesicles(project, d1_name, 0, D, null, vesicleRadius, timer, wpsf);
                break;
            case 2:
                dv = new DiffusantVesicles(project, d1_name, 0, D, null, vesicleRadius, timer, tpsf);
                break;
            case 3:
                dv = new DiffusantVesicles(project, d1_name, 0, D, null, vesicleRadius, timer, epsf);
                break;
        }

        dv.kPhoto = FRAP_kPhoto;
        dv.save.save2TextFile = false;

        d1_num = project.addDiffusant(dv);

        return dv;

    }

    public boolean initDetectors() {

        double sdxy = PSFgauss.computeSTDV(detect_fwhm_xy);
        double sdz = PSFgauss.computeSTDV(detect_fwhm_z);

        PSFgauss dpsf = new PSFgauss(project, null, sdxy, sdxy, sdz);

        //PSFconfocal dpsf = new PSFconfocal(project, null, NA);

        if (igm.eighth_geometry) {
            dpsf.arbitraryVoxelCenter = true;
            dpsf.xySymmetric = false;
        } else {
            dpsf.xySymmetric = true;
        }
        
        dpsf.set("xVoxelCenter", xvcenter);
        dpsf.set("yVoxelCenter", yvcenter);
        dpsf.set("zVoxelCenter", zvcenter);
        
        //dpsf.compute();

        if (d1_num >= 0) {
            DetectorPSF dw1 = Master.addDetectorPSF(d1_num, dpsf);
            dw1.save.saveWhileComputing = saveWhileComputing;
            dw1.save.save2BinaryFile = true;
            dw1.save.save2TextFile = false;
            dw1.save.ydim = "Fluorescence";
            dw1.save.autoDimensions = false;
        }

        if (d2_num >= 0) {
            DetectorPSF dw2 = Master.addDetectorPSF(d2_num, dpsf);
            dw2.save.saveWhileComputing = saveWhileComputing;
            dw2.save.save2BinaryFile = true;
            dw2.save.save2TextFile = false;
            dw2.save.ydim = "Fluorescence";
            dw2.save.autoDimensions = false;
        }

        dPSF = geometry.addPSF(dpsf);

        // dw.psf.checkExists();
        //double sum=dw.psf.sum(geom);

        //System.output.println(sum);

        //c.set(20,20,30,20,20,30);
        //Detector d3 = Master.addDetectorAvg(d1_num);
        //d3.set("saveWhileComputing", saveWhileComputing);

        //c.set(17,17,27,17,17,27);
        //Master.addDetectorAvg(d1_num, c);

        //c.set(14,14,24,14,14,24);
        //Master.addDetectorAvg(d1_num, c);

        //c.set(11,11,21,11,11,21);
        //Master.addDetectorAvg(d1_num, c);

        //c.set(5,5,15,5,5,15);
        //Master.addDetectorAvg(d1_num, c);

        //c.set(18,18,28,23,23,33);
        //c.set(1, 1, 1, 4, 4, 4);
        //Master.addDetectorAvg(d1_num, c);
        //Master.addDetectorDFOF(d1_num, dpsf);

        //dw.set("xVoxelCenter", 0);
        //dw.set("yVoxelCenter", 0);
        //dw.set("zVoxelCenter", 0);

        //dw.psf.checkExists();

        //c.set(0, 0, 0, geom.xVoxels()-1, geom.yVoxels()-1, 0);
        //Master.export3DArrayLaci(directory + "PSF_Gauss.txt", c, dw.psf.p);

        /*
        int istep = (int) ( (geom.xVoxels() - 3) / 5.0);


        for (int ii = 0; ii < geom.xVoxels(); ii += istep) {
        c.set( (int) x0 - 1 - ii, (int) y0 - 1, (int) z0, (int) x0 - ii, (int) y0,
        (int) z0);
        Master.addDetectorAvg(d1_num, c);
        }
         */

        //Master.addDetectorAvg(d1_num, c.set(x0, y0, z0));

        //Master.addDetectorSnapshot(d1_num, "xy", (int) z0, t_onset + t_length);
        //Master.addDetectorSnapshot(d1_num, "yz", (int) x0, t_onset + t_length);

        return false;

    }

    public boolean initSources() {

        Source s0;

        s0 = Master.addSource(d1_num, coordinates.setVoxels(0, 0, 0, geometry.xVoxels - 1, geometry.yVoxels - 1, 0), 1.0);
        s0.clamp = true;

        s0 = Master.addSource(d1_num, coordinates.setVoxels(0, 0, 1, geometry.xVoxels - 1, 0, geometry.zVoxels - 1), 1.0);
        s0.clamp = true;

        s0 = Master.addSource(d1_num, coordinates.setVoxels(0, 1, 1, 0, geometry.yVoxels - 1, geometry.zVoxels - 1), 1.0);
        s0.clamp = true;

        //c.set(25,25,37,25,25,37);
        //c.set(1, 1, 1, 4, 4, 4);
        //Master.addSource(d1_num, c, 156.765, 0, 0.003);

        return false;

    }

    public void initFrap_Zoltan_LongPulse() {

        double sdxy, sdz;

        directory = "/Jason/FrapVGlut1/D3D/";
        folder = "Torok";

        double simTime = 21; // msec
        double saveRate = 400; // hz

        double dx = 0.02;
        double xdim = 2;
        double ydim = 2;
        double zdim = 2;

        //zdim = 4;

        //timeUnits -= dx;
        //concUnits -= dx;
        //zdim -= dx;

        boolean eighth = false;

        double dye_C = 1.0; // dye concentration
        double dye_D = 0.2; // diffusion coefficient
        double kPhoto = 1.0; // photolysis k factor
        //double k = 1.0; // photolysis k factor
        //double I = 1.0; // photolysis k factor
        //kPhoto = k*I;

        double t_onset = 0.1; // ms
        double t_length = 10.0;

        int PSFselect = 0; // (0) Gauss (1) Wilson (2) Torok

        String psf_file = "D3D_Torok_testing_psf.d3d";
        boolean writePSF = false;
        boolean readPSF = false;

        double refrac_water = 1.34; // Wilson PSF
        double lambda_illum = 0.488; // um
        double NA = 1.0;

        double excite_fwhm_xy = 0.343; // Gauss PSF
        double excite_fwhm_z = 1.285;

        double detect_fwhm_xy = 0.26; // Gauss PSF
        double detect_fwhm_z = 0.89;

        int xVoxels, yVoxels, zVoxels;
        double xcenter, ycenter, zcenter;

        //simTimeStr = t_onset + t_length + 0.5;

        project.simTime = simTime;
        project.dx = dx;
        project.saveRate = saveRate;

        // GEOMTRY

        xVoxels = (int) (xdim / dx);
        yVoxels = (int) (ydim / dx);
        zVoxels = (int) (zdim / dx);

        xVoxels += 1;
        zVoxels += 1;

        if (eighth) {

            xVoxels = 1 + (int) (xdim / (dx * 2)); // half
            yVoxels = 1 + (int) (ydim / (dx * 2)); // half

            if (PSFselect < 2) {
                zVoxels = 1 + (int) (zdim / (dx * 2)); // half
            }

        }

        geometry.resize(xVoxels, yVoxels, zVoxels);

        xcenter = geometry.xVoxelCenter;
        ycenter = geometry.yVoxelCenter;
        zcenter = geometry.zVoxelCenter;

        if (eighth) {

            xcenter = xVoxels - 0.5;
            ycenter = yVoxels - 0.5;

            if (PSFselect < 2) {
                zcenter = zVoxels - 0.5;
            }

        }

        // PSFs

        PSFgauss gpsf = null;
        PSFwilson wpsf = null;
        PSFtorok tpsf = null;
        PSFgauss dpsf = null; // detection PSF

        if (readPSF) {
            switch (PSFselect) {
                default:
                    gpsf = (PSFgauss) Master.readObject(directory + psf_file);
                    break;
                case 1:
                    wpsf = (PSFwilson) Master.readObject(directory + psf_file);
                    break;
                case 2:
                    tpsf = (PSFtorok) Master.readObject(directory + psf_file);
                    break;
            }
        } else {

            switch (PSFselect) {
                default:
                    sdxy = PSFgauss.computeSTDV(excite_fwhm_xy);
                    sdz = PSFgauss.computeSTDV(excite_fwhm_z);
                    gpsf = new PSFgauss(project, null, sdxy, sdxy, sdz);
                    if (eighth) {
                        gpsf.setVoxelCenter(xcenter, ycenter, zcenter);
                    }
                    break;
                case 1:
                    wpsf = new PSFwilson(project, null);
                    wpsf.set("numericalAperture", NA);
                    wpsf.set("waveLength", lambda_illum);
                    wpsf.set("refractiveIndex", refrac_water);
                    if (eighth) {
                        wpsf.setVoxelCenter(xcenter, ycenter, zcenter);
                    }
                    break;
                case 2:
                    tpsf = new PSFtorok(project, null);

                    //tpsf.setFitZoltan2();
                    //tpsf.setFitJason();
                    tpsf.setFitDavid();
                    if (eighth) {
                        tpsf.setVoxelCenter(xcenter, ycenter, zcenter);
                    }
                    break;
            }

            if (writePSF) {
                switch (PSFselect) {
                    default:
                        gpsf.compute();
                        gpsf.sum(geometry);
                        Master.writeObject(directory + psf_file, gpsf);
                        break;
                    case 1:
                        wpsf.compute();
                        wpsf.sum(geometry);
                        Master.writeObject(directory + psf_file, wpsf);
                        break;
                    case 2:
                        tpsf.compute();
                        tpsf.sum(geometry);
                        Master.writeObject(directory + psf_file, tpsf);
                        break;
                }
            }

        }

        sdxy = PSFgauss.computeSTDV(detect_fwhm_xy);
        sdz = PSFgauss.computeSTDV(detect_fwhm_z);

        dpsf = new PSFgauss(project, null, sdxy, sdxy, sdz); // detection

        if (eighth) {
            dpsf.setVoxelCenter(xcenter, ycenter, zcenter);
        }

        //c.set(0, tpsf.d.y0, 0, tpsf.d.xVox-1, tpsf.d.y0, tpsf.d.zVox-1);
        //Master.export3DArray(directory + "PSFtorok_xz_densityE.txt", c, tpsf.p);

        //c.set(0, geom.y0(), 0, geom.xVoxels()-1, geom.y0(), geom.zVoxels()-1);
        //Master.export3DArrayLaci(directory + "PSF_Wilson_4x8.txt", c, epsf.p);

        // EXCITATION

        PulseTimer timer = new PulseTimer(project, t_onset, t_length);

        timer.add((t_onset + t_length + 0.5), 0.1);
        timer.add((t_onset + t_length + 2.0), 0.1);
        timer.add((t_onset + t_length + 5.0), 0.1);
        timer.add((t_onset + t_length + 8.5), 0.1);
        timer.add((t_onset + t_length + 18.5), 0.1);
        timer.add((t_onset + t_length + 28.5), 0.1);
        //PulseTimer timer = new PulseTimer("0.1,0.1,1;0.3,0.2,2;");
        //PulseTimer timer = new PulseTimer("0.1,0.1;");

        //int dnum = Master.addDiffusant("DyeBleached", 0, dye_D); // bleached dye
        int dnum = -1;

        DiffusantPhoto frap;

        switch (PSFselect) {
            default:
                frap = new DiffusantPhoto(project, "Dye", dye_C, dye_D, null, dnum, timer, gpsf, kPhoto);
                break;
            case 1:
                frap = new DiffusantPhoto(project, "Dye", dye_C, dye_D, null, dnum, timer, wpsf, kPhoto);
                break;
            case 2:
                frap = new DiffusantPhoto(project, "Dye", dye_C, dye_D, null, dnum, timer, tpsf, kPhoto);
                break;
        }

        int dye_num = project.addDiffusant(frap);

        frap.save.save2TextFile = true;

        // DETECTION

        Master.addDetectorPSF(dye_num, dpsf);

        //int istep = (int) ((shape.xVoxels - 3) / 5.0);

        /*
        for (int ii = 0; ii < geom.xVoxels(); ii += istep) {
        c.set( (int) x0 - 1 - ii, (int) y0 - 1, (int) z0, (int) x0 - ii, (int) y0,
        (int) z0);
        Master.addDetectorAvg(dye_num, c);
        }
         */

        //Master.addDetectorAvg(dye_num, vCoordinates.setVoxels(xcenter, ycenter, zcenter));

        //Master.addDetectorSnapshot(dye_num, "xy", (int) z0, t_onset + t_length);
        //Master.addDetectorSnapshot(dye_num, "yz", (int) x0, t_onset + t_length);

        // SOURCES

        //Master.addSource(dye_num, c.set(0,0,0,geom.xVoxels()-1,geom.yVoxels()-1,0), 1.0);
        //Master.addSource(dye_num, c.set(0,0,0,geom.xVoxels()-1,0,geom.zVoxels()-1), 1.0);
        //Master.addSource(dye_num, c.set(0,0,0,1,geom.yVoxels()-1,geom.zVoxels()-1), 1.0);

        // BATCHES

        //Master.addBatchList("Dye.D", "0.05,0.2,0.3,0.4,0.15,");
        Master.addBatchList("Dye.D", "0.2,0.5,");
        //Master.addBatchList("Dye.kPhoto", "1.0,0.56,");
        //Master.addBatchList("Dye.C0", "0.45,");

        /*
        Master.addBatch("Project.folder", "Torok_D50", false);
        Master.addBatch("Dye.D", 0.5, false);
        Master.addBatch("Dye.pulses", "0.1,0.02;");
        Master.addBatch("Dye.pulses", "0.1,0.25;");

        Master.addBatch("Project.folder", "Torok_D20", false);
        Master.addBatch("Dye.D", 0.2, false);
        Master.addBatch("Dye.pulses", "0.1,0.02;");
        Master.addBatch("Dye.pulses", "0.1,0.25;");
         */

    }

    public void initFrap_Zoltan() {

        double simTime = 2000.0; // msec

        double dx = 0.02;

        double xdim = 1.5;
        double ydim = 1.5;
        double zdim = 5;

        boolean eighth = false;

        int xVoxels, yVoxels, zVoxels;
        double xcenter, ycenter, zcenter;

        double refrac_water = 1.3459;
        double lambda_illum = 0.488; // um
        double NA = 1.0;

        double detect_fwhm_xy = 0.26;
        double detect_fwhm_z = 0.89;

        double vesicle_conc = 1.0;
        double vesicle_D = 1e-4;
        double kPhoto = 1.0; // photolysis k factor

        double t_onset = 0.1; // ms
        double t_length = 0.3;

        //t_onset = -1;

        simTime = t_onset + t_length;

        project.set("simTime", simTime);
        project.set("dx", dx);

        // GEOMTRY

        if (eighth) {
            xVoxels = 1 + (int) (xdim / (dx * 2)); // half
            yVoxels = 1 + (int) (ydim / (dx * 2)); // half
            zVoxels = 1 + (int) (zdim / (dx * 2)); // half
        } else {
            xVoxels = (int) (xdim / dx);
            yVoxels = (int) (ydim / dx);
            zVoxels = (int) (zdim / dx);
        }

        geometry.resize(xVoxels, yVoxels, zVoxels);

        xcenter = geometry.xVoxelCenter;
        ycenter = geometry.yVoxelCenter;
        zcenter = geometry.zVoxelCenter;

        if (eighth) {
            xcenter = xVoxels - 0.5;
            ycenter = yVoxels - 0.5;
            zcenter = zVoxels - 0.5;
        }

        // PSFs

        PSFwilson wpsf = new PSFwilson(project, null);

        wpsf.set("xVoxelCenter", xcenter);
        wpsf.set("yVoxelCenter", ycenter);
        wpsf.set("zVoxelCenter", zcenter);

        wpsf.set("numericalAperture", NA);
        wpsf.set("waveLength", lambda_illum);
        wpsf.set("refractiveIndex", refrac_water);

        double sdxy = PSFgauss.computeSTDV(detect_fwhm_xy);
        double sdz = PSFgauss.computeSTDV(detect_fwhm_z);

        PSFgauss dpsf = new PSFgauss(project, null, sdxy, sdxy, sdz); // detection

        dpsf.set("xVoxelCenter", xcenter);
        dpsf.set("yVoxelCenter", ycenter);
        dpsf.set("zVoxelCenter", zcenter);

        //wpsf.compute();
        //Master.writeObject(directory + "D3D_frap_cube_5um_psf.d3d", wpsf);

        //PSFwilson wpsf = (PSFwilson) Master.readObject(directory+"D3D_bouton_psf.d3d");
        //wpsf = (PSFwilson) Master.readObject(directory+"D3D_frap_cube_5um_psf.d3d");
        //PSFwilson wpsf = (PSFwilson) Master.readObject(directory+"D3D_frap_cube_10um_psf.d3d");


        //c.set(0, geom.y0, geom.z0, geom.xVoxels-1, geom.y0, geom.z0);
        //Master.export3DArray(directory + "PSFgauss.txt", c, wpsf.p);

        // EXCITATION

        PulseTimer timer = new PulseTimer(project, t_onset, t_length);

        //int dnum = Master.addDiffusant("DyeBleached", 0, dye_D); // bleached dye
        int dnum = -1;

        DiffusantPhoto frap = new DiffusantPhoto(project, "Dye", vesicle_conc,
                vesicle_D, null, dnum, timer, dpsf, kPhoto);

        int dye_num = project.addDiffusant(frap);

        //frap.set("saveConc", 1);

        // DETECTION

        Master.addDetectorPSF(dye_num, dpsf);

        int istep = (int) ((geometry.xVoxels - 3) / 5.0);

        for (int ii = 0; ii < geometry.xVoxels; ii += istep) {
            coordinates.setVoxels((int) xcenter - 1 - ii, (int) ycenter - 1, (int) zcenter,
                    (int) xcenter - ii, (int) ycenter, (int) zcenter);
            //Master.addDetectorAvg(dye_num, c);
        }

        // BATCHES

        //Master.addBatchList("Dye.D", "0.1,0.2,");
        //Master.addBatchList("Project.Stability", "0.5,0.4,0.3,0.2");

    }

    @Override
    public boolean addUser(ParamVector pv) {

        super.addUser(pv);

        if (igm != null){
            igm.addUser(pv);
        }

        return true;

    }

    @Override
    public boolean createVector(boolean close) {

        if (!super.createVector(false)) {
            return false;
        }

        if (igm != null){
            addBlankParam();
            igm.createVector(true);
            addVector(igm.getVector());
            igm.addUser(this);
        }

        if (close) {
            closeVector();
        }

        return true;

    }

    @Override
    public void updateVector(ParamObject[] v) {

        super.updateVector(v);

        if (igm != null){
            igm.updateVector(v);
        }

    }

}
