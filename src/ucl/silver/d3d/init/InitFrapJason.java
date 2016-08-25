package ucl.silver.d3d.init;

import ucl.silver.d3d.core.*;

/**
 * <p>Title: D3D</p>
 * <p>Description: 3D Diffusion-Reaction Simulator</p>
 * <p>Copyright: Copyright (c) 2016</p>
 * <p>Company: The Silver Lab at University College London</p>
 * @author Jason S. Rothman
 * @version 1.6
 */

public class InitFrapJason extends InitProject {

    public String directory = "/Jason/D3D/Simulations/";
    public String folder = "Testing";

    public double simTime = 8; // msec
    public double saveRate = 20; // kHz
    public double printRate = 0.002; // kHz

    public boolean saveWhileComputing = false;
    
    public double stabilityFD = 0.00185;

    public double vesicleRadius = 0.0494382 / 2.0; // um
    
    public double mitoVolumePercent = 0.28; // EM data
    public double totalVesicleVolumeFraction = 0.17; // EM data
    public double immobileVesiclePercent = 0.25; // best-match MC FRAP curve
    
    public String d1_name = "Vesicles";
    public double d1_C = 1;
    public double d1_Dshort = 0.060e-3; // best-match FRAP simulation
    public double Dcyto = d1_Dshort / DiffusantVesicle.Dratio_short(totalVesicleVolumeFraction * (1 - immobileVesiclePercent), totalVesicleVolumeFraction * immobileVesiclePercent); // 1.2696940661239775E-4

    public boolean FRAP_on = false;
    public double FRAP_kPhoto = 2; // photolysis k factor
    public double FRAP_onset = 100; // ms
    public double FRAP_bleach_length = 0.5; // ms
    public boolean FRAP_bleach_on = true; // bleach pulse
    public double FRAP_small_length = 2.0; // ms
    public double FRAP_small_amp_ratio = 0.0014; // ratio measured laser power after objective
    public boolean FRAP_small_on = true;

    public int iPSFselect = 2; // (0) Gauss (1) Wilson (2) Torok (3) Impulse
    private int iPSF, cPSF;

    public double refrac_water = 1.338; // at 37C
    public double lambda_illum = 0.488; // um
    public double NA = 0.90; // fit to Zoltan's iPSF

    public double iPSF_fwhm_xy = 0.343; // Gauss iPSF (for testing only)
    public double iPSF_fwhm_z = 1.285; // Gauss iPSF (for testing only)
    
    public double cPSF_fwhm_xy = 0.2551; // measured, bead deconvolved
    public double cPSF_fwhm_z = 0.9157; // measured, bead deconvolved

    private int d1_num = -1;
    private int d2_num = -1;

    private final InitGeometryMito igm;

    private PSFgauss gpsf = null;
    private PSFwilson wpsf = null;
    private PSFtorok tpsf = null;
    private PSF epsf = null;

    private DiffusantPhoto frap = null;
    private DiffusantPhoto frap2 = null;

    private PulseTimer timer = null;

    public InitFrapJason(Project p) {
        super(p);
        igm = new InitGeometryMito(p);
        String[] flist = {"init_FD_FRAP", "init_FD_Axelrod", "init_MC_FRAP", "init_MC_FRAP_Prompt", "init_MC_FRAP_Movie", "init_MC_MSD", "init_MC_MSD_Prompt", "init_MC_Cichocki"};
        initFuncList = flist;
        createVector(true);
    }
    
    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("simTime")) {
            return project.timeUnits;
        }
        if (name.equalsIgnoreCase("saveRate")) {
            return project.freqUnits;
        }
        if (name.equalsIgnoreCase("printRate")) {
            return project.freqUnits;
        }
        if (name.equalsIgnoreCase("vesicleRadius")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("d1_Dshort")) {
            return project.diffusionUnits;
        }
        if (name.equalsIgnoreCase("Dcyto")) {
            return project.diffusionUnits;
        }
        if (name.equalsIgnoreCase("FRAP_kPhoto")) {
            return "1/" + project.timeUnits;
        }                             
        if (name.equalsIgnoreCase("FRAP_onset")) {
            return project.timeUnits;
        }
        if (name.equalsIgnoreCase("FRAP_bleach_length")) {
            return project.timeUnits;
        }
        if (name.equalsIgnoreCase("FRAP_small_length")) {
            return project.timeUnits;
        }
        if (name.equalsIgnoreCase("lambda_illum")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("iPSF_fwhm_xy")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("iPSF_fwhm_z")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("cPSF_fwhm_xy")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("cPSF_fwhm_z")) {
            return project.spaceUnits;
        }   
        return super.units(name);
    }

    @Override
    public boolean initFunction(String initSelect) {

        int i = initFunctionNum(initSelect);

        if (!Master.foundMainStartUpArguments) {
            //i = 0;
        }

        switch (i) {
            
            case 0:
                return init_FD_FRAP();

            case 1:
                return init_FD_Axelrod();
                
            case 2:
                return init_MC_FRAP();

            case 3:
                return init_MC_FRAP_Prompt();

            case 4:
                return init_MC_FRAP_Movie();

            case 5:
                return init_MC_MSD();
                
            case 6:
                return init_MC_MSD_Prompt();

            case 7:
                return init_MC_Cichocki();

            default:
                Master.log("InitFrapJason.initFunction error: failed to find init function " + initSelect);
                return true; // error

        }

    }
    
    private boolean init_MC_FRAP_Prompt() {
        
        double Dshort, IP;
        String DvalueStr, IPstr;
        
        String Dvalues[] = new String[16];
        
        Dvalues[0] = "0.0300";
        Dvalues[1] = "0.0350";
        Dvalues[2] = "0.0400";
        Dvalues[3] = "0.0450";
        Dvalues[4] = "0.0500";
        Dvalues[5] = "0.0525";
        Dvalues[6] = "0.0550";
        Dvalues[7] = "0.0575";
        Dvalues[8] = "0.0600"; // best-match MC simulation
        Dvalues[9] = "0.0625";
        Dvalues[10] = "0.0650";
        Dvalues[11] = "0.0675";
        Dvalues[12] = "0.0700";
        Dvalues[13] = "0.0725";
        Dvalues[14] = "0.0750";
        Dvalues[15] = "0.0800";
        
        DvalueStr = Master.promptForInput(Dvalues, "choose Dshort (um^2/s)", "Init Monte Carlo FRAP", "0.0600");
        
        if ((DvalueStr == null) || (DvalueStr.length() == 0)) {
            return true; // error
        }
        
        Dshort = Double.parseDouble(DvalueStr);
        
        String IPvalues[] = new String[7];
        
        IPvalues[0] = "17.50";
        IPvalues[1] = "20.00";
        IPvalues[2] = "22.50";
        IPvalues[3] = "23.75";
        IPvalues[4] = "25.00"; // best-match MC simulation
        IPvalues[5] = "26.25";
        IPvalues[6] = "27.50";
        
        IPstr = Master.promptForInput(IPvalues, "choose immobile vesicle %", "Init Monte Carlo FRAP", "25.00");
        
        if ((IPstr == null) || (IPstr.length() == 0)) {
            return true; // error
        }
        
        IP = Double.parseDouble(IPstr);
        
        Master.log("Dshort = " + Dshort + "um^2/s, immobile vesicle percent = " + IP + "%");
        
        Dshort /= 1000.0;
        IP /= 100.0;
        
        d1_Dshort = Dshort;
        immobileVesiclePercent = IP;
        
        return init_MC_FRAP();
        
    }

    private boolean init_MC_FRAP() {
        // Figure 3D
        // Figure 4
        // Figure 2-figure supplement 1B (+drift)

        double dx;
        DiffusantVesicles dv;

        project.newMonteCarlo();
        project.name = "Monte Carlo FRAP";

        RunMonteCarlo mc = project.monteCarlo;

        if (mc == null) {
            return true; // error
        }
        
        boolean fixedTissue = false;

        printRate = 0.004;
        
        FRAP_on = true;

        double baselinetime = 500;

        FRAP_onset = baselinetime + 70; // first pulse is 70 ms before bleaching pulse
        simTime = 10000 + FRAP_onset;
        
        //d1_Dshort = 0.060e-3;
        //immobileVesiclePercent = 0.25;
        
        double cubeWidth = 2;
        //double cubeWidth = 3;
        
        igm.mitoVolumeFraction = mitoVolumePercent;
        igm.mito_ijkselect = -2;
        
        if (fixedTissue) {
            mc.mitoCoordinatesOn = true; // mito are coordinate objects like vesicles
            mc.setMitoVolumeFraction = igm.mitoVolumeFraction;
            mc.mitoRadius = igm.mitoRadius;
            mc.mitoAxialRatio = igm.mitoAxialRatio;
            cubeWidth = 3;
            immobileVesiclePercent = 1;
        } else {
            igm.mitoNonSpaceVoxels = true; // mito are non-space voxels
        }

        dx = 2 * vesicleRadius;

        mc.Dcyto = Dcyto;

        igm.xdim = cubeWidth;
        igm.ydim = cubeWidth;
        igm.zdim = cubeWidth;

        if (initProject(dx, false)) {
            return true;
        }
        if (igm.initGeometryMC(iPSFselect)) {
            return true;
        }
        if (initExcitationPSF(igm.xvcenter, igm.yvcenter, igm.zvcenter)) {
            return true;
        }

        project.simTime = simTime;

        dv = initVesicles(d1_Dshort);

        dv.colorReady.setColor("[r=147,g=210,b=21]");
        dv.colorImmobile.setColor("[r=150,g=150,b=150]");

        dv.color.setColor("[r=0,g=0,b=153]");
        dv.color.setGradientColor("[r=255,g=255,b=255]");

        dv.setVolumeFraction = totalVesicleVolumeFraction;
        dv.setImmobilePercent = immobileVesiclePercent;
        
        if (initDetectors(igm.xvcenter, igm.yvcenter, igm.zvcenter)) {
            return true;
        }

        //mc.minVesicleStep = 0.001;
        mc.minVesicleStep = 0.002;

        mc.frapOn = true;
        mc.saveFluorescence = true;

        mc.setPSFSelect(iPSF, cPSF);

        if (fixedTissue) {
            mc.driftOn = true;
            mc.driftOnset = baselinetime;
            mc.PBC = true;
        }

        if (!Master.foundMainStartUpArguments) {
            mc.initAll(); // do not run when running batches
        }

        return false;

    }
    
    private boolean init_MC_FRAP_Movie() { // demo

        double dx;
        DiffusantVesicles dv;

        project.newMonteCarlo();
        project.name = "Monte Carlo FRAP movie";

        RunMonteCarlo mc = project.monteCarlo;

        if (mc == null) {
            return true; // error
        }

        printRate = 0.004;

        FRAP_on = true;

        FRAP_onset = 2;
        simTime = 60 + FRAP_onset;
        FRAP_small_on = false;

        double cubeWidth = 2;

        dx = 2 * vesicleRadius;

        mc.Dcyto = Dcyto;

        igm.xdim = cubeWidth;
        igm.ydim = cubeWidth;
        igm.zdim = cubeWidth;
        
        iPSFselect = 0; // (0) Gauss (1) Wilson (2) Torok (3) Impulse

        if (initProject(dx, false)) {
            return true;
        }
        if (igm.initGeometryMC(iPSFselect)) {
            return true;
        }
        if (initExcitationPSF(igm.xvcenter, igm.yvcenter, igm.zvcenter)) {
            return true;
        }

        dv = initVesicles(d1_Dshort);

        dv.colorReady.setColor("[r=102,g=255,b=0]");
        dv.colorReady.setGradientColor("[r=255,g=0,b=0]");
        
        dv.setVolumeFraction = 0.17;
        dv.setImmobilePercent = 0;

        if (initDetectors(igm.xvcenter, igm.yvcenter, igm.zvcenter)) {
            return true;
        }

        mc.minVesicleStep = 0.002;

        mc.frapOn = true;
        mc.saveFluorescence = false;
        mc.PBC = true;

        mc.setPSFSelect(iPSF, cPSF);
        
        if (!Master.foundMainStartUpArguments) {
            mc.initAll(); // do not run when running batches
        }

        return false;

    }
    
    private boolean init_MC_MSD_Prompt() {
        
        double VP, IP, MP;
        String VPstr, IPstr, MPstr;
        
        String pTitle = "Init Monte Carlo MSD";
        
        String VPvalues[] = new String[7];
        
        VPvalues[0] = "1";
        VPvalues[1] = "10";
        VPvalues[2] = "17";
        VPvalues[3] = "20";
        VPvalues[4] = "30";
        VPvalues[5] = "40";
        VPvalues[6] = "50";
        
        VPstr = Master.promptForInput(VPvalues, "choose vesicle volume % within non-mito volume", pTitle, "17");
        
        if ((VPstr == null) || (VPstr.length() == 0)) {
            return true; // error
        }
        
        VP = Double.parseDouble(VPstr);
        
        String IPvalues[] = new String[4];
        
        IPvalues[0] = "0";
        IPvalues[1] = "12.5";
        IPvalues[2] = "25.0"; // best-match MC simulation
        IPvalues[3] = "50.0";
        
        IPstr = Master.promptForInput(IPvalues, "choose immobile vesicle %", pTitle, "25.0");
        
        if ((IPstr == null) || (IPstr.length() == 0)) {
            return true; // error
        }
        
        IP = Double.parseDouble(IPstr);
        
        String MPvalues[] = new String[2];
        
        MPvalues[0] = "0";
        MPvalues[1] = "28";
        
        MPstr = Master.promptForInput(MPvalues, "choose mitochondria volume %", pTitle, "0");
        
        if ((MPstr == null) || (MPstr.length() == 0)) {
            return true; // error
        }
        
        MP = Double.parseDouble(MPstr);
        
        Master.log("vesicle volume %" + VP + "%, immobile vesicle percent = " + IP + "%, mito volume % = " + MP + "%");
        
        VP /= 100.0;
        IP /= 100.0;
        MP /= 100.0;
        
        totalVesicleVolumeFraction = VP;
        immobileVesiclePercent = IP;
        mitoVolumePercent = MP;
        
        return init_MC_MSD();
        
    }

    private boolean init_MC_MSD() {
        // Figure 4C
        // Figure 5A,C

        double dx;
        DiffusantVesicles dv;

        project.newMonteCarlo();
        project.name = "Monte Carlo MSD";

        RunMonteCarlo mc = project.monteCarlo;

        if (mc == null) {
            return true; // error
        }

        printRate = 0.02;
        saveRate = 20; // kHz

        simTime = 300;

        double cubeWidth = 2.5;

        //totalVesicleVolumeFraction = 0.17;
        //immobileVesiclePercent = 0.25;
        //mitoVolumePercent = 0.28;

        igm.mitoNonSpaceVoxels = true; // mito are non-space voxels
        igm.mitoVolumeFraction = mitoVolumePercent;
        igm.mito_ijkselect = -2;
        
        if (igm.mitoVolumeFraction > 0) {
            
            simTime = 50000;
            
            if (totalVesicleVolumeFraction > 0.01) {
                cubeWidth = 2.0;
                simTime = 30000;
            }
            if (totalVesicleVolumeFraction > 0.10) {
                simTime = 18000;
            }
            if (totalVesicleVolumeFraction > 0.17) {
                simTime = 12000;
            }
            if (totalVesicleVolumeFraction > 0.20) {
                simTime = 6000;
            }
            if (totalVesicleVolumeFraction > 0.30) {
                simTime = 1000;
            }
            if (totalVesicleVolumeFraction > 0.40) {
                simTime = 500;
            }

        } else {
            
            igm.mitoVolumeFraction = 0;

            if (totalVesicleVolumeFraction > 0.11) {
                cubeWidth = 2.0;
            }
            if (totalVesicleVolumeFraction > 0.18) {
                cubeWidth = 1.0;
            }

        }

        dx = 2 * vesicleRadius;
        
        mc.Dcyto = Dcyto;

        igm.xdim = cubeWidth;
        igm.ydim = cubeWidth;
        igm.zdim = cubeWidth;

        if (initProject(dx, false)) {
            return true;
        }
        if (igm.initGeometryMC(iPSFselect)) {
            return true;
        }

        dv = initVesicles(d1_Dshort);

        dv.setVolumeFraction = totalVesicleVolumeFraction;
        dv.setImmobilePercent = immobileVesiclePercent;
        dv.saveMSD = true;
        dv.init();

        //mc.minVesicleStep = 0.0005;
        //mc.minVesicleStep = 0.0010;
        mc.minVesicleStep = 0.0015;

        mc.PBC = true;

        if (dv.setVolumeFraction > 0.3) {
            mc.overlapTrialLimit = 10;
        }

        if (!Master.foundMainStartUpArguments) {
            mc.initAll(); // // do not run when running batches
        }

        return false;

    }
    
    private boolean init_MC_Cichocki() {
        // this code replicates Figure 1 of Cichocki Hinsen (not in paper)
        // Cichocki Hinsen 1990 DYNAMIC COMPUTER SIMULATION OF CONCENTRATED HARD SPHERE SUSPENSIONS

        double dx, lambda, step, t0;

        DiffusantVesicles dv;
        
        project.newMonteCarlo();
        project.name = "Monte Carlo MSD Cichocki";
        
        RunMonteCarlo mc = project.monteCarlo;

        if (mc == null) {
            return true; // error
        }
        
        int lambda_select = 0; // 0, 1 or 2

        switch(lambda_select) {
            case 0:
                lambda = 0.035;
                break;
            case 1:
                lambda = 0.025;
                break;
            case 2:
                lambda = 0.018;
                break;
            default:
                return true;
        }

        saveRate = 500;
        printRate = 0.2;
        
        double cubeWidth = 1.0;

        double vesicleVolumeFraction = 0.20;

        dx = 2 * vesicleRadius;

        d1_Dshort = 0.001;

        igm.xdim = cubeWidth;
        igm.ydim = cubeWidth;
        igm.zdim = cubeWidth;

        step = lambda * vesicleRadius;

        System.out.println("dr = " + step);

        t0 = vesicleRadius * vesicleRadius / ( 6 * d1_Dshort);

        System.out.println("t0 = " + t0);

        simTime = 16 * t0 * 20;

        if (initProject(dx, false)) {
            return true;
        }

        if (igm.initGeometryMC(iPSFselect)) {
            return true;
        }

        dv = initVesicles(d1_Dshort);
        
        immobileVesiclePercent = 0;

        dv.setVolumeFraction = vesicleVolumeFraction;
        dv.setImmobilePercent = immobileVesiclePercent;
        dv.saveXYZ = false;
        dv.saveMSD = true;

        mc.minVesicleStep = step;
        mc.freeDiffusion = false;
        mc.PBC = true;

        if (!Master.foundMainStartUpArguments) {
            mc.initAll(); // // do not run when running batches
        }

        return false;

    }

    public boolean init_FD_FRAP() {
        // Figure 2C, 4D
        // Dlong varied between 0.023 - 0.031e-3
        // immobileVesiclePercent varied between 0.25 - 0.32

        double Dlong = 0.028e-3; // best-match FD simulation
        immobileVesiclePercent = 0.29; // best-match FD simulation
        
        project.newFiniteDifference();
        project.name = "Finite Difference FRAP";

        printRate = 0.001;

        double dx = 0.05;
        double dim = 3;

        igm.xdim = dim;
        igm.ydim = dim;
        igm.zdim = dim;
        
        igm.eighth_geometry = true;

        d1_name = "Mobile";
        d1_C = 1.0 - immobileVesiclePercent;

        FRAP_on = true;
        FRAP_onset = 100;
        simTime = 10000 + FRAP_onset;

        //FRAP_small_on = false;
        //FRAP_bleach_on = false;

        if (initProject(dx, true)) {
            return true;
        }
        if (igm.initGeometryMC(iPSFselect)) {
            return true;
        }
        if (initExcitationPSF(igm.xvcenter, igm.yvcenter, igm.zvcenter)) {
            return true;
        }
        if (initDiffusantsFD(Dlong, immobileVesiclePercent)) {
            return true;
        }
        if (initDetectors(igm.xvcenter, igm.yvcenter, igm.zvcenter)) {
            return true;
        }

        return false;

    }
    
    public boolean init_FD_Axelrod() {
        // Axelrod conditions where iPSF = cPSF = Gaussian beams
        // not in paper, but results show nearly identical FRAP curves to the Axelrod curve fits

        double Dlong = 0.025e-3; // fit to Axelrod equation
        immobileVesiclePercent = 0.29;
        
        project.newFiniteDifference();
        project.name = "Finite Difference FRAP Axelrod";

        printRate = 0.001;

        double dx = 0.05;
        double dim = 3;

        igm.xdim = dim;
        igm.ydim = dim;
        igm.zdim = dim;

        igm.eighth_geometry = true;

        d1_name = "Mobile";
        d1_C = 1.0 - immobileVesiclePercent;

        FRAP_on = true;
        FRAP_onset = 100;
        
        simTime = 10000 + FRAP_onset;

        iPSFselect = 0; // (0) Gauss (1) Wilson (2) Torok (3) Impulse

        iPSF_fwhm_z = -1; // Gaussian beam
        iPSF_fwhm_xy = 0.2708;

        cPSF_fwhm_z = -1;
        cPSF_fwhm_xy = iPSF_fwhm_xy; // cPSF = iPSF

        FRAP_kPhoto = 1.0;

        if (initProject(dx, true)) {
            return true;
        }
        if (igm.initGeometryMC(iPSFselect)) {
            return true;
        }
        if (initExcitationPSF(igm.xvcenter, igm.yvcenter, igm.zvcenter)) {
            return true;
        }
        if (initDiffusantsFD(Dlong, immobileVesiclePercent)) {
            return true;
        }
        if (initDetectors(igm.xvcenter, igm.yvcenter, igm.zvcenter)) {
            return true;
        }

        return false;

    }

    private boolean initProject(double dx, boolean FD) {

        if (!Master.foundMainStartUpArguments) {
           project.simTime = simTime;
           project.directory = directory;
           project.folder = folder;
        }

        project.dx = dx;
        project.set("saveRate", saveRate);
        project.set("printRate", printRate);

        if (FD) {
            project.stability = stabilityFD;
        }

        return false;

    }

    private boolean initExcitationPSF(double xvcenter, double yvcenter, double zvcenter) {

        double sdxy, sdz;

        switch (iPSFselect) {

            case 0:

                sdxy = PSFgauss.computeSTDV(iPSF_fwhm_xy);
                sdz = PSFgauss.computeSTDV(iPSF_fwhm_z);

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
                gpsf.name = "Gaussian iPSF";

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
                wpsf.name = "Wilson iPSF";

                iPSF = geometry.addPSF(wpsf);

                break;

            case 2:

                tpsf = new PSFtorok(project, null);
                tpsf.setFitZoltan3();

                if (igm.eighth_geometry) {
                    tpsf.arbitraryVoxelCenter = true;
                    tpsf.xySymmetric = false;
                } else {
                    tpsf.xySymmetric = true;
                }

                tpsf.set("xVoxelCenter", xvcenter);
                tpsf.set("yVoxelCenter", yvcenter);
                tpsf.set("zVoxelCenter", zvcenter);
                tpsf.name = "Torok iPSF";

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
                epsf.name = "Simple iPSF";

                iPSF = geometry.addPSF(epsf);

                break;

            default:
                return true;
        }

        return false;

    }

    private boolean initDiffusantsFD(double Dlong1, double immobileFraction) {

        int dnum = -1;

        initPulseTimer();

        switch (iPSFselect) {
            default:
                frap = new DiffusantPhoto(project, d1_name, d1_C, Dlong1, null, dnum, timer, gpsf, FRAP_kPhoto);
                break;
            case 1:
                frap = new DiffusantPhoto(project, d1_name, d1_C, Dlong1, null, dnum, timer, wpsf, FRAP_kPhoto);
                break;
            case 2:
                frap = new DiffusantPhoto(project, d1_name, d1_C, Dlong1, null, dnum, timer, tpsf, FRAP_kPhoto);
                break;
            case 3:
                frap = new DiffusantPhoto(project, d1_name, d1_C, Dlong1, null, dnum, timer, epsf, FRAP_kPhoto);
                break;
        }

        d1_num = project.addDiffusant(frap);

        if ((immobileFraction >= 0) && (immobileFraction < 1)) {

            switch (iPSFselect) {
                default:
                    frap2 = new DiffusantPhoto(project, "Immobile", immobileFraction, 0, null, dnum, timer, gpsf, FRAP_kPhoto);
                    break;
                case 1:
                    frap2 = new DiffusantPhoto(project, "Immobile", immobileFraction, 0, null, dnum, timer, wpsf, FRAP_kPhoto);
                    break;
                case 2:
                    frap2 = new DiffusantPhoto(project, "Immobile", immobileFraction, 0, null, dnum, timer, tpsf, FRAP_kPhoto);
                    break;
                case 3:
                    frap2 = new DiffusantPhoto(project, "Immobile", immobileFraction, 0, null, dnum, timer, epsf, FRAP_kPhoto);
                    break;
            }

            d2_num = project.addDiffusant(frap2);

        }

        return false;

    }

    private void initPulseTimer() {

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

    private DiffusantVesicles initVesicles(double D) {

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

    private boolean initDetectors(double xvcenter, double yvcenter, double zvcenter) {

        double sdxy = PSFgauss.computeSTDV(cPSF_fwhm_xy);
        double sdz = PSFgauss.computeSTDV(cPSF_fwhm_z);

        PSFgauss dpsf = new PSFgauss(project, null, sdxy, sdxy, sdz);

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

        cPSF = geometry.addPSF(dpsf);

        return false;

    }
    
    private void writePSF(double[][][] psf) {

        for (int i = 0; i < psf.length; i++) {
            for (int j = 0; j < psf[0].length; j++) {
                for (int k = 0; k < psf[0][0].length; k++) {
                    Master.log("xyz[" + i + "][" + j + "][" + k + "] = " + psf[i][j][k] + ";");
                }
            }
        }

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
