package ucl.silver.d3d.init;

import ucl.silver.d3d.core.*;
import ucl.silver.d3d.utils.*;
import java.io.*;
import java.util.Scanner;

/**
 * <p>Title: D3D</p>
 * <p>Description: 3D Diffusion-Reaction Simulator</p>
 * <p>Copyright: Copyright (c) 2016</p>
 * <p>Company: The Silver Lab at University College London</p>
 * @author Jason S. Rothman
 * @version 1.6
 */

public class InitAZEMJason extends InitProject {

    public String directory = "/Jason/D3D/Simulations/";
    public String folder = "Testing";
    public String EMdirectory = ""; // "/Jason/FrapVGlut1/Zoltan/EM/";

    public double simTime = 100000; // msec
    public double saveRate = 20; // kHz
    public double printRate = 0.004; // kHz

    public double vesicleRadius = 0.044/2.0; // um

    public double Dshort = 0.060e-3;
    
    public int EMseries = 3;
    public int EMseriesAZ = 0;
    
    //public boolean openGeometry = true;
    public boolean openGeometry = false;
    
    //public boolean connectVesicles = true;
    public boolean connectVesicles = false;
    
    public boolean stim100Hz = false;
    //public boolean stim100Hz = true;

    private final InitGeometryMito igm;

    public InitAZEMJason(Project p) {
        super(p);
        igm = new InitGeometryMito(p);
        String[] flist = {"init_MC_AZ", "init_MC_AZ_Prompt"};
        initFuncList = flist;
        createVector(true);
    }

    @Override
    public boolean initFunction(String initSelect) {

        int i = initFunctionNum(initSelect);

        if (!Master.foundMainStartUpArguments) {
            //i = 0;
        }

        switch (i) {
            
            case 0:
                return init_MC_AZ();

            case 1:
                return init_MC_AZ_Prompt();

            default:
                Master.log("InitAZEMJason.initFunction error: failed to find init function " + initSelect);
                return true; // error

        }

    }
    
    private boolean init_MC_AZ_Prompt() {
        
        int EMvalue, AZvalue;
        String EMstr, AZstr, YNstr;
        
        String EMvalues[] = new String[7];
        
        EMvalues[0] = "1";
        EMvalues[1] = "3";
        EMvalues[2] = "4";
        EMvalues[3] = "5";
        EMvalues[4] = "10";
        EMvalues[5] = "16";
        EMvalues[6] = "18";
        
        EMstr = Master.promptForInput(EMvalues, "choose EM series #", "Init Monte Carlo AZ", "3");
        
        if ((EMstr == null) || (EMstr.length() == 0)) {
            return true; // cancel
        }
        
        EMvalue = Integer.parseInt(EMstr);
        
        String AZvalues[] = null;
        
        switch (EMvalue) {
            case (4):
                AZvalues = new String[4];
                AZvalues[0] = "0";
                AZvalues[1] = "1";
                AZvalues[2] = "2";
                AZvalues[3] = "3";
                break;
            case (10):
                AZvalues = new String[5];
                AZvalues[0] = "0";
                AZvalues[1] = "1";
                AZvalues[2] = "2";
                AZvalues[3] = "3";
                AZvalues[4] = "4";
                break;
        }
        
        if (AZvalues == null) {
            
            AZvalue = 0;
            
        } else {

            AZstr = Master.promptForInput(AZvalues, "choose AZ #", "Init Monte Carlo AZ", "0");

            if ((AZstr == null) || (AZstr.length() == 0)) {
                return true; // cancel
            }

            AZvalue = Integer.parseInt(AZstr);

        }
        
        EMseries = EMvalue;
        EMseriesAZ = AZvalue;
        
        String YesNo[] = new String[2];
        YesNo[0] = "Yes";
        YesNo[1] = "No";
        
        YNstr = Master.promptForInput(YesNo, "simulate with open geometry configuration?", "Init Monte Carlo AZ", "No");
        
        if ((YNstr == null) || (YNstr.length() == 0)) {
            return true; // cancel
        }
        
        openGeometry = YNstr.equalsIgnoreCase("Yes");
        
        YNstr = Master.promptForInput(YesNo, "simulate with connectors/tethers?", "Init Monte Carlo AZ", "No");
        
        if ((YNstr == null) || (YNstr.length() == 0)) {
            return true; // cancel
        }
        
        connectVesicles = YNstr.equalsIgnoreCase("Yes");
        
        YNstr = Master.promptForInput(YesNo, "simulate with docking/priming and 100 Hz stimulation?", "Init Monte Carlo AZ", "No");
        
        if ((YNstr == null) || (YNstr.length() == 0)) {
            return true; // cancel
        }
        
        stim100Hz = YNstr.equalsIgnoreCase("Yes");
        
        if (openGeometry) {
            Master.log("simulating with open geometry");
        }
        
        if (connectVesicles) {
            Master.log("simulating with connectors/tethers");
        }
        
        if (stim100Hz) {
            Master.log("simulating with docking/priming and 100 Hz stimulation");
        }
        
        return init_MC_AZ();
        
    }

    private boolean init_MC_AZ() {
        // Figures 6, 7 and 8

        int nvesicles;
        DiffusantVesiclesAZ dv;

        project.newMonteCarloAZEM();
        project.name = "Monte Carlo AZ EM";

        RunMonteCarloAZEM mc = (RunMonteCarloAZEM) project.monteCarlo;

        if (mc == null) {
            return true; // error
        }

        if (!Master.foundMainStartUpArguments) {
            simTime = 100000; // 100;
        }

        vesicleRadius = 0.044/2.0;
        double dx = 0.044;

        //mc.minVesicleStep = 0.001;
        //mc.minVesicleStep = 0.002;
        mc.minVesicleStep = 0.005;

        boolean EManalysis = false; // Figure 7, 8
        //boolean EManalysis = true; // Figure 6

        boolean extraReserve = true;
        //boolean extraReserve = false;

        //boolean emptyGeometry = true;
        boolean emptyGeometry = false;

        boolean hydrodynamics = true;
        //boolean hydrodynamics = false;
       
        String wrlFileName;

        if (EManalysis) {
            extraReserve = false;
            openGeometry = false;
            connectVesicles = false;
            hydrodynamics = false;
            stim100Hz = false;
        }

        if (emptyGeometry) {
            EManalysis = false;
            openGeometry = false;
            hydrodynamics = false;
        }
        
        //EMseries = 3; // 1, 3, 4, 5, 10, 16, 18
        //EMseriesAZ = 0;
        
        Master.log("initializing EM series #" + EMseries + ", AZ #" + EMseriesAZ);

        if (openGeometry) {
            mc.importVesicleXYZ_dr = "/Jason/FrapVGlut1/Simulations/Sims24_AZ_SupplyRate_EMs/VesiclePositionsOpen/";
        } else {
            mc.importVesicleXYZ_dr = "/Jason/FrapVGlut1/Simulations/Sims24_AZ_SupplyRate_EMs/VesiclePositions/";
        }

        project.batchNum = 0;

        wrlFileName = "ser " + Integer.toString(EMseries) + ".wrl";
        
        boolean importEMdata = false;
        
        double[][][] azXYZtemp, azXYZ;
        double[][] vesicleXYZ, vxyz;
        
        if (importEMdata) {
            Master.log("importing EM file " + EMdirectory + wrlFileName);
            azXYZtemp = importAZs(EMdirectory + wrlFileName);
            //writeActiveZoneArray(azXYZtemp);
        } else {
            azXYZtemp = InitAZEMJason2.getAZs(EMseries);
        }

        if (azXYZtemp == null) {
            return false;
        }

        azXYZ = azInterpolate(azXYZtemp, 0.01, 2); // dx = 0.015 // USED IN SIMS
        //azXYZ = azXYZtemp;

        if (azXYZ == null) {
            return true;
        }
        
        if (importEMdata) {
            
            vesicleXYZ = importVesicles(EMdirectory + wrlFileName);

            if (vesicleXYZ == null) {
                return true;
            }

            double xyTolerance = 0.010;

            removeVesicleDuplicates(vesicleXYZ, xyTolerance);

            vxyz = removeNANs(vesicleXYZ);

            //writeVesicleArray(vxyz);

        } else {
            
            vxyz = InitAZEMJason2.getVesicles(EMseries);

        }
        
        nvesicles = vxyz.length;
        
        Master.log("number of EM vesicles = " + nvesicles);

        mc.EMseries = EMseries;
        mc.EMseriesAZ = EMseriesAZ;
        mc.openGeometry = openGeometry;
        mc.extraReserve = extraReserve;

        mc.emptyGeometry = emptyGeometry;

        if (emptyGeometry) {
            //mc.PBC = true;
        }

        configureGeometry(mc, azXYZ, vxyz, dx);

        if (initProject(dx)) {
            return true;
        }

        if (igm.initGeometryMC(-1)) {
            return true;
        }

        double totalVesicleVolumeFraction = 0.17;
        double setImmobileVesiclePercent = 0.25;

        double mobileVolumeFraction = totalVesicleVolumeFraction * (1 - setImmobileVesiclePercent);
        double immobileVolumeFraction = totalVesicleVolumeFraction * setImmobileVesiclePercent;

        mc.Dcyto = 0.060e-3 / DiffusantVesicle.Dratio_short(0.17 * 0.75, 0.17 * 0.25); // 1.2696940661239775E-4

        if (hydrodynamics) {
            Dshort = mc.Dcyto * DiffusantVesicle.Dratio_short(mobileVolumeFraction, immobileVolumeFraction);
        } else {
            Dshort = 0.060e-3; // mc.Dcyto;
        }

        dv = initVesicles(Dshort);
        dv.setVolumeFraction = 0;
        dv.setImmobilePercent = setImmobileVesiclePercent;
        dv.setDensity = 0;
        dv.numVesicles = vxyz.length;
        dv.setNumVesiclesReady = vxyz.length;
        dv.xyz = vxyz;

        dv.colorReady.setColor("[r=147,g=210,b=21]");
        dv.colorImmobile.setColor("[r=150,g=150,b=150]");

        dv.saveXYZ = true;
        
        mc.azXYZ = azXYZ;

        //mc.saveMSDspatial = true;
        mc.MSDspatialTbgn = 200;

        mc.extendAZ = 0.01; // um, extend AZ boundary a little for "inside az" macro

        if (EManalysis) {
            mc.importVesicleXYZ = false;
            mc.fixNonSpaceVoxels = false;
            mc.removeShrinkage = false;
            mc.removeVesicleOverlap = false;
        } else {
            //mc.importVesicleXYZ = true;
            mc.fixNonSpaceVoxels = true;
            mc.removeShrinkage = true;
            mc.removeVesicleOverlap = true;
            //mc.removeVesicleOverlap = false;
        }

        mc.reserveImmobileFraction = setImmobileVesiclePercent;

        mc.initVesicleRandom = false;
        
        if (stim100Hz) {
            mc.releaseRate = 0.1; // kHz
            mc.releaseProb = 0.5;
            mc.dockRefractoryPeriod = 1000.0/60.0;
        } else {
            mc.releaseRate = Double.POSITIVE_INFINITY; // kHz
            mc.releaseProb = 1.0;
        }
        
        mc.releaseLimit = 0; // limit of number of release events

        if (openGeometry) {
            mc.replenishRate = Double.POSITIVE_INFINITY; // kHz
        }

        mc.dockingOn = !emptyGeometry;
        
        mc.dockedVesiclesMax = 1;
        //mc.dockedVesiclesMax = 2;
        
        mc.saveRelease = mc.releaseRate > 0;

        //mc.timer = false;

        if (hydrodynamics) {
            mc.hydroWallZ = true;
            mc.hydrodynamicsLocalD = true;
            mc.hydrodynamicsLocalDVoxels = true;
        }

        //mc.freeDiffusion = true; // no steric interactions

        if (connectVesicles) {

            mc.connectVesicles = true;
            mc.maxNumConnectors = 3;
            mc.connectRate = 10; // P/ms
            mc.unconnectRate = 0.01;
            //mc.unconnectRate = 1.0;
            mc.azPermitConnectorRadius = 0.15; // um

            mc.tetherVesicles = true;
            mc.azNumTethers = 2;

        }

        if (!Master.foundMainStartUpArguments) {
            mc.initAll(); // do not run when running batches
        }

        return false;

    }

    private double[][][] importAZs(String WRL_fileName) {

        String s;
        String tagAZ = "# AZ";
        String tagPoint = "point [";
        String ss;

        int iAZ = 0, ipoint = 0, ilimit, k1, k2;

        boolean finished;

        int numAZs = 5;
        int numPoints = 200;

        double[][][] xyz = new double[numAZs][numPoints][3];

        int tagAZLength = tagAZ.length();
        int tagPointLength = tagPoint.length();

        for (int i = 0; i < xyz.length; i++) {
            for (int j = 0; j < xyz[0].length; j++) {
                xyz[i][j][0] = Double.NaN;
                xyz[i][j][1] = Double.NaN;
                xyz[i][j][2] = Double.NaN;
            }
        }

        try {

            s = new Scanner(new File(WRL_fileName)).useDelimiter("\\A").next();
            ilimit = s.length();

            for (int i = 0; i < ilimit-tagAZLength; i++) {

                if (s.substring(i, i + tagAZLength).equalsIgnoreCase(tagAZ)) {

                    finished = false;

                    //Master.log("importing AZ #" + iAZ + " points...");

                    for (int j = i + tagAZLength; j < ilimit - tagPointLength; j++) {

                        if (s.substring(j, j + tagPointLength).equalsIgnoreCase(tagPoint)) {

                            k1 = j + tagPointLength + 4;

                            for (int k = k1; k < ilimit; k++) {

                                if (s.substring(k, k+1).equalsIgnoreCase("]")) {
                                    finished = true;
                                    break;
                                }

                                if (s.substring(k, k + 1).equalsIgnoreCase(",")) {

                                    k2 = k + 1;
                                    ss = s.substring(k1, k2);
                                    ss = ss.replaceAll(" ", ",");

                                    xyz[iAZ][ipoint][0] = Double.parseDouble(D3DUtility.itemFromList(ss, 0));
                                    xyz[iAZ][ipoint][1] = Double.parseDouble(D3DUtility.itemFromList(ss, 1));
                                    xyz[iAZ][ipoint][2] = Double.parseDouble(D3DUtility.itemFromList(ss, 2));
                                    //Master.log("x=" + xyz[ipoint][0] + ",y=" + xyz[ipoint][1] + ",z=" + xyz[ipoint][2]);
                                    //Master.log("" + xyz[iAZ][ipoint][2]);

                                    k1 = k2 + 4;

                                    ipoint++;

                                }

                            }

                            if (ipoint >= numPoints) {
                                Master.log("need bigger vesicle array");
                                return null;
                            }

                        }

                        if (finished) {
                            break;
                        }

                    }

                    Master.log("imported AZ #" + iAZ + " positions n = " + ipoint);

                    iAZ++;
                    ipoint = 0;

                }

            }

        } catch (IOException e) {
            System.err.println("unable to read file: " + WRL_fileName);
            return null;
        }

        return xyz;

    }

    private double azMaxSegment(double[][][] azXYZ) {

        double d, maxD = 0;
        double ax1, ay1, az1;
        double ax2, ay2, az2;

        if (azXYZ == null) {
            return Double.NaN;
        }

        for (int iAZ = 0; iAZ < azXYZ.length; iAZ++) {
            for (int ipnt = 1; ipnt < azXYZ[0].length; ipnt++) {

                ax1 = azXYZ[iAZ][ipnt - 1][0];
                ay1 = azXYZ[iAZ][ipnt - 1][1];
                az1 = azXYZ[iAZ][ipnt - 1][2];

                ax2 = azXYZ[iAZ][ipnt][0];
                ay2 = azXYZ[iAZ][ipnt][1];
                az2 = azXYZ[iAZ][ipnt][2];

                if (Double.isNaN(ax1 * ay1 * az1 * ax2 * ay2 * az2)) {
                    continue;
                }

                if (az1 == az2) {

                    d = Math.sqrt((ax2 - ax1) * (ax2 - ax1) + (ay2 - ay1) * (ay2 - ay1) + (az2 - az1) * (az2 - az1));

                    if (d > maxD) {
                        maxD = d;
                    }

                }

            }

        }

        return maxD;

    }

    private double[][][] azInterpolate(double[][][] azXYZ, double maxSegmentLength, int numZstacks) {

        int i, ii, tsteps;
        double d, t;
        double ax1, ay1, az1;
        double ax2, ay2, az2;
        double mx, my, mz;

        if (azXYZ == null) {
            return null;
        }

        double zstackincrement = 0.03 / (1.0 * numZstacks);

        double largestSegmentLength = azMaxSegment(azXYZ);

        Master.log("max AZ segment length = " + largestSegmentLength);

        int npnts = (int) Math.ceil( largestSegmentLength / maxSegmentLength);

        npnts = azXYZ[0].length * npnts * numZstacks;

        double[][][] az = new double[azXYZ.length][npnts][3];

        double[][] aztemp1;
        double[][] aztemp2 = new double[npnts][3];

        for (int iAZ = 0; iAZ < az.length; iAZ++) {
            for (int ipnt = 0; ipnt < az[0].length; ipnt++) {
                az[iAZ][ipnt][0] = Double.NaN;
                az[iAZ][ipnt][1] = Double.NaN;
                az[iAZ][ipnt][2] = Double.NaN;
            }
        }

        for (int iAZ = 0; iAZ < azXYZ.length; iAZ++) {

            i = 0;

            for (int k = 1; k <= 10; k++) {

                aztemp1 = getZstackAZ(azXYZ, iAZ, k);

                if (aztemp1 == null) {
                    break;
                }

                for (int ipnt = 1; ipnt < aztemp2.length; ipnt++) {
                    aztemp2[ipnt][0] = Double.NaN;
                    aztemp2[ipnt][1] = Double.NaN;
                    aztemp2[ipnt][2] = Double.NaN;
                }

                aztemp2[0][0] = aztemp1[0][0];
                aztemp2[0][1] = aztemp1[0][1];
                aztemp2[0][2] = aztemp1[0][2];

                ii = 1;

                for (int ipnt = 1; ipnt < aztemp1.length; ipnt++) {

                    ax1 = aztemp1[ipnt - 1][0];
                    ay1 = aztemp1[ipnt - 1][1];
                    az1 = aztemp1[ipnt - 1][2];

                    ax2 = aztemp1[ipnt][0];
                    ay2 = aztemp1[ipnt][1];
                    az2 = aztemp1[ipnt][2];

                    if (Double.isNaN(ax1 * ay1 * az1 * ax2 * ay2 * az2)) {
                        break;
                    }

                    d = Math.sqrt((ax2 - ax1) * (ax2 - ax1) + (ay2 - ay1) * (ay2 - ay1) + (az2 - az1) * (az2 - az1));

                    if (d < maxSegmentLength) {

                        aztemp2[ii][0] = ax2;
                        aztemp2[ii][1] = ay2;
                        aztemp2[ii][2] = az2;

                        ii++;

                    } else {

                        mx = ax2 - ax1;
                        my = ay2 - ay1;
                        mz = az2 - az1;

                        tsteps = (int) (d / maxSegmentLength);

                        if (tsteps < 2) {
                            tsteps = 2;
                        }

                        for (int it = 1; it <= tsteps; it++) {

                            t = it / (tsteps * 1.0);

                            aztemp2[ii][0] = ax1 + t * mx;
                            aztemp2[ii][1] = ay1 + t * my;
                            aztemp2[ii][2] = az1 + t * mz;

                            ii++;

                        }

                    }

                }

                // copy temp array to new az array

                for (int iz = 1; iz <= numZstacks; iz++) {

                    //Master.log("" + iAZ + "," + iz + "," + i);
                    for (double[] xyz : aztemp2) {
                        if (Double.isNaN(xyz[0])) {
                            break;
                        }
                        az[iAZ][i][0] = xyz[0];
                        az[iAZ][i][1] = xyz[1];
                        az[iAZ][i][2] = xyz[2] + zstackincrement * (iz - 1);
                        i++;
                    }

                }

            }

        }

        Master.log("new max AZ segment length = " + azMaxSegment(az));

        return az;

    }

    private double[][] getZstackAZ(double[][][] azXYZ, int azSelect, int zSelect) {

        double z, lastZ = -1;
        int ii = 0, counter = 0;
        boolean found = false;

        double[][] aztemp = new double[azXYZ[0].length][3];

        for (double[] xyz : aztemp) {
            xyz[0] = Double.NaN;
            xyz[1] = Double.NaN;
            xyz[2] = Double.NaN;
        }

        for (int ipnt = 0; ipnt < azXYZ[0].length; ipnt++) {

            z = azXYZ[azSelect][ipnt][2];

            if (Double.isNaN(z)) {
                break;
            }

            if (z != lastZ) {
                counter++;
            }

            if (counter > zSelect) {
                break;
            }

            if (counter == zSelect) {
                aztemp[ii][0] = azXYZ[azSelect][ipnt][0];
                aztemp[ii][1] = azXYZ[azSelect][ipnt][1];
                aztemp[ii][2] = azXYZ[azSelect][ipnt][2];
                ii++;
                found = true;
            }

            lastZ = z;

        }

        if (found) {
            return aztemp;
        } else {
            return null;
        }

    }

    private double[][] importVesicles(String WRL_fileName) {

        String s;
        int kk, ivesicle, imax = 0;

        boolean finished;
        boolean foundx, foundy, foundz;

        int numVesicles = 10000;

        double[][] xyz = new double[numVesicles][3];

        for (double[] d : xyz) {
            d[0] = Double.NaN;
            d[1] = Double.NaN;
            d[2] = Double.NaN;
        }

        try {

            s = new Scanner(new File(WRL_fileName)).useDelimiter("\\A").next();

            for (int i = 0; i < s.length()-10; i++) {

                if (s.substring(i, i + 10).equalsIgnoreCase("# vesicles")) {

                    ivesicle = -1;

                    for (int j = i + 10; j < i + 10 + 10; j++) {
                        if (s.substring(j, j + 1).equalsIgnoreCase(" ")) {
                            ivesicle = Integer.parseInt(s.substring(i + 10, j));
                            break;
                        }
                    }

                    if (ivesicle <= 0) {
                        return null;
                    }

                    if (ivesicle >= numVesicles) {
                        Master.log("need bigger vesicle array");
                        return null;
                    }

                    if (ivesicle > imax) {
                        imax = ivesicle;
                    }

                    finished = false;
                    foundx = false;
                    foundy = false;
                    foundz = false;

                    kk = 0;

                    for (int j = i + 10; j < i + 10 + 50; j++) {

                        if (s.substring(j, j + 8).equalsIgnoreCase("center: ")) {

                            for (int k = j + 8; k < j + 8 + 12; k++) {
                                if (s.substring(k, k+1).equalsIgnoreCase(" ")) {
                                    xyz[ivesicle][0] = Double.parseDouble(s.substring(j + 8, k));
                                    foundx = true;
                                    kk = k+1;
                                    break;
                                }
                            }

                            if (!foundx) {
                                return null;
                            }

                            for (int k = kk; k < kk + 12; k++) {
                                if (s.substring(k, k+1).equalsIgnoreCase(" ")) {
                                    xyz[ivesicle][1] = Double.parseDouble(s.substring(kk, k));
                                    foundy = true;
                                    kk = k+1;
                                    break;
                                }
                            }

                            if (!foundy) {
                                return null;
                            }

                            for (int k = kk; k < kk + 12; k++) {
                                if (s.substring(k, k+1).equalsIgnoreCase(" ")) {
                                    xyz[ivesicle][2] = Double.parseDouble(s.substring(kk, k));
                                    foundz = true;
                                    kk = k+1;
                                    break;
                                }
                            }

                            if (!foundz) {
                                return null;
                            }

                            if ((xyz[ivesicle][0] < 0) || (xyz[ivesicle][1] < 0) || (xyz[ivesicle][2] < 0)) {
                                xyz[ivesicle][0] = Double.NaN;
                                xyz[ivesicle][1] = Double.NaN;
                                xyz[ivesicle][2] = Double.NaN;
                            }

                            //Master.log("x=" + xyz[ivesicle][0] + ", y=" + xyz[ivesicle][1] + ", z=" + xyz[ivesicle][2]);
                            //Master.log("" + xyz[ivesicle][2]);

                            //ivesicle++;

                            finished = true;


                        }

                        if (finished) {
                            break;
                        }

                    }

                }

            }

        } catch (IOException e) {
            System.err.println("unable to read file: " + WRL_fileName);
            return null;
        }

        Master.log("imported vesicles positions n = " + imax);

        return xyz;

    }

    private double[][] importVesicles2(String WRL_fileName) {

        String s, ss;
        int k1, k2, ilimit, ivesicle, imax = 0;
        double xvalue, yvalue, zvalue;
        double xsum, ysum, zsum, ipoints;

        boolean finished;

        int numVesicles = 10000;

        double[][] vxyz = new double[numVesicles][3];

        String tagVesicle = "# vesicles";
        String tagPoint = "point [";

        int tagVesicleLength = tagVesicle.length();
        int tagPointLength = tagPoint.length();

        for (double[] xyz : vxyz) {
            xyz[0] = Double.NaN;
            xyz[1] = Double.NaN;
            xyz[2] = Double.NaN;
        }

        try {

            s = new Scanner(new File(WRL_fileName)).useDelimiter("\\A").next();
            ilimit = s.length();

            for (int i = 0; i < s.length()-10; i++) {

                if (s.substring(i, i + tagVesicleLength).equalsIgnoreCase(tagVesicle)) {

                    ivesicle = -1;

                    for (int j = i + tagVesicleLength; j < i + tagVesicleLength + 10; j++) {
                        if (s.substring(j, j + 1).equalsIgnoreCase(" ")) {
                            ivesicle = Integer.parseInt(s.substring(i + tagVesicleLength, j));
                            break;
                        }
                    }

                    if (ivesicle <= 0) {
                        return null;
                    }

                    if (ivesicle >= numVesicles) {
                        Master.log("need bigger vesicle array");
                        return null;
                    }

                    if (ivesicle > imax) {
                        imax = ivesicle;
                    }

                    finished = false;

                    xsum = 0;
                    ysum = 0;
                    zsum = 0;
                    ipoints = 0;

                    for (int j = i + tagVesicleLength; j < i + tagVesicleLength + 2000; j++) {

                        if (s.substring(j, j + tagPointLength).equalsIgnoreCase(tagPoint)) {

                            k1 = j + tagPointLength + 4;

                            for (int k = k1; k < ilimit; k++) {

                                if (s.substring(k, k+1).equalsIgnoreCase("]")) {
                                    finished = true;
                                    break;
                                }

                                if (s.substring(k, k + 1).equalsIgnoreCase(",")) {

                                    k2 = k + 1;
                                    ss = s.substring(k1, k2);
                                    ss = ss.replaceAll(" ", ",");

                                    xvalue = Double.parseDouble(D3DUtility.itemFromList(ss, 0));
                                    yvalue = Double.parseDouble(D3DUtility.itemFromList(ss, 1));
                                    zvalue = Double.parseDouble(D3DUtility.itemFromList(ss, 2));

                                    xsum += xvalue;
                                    ysum += yvalue;
                                    zsum += zvalue;
                                    ipoints++;

                                    //Master.log("x=" + xvalue + ",y=" + yvalue + ",z=" + zvalue);
                                    //Master.log("" + xyz[ipoint][2]);

                                    k1 = k2 + 4;

                                }

                            }

                        }

                        if (finished) {

                            if (ipoints == 0) {
                                Master.log("failed to find vertices for vesicle #" + ivesicle);
                                return null;
                            } else {
                                break;
                            }

                        }

                    }

                    vxyz[ivesicle][0] = xsum / ipoints;
                    vxyz[ivesicle][1] = ysum / ipoints;
                    vxyz[ivesicle][2] = zsum / ipoints;

                    //Master.log("x=" + xsum + ",y=" + ysum + ",z=" + zsum + ", count=" + ipoints);

                    //if (true) {
                    //    return null;
                    //}

                    //Master.log("x=" + xyz[ivesicle][0] + ",y=" + xyz[ivesicle][1] + ",z=" + xyz[ivesicle][2]);



                }

            }

        } catch (IOException e) {
            System.err.println("unable to read file: " + WRL_fileName);
            return null;
        }

        Master.log("imported vesicles positions n = " + imax);

        return vxyz;

    }

    private void removeVesicleDuplicates(double[][] vesicleXYZ, double xyTolerance) {

        double d, vx1, vy1, vz1, vx2, vy2, vz2, minD = Double.POSITIVE_INFINITY;
        int ivz1, ivz2, count = 0;
        double vdx, vdy;

        for (int i = 0; i < vesicleXYZ.length; i++) {

            if (Double.isNaN(vesicleXYZ[i][0])) {
                continue;
            }

            vx1 = vesicleXYZ[i][0];
            vy1 = vesicleXYZ[i][1];
            vz1 = vesicleXYZ[i][2];

            for (int j = i + 1; j < vesicleXYZ.length; j++) {

                if (Double.isNaN(vesicleXYZ[i][0])) {
                    continue;
                }

                vx2 = vesicleXYZ[j][0];
                vy2 = vesicleXYZ[j][1];
                vz2 = vesicleXYZ[j][2];

                d = Math.sqrt((vx1 - vx2) * (vx1 - vx2) + (vy1 - vy2) * (vy1 - vy2) + (vz1 - vz2) * (vz1 - vz2));

                if (Double.isNaN(d)) {
                    continue;
                }

                if (d < minD) {
                    minD = d;
                }

                if (d < 0.04345) {

                    ivz1 = (int) (vz1 * 1000.0);
                    ivz2 = (int) (vz2 * 1000.0);

                    if (ivz2 == ivz1) {
                        continue; // in the same z-plane
                    }

                    if ((ivz2 > ivz1 + 30) || (ivz2 < ivz1 - 30)) {
                        continue;
                    }

                    vdx = Math.abs(vx2 - vx1);
                    vdy = Math.abs(vy2 - vy1);

                    if ((vdx > xyTolerance) || (vdy > xyTolerance)) {
                        continue;
                    }

                    vesicleXYZ[j][0] = Double.NaN;
                    vesicleXYZ[j][1] = Double.NaN;
                    vesicleXYZ[j][2] = Double.NaN;
                    count++;

                    // looks like the same vesicle in adjacent planes

                }

            }

        }

        Master.log("removed vesicle duplicates, n = " + count);

    }

    private double[][] removeNANs(double[][] vesicleXYZ) {

        int numVesicles = 0;

        for (double[] xyz : vesicleXYZ) {
            if (!Double.isNaN(xyz[0])) {
                numVesicles++;
            }
        }

        double[][] vxyz = new double[numVesicles][3];

        int ii = 0;

        for (double[] xyz : vesicleXYZ) {
            if (!Double.isNaN(xyz[0])) {
                vxyz[ii][0] = xyz[0];
                vxyz[ii][1] = xyz[1];
                vxyz[ii][2] = xyz[2];
                ii++;
            }
        }

        Master.log("created vesicle array n = " + vxyz.length);

        return vxyz;

    }
    
    private void writeActiveZoneArray(double[][][] azXYZ) {

        for (int i = 0; i < azXYZ.length; i++) {
            for (int j = 0; j < azXYZ[0].length; j++) {
                if (!Double.isNaN(azXYZ[i][j][0])) {
                    Master.log("xyz[" + i + "][" + j + "][0] = " + azXYZ[i][j][0] + ";");
                    Master.log("xyz[" + i + "][" + j + "][1] = " + azXYZ[i][j][1] + ";");
                    Master.log("xyz[" + i + "][" + j + "][2] = " + azXYZ[i][j][2] + ";");
                }
            }
        }

    }
    
    private void writeVesicleArray(double[][] vesicleXYZ) {

        for (int i = 0 ; i < vesicleXYZ.length; i++) {
            Master.log("xyz[" + i + "][0] = " + vesicleXYZ[i][0] + ";");
            Master.log("xyz[" + i + "][1] = " + vesicleXYZ[i][1] + ";");
            Master.log("xyz[" + i + "][2] = " + vesicleXYZ[i][2] + ";");
        }
        
    }

    private void configureGeometry(RunMonteCarloAZEM mc, double[][][] azXYZ, double[][] vxyz, double dx) {

        double xwidth, ywidth, zwidth;
        double x0, y0, z0;

        double xmax = Double.NEGATIVE_INFINITY, xmin = Double.POSITIVE_INFINITY;
        double ymax = Double.NEGATIVE_INFINITY, ymin = Double.POSITIVE_INFINITY;
        double zmax = Double.NEGATIVE_INFINITY, zmin = Double.POSITIVE_INFINITY;

        for (double[] xyz : vxyz) {
            xmax = Math.max(xmax, xyz[0]);
            xmin = Math.min(xmin, xyz[0]);
            ymax = Math.max(ymax, xyz[1]);
            ymin = Math.min(ymin, xyz[1]);
            zmax = Math.max(zmax, xyz[2]);
            zmin = Math.min(zmin, xyz[2]);
        }

        xwidth = xmax - xmin;
        ywidth = ymax - ymin;
        zwidth = zmax - zmin;

        x0 = (xmin + xmax) / 2.0;
        y0 = (ymin + ymax) / 2.0;
        z0 = (zmin + zmax) / 2.0;

        z0 += 0.015/2.0;

        // translate vesicles to center of geometry

        Master.log("transaling x, y, z..." + x0 + "," + y0 + "," + z0);

        for (double[] xyz : vxyz) {
            xyz[0] -= x0;
            xyz[1] -= y0;
            xyz[2] -= z0;
        }

        mc.azx1 = new double[azXYZ.length];
        mc.azx2 = new double[azXYZ.length];
        mc.azy1 = new double[azXYZ.length];
        mc.azy2 = new double[azXYZ.length];
        mc.azz1 = new double[azXYZ.length];
        mc.azz2 = new double[azXYZ.length];

        for (int iAZ = 0; iAZ < azXYZ.length; iAZ++) {

            mc.azx1[iAZ] = Double.POSITIVE_INFINITY;
            mc.azx2[iAZ] = Double.NEGATIVE_INFINITY;
            mc.azy1[iAZ] = Double.POSITIVE_INFINITY;
            mc.azy2[iAZ] = Double.NEGATIVE_INFINITY;
            mc.azz1[iAZ] = Double.POSITIVE_INFINITY;
            mc.azz2[iAZ] = Double.NEGATIVE_INFINITY;

            // compute AZ coordinates

            for (int ipnt = 0; ipnt < azXYZ[0].length; ipnt++) {

                if (Double.isNaN(azXYZ[iAZ][ipnt][0])) {
                    continue;
                }

                //Master.log("" + azXYZ[iAZ][ipnt][2]);

                azXYZ[iAZ][ipnt][0] -= x0;
                azXYZ[iAZ][ipnt][1] -= y0;
                azXYZ[iAZ][ipnt][2] -= z0;

                mc.azx1[iAZ] = Math.min(mc.azx1[iAZ], azXYZ[iAZ][ipnt][0]);
                mc.azx2[iAZ] = Math.max(mc.azx2[iAZ], azXYZ[iAZ][ipnt][0]);

                mc.azy1[iAZ] = Math.min(mc.azy1[iAZ], azXYZ[iAZ][ipnt][1]);
                mc.azy2[iAZ] = Math.max(mc.azy2[iAZ], azXYZ[iAZ][ipnt][1]);

                mc.azz1[iAZ] = Math.min(mc.azz1[iAZ], azXYZ[iAZ][ipnt][2]);
                mc.azz2[iAZ] = Math.max(mc.azz2[iAZ], azXYZ[iAZ][ipnt][2]);

            }

        }

        double extraForLongerSims = 2; // um

        double reserveWidthX = 1.0 * extraForLongerSims; // um
        double reserveWidthY = 1.0 * extraForLongerSims; // um
        double reserveWidthZ = 4 * dx;

        switch(mc.EMseries) {
            case 1:
                if (mc.openGeometry) {
                    extraForLongerSims = 0.96;
                    mc.openGeometryForwardX = -1;
                    mc.openGeometryForwardY = -1;
                    reserveWidthZ = 1.2 * extraForLongerSims;
                } else {
                    extraForLongerSims = 2.0;
                }
                reserveWidthX = 1.2 * extraForLongerSims;
                reserveWidthY = 1.2 * extraForLongerSims;
                mc.xStep_T = 0.5;
                break;
            case 3:
                if (mc.openGeometry) {
                    extraForLongerSims = 1.1;
                    mc.openGeometryForwardX = 1;
                    mc.openGeometryForwardY = -1;
                    reserveWidthZ = 1.0 * extraForLongerSims; // um
                } else {
                    extraForLongerSims = 2.1;
                }
                reserveWidthX = 1.0 * extraForLongerSims; // um
                reserveWidthY = 1.0 * extraForLongerSims; // um
                mc.xStep_T = 0.5;
                break;
            case 4:
                mc.reserveY2 = false;
                switch(mc.EMseriesAZ) {
                    case 0:
                        if (mc.openGeometry) {
                            extraForLongerSims = 1.5;
                            mc.openGeometryForwardY = 1;
                        } else {
                            extraForLongerSims = 1.5;
                        }
                        reserveWidthX = 0.7 * extraForLongerSims;
                        reserveWidthY = 0.7 * extraForLongerSims;
                        reserveWidthZ = 0.7 * extraForLongerSims;
                        mc.reserveZ1 = true;
                        mc.xStep_T = 0.4;
                        break;
                    case 1:
                        if (mc.openGeometry) {
                            extraForLongerSims = 1.1;
                            mc.openGeometryForwardX = -1;
                            mc.openGeometryForwardY = -1;
                            reserveWidthZ = 1.0 * extraForLongerSims; // um
                        } else {
                            extraForLongerSims = 2.0;
                        }
                        reserveWidthX = 1.0 * extraForLongerSims;
                        reserveWidthY = 1.0 * extraForLongerSims;
                        mc.xStep_T = 0.8;
                        break;
                    case 2:
                        if (mc.openGeometry) {
                            extraForLongerSims = 1.5;
                            mc.openGeometryForwardX = -1;
                            mc.openGeometryForwardY = 1;
                        } else {
                            extraForLongerSims = 1.6;
                        }
                        reserveWidthX = 0.7 * extraForLongerSims;
                        reserveWidthY = 0.7 * extraForLongerSims;
                        reserveWidthZ = 0.7 * extraForLongerSims;
                        mc.reserveX1 = false;
                        mc.reserveZ1 = true;
                        mc.xStep_T = 0.5;
                        break;
                    case 3:
                        if (mc.openGeometry) {
                            extraForLongerSims = 1.5;
                            mc.openGeometryForwardY = 1;
                        } else {
                            extraForLongerSims = 1.5;
                        }
                        reserveWidthX = 0.7 * extraForLongerSims;
                        reserveWidthY = 0.7 * extraForLongerSims;
                        reserveWidthZ = 0.7 * extraForLongerSims;
                        mc.reserveZ1 = true;
                        mc.xStep_T = 0;
                        break;
                }
                break;
            case 5:
                if (mc.openGeometry) {
                    extraForLongerSims = 1.0;
                    //mc.openGeometryForwardX = 1;
                    mc.openGeometryForwardY = 1;
                    reserveWidthZ = 1.1 * extraForLongerSims;
                } else {
                    extraForLongerSims = 2.0;
                }
                reserveWidthX = 1.1 * extraForLongerSims;
                reserveWidthY = 1.1 * extraForLongerSims;
                mc.reserveY2 = false;
                mc.xStep_T = 0;
                break;
            case 10:
                extraForLongerSims = 1.7;
                reserveWidthX = 4 * dx;
                reserveWidthY = 0.9 * extraForLongerSims;
                reserveWidthZ = 0.8 * extraForLongerSims;
                mc.reserveX1 = false;
                mc.reserveX2 = false;
                mc.reserveY1 = false;
                mc.reserveZ2 = true;
                switch(mc.EMseriesAZ) {
                    case 0:
                        //reserveWidthX = 0.9 * extraForLongerSims;
                        //mc.reserveX1 = true;
                        //mc.reserveX2 = false;
                        if (mc.openGeometry) {
                            extraForLongerSims = 0.6;
                            reserveWidthX = 0.9 * extraForLongerSims;
                            mc.openGeometryForwardX = -1;
                            mc.openGeometryForwardZ = 1;
                        }
                        mc.reserveY2 = false;
                        reserveWidthY = 4 * dx;
                        reserveWidthZ = 2 * extraForLongerSims;
                        mc.xStep_T = 1;
                        break;
                    case 1:
                        if (mc.openGeometry) {
                            extraForLongerSims = 0.8;
                            reserveWidthX = 0.9 * extraForLongerSims;
                            mc.openGeometryForwardY = -1;
                        }
                        reserveWidthY = 0.9 * extraForLongerSims;
                        reserveWidthZ = 0.8 * extraForLongerSims;
                        mc.xStep_T = 0;
                        break;
                    case 2:
                        if (mc.openGeometry) {
                            extraForLongerSims = 0.8;
                            reserveWidthX = 0.9 * extraForLongerSims;
                            mc.openGeometryForwardX = 1;
                            mc.openGeometryForwardZ = -1;
                        }
                        reserveWidthY = 0.9 * extraForLongerSims;
                        reserveWidthZ = 0.8 * extraForLongerSims;
                        mc.xStep_T = 1;
                        break;
                    case 3:
                        if (mc.openGeometry) {
                            extraForLongerSims = 0.8;
                            reserveWidthX = 0.9 * extraForLongerSims;
                            mc.openGeometryForwardX = -1;
                            mc.openGeometryForwardY = -1;
                        }
                        reserveWidthY = 0.9 * extraForLongerSims;
                        reserveWidthZ = 0.8 * extraForLongerSims;
                        mc.xStep_T = 0.5;
                        break;
                    case 4:
                        mc.reserveZ1 = true;
                        mc.reserveZ2 = false;
                        if (mc.openGeometry) {
                            extraForLongerSims = 0.8;
                            reserveWidthX = 0.9 * extraForLongerSims;
                            mc.openGeometryForwardX = -1;
                            mc.openGeometryForwardY = -1;
                        }
                        reserveWidthY = 0.9 * extraForLongerSims;
                        reserveWidthZ = 0.8 * extraForLongerSims;
                        mc.xStep_T = 1;
                        break;
                }
                break;
            case 16:
                if (mc.openGeometry) {
                    extraForLongerSims = 1.2;
                    reserveWidthZ = 0.8 * extraForLongerSims;
                    mc.openGeometryForwardY = 1;
                } else {
                    extraForLongerSims = 2.0;
                }
                reserveWidthX = 0.8 * extraForLongerSims;
                reserveWidthY = 0.8 * extraForLongerSims;
                mc.xStep_T = 0.2;
                break;
            case 18:
                if (mc.openGeometry) {
                    extraForLongerSims = 1.0;
                    reserveWidthZ = 1.3 * extraForLongerSims;
                    mc.openGeometryForwardY = 1;
                } else {
                    extraForLongerSims = 1.9;
                }
                reserveWidthX = 1.3 * extraForLongerSims;
                reserveWidthY = 1.3 * extraForLongerSims;
                mc.reserveY1 = false;
                mc.xStep_T = 0.4;
                break;
        }

        if(mc.xStep_T < 0) {
            Master.exit("bad xStep_T value");
        }

        if (mc.openGeometry) {
            mc.reserveX1 = true;
            mc.reserveX2 = true;
            mc.reserveY1 = true;
            mc.reserveY2 = true;
            mc.reserveZ1 = true;
            mc.reserveZ2 = true;
        }

        if (mc.extraReserve) {
            igm.xdim = xwidth + reserveWidthX;
            igm.ydim = ywidth + reserveWidthY;
            igm.zdim = zwidth + reserveWidthZ;
            Master.log("reserveWidthX = " + reserveWidthX);
            Master.log("reserveWidthY = " + reserveWidthY);
            Master.log("reserveWidthZ = " + reserveWidthZ);
        } else {
            igm.xdim = xwidth + 4 * dx;
            igm.ydim = ywidth + 4 * dx;
            igm.zdim = zwidth + 4 * dx;
        }

    }

    private boolean initProject(double dx) {

        if (!Master.foundMainStartUpArguments) {
           project.simTime = simTime;
           project.directory = directory;
           project.folder = folder;
        }

        project.dx = dx;
        project.set("saveRate", saveRate);
        project.set("printRate", printRate);

        return false;

    }

    private DiffusantVesiclesAZ initVesicles(double D) {

        DiffusantVesiclesAZ dvs;

        String d1_name = "Vesicles";

        dvs = new DiffusantVesiclesAZ(project, d1_name, 0, D, null, vesicleRadius, null);

        dvs.save.save2TextFile = false;

        project.addDiffusant(dvs);

        return dvs;

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
