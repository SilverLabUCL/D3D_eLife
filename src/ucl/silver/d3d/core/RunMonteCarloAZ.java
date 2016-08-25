package ucl.silver.d3d.core;

/**
 * <p>Title: D3D</p>
 * <p>Description: 3D Diffusion-Reaction Simulator</p>
 * <p>Copyright: Copyright (c) 2016</p>
 * <p>Company: The Silver Lab at University College London</p>
 * @author Jason S. Rothman
 * @version 1.6
 */

public class RunMonteCarloAZ
        extends RunMonteCarlo {

    public Coordinates activeZone = null; // active zone (az) coordinates

    public double[][][] azXYZ = null;

    public double azWidth = 0.16; // active zone xy-width
    //public double azHeight = 0.025; // active zone z-height
    public double azHeight = 0.05; // active zone z-height

    public double azx1[] = null;
    public double azx2[] = null;
    public double azy1[] = null;
    public double azy2[] = null;
    public double azz1[] = null;
    public double azz2[] = null;
    public String azPlane = "xy";

    public double azKeepClearExtra = 0; // 0.05; // um (keeps AZ clear of immobile vesicles)
    public boolean azExcludeReadyInit = true; // keep AZ clear when initiating ready vesicles
    
    public boolean dockingOn = false;
    public boolean dockingOffAfterInit = false;
    public double dockingOnTime = 0; // ms

    public int dockedVesicles = 0; // number of docked vesicles at active zone
    public int dockedVesiclesMax = -1; // use -1 for no limit
    private DiffusantVesicleAZ firstDocked = null;

    public int reserveVesicles = 0; // number of reserve vesicles
    public int reserveVesiclesMax = -1; // use -1 for no limit
    private DiffusantVesicleAZ firstReserve = null;

    public double releaseRate = 0; // release docked, move to reserve (kHz)
    public double releaseStartTime = 0; // ms, min time to allow release
    public int releaseCounter = 0;
    public int releaseLimit = 0; // use value > 0 to end simulation
    public double releaseTimeOutLimit = 0; // abort simulation if no releases occur before this time
    public double releaseProb = 1.0; // release probability

    public double replenishRate = 0; // move reserve to ready (kHz)
    public double replenishStartTime = 0; // ms, min time to allow replacement
    private int replenishCounter = 0;

    public Coordinates replenishCoordinates = null; // where reserve vesicles are placed
    public Coordinates replenishCoordinatesExclude = null;
    public double replenishFromAZdistance = -1; // number of voxels from AZ for setting replenish coordinates

    public boolean tetherVesicles = false;
    
    public double azPermitConnectorRadius = 0.2; // um
    public int azNumTethers = 1;
    public double azTetherLength = 0.008; // um
    public Coordinates activeZoneTethering = null;
    public DiffusantVesicle azDV = null;

    public boolean saveRelease = false;

    private long iRelease, iReleaseMax;
    private long iReplenish, iReplenishMax;

    Save save_Release = null; // saving data to file and/or internal array

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("azWidth")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("azHeight")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("azKeepClearExtra")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("dockingOnTime")) {
            return project.timeUnits;
        }
        if (name.equalsIgnoreCase("releaseRate")) {
            return "1/" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("releaseStartTime")) {
            return project.timeUnits;
        }
        if (name.equalsIgnoreCase("releaseTimeOutLimit")) {
            return project.timeUnits;
        }
        if (name.equalsIgnoreCase("replenishRate")) {
            return "1/" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("replenishStartTime")) {
            return project.timeUnits;
        }
        if (name.equalsIgnoreCase("azPermitConnectorRadius")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("azTetherLength")) {
            return project.spaceUnits;
        }     
        return super.units(name);
    }

    @Override
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("dockedVesicles")) {
            return false;
        }
        if (name.equalsIgnoreCase("tetheredVesicles")) {
            return false;
        }
        if (name.equalsIgnoreCase("reserveVesicles")) {
            return false;
        }
        return super.canEdit(name);
    }

    public RunMonteCarloAZ(Project p) {

        super(p);

        save_Release = new Save(p);
        save_Release.name = "Release";
        
        activeZone = new Coordinates(project, 0, 0, 0, 0, 0, 0);
        activeZone.name = "Active Zone";
        activeZone.setShape("cylinderz");
        activeZone.color.name = "AZcolor";
        activeZone.color.setColor("[r=204,g=0,b=0]");

        createVector(true);

    }

    @Override
    public void init() {

        super.init();

        if (saveRelease) {
            save_Release.init();
            save_Release.save2TextFile = true;
            save_Release.fileName("Release", "");
            save_Release.ydim = project.timeUnits;
            save_Release.saveWhileComputing = true;
            save_Release.saveWhileComputingAppend = true;
        }

        if (PBC) {
            //Master.exit("RunMonteCarloAZ: init: PBC not allowed in active zone simulations");
        }

    }

    @Override
    public void setOutputRate(double newRate) {

        super.setOutputRate(newRate);

        if (save_Release != null) {
            save_Release.setOutputRate(newRate);
        }

    }

    @Override
    public boolean initSimulation() {

        if (checkVariables()) {
            return true;
        }

        if (autoInit) {
            if (initAll()) {
                return true;
            }
        }

        iRelease = 0;
        iReplenish = 0;

        if (releaseRate > 0) {
            iReleaseMax = (int) ((1.0 / releaseRate) / project.dt);
            iReleaseMax = Math.max(iReleaseMax, 1); // at least one time rstep
            //Master.log("iReleaseMax: " + iReleaseMax);
            iRelease = iReleaseMax; // start at time = 0
        }

        if (replenishRate > 0) {
            iReplenishMax = (int) (1.0 / (replenishRate * project.dt));
            iReplenishMax = Math.max(iReplenishMax, 1); // at least one time rstep
            replenishCounter = 0;
        }

        time = 0;

        timer1.start();

        //saveXYZincrement.impulse = true;
        //saveXYZincrement.name = "saveXYZincrement";
        //saveXYZincrement.initTimer();

        initialized = true;

        printParameters();

        Master.log("initialized Monte Carlo simulation");

        return false;

    }

    @Override
    public void finishSimulation(){

        finishSave();

        if (hydrodynamicsLocalD) {
            localDensityAll(true);
        }

        vesicleStats();

        if (connectVesicles) {
            avgConnectorLifeTime /= connectorLifeTimeCounter;
            Master.log("average connector life time (ms): " + avgConnectorLifeTime + " (n=" + Integer.toString(connectorLifeTimeCounter) + ")");
        }

        Master.log("final docked vesicles: " + dockedVesicles);
        Master.log("final reserve vesicles: " + reserveVesicles);
        Master.log("final vesicle volume fraction: " + vesicleVolumeFraction);

        if (replenishRate > 0) {
            Master.log("replenished vesicles: " + replenishCounter);
        }

    }

    @Override
    public double spaceVolume() {
        return project.geometry.spaceVolume - mitoVolume - activeZone.geometryVolume();
    }

    @Override
    public void run() {

        int icount, numReturned;

        DiffusantVesicleAZ dv;

        while (runSimulation) { // runSimulation thru batches

            if (!initialized) {

                if (project.simulationInit(preview)) {
                    cancelSimulation(); // ERROR
                }

                saveVesiclePositions(0);
                saveVesicleDensity("start");

            }

            if (time < project.simTime) {
                Master.log("starting Monte Carlo simulation at time " + time + " " + project.timeUnits);
            }

            while (runSimulation && !cancel && (time < project.simTime)) {

                for (int s = 0; s < diffusant.length; s++) {
                    diffusant[s].save(); // save diffusant variables
                }

                runMSDspatial();

                moveVesiclesActiveZone();

                if (hydrodynamicsLocalD) {
                    localDensityAll(false);
                }

                if (connectVesicles) {
                    connectVesicles();
                    unconnectVesicles();
                }

                if ((releaseRate > 0) && (iRelease == iReleaseMax)) {

                    //if (releaseRate == Double.POSITIVE_INFINITY) {
                    //    icount = 9999;
                    //} else {
                    //    icount = dockedVesiclesMax;
                    //}

                    icount = dockedVesicles;

                    if (time < releaseStartTime) {
                        icount = 0;
                    }

                    for (int i = 0; i < icount; i++) {
                        
                        dv = lastDockedVesicle();
                        
                        if (dv == null) {
                            break;
                        }

                        releaseDockedVesicle(dv, time);

                    }

                    iRelease = 0; // reset counter

                }

                if ((replenishRate > 0) && (iReplenish == iReplenishMax)) {

                    if (replenishRate == Double.POSITIVE_INFINITY) {
                        icount = 999;
                    } else {
                        icount = 1;
                    }

                    if (time < replenishStartTime) {
                        icount = 0;
                    }

                    numReturned = moveReserveToReady(icount);

                    //if (numReturned > 0) {
                    //    Master.log("added reserve vesicles: " + numReturned);
                    //}

                    replenishCounter += numReturned;

                    iReplenish = 0; // reset counter

                }

                if ((grid != null) && preview) {
                    grid.repaint();
                }

                iRelease++;
                iReplenish++;

                itime += 1;
                time += project.dt;

                if (timer) {
                    timer2.timer(time);
                }

                if ((releaseLimit > 0) && (releaseCounter >= releaseLimit)) {
                    cancel = true;
                }

                if ((releaseTimeOutLimit > 0) && (time >= releaseTimeOutLimit) && (releaseCounter == 0)) {
                    Master.log("cancelled MC AZ simulation: no release events within " + releaseTimeOutLimit + " ms");
                    cancel = true;
                }

            }

            if (cancel || (time >= project.simTime)) {
                finish();
            }

        }

    }

    @Override
    public boolean initAll(){

        releaseCounter = 0;

        if (checkVariables()) {
            return true;
        }

        if (initDiffusantVesiclesArray()) {
            return true;
        }

        if (initVoxels()) {
            return true;
        }

        if (initActiveZone()) {
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

        if (initVesiclesDocked()){
            return true;
        }

        if (initVesicleRandom && initVesiclesRandom()) {
            return true;
        }

        if (initDockedList()) {
            return true;
        }

        if (initReserveList()) {
            return true;
        }

        if (removeVesicleOverlap && removeVesicleOverlap()) {
            return true;
        }

        if (initVesiclesImmobile()) {
            return true;
        }

        initConnectors();

        vesicleStats();

        replenishCoordinates = new Coordinates(project, geometry);

        if (!freeDiffusion && (vesicleVolumeFraction > 0.44)) {
            replenishRate = 0; // TAKES TOO LONG
            Master.log("TURNED OFF REPLENISHING: replenishRate = 0");
        }

        if (hydrodynamicsLocalD) {
            localDensityAll(true);
        }

        if (dockingOffAfterInit) {
            dockingOn = false; // TEMPORARY removal of docking
        }

        return false;

    }

    public boolean initActiveZone() {

        double radius, height, SA, diameter;
        double x0, x1, y0, y1, z0, z1;

        if ((geometry == null) || (activeZone == null)) {
            return true;
        }

        radius = azWidth / 2.0;
        height = azHeight / 2.0;

        if (azPlane.equalsIgnoreCase("xy")) {
            x0 = 0 - radius;
            y0 = 0 - radius;
            z0 = geometry.z2 - height;
            x1 = 0 + radius;
            y1 = 0 + radius;
            z1 = geometry.z2 + height;
            activeZone.setShape("cylinderz");
        } else if (azPlane.equalsIgnoreCase("yz")) {
            x0 = geometry.x2 - height;
            y0 = 0 - radius;
            z0 = 0 - radius;
            x1 = geometry.x2 + height;
            y1 = 0 + radius;
            z1 = 0 + radius;
            activeZone.setShape("cylinderx");
        } else if (azPlane.equalsIgnoreCase("zx")) {
            x0 = 0 - radius;
            y0 = geometry.y2 - height;
            z0 = 0 - radius;
            x1 = 0 + radius;
            y1 = geometry.y2 + height;
            z1 = 0 + radius;
            activeZone.setShape("cylindery");
        } else {
            return true;
        }

        activeZone.setCoordinates(x0, y0, z0, x1, y1, z1);

        if (tetherVesicles) {

            radius += azTetherLength;
            height += azTetherLength;

            activeZoneTethering = new Coordinates(project, 0, 0, 0, 0, 0, 0);
            activeZoneTethering.setCoordinates(x0, y0, z0, x1, y1, z1);
            activeZoneTethering.setShape(activeZone.shape);

            if (azNumTethers > 0) {
                azDV = new DiffusantVesicle(project, "AZ", 0, 0, 0, 0, 0);
                azDV.connectTo = new DiffusantVesicle[azNumTethers];
                azDV.connectorOffTime = new double[azNumTethers];
            }

        }

        return false;

    }

    @Override
    public boolean initVesiclesImmobile() {

        double xx1, xx2, yy1, yy2, zz1, zz2;

        Coordinates keepClear = null;

        if (diffusant == null) {
            return true;
        }

        if (azKeepClearExtra >= 0) {

            xx1 = activeZone.x1 - 1 * azKeepClearExtra;
            xx2 = activeZone.x2 + 1 * azKeepClearExtra;
            yy1 = activeZone.y1 - 1 * azKeepClearExtra;
            yy2 = activeZone.y2 + 1 * azKeepClearExtra;
            zz1 = activeZone.z1 - 1 * azKeepClearExtra;
            zz2 = activeZone.z2 + 1 * azKeepClearExtra;

            keepClear = new Coordinates(project, xx1, yy1, zz1, xx2, yy2, zz2);

        }

        for (int i = 0; i < diffusant.length; i++) {

            if (diffusant[i] == null) {
                continue;
            }

            //Master.log(keepClear.printDimensions());

            if (diffusant[i].initImmobileVesicles(keepClear)){
                return true;
            }

        }

        return false;

    }

    public boolean initVesicleRandom(DiffusantVesicle dv, Coordinates c, Coordinates exclude, boolean allowOverlaps, boolean excludeAZ, double distanceFromAZ, int overlapTrialLimit, int absLimit, double az_z) {

        int trial = 0;
        double x, y, z;
        double dx, dy, dz, d2AZ;
        boolean overlap;

        if (dv == null) {
            return true;
        }

        if (c == null) {
            return true;
        }

        while (true) {

            trial++;

            if (trial >= absLimit) {
                Master.log("warning: initVesicleRandom: reached trial limit");
                return true;
            }

            x = (mt.nextDouble() * (c.x2 - c.x1)) + c.x1;
            y = (mt.nextDouble() * (c.y2 - c.y1)) + c.y1;

            if (Double.isNaN(az_z)) {
                z = (mt.nextDouble() * (c.z2 - c.z1)) + c.z1;
            } else {
                z = az_z;
            }

            if (!c.isInside(x, y, z, dv.radius)) {
                continue;
            }

            if ((exclude != null) && exclude.isInside(x, y, z, dv.radius)) {
                continue;
            }

            if (distanceFromAZ > 0) {

                dx = x - activeZone.xCenter;
                dy = y - activeZone.yCenter;
                dz = z - activeZone.zCenter;

                d2AZ = Math.sqrt(dx * dx + dy * dy + dz * dz);

                if (d2AZ < distanceFromAZ) {
                    continue;
                }

            }

            if (excludeAZ && (activeZone != null)) {
                if (activeZone.isInside(x, y, z)) {
                    continue;
                }
            }

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
                overlap = testVesicleOverlap(dv, null) != null;
            }

            if (!overlap || (allowOverlaps && (trial > overlapTrialLimit))) {

                if (!addToVoxelList(dv)) {
                    return true;
                }

                return false;

            }

        }

    }

    @Override
    public boolean setVesicleLocation(DiffusantVesicle dv, double x, double y, double z, boolean initStartLocation) {

        int xVoxel, yVoxel, zVoxel;

        if (PBC) {

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

            dv.x = x;
            dv.y = y;
            dv.z = z;

            if (initStartLocation) {
                dv.x0 = x;
                dv.y0 = y;
                dv.z0 = z;
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

    @Override
    public boolean initVesiclesRandom() {

        boolean allowOverlaps = true;

        Coordinates c;
        DiffusantVesicleAZ dv;

        if (diffusant == null) {
            return true;
        }

        Master.log("randomly placing ready vesicles...");

        for (int i = 0; i < diffusant.length; i++) {

            if ((diffusant[i] == null) || (diffusant[i].vesicle == null)) {
                continue;
            }

            if (diffusant[i].coordinates == null) {
                c = new Coordinates(project, geometry);
            } else {
                c = new Coordinates(project, diffusant[i].coordinates);
            }

            for (int j = 0; j < diffusant[i].vesicle.length; j++) {

                if (!(diffusant[i].vesicle[j] instanceof DiffusantVesicleAZ)) {
                    continue;
                }

                dv = (DiffusantVesicleAZ) diffusant[i].vesicle[j];

                if (!dv.isReady) {
                    continue;
                }

                if (initVesicleRandom(dv, c, null, allowOverlaps, azExcludeReadyInit, -1, overlapTrialLimit, absLimit, Double.NaN)) {
                    return true;
                }

            }

        }

        return false;

    }

    @Override
    public boolean removeVesicleOverlap() {

        int radiusNM, lastRadiusNM = 0, count = 0;
        double halfDistance = 0;
        boolean saveDockingOn = dockingOn;

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
        dockingOn = false;

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

                //if ((diffusant[j] == null) || (diffusant[j].vesicle == null)) {
                //    continue;
                //}

                for (int j = 0; j < diffusant[i].vesicle.length; j++) {
                    //if (diffusant[i].vesicle[j].insideGeometry) {
                        //diffusant[j].vesicle[j].radius = Math.min(halfDistance, diffusant[j].meanRadius);
                        diffusant[i].vesicle[j].radius = halfDistance;
                    //}
                }

            }

            moveVesiclesActiveZone();

            count++;

            if (count > absLimit) {
                break;
            }

        }

        dockingOn = saveDockingOn;
        removingVesicleOverlap = false;

        if (testVesicleOverlap()) {

            error("removeVesicleOverlap: failed to remove vesicle overlap.");

            for (int i = 0; i < diffusant.length; i++) {

                if ((diffusant[i] == null) || (diffusant[i].vesicle == null)) {
                    continue;
                }

                for (int j = 0; j < diffusant[i].vesicle.length; j++) {
                    diffusant[i].vesicle[j].radius = diffusant[i].meanRadius;
                    diffusant[i].vesicle[j].x = diffusant[i].vesicle[j].x0;
                    diffusant[i].vesicle[j].y = diffusant[i].vesicle[j].y0;
                    diffusant[i].vesicle[j].z = diffusant[i].vesicle[j].z0;
                }

            }

            return true;

        }

        for (int i = 0; i < diffusant.length; i++) {

            if ((diffusant[i] == null) || (diffusant[i].vesicle == null)) {
                continue;
            }

            for (int j = 0; j < diffusant[i].vesicle.length - 1; j++) {
                diffusant[i].vesicle[j].radius = diffusant[i].meanRadius;
                diffusant[i].vesicle[j].x0 = diffusant[i].vesicle[j].x;
                diffusant[i].vesicle[j].y0 = diffusant[i].vesicle[j].y;
                diffusant[i].vesicle[j].z0 = diffusant[i].vesicle[j].z;
            }

        }

        return false;

    }

    public boolean initVesiclesDocked() {

        double az_z;
        double zwidth;
        double vesVolume;

        int count = 0;
        int initNum = 0;

        boolean allowOverlaps = true;
        boolean excludeAZ = false;

        DiffusantVesiclesAZ dvs;
        DiffusantVesicleAZ dv;

        if (!dockingOn) {
            return false;
        }

        if (diffusant == null) {
            return true;
        }

        zwidth = 0.5 * azHeight;

        //Master.log("randomly placing docked vesicles...");

        for (int i = 0; i < diffusant.length; i++) {

            if ((diffusant[i] == null) || (diffusant[i].vesicle == null)) {
                continue;
            }

            if (!(diffusant[i] instanceof DiffusantVesiclesAZ)) {
                continue;
            }

            dvs = (DiffusantVesiclesAZ) diffusant[i];

            if (dvs.dockedVesiclesInit > 0) {
                initNum = dvs.dockedVesiclesInit;
            } else if (dvs.dockedVesiclesInitFraction > 0) {
                vesVolume = 4.0 * Math.PI * Math.pow(dvs.setMeanRadius, 3) / 3.0;
                initNum = (int) (dvs.dockedVesiclesInitFraction * activeZone.geometryVolume() / vesVolume) ;
            } else {
                continue;
            }

            count = 0;

            for (int j = 0; j < diffusant[i].vesicle.length; j++) {

                if (!(diffusant[i].vesicle[j] instanceof DiffusantVesicleAZ)) {
                    continue;
                }

                dv = (DiffusantVesicleAZ) diffusant[i].vesicle[j];

                if (!dv.isReady) {
                    continue;
                }

                az_z = geometry.z2 - zwidth;

                if (initVesicleRandom(dv, activeZone, null, allowOverlaps, excludeAZ, -1, overlapTrialLimit, absLimit, az_z)) {
                    Master.log("failed to place docked vesicle #" + j);
                    return true;
                }

                dv.setType("docked");

                count++;

                if (count >= initNum) {
                    break; // finished
                }

            }

        }

        Master.log("placed docked vesicles: " + count);

        return false;

    }

    public boolean initDockedList() {

        DiffusantVesicleAZ dv;

        dockedVesicles = 0;

        firstDocked = null;

        if (diffusant == null) {
            return true;
        }

        for (int i = 0; i < diffusant.length; i++) {

            if ((diffusant[i] == null) || (diffusant[i].vesicle == null)) {
                continue;
            }

            for (int j = 0; j < diffusant[i].vesicle.length; j++) {

                if (diffusant[i].vesicle[j] == null) {
                    continue;
                }

                if (!(diffusant[i].vesicle[j] instanceof DiffusantVesicleAZ)) {
                    continue;
                }

                dv = (DiffusantVesicleAZ) diffusant[i].vesicle[j];

                if (dv.isDocked) {
                    addToDockedList(dv, -99999);
                }

            }

        }

        Master.log("initial docked vesicles: " + dockedVesicles);

        return false;

    }

    public boolean addToDockedList(DiffusantVesicleAZ dv, double time) {

        if ((dockedVesiclesMax >= 0) && (dockedVesicles >= dockedVesiclesMax)) {
            return false;
        }

        dv.setType("docked");
        dv.dockTime = time;
        dv.nextDocked = firstDocked;
        firstDocked = dv;
        dockedVesicles += 1;

        //System.out.println("docked vesicle " + dv.dockTime + ", type " + dv.name);

        return true;

    }

    public boolean removeFromDockedList(DiffusantVesicleAZ dv) {

        DiffusantVesicleAZ idv = firstDocked;
        DiffusantVesicleAZ idvnext;
        DiffusantVesicleAZ idvhold = null;

        boolean found = false;

        while (idv != null) {

            idvnext = idv.nextDocked;

            if (idv == dv) {
                if (idvhold == null) {
                    firstDocked = idvnext; // dv was first list
                } else {
                    idvhold.nextDocked = idvnext;
                }
                found = true;
            }

            idvhold = idv;
            idv = idvnext;

        }

        if (found) {
            dockedVesicles -= 1;
            dv.dockTime = -1;
            dv.name = "";
        }

        return found;

    }

    public DiffusantVesicleAZ lastDockedVesicle() {

        DiffusantVesicleAZ idv = firstDocked;

        while (idv != null) {

            if (idv.nextDocked == null) {
                return idv;
            }

            idv = idv.nextDocked;

        }

        return null;

    }

    public boolean releaseDockedVesicle(DiffusantVesicleAZ dv, double time) {

        double d;
        double dlimit;

        if (dv == null) {
            return false;
        }

        if (dv.dockRefractoryPeriod > 0) {
            if (time - dv.dockTime < dv.dockRefractoryPeriod) {
                return false; // still in refractory period
            }
        }

        if ((releaseProb < 1.0) && (mt.nextDouble() > releaseProb)) {
            return false;
        }

        d = dv.squareDisplacement();
        d = Math.sqrt(d);

        if (connectVesicles) {

            dv.unconnectAll();

            if (dv.isTetheredToAZ) {
                azDV.unconnectFrom(dv);
                dv.isTetheredToAZ = false;
            }

        }

        //dv.z = geometry.z2 + project.dx * 0.5; // move out of the way
        dv.x = geometry.x2 + project.dx * 1.0; // move out of the way
        dv.y = geometry.y1 + releaseCounter * project.dx;

        //Master.log("released docked vesicle " + releaseCounter + " at time " + time);

        removeFromVoxelList(dv, dv.voxel);
        removeFromDockedList(dv);
        addToReserveList(dv);
        releaseCounter++;

        //Toolkit.getDefaultToolkit().beep();

        //d2az = this.activeZone.distanceToCenter(d, d, d);

        
        //Master.log("released docked vesicle " + releaseCounter + ", docking time " + (time - dv.dockTime));
        //Master.log("" + d + ", " + dv.d2AZ);

        if (!preview && saveRelease) {
            if (dv.d2AZ == 0) {
                saveReleaseData(time, 0); // reserve vesicle
            } else {
                saveReleaseData(time, d);
            }
        }

        return true;

    }

    @Override
    public void connectVesicles() {

        int neighbors;
        double sqrDistance, d;
        boolean inside;

        DiffusantVesicleAZ dv, kvesicle;

        if (tetherVesicles && (azDV == null)) {
            Master.exit("azDV is null");
        }

        for (int i = 0; i < diffusant.length; i++) {

            if ((diffusant[i] == null) || (diffusant[i].vesicle == null)) {
                continue;
            }

            for (int j = 0; j < diffusant[i].vesicle.length; j++) {

                dv = (DiffusantVesicleAZ) diffusant[i].vesicle[j];

                if (dv.isReserve) {
                    continue;
                }

                if (tetherVesicles && !dv.isTetheredToAZ && !azDV.noMoreConnectors()) {

                    inside = isInsideActiveZoneTethering(dv);

                    if (inside) {
                        if (azDV.connectToNew(dv, mt, Double.POSITIVE_INFINITY, time)) {
                            dv.isTetheredToAZ = true;
                        } else {
                            Master.exit("tether to AZ failure ");
                        }
                    }

                }

                if (dv.noMoreConnectors()) {
                    continue;
                }

                sqrDistance = dv.squareDistanceFromActiveZone(activeZone);

                if (sqrDistance >= azPermitConnectorRadius * azPermitConnectorRadius) {
                    continue;
                }

                neighbors = dv.voxel.numNeighbors;

                for (int k = 0; k < neighbors; k++) {

                    kvesicle = (DiffusantVesicleAZ) dv.voxel.neighbors[k].firstReady;

                    while (kvesicle != null) {

                        if (dv.noMoreConnectors()) {
                            break;
                        }

                        if (kvesicle.isReady && (kvesicle != dv) && !dv.isConnectedTo(kvesicle) && !kvesicle.noMoreConnectors()) {

                            if (dv.overlap(kvesicle, connectorLength)) {

                                sqrDistance = kvesicle.squareDistanceFromActiveZone(activeZone);

                                if (sqrDistance < azPermitConnectorRadius * azPermitConnectorRadius) {

                                    if (mt.nextDouble() < connectRate * project.dt) {

                                        if (kvesicle.connectToNew(dv, mt, unconnectRate, time)) {

                                            avgConnectorLifeTime += dv.connectorLifeTime;
                                            connectorLifeTimeCounter += 1;

                                            if (!dv.connectTo(kvesicle)) {
                                                Master.exit("connection failure");
                                            } else {
                                                //Master.log("connected vesicles");
                                            }

                                        } else {
                                            Master.exit("connection failure");
                                        }

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

    @Override
    public double[] computeNumberOfConnectors() {

        double sqrDistance, avg = 0, connected = 0;
        int n, ntotal = 0;

        double[] results = new double[2];

        DiffusantVesicleAZ dv;

        for (int i = 0; i < diffusant.length; i++) {

            if ((diffusant[i] == null) || (diffusant[i].vesicle == null)) {
                continue;
            }

            for (int j = 0; j < diffusant[i].vesicle.length; j++) {

                dv = (DiffusantVesicleAZ) diffusant[i].vesicle[j];

                if (dv.isReserve) {
                    continue;
                }

                sqrDistance = dv.squareDistanceFromActiveZone(activeZone);

                if (sqrDistance >= azPermitConnectorRadius * azPermitConnectorRadius) {
                    continue;
                }

                n = dv.numberOfConnections(true);
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

    public double minDistanceToAZ(double x, double y, double z) {

        double dx2, dy2;

        double AZx0 = 0;
        double AZy0 = 0;
        double AZradius = azWidth * 0.5;

        double dx = AZx0 - x;
        double dy = AZy0 - y;
        double dz = geometry.z2 - z;

        double d = Math.sqrt(dx * dx + dy * dy + dz * dz);

        double dxy = Math.sqrt(dx * dx + dy * dy);

        if (dxy < AZradius) {

            d = dz;

        } else {

            dx2 = AZradius * dx / dxy;
            dy2 = AZradius * dy / dxy;
            dx2 = dx - dx2;
            dy2 = dy - dy2;

            d = Math.sqrt(dx2 * dx2 + dy2 * dy2 + dz * dz);

        }

        return d;

    }

    public boolean initReserveList() {

        DiffusantVesicleAZ dv;

        reserveVesicles = 0;

        firstReserve = null;

        if (diffusant == null) {
            return true;
        }

        for (int i = 0; i < diffusant.length; i++) {

            if ((diffusant[i] == null) || (diffusant[i].vesicle == null)) {
                continue;
            }

            for (int j = 0; j < diffusant[i].vesicle.length; ++j) {

                if (diffusant[i].vesicle[j] == null) {
                    continue;
                }

                if (!(diffusant[i].vesicle[j] instanceof DiffusantVesicleAZ)) {
                    continue;
                }

                dv = (DiffusantVesicleAZ) diffusant[i].vesicle[j];

                if (dv.isReserve) {
                    addToReserveList(dv);
                }

            }

        }

        Master.log("initial reserve vesicles: " + reserveVesicles);

        return false;

    }

    public boolean addToReserveList(DiffusantVesicleAZ dv) {

        if ((reserveVesiclesMax >= 0) && (reserveVesicles >= reserveVesiclesMax)) {
            return false;
        }

        dv.setType("reserve");
        dv.nextReserve = firstReserve;
        firstReserve = dv;
        reserveVesicles += 1;

        return true;

    }

    public boolean removeFromReserveList(DiffusantVesicleAZ dv) {

        DiffusantVesicleAZ idv = firstReserve;
        DiffusantVesicleAZ idvnext;
        DiffusantVesicleAZ idvhold = null;

        boolean found = false;

        while (idv != null) {

            idvnext = idv.nextReserve;

            if (idv == dv) {
                if (idvhold == null) {
                    firstReserve = idvnext; // dv was first list
                } else {
                    idvhold.nextReserve = idvnext;
                }
                found = true;
            }

            idvhold = idv;
            idv = idvnext;

        }

        if (found) {
            reserveVesicles -= 1;
        }

        if (dv.isReserve){
            dv.name = "";
        }

        return found;

    }

     private int moveReserveToReady(int nvesicles) {

        int count = 0;
        double x, y, z, dx, dy, dz;
        double d2AZ;
        boolean error;

        boolean allowOverlaps = false;
        boolean excludeAZ = true;

        DiffusantVesicleAZ dv;

        if (nvesicles <= 0) {
            return 0;
        }

        for (int i = 0; i < nvesicles; i++) {

            dv = firstReserve;

            if (dv == null) {
                return count;
            }

            x = dv.x;
            y = dv.y;
            z = dv.z;

            error = initVesicleRandom(dv, replenishCoordinates, replenishCoordinatesExclude, allowOverlaps, excludeAZ, replenishFromAZdistance, overlapTrialLimit, absLimit, Double.NaN);
            
            if (error) {

                dv.x = x;
                dv.y = y;
                dv.z = z;

                return count;

            } else {

                dv.setType("ready");
                removeFromReserveList(dv);
                dv.isTetheredToAZ = false;
                dv.isDocked = false;
                dv.isReserve = false;
                
                count++;

            }

        }

        return count;

    }

    public boolean setVesicleLocation(DiffusantVesicleAZ dv, double x, double y, double z, boolean initStartLocation) {
        if (super.setVesicleLocation(dv, x, y, z, initStartLocation)) {
            return true;
        } else {
            return false;
        }
    }

    public boolean moveVesicleGauss(DiffusantVesicleAZ dv) {

        double step3;

        if (hydrodynamicsLocalD) {
            step3 = dv.localStep3;
        } else {
            step3 = dv.step3;
        }

        if (removingVesicleOverlap) {
            step3 = removeVesicleOverlapStep3;
        }

        stepx = ranGauss() * step3;
        stepy = ranGauss() * step3;
        stepz = ranGauss() * step3;

        if (dv.isDocked) {
            stepz = 0;
        }

        return setVesicleLocation(dv, dv.x + stepx, dv.y + stepy, dv.z + stepz, false);

    }

    public boolean moveVesicleGaussHydroWallz(DiffusantVesicleAZ dv) {

        double step3;
        double dz, b_ll, b_T;

        if (hydrodynamicsLocalD) {
            step3 = dv.localStep3;
        } else {
            step3 = dv.step3;
        }

        if (removingVesicleOverlap) {

            step3 = removeVesicleOverlapStep3;

            stepx = ranGauss() * step3;
            stepy = ranGauss() * step3;
            stepz = ranGauss() * step3;

        } else {

            dz = geometry.z2 - dv.z;

            b_ll = hydroWall_ll(dv.radius, dz);
            b_T = hydroWall_T(dv.radius, Math.abs(dz - dv.radius));

            b_ll = 1 / (1 + dv.DsDff * ((1 / b_ll) - 1)); // Michailidou et al. 2009
            b_T = 1 / (1 + dv.DsDff * ((1 / b_T) - 1)); // Michailidou et al. 2009

            b_ll = Math.sqrt(b_ll);
            b_T = Math.sqrt(b_T);

            stepx = ranGauss() * step3 * b_ll;
            stepy = ranGauss() * step3 * b_ll;
            stepz = ranGauss() * step3 * b_T;

        }

        return setVesicleLocation(dv, dv.x + stepx, dv.y + stepy, dv.z + stepz, false);

    }

    public int moveVesiclesActiveZone() {

        int moved = 0;
        double d, d2AZ;
        boolean outOfBounds, overlap, sameVoxel, insideAZ = false, outsideAZ = true, moveTethered;
        boolean ok;

        DiffusantVesicleAZ dv;
        DiffusantVesicleAZ testVesicle = new DiffusantVesicleAZ(project, "ready", 0, 0, 0, 0, 0);

        if (diffusant == null) {
            return 0;
        }

        for (int i = 0; i < diffusant.length; i++) {

            if ((diffusant[i] == null) || (diffusant[i].vesicle == null)) {
                continue;
            }

            diffusant[i].shuffleVesicles();

            for (int j = 0; j < diffusant[i].vesicle.length; j++) {

                if (!(diffusant[i].vesicle[j] instanceof DiffusantVesicleAZ)) {
                    continue;
                }

                dv = (DiffusantVesicleAZ) diffusant[i].vesicle[j];

                if (!dv.mobile || dv.isDocked || dv.isReserve || dv.isBound) {
                    continue;
                }

                testVesicle.copy(dv);

                if (hydroWallZ) {
                    ok = moveVesicleGaussHydroWallz(testVesicle);
                    if (!ok) {
                        //    continue;
                    }
                } else {
                    ok = moveVesicleGauss(testVesicle);
                    if (!ok) {
                    //    continue;
                    }
                }

                if (ok && !removingVesicleOverlap) {

                    insideAZ = isInsideActiveZone(testVesicle);

                    if (dockingOn && dv.isDocked && !insideAZ) {
                        ok = false;
                        //continue; // dont let docked vesicles move off AZ
                    }

                }

                if (ok && tetherVesicles && dv.isTetheredToAZ) {
                    if (!isInsideActiveZoneTethering(testVesicle)) {
                        ok = false;
                        //continue; // dont let tethered vesicles move outside AZ tethering region
                    }
                }

                if (ok) {
                    outOfBounds = outOfBounds(testVesicle);
                } else {
                    outOfBounds = false;
                }

                if (ok && outOfBounds) {
                    if (dockingOn && insideAZ) {

                        testVesicle.z = activeZone.z1;

                    } else if (PBC) {
                        outOfBounds = wrapAtBorder(testVesicle);
                    } else {
                        ok = false;
                        //continue;
                    }
                }

                if (outOfBounds) {
                    ok = false;
                    //continue;
                }

                if (ok) {
                    if (testMitochondriaOverlap(testVesicle)) {
                        ok = false;
                        //continue;
                    }
                }

                if (ok && !freeDiffusion) {

                    overlap = testVesicleOverlap(testVesicle, dv) != null;

                    if (overlap) {
                        ok = false;
                        //continue;
                    }

                }

                if (ok && connectVesicles && (dv.connectTo != null)) {

                    moveTethered = true;

                    for (int k = 0; k < dv.connectTo.length; k++) {

                        if (dv.connectTo[k] == null) {
                            continue;
                        }

                        if (!testVesicle.overlap(dv.connectTo[k], connectorLength)) {
                            moveTethered = false;
                            break;
                        }

                    }

                    if (!moveTethered) {
                        ok = false;
                        //continue;
                    }

                }

                if (saveMSDspatial && !removingVesicleOverlap) {
                    updateDspatial(dv, ok);
                }

                if (!ok) {
                    continue;
                }

                moved++;

                if (dockingOn && (time >= dockingOnTime) && !dv.isDocked) {

                    if (insideAZ) {
                        if (!addToDockedList(dv, time)) {
                            continue;
                        }
                    }

                }

                sameVoxel = false;

                if (testVesicle.voxel == dv.voxel) {
                    sameVoxel = true;
                }

                if (!sameVoxel) {
                    removeFromVoxelList(dv, dv.voxel);
                }

                setVesicleLocation(dv, testVesicle.x, testVesicle.y, testVesicle.z, false);

                if (!sameVoxel) {
                    addToVoxelList(dv);
                }

            }

        }

        return moved;

    }

    public boolean isInsideActiveZone(DiffusantVesicleAZ dv) {
        return activeZone.isInsideCylinderZ(dv.x, dv.y, dv.z);
    }

    public boolean isInsideActiveZoneTethering(DiffusantVesicleAZ dv) {
        return activeZoneTethering.isInsideCylinderZ(dv.x, dv.y, dv.z);
    }

    public void updateDspatial(DiffusantVesicleAZ dv, boolean ok) {
    }

    @Override
    public boolean initSave() {

        super.initSave();

        int dataPoints = 1;

        if (save_Release != null) {
            save_Release.init("Monte Carlo AZ Release Times", null, -1, dataPoints);
        }

        return true;

    }

    @Override
    public boolean finishSave() {

        super.finishSave();

        if (save_Release != null) {
            save_Release.finish("Monte Carlo AZ Release Times", null, -1);
        }

        return true;

    }

    public boolean saveReleaseData(double t, double d) {
        if (d == Double.NaN) {
            return save_Release.writeData(t);
        } else {
            return save_Release.writeString(Double.toString(t) + '\t' + Double.toString(d));
        }
    }

    @Override
    public boolean addUser(ParamVector pv) {

        if (pv == null) {
            return false;
        }

        super.addUser(pv);

        if (activeZone != null) {
            activeZone.addUser(pv);
        }

        if (save_Release != null) {
            save_Release.addUser(pv);
        }

        return true;

    }

    @Override
    public boolean createVector(boolean close) {

        if (!super.createVector(false)) {
            return false;
        }

        if (activeZone != null) {
            addBlankParam();
            activeZone.createVector(true);
            addVector(activeZone.getVector());
            activeZone.addUser(this);
        }

        if (save_Release != null) {
            addBlankParam();
            save_Release.createVector(true);
            addVector(save_Release.getVector());
            save_Release.addUser(this);
        }

        if (close) {
            closeVector();
        }

        return true;

    }

    @Override
    public void updateVector(ParamObject[] v) {

        super.updateVector(v);

        if (activeZone != null) {
            activeZone.updateVector(v);
        }

        if (save_Release != null) {
            save_Release.updateVector(v);
        }

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        String n = o.getName();

        if (!(o.paramVector instanceof RunMonteCarloAZ)) {
            return false;
        }

        if (super.setMyParams(o, v)) {
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
        if (n.equalsIgnoreCase("releaseRate")) {
            releaseRate = v;
            return true;
        }
        if (n.equalsIgnoreCase("releaseStartTime")) {
            releaseStartTime = v;
            return true;
        }
        if (n.equalsIgnoreCase("releaseLimit")) {
            releaseLimit = (int) v;
            return true;
        }
        if (n.equalsIgnoreCase("replenishRate")) {
            replenishRate = v;
            return true;
        }
        if (n.equalsIgnoreCase("azWidth")) {
            azWidth = Math.abs(v);
            return true;
        }
        if (n.equalsIgnoreCase("azHeight")) {
            azHeight = Math.abs(v);
            return true;
        }
        if (n.equalsIgnoreCase("tetherVesicles")) {
            if (v == 1) {
                tetherVesicles = true;
            } else {
                tetherVesicles = false;
            }
            return true;
        }
        if (n.equalsIgnoreCase("azNumTethers")) {
            azNumTethers = (int) v;
            return true;
        }
        return false;
    }

}
