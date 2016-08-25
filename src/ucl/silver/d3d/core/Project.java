package ucl.silver.d3d.core;

import java.io.*;
import ucl.silver.d3d.init.*;
import ucl.silver.d3d.utils.*;

/**
 * <p>Title: D3D</p>
 * <p>Description: 3D Diffusion-Reaction Simulator</p>
 * <p>Copyright: Copyright (c) 2016</p>
 * <p>Company: The Silver Lab at University College London</p>
 * @author Jason S. Rothman
 * @version 1.6
 */

public class Project extends ParamVector {

    public String initClassAndFunction = null; // (e.g. "InitProject.initCube")
    public String file = null; // file name, if this Project was opened from previously saved class
    public String psfFile = null; // psf file to open
    public String geomFile = null; // geometry file to open

    public double simTime = 1; // simulation time (ms)
    public double dt = 1; // time step (ms), computed below
    public double dx = 1; // voxel width (um)
    public double stability = 0.4; // stability factor for dt computation ( see below )

    public double maxD = 0; // max diffusion coefficient of all diffusants ( see maxD() )
    public double maxC0 = 0; // maximum concentration (mM) of diffusants and sources ( see maxC0() )

    public double saveRate = 100; // sample rate of output files (kHz)
    public double saveDT = 1 / saveRate; // sample rate of output files (ms)

    public double printRate = 1; // rate simulation results are printed (kHz)
    public double printDT = 1 / printRate; // rate simulation results are printed (ms)

    public int batchNum = -1; // file batch number, use negative number (-1) for no batch number
    private int batchCounter = -1; // Batch array counter

    public String seed = ""; // seed for random number generator

    public String timeUnits = "ms"; // time dimension
    public String freqUnits = "kHz"; // inverse time dimension
    public String spaceUnits = "um"; // space dimension
    public String volumeUnits = spaceUnits + "^3"; // volume dimension
    public String concUnits = "mM"; // concentration dimension
    public String diffusionUnits = spaceUnits + "^2/" + timeUnits; // volume dimension

    public String directory = ""; // diretory where folder/files are saved
    public String folder = "D3Doutput";
    private String subfolder = ""; // subdirectory for batches

    public String date = "";

    public InitProject initProject = null;
    public Geometry geometry = null;
    public Diffusant[] diffusants = null;
    public Source[] sources = null;
    public Detector[] detectors = null;
    public Batch[] batches = null;
    private Batch saveBatch = null;
    public D3Derror[] errors = null;

    public double[] sourceArray = null;

    public RunFiniteDifference finiteDifference = null;
    public RunMonteCarlo monteCarlo = null;

    public transient StopWatch timer = new StopWatch();

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("simTime")) {
            return timeUnits;
        }
        if (name.equalsIgnoreCase("dt")) {
            return timeUnits;
        }
        if (name.equalsIgnoreCase("dx")) {
            return spaceUnits;
        }
        if (name.equalsIgnoreCase("maxD")) {
            return diffusionUnits;
        }
        if (name.equalsIgnoreCase("maxC0")) {
            return concUnits;
        }
        if (name.equalsIgnoreCase("saveRate")) {
            return freqUnits;
        }
        if (name.equalsIgnoreCase("saveDT")) {
            return timeUnits;
        }
        if (name.equalsIgnoreCase("printRate")) {
            return freqUnits;
        }
        if (name.equalsIgnoreCase("printDT")) {
            return timeUnits;
        }
        if (name.equalsIgnoreCase("directory")) {
            return "DIR";
        }
        return super.units(name);
    }

    @Override
    public String help(String name) {
        if (name.equalsIgnoreCase("name")) {
            return "project name";
        }
        if (name.equalsIgnoreCase("initClassAndFunction")) {
            return "initialization Class and function";
        }
        if (name.equalsIgnoreCase("simTime")) {
            return "simulation time";
        }
        if (name.equalsIgnoreCase("dt")) {
            return "time step";
        }
        if (name.equalsIgnoreCase("dx")) {
            return "voxel cube dimensions";
        }
        if (name.equalsIgnoreCase("stability")) {
            return "used to compute dt, dt=(stability*dx*dx/(3*maxD);";
        }
        if (name.equalsIgnoreCase("maxD")) {
            return "maximum diffusion constant, used to compute dt";
        }
        if (name.equalsIgnoreCase("maxC0")) {
            return "maximum initial concentration";
        }
        if (name.equalsIgnoreCase("timeUnits")) {
            return "time units";
        }
        if (name.equalsIgnoreCase("freqUnits")) {
            return "inverse time units";
        }
        if (name.equalsIgnoreCase("spaceUnits")) {
            return "space dimension units";
        }
        if (name.equalsIgnoreCase("concUnits")) {
            return "concentration units";
        }
        if (name.equalsIgnoreCase("saveRate")) {
            return "default rate for saving data";
        }
        if (name.equalsIgnoreCase("saveDT")) {
            return "default time step for saving data";
        }
        if (name.equalsIgnoreCase("printRate")) {
            return "rate of simulation progress display";
        }
        if (name.equalsIgnoreCase("printDT")) {
            return "time step of simulation progress display";
        }
        return super.units(name);
    }

    @Override
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("file")) {
            return false;
        }
        if (name.equalsIgnoreCase("initClassAndFunction")) {
            return false;
        }
        if (name.equalsIgnoreCase("psfFile")) {
            return false;
        }
        if (name.equalsIgnoreCase("geomFile")) {
            return false;
        }
        if (name.equalsIgnoreCase("maxD")) {
            return false;
        }
        if (name.equalsIgnoreCase("maxC0")) {
            return false;
        }
        return true;
    }

    public Project() {
        super(null);
        geometry = new Geometry(this, 10, 10, 10);
        finiteDifference = new RunFiniteDifference(this);
        name = "D3D Finite Difference";
        createVector(true);
    }

    public Project(String type) {

        super(null);

        geometry = new Geometry(this, 10, 10, 10);

        if (type.equalsIgnoreCase("FD")) {
            newFiniteDifference();
        } else if (type.equalsIgnoreCase("MC")) {
            newMonteCarlo();
        } else if (type.equalsIgnoreCase("MCAZ")) {
            newMonteCarloAZ();
        }

        createVector(true);

    }

    public void newFiniteDifference() {
        monteCarlo = null;
        finiteDifference = new RunFiniteDifference(this);
        name = "D3D Finite Difference";
        Master.log("started Finite Difference project");
    }

    public void newMonteCarlo() {
        finiteDifference = null;
        monteCarlo = new RunMonteCarlo(this);
        name = "D3D Monte Carlo";
        Master.log("started Monte Carlo project");
    }

    public void newMonteCarloAZ() {
        finiteDifference = null;
        monteCarlo = new RunMonteCarloAZ(this);
        name = "D3D Monte Carlo Active Zone";
        Master.log("started Monte Carlo active zone project");
    }

    public void newMonteCarloAZEM() {
        finiteDifference = null;
        monteCarlo = new RunMonteCarloAZEM(this);
        name = "D3D Monte Carlo Active Zone";
        Master.log("started Monte Carlo active zone project");
    }

    public void setDate(){
        date = Master.currentDate();
        updateVectors();
    }

    public void nullArrays() {
        diffusants = null;
        sources = null;
        detectors = null;
        batches = null;
        errors = null;
    }

    public void init() {

        setDate();

        maxD = maxD();

        if (maxD > 0){
            dt = (stability * Math.pow(dx, 2.0)) / (3 * maxD); // "Mathematics of Diffusion", J. Crank, p. 151, Eq. 8.50
        }

        setParamError("dt", null);
        
        if ((dt <= 0 ) || (Double.isNaN(dt)) || (Double.isInfinite(dt))) {
            error("Project.init", "dt", "bad value");
        }

        geometry.init();

        if (finiteDifference != null ) {
            finiteDifference.init();
        } else if (monteCarlo != null ) {
            monteCarlo.init();
        }

        initDiffusants(-1);
        initSources(-1);
        initDetectors(-1);
        initBatches();

        maxC0 = maxC0();

        updateVectors();

        Master.updateMainFrameTitle();
        

        //Master.log("initialized project: " + name);

    }

    public double maxC0() {

        double max = 0;

        if (diffusants != null) {
            for (int i = 0; i < diffusants.length; i++) {
                max = Math.max(max, diffusants[i].C0);
            }
        }

        if (sources != null) {
            for (int i = 0; i < sources.length; i++) {
                max = Math.max(max, sources[i].C0);
            }
        }

        return max;

    }

    public void setOutputRate(double newRate) {

        if (newRate > 0) {

            saveRate = newRate;
            saveDT = 1.0 / newRate;

            if (diffusants != null) {
                for (int i = 0; i < diffusants.length; i += 1) {
                    diffusants[i].save.setOutputRate(newRate);
                }
            }

            if (detectors != null) {
                for (int i = 0; i < detectors.length; i += 1) {
                    detectors[i].save.setOutputRate(newRate);
                }
            }

            if (sources != null) {
                for (int i = 0; i < sources.length; i += 1) {
                    sources[i].save.setOutputRate(newRate);
                }
            }

            if (monteCarlo != null) {
                monteCarlo.setOutputRate(newRate);
            }

            updateVectors();

        }

    }

    public void setPrintRate(double newRate){

        if (newRate > 0) {
            printRate = newRate;
            printDT = 1.0 / newRate;
            updateVectors();
        }

    }

    public void setSaveRate(double newRate){

        if (newRate > 0) {
            saveRate = newRate;
            saveDT = 1.0 / newRate;
            updateVectors();
        }

    }

    public String fullDirectory() {

        String dir = directory;

        if (folder.length() > 0) {
            dir += folder + "/";
        }

        if (subfolder.length() > 0) {
            dir += subfolder + "/";
        }

        return dir;

    }

    public boolean checkDirectory() {

        File f = new File(fullDirectory());

        if (!f.isDirectory()) {
            f.mkdirs();
        }

        return f.isDirectory();

    }

    //
    // DIFFUSANT FUNCTIONS
    //

    public double maxD() {

        double max = 0;

        if (diffusants != null) {
            for (int i = 0; i < diffusants.length; i++) {
                max = Math.max(max, diffusants[i].maxD());
            }
        }

        return max;

    }

    public boolean reactionExists() {
        if (diffusants != null) {
            for (int i = 0; i < diffusants.length; i++) {
                if ((diffusants[i].reaction)) {
                    return true;
                }
            }
        }
        return false;
    }

    public int numDiffusants() {
        if (diffusants == null) {
            return 0;
        }
        return diffusants.length;
    }

    public boolean checkDiffusantNum(int arrayNum) {
        if ((diffusants != null) && (arrayNum >= 0) && (arrayNum < diffusants.length)) {
            return true;
        }
        return false;
    }

    public String diffusantName(int arrayNum) {
        if ((diffusants != null) && (arrayNum >= 0) && (arrayNum < diffusants.length)) {
            return diffusants[arrayNum].name;
        }
        return null;
    }

    public void diffusantsUpdateVectors() {
        if (diffusants != null) {
            for (int i = 0; i < diffusants.length; i++) {
                diffusants[i].updateVectors();
            }
        }
    }

    public void checkDiffusants() {
        if (diffusants != null) {
            for (int i = 0; i < diffusants.length; i++) {
                diffusants[i].check();
            }
        }
    }

    public void initDiffusantsPulseTimer() {
        if (diffusants != null) {
            for (int i = 0; i < diffusants.length; i++) {
                diffusants[i].initPulseTimer();
            }
        }
    }

    public void openDiffusants() {

        checkDirectory();

        if (diffusants != null) {
            for (int i = 0; i < diffusants.length; i++) {
                diffusants[i].saveInit();
            }
        }

    }

    public void closeDiffusants() {
        if (diffusants != null) {
            for (int i = 0; i < diffusants.length; i++) {
                diffusants[i].saveFinish();
            }
        }
    }

    public void initDiffusants(int dNum) {

        int ibgn = dNum;
        int iend = dNum;

        if (diffusants == null) {
            return;
        }
        
        if ((ibgn < 0) || (iend >= diffusants.length)) {
            ibgn = 0;
            iend = diffusants.length - 1;
        }

        for (int i = ibgn; i <= iend; i++) {

            if (diffusants[i].name.length() == 0) {
                diffusants[i].name = "Diffusant" + Integer.toString(i);
            }

            diffusants[i].init();

        }

    }

    public Diffusant getDiffusant(int diffusantNum) {

        if (diffusants == null) {
            return null;
        }

        if ((diffusantNum >= 0) && (diffusantNum < diffusants.length)) {
            return diffusants[diffusantNum];
        }

        return null;

    }

    public Diffusant getDiffusant(String diffusantName) {

        if ((diffusantName == null) || (diffusantName.length() == 0)) {
            return null;
        }

        if (diffusants == null) {
            return null;
        }

        for (int i = 0; i < diffusants.length; i++) {

            if (diffusants[i] == null) {
                continue;
            }

            if (diffusants[i].name.equalsIgnoreCase(diffusantName)) {
                return diffusants[i];
            }

        }

        return null;

    }

    public int getDiffusantNum(String diffusantName) {

        if ((diffusantName == null) || (diffusantName.length() == 0)) {
            return -1;
        }

        if (diffusants == null) {
            return -1;
        }

        for (int i = 0; i < diffusants.length; i++) {

            if (diffusants[i] == null) {
                continue;
            }

            if (diffusants[i].name.equalsIgnoreCase(diffusantName)) {
                return i;
            }

        }

        return -1;

    }

    public int addDiffusant(Diffusant newDiffusant) {

        int i = 0;

        if (diffusants != null) {
            i = diffusants.length;
        }

        Diffusant[] newArray = new Diffusant[i+1];

        if (i > 0) {
            System.arraycopy(diffusants, 0, newArray, 0, i);
        }

        newArray[i] = newDiffusant;

        diffusants = newArray; // replace old array with new one

        Master.log("added Diffusant #" + i + " : " + newDiffusant.name);

        initDiffusants(i);

        return i;

    }

    public boolean killDiffusant(int i) {

        int k = 0;

        if (!checkDiffusantNum(i)) {
            return false;
        }

        if (diffusants.length == 1) {
            diffusants = null;
            Master.log("killed Diffusant #0");
            return true;
        }

        Diffusant[] part = new Diffusant[diffusants.length - 1]; // new array

        for (int j = 0; j < diffusants.length; j++) {

            if (j == i) {
                continue;
            }

            part[k] = diffusants[j];
            k++;

        }

        diffusants = part; // replace old array with new one

        Master.log("killed Diffusant #" + i);

        return true;

    }

    public boolean setDiffusantParam(int arrayNum, String varName, double v) {

        if ((diffusants == null) || (arrayNum < 0) || (arrayNum >= diffusants.length)) {
            return false;
        }

        return diffusants[arrayNum].set(varName, v);

    }

    public boolean setDiffusantParam(int arrayNum, String varName, String s) {

        if ((diffusants == null) || (arrayNum < 0) || (arrayNum >= diffusants.length)) {
            return false;
        }

        return diffusants[arrayNum].set(varName, s);

    }

    //
    // SOURCE FUNCTIONS
    //

    public boolean sourceExists() {

        if (sources != null) {
            return true;
        }

        return false;

    }

    public double maxSource() {

        double max = 0;

        if (sources != null) {
            for (int i = 0; i < sources.length; i++) {
                max = Math.max(max, sources[i].C0);
            }
        }

        return max;

    }

    public int numSources() {

        if (sources == null) {
            return 0;
        }

        return sources.length;

    }

    public boolean checkSourceNum(int i) {

        if ((sources != null) && (i >= 0) && (i < sources.length)) {
            return true;
        }

        return false;

    }

    public void sourcesUpdateVectors() {
        if (sources != null) {
            for (int i = 0; i < sources.length; i++) {
                sources[i].updateVectors();
            }
        }
    }

    public void checkSources() {
        if (sources != null) {
            for (int i = 0; i < sources.length; i++) {
                sources[i].check();
            }
        }
    }

    public void initSourcesPulseTimer() {
        if (sources != null) {
            for (int i = 0; i < sources.length; i++) {
                sources[i].initPulseTimer();
            }
        }
    }

    public void openSources() {

        checkDirectory();

        if (sources != null) {
            for (int i = 0; i < sources.length; i++) {
                sources[i].saveInit();
            }
        }

    }

    public void closeSources() {
        if (sources != null) {
            for (int i = 0; i < sources.length; i++) {
                sources[i].saveFinish();
            }
        }
    }

    public void initSources(int sNum) {
        
        int ibgn = sNum;
        int iend = sNum;

        if (sources == null) {
            return;
        }

        if ((sNum < 0) || (sNum >= sources.length)) {
            ibgn = 0;
            iend = sources.length - 1;
        }

        for (int i = ibgn; i <= iend; i++) {

            if (sources[i].name.length() == 0) {
                sources[i].name = "Source" + Integer.toString(i);
            }

            sources[i].init();

        }

    }

    public Source getSource(int sourceNum) {

        if (sources == null) {
            return null;
        }

        if ((sourceNum >= 0) && (sourceNum < sources.length)) {
            return sources[sourceNum];
        }

        return null;

    }

    public Source getSource(String sourceName) {

        if ((sourceName == null) || (sourceName.length() == 0)) {
            return null;
        }

        if (sources == null) {
            return null;
        }

        for (int i = 0; i < sources.length; i++) {

            if (sources[i] == null) {
                continue;
            }

            if (sources[i].name.equalsIgnoreCase(sourceName)) {
                return sources[i];
            }

        }

        return null;

    }

    public int getSourceNum(String sourceName) {

        if ((sourceName == null) || (sourceName.length() == 0)) {
            return -1;
        }

        if (sources == null) {
            return -1;
        }

        for (int i = 0; i < sources.length; i++) {

            if (sources[i] == null) {
                continue;
            }

            if (sources[i].name.equalsIgnoreCase(sourceName)) {
                return i;
            }

        }

        return -1;

    }

    public int addSource(Source newSource) {

        int i = 0;

        if (sources != null) {
            i = sources.length;
        }

        Source[] newArray = new Source[i+1];
        
        if (i > 0) {
            System.arraycopy(sources, 0, newArray, 0, i);
        }

        newArray[i] = newSource;

        sources = newArray; // replace old array with new one

        Master.log("added Source #" + i + " : " + newSource.name);

        initSources(i);

        return i;

    }

    public boolean killSource(int i) {
        int k = 0;

        if (!checkSourceNum(i)) {
            return false;
        }

        if (sources.length == 1) {
            sources = null;
            Master.log("killed Source #0");
            return true;
        }

        Source[] s = new Source[sources.length - 1]; // new array

        for (int j = 0; j < sources.length; j++) {
            if (j == i) {
                continue;
            }
            s[k] = sources[j];
            k++;
        }

        sources = s; // replace old array with new one

        Master.log("killed Source #" + i);

        return true;

    }

    public boolean setSourceParam(int arrayNum, String varName, double v) {
        if ((sources == null) || (arrayNum < 0) || (arrayNum >= sources.length)) {
            return false;
        }
        return sources[arrayNum].set(varName, v);
    }

    public boolean setSourceParam(int arrayNum, String varName, String s) {
        if ((sources == null) || (arrayNum < 0) || (arrayNum >= sources.length)) {
            return false;
        }
        return sources[arrayNum].set(varName, s);
    }

    //
    // DETECTOR FUNCTIONS
    //
    public boolean detectorExists() {
        if (detectors != null) {
            return true;
        }
        return false;
    }

    public int numDetectors() {
        if (detectors == null) {
            return 0;
        }
        return detectors.length;
    }

    public boolean checkDetectorNum(int i) {
        if ((detectors != null) && (i >= 0) && (i < detectors.length)) {
            return true;
        }
        return false;
    }

    public void detectorsUpdateVectors() {
        if (detectors != null) {
            for (int i = 0; i < detectors.length; i++) {
                detectors[i].updateVectors();
            }
        }
    }

    public void checkDetectors() {
        if (detectors != null) {
            for (int i = 0; i < detectors.length; i++) {
                detectors[i].check();
            }
        }
    }

    public void initDetectorsPulseTimer() {
        if (detectors != null) {
            for (int i = 0; i < detectors.length; i++) {
                detectors[i].initPulseTimer();
            }
        }
    }

    public void openDetectors() {

        checkDirectory();

        if (detectors != null) {
            for (int i = 0; i < detectors.length; i++) {
                detectors[i].saveInit();
            }
        }

    }

    public void closeDetectors() {
        if (detectors != null) {
            for (int i = 0; i < detectors.length; i++) {
                detectors[i].saveFinish();
            }
        }
    }

    public void initDetectors(int dNum) {

        int ibgn = dNum;
        int iend = dNum;

        if (detectors == null) {
            return;
        }

        if ((dNum < 0) || (dNum >= detectors.length)) {
            ibgn = 0;
            iend = detectors.length - 1;
        }

        for (int i = ibgn; i <= iend; i++) {

            if (detectors[i].name.length() == 0) {
                detectors[i].name = "Detector" + Integer.toString(i);
            }

            detectors[i].init();

        }

    }

    public Detector getDetector(int detectorNum) {

        if (detectors == null) {
            return null;
        }

        if ((detectorNum >= 0) && (detectorNum < detectors.length)) {
            return detectors[detectorNum];
        }

        return null;

    }

    public Detector getDetector(String detectorName) {

        if ((detectorName == null) || (detectorName.length() == 0)) {
            return null;
        }

        if (detectors == null) {
            return null;
        }

        for (int i = 0; i < detectors.length; i++) {

            if (detectors[i] == null) {
                continue;
            }

            if (detectors[i].name.equalsIgnoreCase(detectorName)) {
                return detectors[i];
            }

        }

        return null;

    }

    public int getDetectorNum(String detectorName) {

        if ((detectorName == null) || (detectorName.length() == 0)) {
            return -1;
        }

        if (detectors == null) {
            return -1;
        }

        for (int i = 0; i < detectors.length; i++) {

            if (detectors[i] == null) {
                continue;
            }

            if (detectors[i].name.equalsIgnoreCase(detectorName)) {
                return i;
            }

        }

        return -1;

    }

    public int addDetector(Detector newDetector) {

        int i = 0;

        if (detectors != null) {
            i = detectors.length;
        }

        Detector[] newArray = new Detector[i+1];
        
        if (i > 0) {
            System.arraycopy(detectors, 0, newArray, 0, i);
        }

        newArray[i] = newDetector;

        detectors = newArray; // replace old array with new one

        Master.log("added detector #" + i + " : " + newDetector.name);

        initDetectors(i);

        return i;

    }

    public boolean killDetector(int i) {
        int k = 0;

        if (!checkDetectorNum(i)) {
            return false;
        }

        if (detectors.length == 1) {
            detectors = null;
            Master.log("killed Detector #0");
            return true;
        }

        Detector[] detect = new Detector[detectors.length - 1]; // new array

        for (int j = 0; j < detectors.length; j++) {
            if (j == i) {
                continue;
            }
            detect[k] = detectors[j];
            k++;
        }

        detectors = detect; // replace old array with new one

        Master.log("killed Detector #" + i);

        return true;

    }

    public boolean setDetectorParam(int arrayNum, String varName, double v) {
        if ((detectors == null) || (arrayNum < 0) || (arrayNum >= detectors.length)) {
            return false;
        }
        return detectors[arrayNum].set(varName, v);
    }

    public boolean setDetectorParam(int arrayNum, String varName, String s) {
        if ((detectors == null) || (arrayNum < 0) || (arrayNum >= detectors.length)) {
            return false;
        }
        return detectors[arrayNum].set(varName, s);
    }

    //
    // BATCH FUNCTIONS
    //

    public boolean setBatchNum(int i) {

        if (i >= 0) {
            batchNum = i;
        } else {
            batchNum = -1;
        }

        init();

        return true;
    }

    public int numBatches() {
        if (batches == null) {
            return 0;
        }
        return batches.length;
    }

    public boolean checkBatchNum(int i) {
        if ((batches != null) && (i >= 0) && (i < batches.length)) {
            return true;
        }
        return false;
    }

    public void initBatches() {

        String parameter, strValue, fname;

        if (batches != null) {

            for (int i = 0; i < batches.length; i++) {

                parameter = batches[i].parameter;

                if (batches[i].isString()) {
                    strValue = "";
                } else {
                    strValue = Float.toString((float) batches[i].value);
                }

                strValue = strValue.replace('.', '_');
                strValue = strValue.replace(';', '_');
                strValue = strValue.replace(',', '_');

                batches[i].init();

                if (!paramVectorVariableExists(parameter)) {
                    batches[i].error("Project.initBatches", "parameter", "does not exist");
                }

                fname = "Batch" + Integer.toString(batches[i].batchNum) + "_" + parameter;

                if (strValue.length() > 0) {
                    fname += "_" + strValue;
                }

                fname = fname.replace('.', '_');

                if (batches[i].execute) {
                    if (batches[i].folder.length() == 0) {
                        batches[i].folder = fname;
                    }
                }

            }

        }

    }

    public boolean nextBatch() {

        double saveVar;
        String saveStr = "", printStr  = "";

        if (batches == null) {
            return false;
        }

        if (saveBatch != null) { // restore original parameter value if it exists

            if (saveBatch.isString()) {
                setObjectParam(saveBatch.parameter, saveBatch.strValue, false, false);
            } else {
                setObjectParam(saveBatch.parameter, saveBatch.value);
            }

        }

        for (int i = 0; i < batches.length; i += 1) {

            if (!batches[i].finished) {

                batchCounter = i;
                batchNum = batches[i].batchNum;
                subfolder = batches[i].folder;

                if (batches[i].isString()) {
                    saveStr = getObjectParamStr(batches[i].parameter);
                    saveBatch = new Batch(batchNum, batches[i].parameter, saveStr, batches[i].execute);
                    setObjectParam(batches[i].parameter, batches[i].strValue, false, false);
                    printStr = "New Batch: " + batches[i].parameter + ", " + batches[i].strValue;
                } else {
                    saveVar = getObjectParamVar(batches[i].parameter);
                    saveBatch = new Batch(batchNum, batches[i].parameter, saveVar, batches[i].execute);
                    setObjectParam(batches[i].parameter, batches[i].value);
                    printStr = "New Batch: " + batches[i].parameter + ", " + batches[i].value;
                }

                if (batches[i].execute) {
                    Master.log(printStr);
                    return true;
                }

            }

        }

        return false;

    }

    public void closeBatch() {
        if (checkBatchNum(batchCounter)) {
            batches[batchCounter].finished = true;
            subfolder = "";
        }
    }

    public int addBatch(Batch newBatch) {

        int i = 0;

        if (batches != null) {
            i = batches.length;
        }

        Batch[] newArray = new Batch[i+1];

        if (i > 0) {
            System.arraycopy(batches, 0, newArray, 0, i);
        }

        newArray[i] = newBatch;

        batches = newArray; // replace old array with new one

        //Master.log("added batch #" + i);

        return i;

    }

    public boolean killBatch(int i) {

        int k = 0;

        if (!checkBatchNum(i)) {
            return false;
        }

        if (batches.length == 1) {
            batches = null;
            Master.log("killed Batch #0");
            return true;
        }

        Batch[] b = new Batch[batches.length - 1]; // new array

        for (int j = 0; j < batches.length; j++) {
            if (j == i) {
                continue;
            }
            b[k] = batches[j];
            k++;
        }

        batches = b; // replace old array with new one

        Master.log("killed Batch #" + i);

        return true;

    }

    //
    // ERROR FUNCTIONS
    //

    public int addError(String errorStr) {
        return addError(new D3Derror(toString(), "", "", errorStr));
    }

    public int addError(D3Derror newError) {

        int i = 0;

        if (errors != null) {

            for (i = 0; i < errors.length; i++) {

                if (newError.object.equalsIgnoreCase(errors[i].object)) {
                    if (newError.function.equalsIgnoreCase(errors[i].function)) {
                        if (newError.parameter.equalsIgnoreCase(errors[i].parameter)) {
                            if (newError.error.equalsIgnoreCase(errors[i].error)) {
                                return -1; // already exists
                            }
                        }
                    }
                }

            }

            i = errors.length;

        }

        D3Derror[] newArray = new D3Derror[i+1];

        if (i > 0) {
            System.arraycopy(errors, 0, newArray, 0, i);
        }

        newArray[i] = newError;

        errors = newArray; // replace old array with new one

        return i;

    }

    public int numErrors() {

        if (errors == null) {
            return 0;
        }

        return errors.length;

    }
    
    public boolean checkErrorNum(int i) {
        if ((errors != null) && (i >= 0) && (i < errors.length)) {
            return true;
        }
        return false;
    }

    public boolean killError(int i) {
        int k = 0;

        if (!checkErrorNum(i)) {
            return false;
        }

        if (errors.length == 1) {
            errors = null;
            Master.log("killed Error #0");
            return true;
        }

        D3Derror[] e = new D3Derror[errors.length - 1]; // new array

        for (int j = 0; j < errors.length; j++) {
            if (j == i) {
                continue;
            }
            e[k] = errors[j];
            k++;
        }

        errors = e; // replace old array with new one

        Master.log("killed Error #" + i);

        return true;

    }

    //
    // Finite Difference and Monte Carlo Simulation Functions
    //

    public boolean simulationInit(boolean preview) {

        Master.writeToFiles = false;

        if (batches != null) {
            if (!nextBatch()) {
                return true; // finished batches
            }
        }

        if (monteCarlo != null) {
            if (monteCarlo.initSimulation()){
                error("ABORTED D3D init simulation");
                return true;
            }
        }

        init();

        if (numErrors() > 0) {
            error("ABORTED D3D init simulation");
            return true;
        }

        if (preview) {
            Master.log("preview D3D simulation, " + Master.currentDate());
        } else {
            Master.writeToFiles = true;
            Master.log("run D3D simulation, " + Master.currentDate());
        }

        if (monteCarlo != null) {
            monteCarlo.initSave();
        }

        geometry.checkSpace();
        geometry.checkPSFs();

        checkDiffusants();
        checkSources();
        checkDetectors();

        initDiffusantsPulseTimer();
        initDetectorsPulseTimer();
        initSourcesPulseTimer();

        openDiffusants();
        openSources();
        openDetectors();

        writeParamFile();

        return false;

    }

    public void simulationStart(boolean preview) {

        timer.start();

        init();

        if (numErrors() > 0) {
            Master.log("cannot start simulation due to registered error(s)");
            return;
        }

        if (monteCarlo != null) {

            monteCarlo.startSimulation(preview);

        } else {

            if ((finiteDifference == null) || finiteDifference.finished) {
                finiteDifference = new RunFiniteDifference(this);
            }

            finiteDifference.startSimulation(preview);

        }

    }

    public void simulationPause() {

        if (numErrors() > 0) {
            Master.log("cannot pause simulation due to registered error(s)");
            return;
        }

        if (monteCarlo != null) {
            monteCarlo.pauseSimulation();
        }

        if (finiteDifference != null) {
            finiteDifference.pauseSimulation();
        }

    }

    public void simulationFinish() {

        if (monteCarlo != null) {
            monteCarlo.finishSimulation();
        }

        closeDiffusants();
        closeSources();
        closeDetectors();
        closeBatch();

        timer.stop();

        Master.log(Master.currentDate());
        Master.log("total simulation session time = " + timer.toString());

    }

    public void simulationCancel() {

        if (numErrors() > 0) {
            Master.log("cannot cancel simulation due to registered error(s)");
            return;
        }

        if (monteCarlo != null) {
            monteCarlo.cancelSimulation();
        }

        if (finiteDifference != null) {
            finiteDifference.stopSimulation();
        }

    }

    //
    // Log File Functions
    //

    public String logFileName(){

        checkDirectory();

        String logFile = "D3Dlog";

        if (batchNum >= 0) {
            logFile += "_" + Integer.toString(batchNum);
        }

        logFile += ".dat";

        return logFile;

    }

    //
    // PARAMETER VECTOR FUNCTIONS
    //
    public void writeParamFile() {

        if (!Master.writeToFiles) {
            return;
        }

        String fileName = "D3Dconfig";
        
        if (batchNum >= 0) {
            fileName += "_" + Integer.toString(batchNum);
        }
        
        checkDirectory();

        geometry.updateVectors();
        diffusantsUpdateVectors();
        sourcesUpdateVectors();
        detectorsUpdateVectors();
        initProject.updateVectors();

        writeParamVectors(fileName);
        writeParamVectorsXML(fileName);

    }

    private boolean writeParamVectors(String fileName) {

        if (!Master.writeParamVectorOpen(fullDirectory(), fileName)) {
            return false;
        }
        
        Master.writeParamVector(this);
        Master.writeParamVector(geometry);
        Master.writeParamVectorArray(diffusants);
        Master.writeParamVectorArray(sources);
        Master.writeParamVectorArray(detectors);
        Master.writeParamVector(initProject);
        Master.writeParamVectorClose();

        return true;

    }

    private boolean writeParamVectorsXML(String fileName) {

        if (!Master.writeParamVectorXMLopen(fullDirectory(), fileName)) {
            return false;
        }
        
        Master.writeParamVectorXML(this);
        Master.writeParamVectorXML(geometry);
        Master.writeParamVectorArrayXML(diffusants);
        Master.writeParamVectorArrayXML(sources);
        Master.writeParamVectorArrayXML(detectors);
        Master.writeParamVectorXML(initProject);
        Master.writeParamVectorXMLclose();

        return true;

    }

    public boolean paramVectorVariableExists(String longName) {

        if (getObjectParamStr(longName) != null) {
            return true;
        }

        return false;

    }

    public double getObjectParamVar(String longName) {

        String parent = D3DUtility.parentClass(longName);
        String child = D3DUtility.childClass(longName);
        String parameter = D3DUtility.parameterName(longName);

        if (parent.equalsIgnoreCase("initProject")) {
            return initProject.getVar(parameter);
        }

        if (parent.equalsIgnoreCase("project")) {
            return getVar(parameter);
        }

        if (parent.equalsIgnoreCase("geometry")) {
            if (geometry == null) {
                return Double.NaN;
            }
            return geometry.getVar(parameter);
        }

        if (parent.equalsIgnoreCase("finiteDifferences")) {
            if (finiteDifference == null) {
                return Double.NaN;
            }
            return finiteDifference.getVar(parameter);
        }

        if (parent.equalsIgnoreCase("monteCarlo")) {

            if (monteCarlo == null) {
                return Double.NaN;
            }

            if (child == null) {
                return monteCarlo.getVar(parameter);
            }

            return Double.NaN;

        }

        if (parent.equalsIgnoreCase("monteCarloActiveZone")) {

            if (monteCarlo == null) {
                return Double.NaN;
            }

            RunMonteCarloAZ mcaz = (RunMonteCarloAZ) monteCarlo;

            if (child == null) {
                return mcaz.getVar(parameter);
            }

            if (child.equalsIgnoreCase("ActiveZone")) {
                if (mcaz.activeZone == null) {
                    return Double.NaN;
                }
                return mcaz.activeZone.getVar(parameter);
            }

            if (child.equalsIgnoreCase("ActiveZone")) {
                if (mcaz.activeZone == null) {
                    return Double.NaN;
                }
                return mcaz.activeZone.getVar(parameter);
            }

            if (mcaz.save_Release != null) {
                if (mcaz.save_Release.name.equalsIgnoreCase(child)) {
                    return mcaz.save_Release.getVar(parameter);
                }
            }

            return Double.NaN;

        }

        if (diffusants != null) {
            for (int i = 0; i < diffusants.length; i++) {
                if (diffusants[i] == null) {
                    continue;
                }
                if (diffusants[i].name.equalsIgnoreCase(parent)) {
                    if (child == null) {
                        return diffusants[i].getVar(parameter);
                    }
                    if (child.equalsIgnoreCase("Save")) {
                        if (diffusants[i].save == null) {
                            return Double.NaN;
                        }
                        return diffusants[i].save.getVar(parameter);
                    }
                    if (child.equalsIgnoreCase("PulseTimer")) {
                        if (diffusants[i].pulseTimer == null) {
                            return Double.NaN;
                        }
                        return diffusants[i].pulseTimer.getVar(parameter);
                    }
                    if (child.equalsIgnoreCase("PSF")) {
                        if (diffusants[i].psf == null) {
                            return Double.NaN;
                        }
                        return diffusants[i].psf.getVar(parameter);
                    }
                    return Double.NaN;
                }
            }
        }

        if (sources != null) {
            for (int i = 0; i < sources.length; i++) {
                if (sources[i] == null) {
                    continue;
                }
                if (sources[i].name.equalsIgnoreCase(parent)) {
                    if (child == null) {
                        return sources[i].getVar(parameter);
                    }
                    if (child.equalsIgnoreCase("Save")) {
                        if (sources[i].save == null) {
                            return Double.NaN;
                        }
                        return sources[i].save.getVar(parameter);
                    }
                    if (child.equalsIgnoreCase("PulseTimer")) {
                        if (sources[i].pulseTimer == null) {
                            return Double.NaN;
                        }
                        return sources[i].pulseTimer.getVar(parameter);
                    }
                    return Double.NaN;
                }
            }
        }

        if (detectors != null) {
            for (int i = 0; i < detectors.length; i++) {
                if (detectors[i] == null) {
                    continue;
                }
                if (detectors[i].name.equalsIgnoreCase(parent)) {
                    if (child == null) {
                        return detectors[i].getVar(parameter);
                    }
                    if (child.equalsIgnoreCase("Save")) {
                        if (detectors[i].save == null) {
                            return Double.NaN;
                        }
                        return detectors[i].save.getVar(parameter);
                    }
                    if (child.equalsIgnoreCase("PulseTimer")) {
                        if (detectors[i].pulseTimer == null) {
                            return Double.NaN;
                        }
                        return detectors[i].pulseTimer.getVar(parameter);
                    }
                    if (child.equalsIgnoreCase("PSF")) {
                        if (detectors[i].psf == null) {
                            return Double.NaN;
                        }
                        return detectors[i].psf.getVar(parameter);
                    }
                    return Double.NaN;
                }
            }
        }

        return Double.NaN;

    }

    public String getObjectParamStr(String longName) {

        String parent = D3DUtility.parentClass(longName);
        String child = D3DUtility.childClass(longName);
        String parameter = D3DUtility.parameterName(longName);

        //Master.log(parent + ", " + child + ", " + parameter);

        if (parent.equalsIgnoreCase("initProject")) {
            return initProject.getStr(parameter);
        }

        if (parent.equalsIgnoreCase("project")) {
            return getStr(parameter);
        }

        if (parent.equalsIgnoreCase("geometry")) {
            if (geometry == null) {
                return null;
            }
            return geometry.getStr(parameter);
        }

        if (parent.equalsIgnoreCase("FiniteDifferences")) {
            if (finiteDifference == null) {
                return null;
            }
            return finiteDifference.getStr(parameter);
        }

        if (parent.equalsIgnoreCase("MonteCarlo")) {

            if (monteCarlo == null) {
                return null;
            }

            if (child == null) {
                return monteCarlo.getStr(parameter);
            }

            return null;

        }

        if (parent.equalsIgnoreCase("MonteCarloActiveZone")) {

            if (monteCarlo == null) {
                return null;
            }

            RunMonteCarloAZ mcaz = (RunMonteCarloAZ) monteCarlo;

            if (child == null) {
                return mcaz.getStr(parameter);
            }

            if (child.equalsIgnoreCase("ActiveZone")) {
                if (mcaz.activeZone == null) {
                    return null;
                }
                return mcaz.activeZone.getStr(parameter);
            }

            if (mcaz.save_Release != null) {
                if (mcaz.save_Release.name.equalsIgnoreCase(child)) {
                    return mcaz.save_Release.getStr(parameter);
                }
            }

            return null;

        }

        if (diffusants != null) {
            for (int i = 0; i < diffusants.length; i++) {
                if (diffusants[i] == null) {
                    continue;
                }
                if (diffusants[i].name.equalsIgnoreCase(parent)) {
                    if (child == null) {
                        return diffusants[i].getStr(parameter);
                    }
                    if (child.equalsIgnoreCase("Save")) {
                        if (diffusants[i].save == null) {
                            return null;
                        }
                        return diffusants[i].save.getStr(parameter);
                    }
                    if (child.equalsIgnoreCase("PulseTimer")) {
                        if (diffusants[i].pulseTimer == null) {
                            return null;
                        }
                        return diffusants[i].pulseTimer.getStr(parameter);
                    }
                    if (child.equalsIgnoreCase("PSF")) {
                        if (diffusants[i].psf == null) {
                            return null;
                        }
                        return diffusants[i].psf.getStr(parameter);
                    }
                    return null;
                }
            }
        }

        if (sources != null) {
            for (int i = 0; i < sources.length; i++) {
                if (sources[i] == null) {
                    continue;
                }
                if (sources[i].name.equalsIgnoreCase(parent)) {
                    if (child == null) {
                        return sources[i].getStr(parameter);
                    }
                    if (child.equalsIgnoreCase("Save")) {
                        if (sources[i].save == null) {
                            return null;
                        }
                        return sources[i].save.getStr(parameter);
                    }
                    if (child.equalsIgnoreCase("PulseTimer")) {
                        if (sources[i].pulseTimer == null) {
                            return null;
                        }
                        return sources[i].pulseTimer.getStr(parameter);
                    }
                    return null;
                }
            }
        }

        if (detectors != null) {
            for (int i = 0; i < detectors.length; i++) {
                if (detectors[i] == null) {
                    continue;
                }
                if (detectors[i].name.equalsIgnoreCase(parent)) {
                    if (child == null) {
                        return detectors[i].getStr(parameter);
                    }
                    if (child.equalsIgnoreCase("Save")) {
                        if (detectors[i].save == null) {
                            return null;
                        }
                        return detectors[i].save.getStr(parameter);
                    }
                    if (child.equalsIgnoreCase("PulseTimer")) {
                        if (detectors[i].pulseTimer == null) {
                            return null;
                        }
                        return detectors[i].pulseTimer.getStr(parameter);
                    }
                    if (child.equalsIgnoreCase("PSF")) {
                        if (detectors[i].psf == null) {
                            return null;
                        }
                        return detectors[i].psf.getStr(parameter);
                    }
                    return null;
                }
            }
        }

        return null;

    }

    public boolean setObjectParam(String longName, double value) {

        String parent = D3DUtility.parentClass(longName);
        String child = D3DUtility.childClass(longName);
        String parameter = D3DUtility.parameterName(longName);

        if (parent.equalsIgnoreCase("initProject")) {
            if (initProject == null) {
                return false;
            }
            return initProject.set(parameter, value);
        }

        if (parent.equalsIgnoreCase("project")) {
            return set(parameter, value);
        }

        if (parent.equalsIgnoreCase("geometry")) {
            if (geometry == null) {
                return false;
            }
            return geometry.set(parameter, value);
        }

        if (parent.equalsIgnoreCase("FiniteDifferences")) {
            if (finiteDifference == null) {
                return false;
            }
            return finiteDifference.set(parameter, value);
        }

        if (parent.equalsIgnoreCase("MonteCarlo")) {

            if (monteCarlo == null) {
                return false;
            }

            if (child == null) {
                return monteCarlo.set(parameter, value);
            }

            return false;

        }

        if (parent.equalsIgnoreCase("MonteCarloActiveZone")) {

            if (monteCarlo == null) {
                return false;
            }

            RunMonteCarloAZ mcaz = (RunMonteCarloAZ) monteCarlo;

            if (child == null) {
                return mcaz.set(parameter, value);
            }

            if (child.equalsIgnoreCase("ActiveZone")) {
                if (mcaz.activeZone == null) {
                    return false;
                }
                return mcaz.activeZone.set(parameter, value);
            }

            if (mcaz.save_Release != null) {
                if (mcaz.save_Release.name.equalsIgnoreCase(child)) {
                    return mcaz.save_Release.set(parameter, value);
                }
            }

            return false;

        }

        if (diffusants != null) {
            for (int i = 0; i < diffusants.length; i++) {
                if (diffusants[i] == null) {
                    continue;
                }
                if (diffusants[i].name.equalsIgnoreCase(parent)) {
                    if (child == null) {
                        return diffusants[i].set(parameter, value);
                    }
                    if (child.equalsIgnoreCase("Save")) {
                        if (diffusants[i].save == null) {
                            return false;
                        }
                        return diffusants[i].save.set(parameter, value);
                    }
                    if (child.equalsIgnoreCase("PulseTimer")) {
                        if (diffusants[i].pulseTimer == null) {
                            return false;
                        }
                        return diffusants[i].pulseTimer.set(parameter, value);
                    }
                    if (child.equalsIgnoreCase("PSF")) {
                        if (diffusants[i].psf == null) {
                            return false;
                        }
                        return diffusants[i].psf.set(parameter, value);
                    }
                    return false;
                }
            }
        }

        if (sources != null) {
            for (int i = 0; i < sources.length; i++) {
                if (sources[i] == null) {
                    continue;
                }
                if (sources[i].name.equalsIgnoreCase(parent)) {
                    if (child == null) {
                        return sources[i].set(parameter, value);
                    }
                    if (child.equalsIgnoreCase("Save")) {
                        if (sources[i].save == null) {
                            return false;
                        }
                        return sources[i].save.set(parameter, value);
                    }
                    if (child.equalsIgnoreCase("PulseTimer")) {
                        if (sources[i].pulseTimer == null) {
                            return false;
                        }
                        return sources[i].pulseTimer.set(parameter, value);
                    }
                    return false;
                }
            }
        }

        if (detectors != null) {
            for (int i = 0; i < detectors.length; i++) {
                if (detectors[i] == null) {
                    continue;
                }
                if (detectors[i].name.equalsIgnoreCase(parent)) {
                    if (child == null) {
                        return detectors[i].set(parameter, value);
                    }
                    if (child.equalsIgnoreCase("Save")) {
                        if (detectors[i].save == null) {
                            return false;
                        }
                        return detectors[i].save.set(parameter, value);
                    }
                    if (child.equalsIgnoreCase("PulseTimer")) {
                        if (detectors[i].pulseTimer == null) {
                            return false;
                        }
                        return detectors[i].pulseTimer.set(parameter, value);
                    }
                    if (child.equalsIgnoreCase("PSF")) {
                        if (detectors[i].psf == null){
                            return false;
                        }
                        return detectors[i].psf.set(parameter, value);
                    }
                    return false;
                }
            }
        }

        return false;

    }
    
    public boolean setObjectParamProject(String longName, String strValue) {

        String parent = D3DUtility.parentClass(longName);
        String child = D3DUtility.childClass(longName);
        String parameter = D3DUtility.parameterName(longName);

        //Master.log(parent + "," + child + "," + parameter + "," + strValue);

        if (parent.equalsIgnoreCase("initProject")) {
            return initProject.set(parameter, strValue);
        }

        return false;

    }

    public boolean setObjectParamInitProject(String longName, String strValue) {

        String parent = D3DUtility.parentClass(longName);
        String child = D3DUtility.childClass(longName);
        String parameter = D3DUtility.parameterName(longName);

        //Master.log(parent + "," + child + "," + parameter + "," + strValue);

        if (parent.equalsIgnoreCase("initProject")) {
            return initProject.set(parameter, strValue);
        }

        return false;

    }

    public boolean setObjectParam(String longName, String strValue, boolean onlyProject, boolean onlyInitProject) {

        boolean set;

        String parent = D3DUtility.parentClass(longName);
        String child = D3DUtility.childClass(longName);
        String parameter = D3DUtility.parameterName(longName);
        
        if (onlyProject) {

            if (parent.equalsIgnoreCase("project")) {
                return set(parameter, strValue);
            }

            return false; // nothing

        }

        if (onlyInitProject) {

            if (parent.equalsIgnoreCase("initProject")) {
                return initProject.set(parameter, strValue);
            }

            return false; // nothing

        }

        if (parent.equalsIgnoreCase("geometry")) {

            if (geometry == null) {
                return false; // error
            }

            set = geometry.set(parameter, strValue);

            if (set) {
                if (initProject.initGeometry()) {
                    return false; // error
                } else {
                    return true; // OK
                }
            }

            return set;

        }

        if (parent.equalsIgnoreCase("finitedifferences")) {
            if (finiteDifference == null) {
                return false;
            }
            return finiteDifference.set(parameter, strValue);
        }

        if (parent.equalsIgnoreCase("MonteCarlo")) {

            if (monteCarlo == null) {
                return false;
            }

            if (child == null) {
                return monteCarlo.set(parameter, strValue);
            }

            return false;

        }

        if (parent.equalsIgnoreCase("MonteCarloActiveZone")) {

            if (monteCarlo == null) {
                return false;
            }

            RunMonteCarloAZ mcaz = (RunMonteCarloAZ) monteCarlo;

            if (child == null) {
                return mcaz.set(parameter, strValue);
            }

            if (child.equalsIgnoreCase("ActiveZone")) {
                if (mcaz.activeZone == null) {
                    return false;
                }
                return mcaz.activeZone.set(parameter, strValue);
            }

            if (mcaz.save_Release != null) {
                if (mcaz.save_Release.name.equalsIgnoreCase(child)) {
                    return mcaz.save_Release.set(parameter, strValue);
                }
            }

            return false;

        }

        if (diffusants != null) {
            for (int i = 0; i < diffusants.length; i++) {
                if (diffusants[i] == null) {
                    continue;
                }
                if (diffusants[i].name.equalsIgnoreCase(parent)) {
                    if (child == null) {
                        return diffusants[i].set(parameter, strValue);
                    }
                    if (child.equalsIgnoreCase("Save")) {
                        if (diffusants[i].save == null) {
                            return false;
                        }
                        return diffusants[i].save.set(parameter, strValue);
                    }
                    if (child.equalsIgnoreCase("PulseTimer")) {
                        if (diffusants[i].pulseTimer == null) {
                            return false;
                        }
                        return diffusants[i].pulseTimer.set(parameter, strValue);
                    }
                    if (child.equalsIgnoreCase("PSF")) {
                        if (diffusants[i].psf == null) {
                            return false;
                        }
                        return diffusants[i].psf.set(parameter, strValue);
                    }
                    return false;
                }
            }
        }

        if (sources != null) {
            for (int i = 0; i < sources.length; i++) {
                if (sources[i] == null) {
                    continue;
                }
                if (sources[i].name.equalsIgnoreCase(parent)) {
                    if (child == null) {
                        return sources[i].set(parameter, strValue);
                    }
                    if (child.equalsIgnoreCase("Save")) {
                        if (sources[i].save == null) {
                            return false;
                        }
                        return sources[i].save.set(parameter, strValue);
                    }
                    if (child.equalsIgnoreCase("PulseTimer")) {
                        if (sources[i].pulseTimer == null) {
                            return false;
                        }
                        return sources[i].pulseTimer.set(parameter, strValue);
                    }
                    return false;
                }
            }
        }

        if (detectors != null) {
            for (int i = 0; i < detectors.length; i++) {
                if (detectors[i] == null) {
                    continue;
                }
                if (detectors[i].name.equalsIgnoreCase(parent)) {
                    if (child == null) {
                        return detectors[i].set(parameter, strValue);
                    }
                    if (child.equalsIgnoreCase("Save")) {
                        if (detectors[i].save == null) {
                            return false;
                        }
                        return detectors[i].save.set(parameter, strValue);
                    }
                    if (child.equalsIgnoreCase("PulseTimer")) {
                        if (detectors[i].pulseTimer == null) {
                            return false;
                        }
                        return detectors[i].pulseTimer.set(parameter, strValue);
                    }
                    if (child.equalsIgnoreCase("PSF")) {
                        if (detectors[i].psf == null) {
                            return false;
                        }
                        return detectors[i].psf.set(parameter, strValue);
                    }
                    return false;
                }
            }
        }

        return false;

    }

    @Override
    public boolean addUser(ParamVector pv) {

        super.addUser(pv);

        if (finiteDifference != null){
            finiteDifference.addUser(pv);
        }

        if (monteCarlo != null){
            monteCarlo.addUser(pv);
        }

        return true;

    }

    @Override
    public boolean createVector(boolean close) {

        if (!super.createVector(false)) {
            return false;
        }

        if (finiteDifference != null){
            addBlankParam();
            finiteDifference.createVector(true);
            addVector(finiteDifference.getVector());
            finiteDifference.addUser(this);
        }

        if (monteCarlo != null){
            addBlankParam();
            monteCarlo.createVector(true);
            addVector(monteCarlo.getVector());
            monteCarlo.addUser(this);
        }

        if (close) {
            closeVector();
        }

        return true;

    }

    @Override
    public void updateVector(ParamObject[] v) {

        super.updateVector(v);

        if (finiteDifference != null){
            finiteDifference.updateVector(v);
        }

        if (monteCarlo != null){
            monteCarlo.updateVector(v);
        }

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {
        
        if (o == null) {
            return false;
        }
        
        String n = o.getName();

        if (!(o.paramVector instanceof Project)) {
            return false;
        }

        if (n.equalsIgnoreCase("dx")) {
            if (v <= 0) {
                return false;
            }
            dx = v;
            init();
            return true;
        }
        if (n.equalsIgnoreCase("dt")) {
            if (v <= 0) {
                return false;
            }
            stability = v * 3 * maxD / Math.pow(dx, 2.0);
            init();
            return true;
        }
        if (n.equalsIgnoreCase("simTime")) {
            if (v <= 0) {
                return false;
            }
            simTime = v;
            return true;
        }
        if (n.equalsIgnoreCase("maxT")) {
            if (v <= 0) {
                return false;
            }
            simTime = v;
            return true;
        }
        if (n.equalsIgnoreCase("stability")) {
            if (v <= 0) {
                return false;
            }
            stability = v;
            init();
            return true;
        }
        if (n.equalsIgnoreCase("saveRate")) {
            if (v <= 0) {
                return false;
            }
            setOutputRate(v);
            return true;
        }
        if (n.equalsIgnoreCase("saveDT")) {
            if (v <= 0) {
                return false;
            }
            setOutputRate(1/v);
            return true;
        }
        if (n.equalsIgnoreCase("printRate")) {
            if (v < 0) {
                return false;
            }
            setPrintRate(v);
            return true;
        }
        if (n.equalsIgnoreCase("printDT")) {
            if (v < 0) {
                return false;
            }
            setPrintRate(1/v);
            return true;
        }
        if (n.equalsIgnoreCase("batchNum")) {
            return setBatchNum((int) v);
        }
        return false;
    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, String s) {
        
        if ((o == null) || (s == null)) {
            return false;
        }

        if (!(o.paramVector instanceof Project)) {
            return false;
        }

        String n = o.getName();

        if (n.equalsIgnoreCase("name")) {
            name = s;
            Master.updateMainFrameTitle();
            return true;
        }
        if (n.equalsIgnoreCase("directory")) {
            directory = s;
            init(); // reset output file names
            return true;
        }
        if (n.equalsIgnoreCase("folder")) {
            folder = s;
            init(); // reset output file names
            return true;
        }
        if (n.equalsIgnoreCase("timeUnits")) {
            timeUnits = s;
            init();
            return true;
        }
        if (n.equalsIgnoreCase("spaceUnits")) {
            spaceUnits = s;
            init();
            return true;
        }
        if (n.equalsIgnoreCase("concUnits")) {
            concUnits = s;
            init();
            return true;
        }
        if (n.equalsIgnoreCase("date")) {
            if (s.length() == 0) {
                setDate();
            } else {
                date = s;
            }
            return true;
        }
        if (n.equalsIgnoreCase("seed")) {
            seed = s;
            return true;
        }
        return false;
    }
}
