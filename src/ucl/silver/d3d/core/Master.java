package ucl.silver.d3d.core;

import ucl.silver.d3d.gui.*;
import ucl.silver.d3d.init.*;
import ucl.silver.d3d.utils.*;
import java.awt.*;
import java.text.SimpleDateFormat;
import java.util.Date;
import javax.swing.*;
import java.io.*;

/**
 * <p>Title: D3D</p>
 * <p>Description: 3D Diffusion-Reaction Simulator</p>
 * <p>Copyright: Copyright (c) 2016</p>
 * <p>Company: The Silver Lab at University College London</p>
 * @author Jason S. Rothman
 * @version 1.6
 */

public final class Master {

    public static final String D3Dversion = "1.6";

    public static Date startDate = null;

    public static Project project = null;
    public static MainFrame mainframe = null;

    private static Project[] projects = null;
    
    //private static String initClassAndFunction = "InitProject.initCube";
    //private static String initClassAndFunction = "InitProject.initCylinderCrank";
    //private static String initClassAndFunction = "InitProject.initSimpleMovie";

    //private static String initClassAndFunction = "InitFrapJason.init_FD_FRAP";
    //private static String initClassAndFunction = "InitFrapJason.init_FD_Axelrod";
    //private static String initClassAndFunction = "InitFrapJason.init_MC_FRAP";
    private static String initClassAndFunction = "InitFrapJason.init_MC_FRAP_Movie";
    //private static String initClassAndFunction = "InitFrapJason.init_MC_MSD";
    //private static String initClassAndFunction = "InitFrapJason.init_MC_Cichocki";

    //private static String initClassAndFunction = "InitAZEMJason.init_MC_AZ";

    public static String[] initProjectList = {"InitProject", "InitFrapJason", "InitAZEMJason"};
    public static String[] diffusantList = {"Diffusant", "DiffusantPhoto", "DiffusantVesicle", "DiffusantVesicles"};
    public static String[] detectorList = {"Detector", "DetectorWeighted"};
    public static String[] sourceList = {"Source"};
    
    public static String[] psfList = {"PSF", "PSFgauss", "PSFtorok", "PSFwilson"};

    public static String[] colorScaleList = {"BTC", "BTY", "Gray", "LinGray", "Heat", "Optimal", "LinOptimal", "Magenta", "Rainbow"};

    public static ColorScaleBTC colorScaleBTC = null;
    public static ColorScaleBTY colorScaleBTY = null;
    public static ColorScaleGray colorScaleGray = null;
    public static ColorScaleLinGray colorScaleLinGray = null;
    public static ColorScaleHeat colorScaleHeat = null;
    public static ColorScaleOptimal colorScaleOptimal = null;
    public static ColorScaleLinOptimal colorScaleLinOptimal = null;
    public static ColorScaleMagenta colorScaleMagenta = null;
    public static ColorScaleRainbow colorScaleRainbow = null;

    public static BufferedWriter bw1 = null;
    private static BufferedWriter bw2 = null;

    private static String[] logList = null;

    public static boolean writeToFiles = false; // allows writing to external file

    public static boolean foundMainStartUpArguments = false;

    private Master(){
        // cannot be instantiated
    }

    public static boolean initMainFrame(boolean packFrame) {

        mainframe = new MainFrame();

        if (packFrame) {
            mainframe.pack();
        } else {
            mainframe.validate();
        }

        // center the window
        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
        Dimension frameSize = mainframe.getSize();

        if (frameSize.height > screenSize.height) {
            frameSize.height = screenSize.height;
        }

        if (frameSize.width > screenSize.width) {
            frameSize.width = screenSize.width;
        }

        mainframe.setLocation((screenSize.width - frameSize.width) / 2,
                (screenSize.height - frameSize.height) / 2);

        mainframe.initTabs();

        updateMainFrameTitle();

        return true;

    }

    public static void updatePanel2D() {

        if ((mainframe == null) || (mainframe.panel2D == null)) {
            return;
        }

        mainframe.panel2D.updateControls();

    }

    public static Grid grid() {

        if ((mainframe == null) || (mainframe.panel2D == null)) {
            return null;
        }

        return mainframe.panel2D.grid2D;

    }

    public static void exit(String errorStr) {
        log("encountered fatal error: " + errorStr);
        System.exit(0);
    }

    public static void log(String message) {

        if (message == null) {
            return;
        }

        System.out.println(message);

        if (project != null) {
            writeToLogFile(message, false);
        }

        if (mainframe != null) {
            mainframe.panelLog.append(message, true);
        }

    }

    public void fatalError(String message, boolean exit) {

        if (message != null) {
            System.err.println(message + " (" + toString() + ")");
        }

        if (project != null) {
            writeToLogFile(message, false);
        }

        if (mainframe != null) {
            mainframe.panelLog.append(message, true);
        }

        if (exit) {

            if ((project != null) && (mainframe != null)) {
                mainframe.dispose();
            }

            System.exit(0);

        }

    }

    ///////////////////////////////////////////////////////////////////
    //
    //     Log Functions
    //
    ///////////////////////////////////////////////////////////////////

    public static boolean writeToLogFile(String s, boolean newFile) {

        String logFile;
        BufferedWriter bw;
        
        if (!writeToFiles) {
            return writeToLogList(s);
        }

        logFile = project.fullDirectory() + project.logFileName();

        try {

            if (newFile){
                bw = new BufferedWriter(new FileWriter(logFile)); // new file
            } else {
                bw = new BufferedWriter(new FileWriter(logFile, true)); // append
            }

            if (logList != null) {

                for (int i = 0; i < logList.length; i++) {

                    if ((logList[i] == null) || logList[i].equalsIgnoreCase("")) {
                        break;
                    }

                    bw.write(logList[i]); // write saved logs first
                    bw.newLine();

                }

                logList = null;

            }

            bw.write(s);
            bw.newLine();
            bw.close();
            bw = null;

        } catch (IOException e) {
            Master.exit("fatal error: unable to write to log file: " + logFile);
        }

        return true;

    }

    public static boolean writeToLogList(String s) {

        if (logList == null) {
            logList = new String[1000];
        }

        for (int i = 0; i < logList.length; i++) {
            if ((logList[i] == null) || logList[i].equalsIgnoreCase("")) {
                logList[i] = s;
                break;
            }
        }

        return true;

    }

    ///////////////////////////////////////////////////////////////////
    //
    //     Project Functions
    //
    ///////////////////////////////////////////////////////////////////

    public static int addProject(Project newProject) {

        int i = 0;

        if (projects != null) {
            i = projects.length;
        }

        Project[] newArray = new Project[i+1];

        if (i > 0) {
            System.arraycopy(projects, 0, newArray, 0, i);
        }

        newArray[i] = newProject;

        projects = newArray; // replace old array with new one

        Master.log("created Project #" + i + " : " + newProject.name);

        return i;

    }

    public static boolean killProject(int i) {

        int k = 0;
        String name = "";

        if ((i < 0) || (i >= projects.length)) {
            return false;
        }

        if (projects.length == 1) {
            projects = null;
            name = project.name;
            project = null;
            Master.log("killed Project #0: " + name);
            return true;
        }

        Project[] newArray = new Project[projects.length - 1]; // new array

        for (int j = 0; j < projects.length; j++) {

            if (j == i) {
                name = projects[j].name;
                continue;
            }

            newArray[k] = projects[j];
            k++;

        }

        projects = newArray; // replace old array with new one

        if (name.length() > 0) {
            Master.log("killed Project #" + i + ": " + name);
        }

        return true;

    }

    public static String[] projectList() {

        int numProjects = projects.length;

        String projectList[] = new String[numProjects];

        for (int i = 0; i < numProjects; i++) {
            if (projects[i] != null) {
                projectList[i] = projects[i].name;
            }
        }

        return projectList;

    }

    public static boolean projectSelect(String projectName) {

        if (projects == null) {
            return false;
        }

        if (project.name.equalsIgnoreCase(projectName)) {
            return false; // already selected
        }

        for (int i = 0; i < projects.length; i++) {
            if (projects[i] == null) {
                continue;
            }
            if (projects[i].name.equalsIgnoreCase(projectName)) {
                project = projects[i];
                mainframe.panelParams.setParamSelect(-1, 0);
                Master.log("changed to Project #" + i + " : " + projectName);
                return true;
            }
        }

        return false;

    }

    public static boolean projectSelect(int projectNum) {

        if ((projects == null) || (projectNum < 0) || (projectNum >= projects.length)) {
            return false;
        }

        if (projects[projectNum] == null) {
            return false;
        }

        project = projects[projectNum];

        mainframe.panelParams.setParamSelect(-1, 0);

        return true;

    }
    
    public static void createProject(String projectType, String projectName, String[][] startupArgsProject) {

        project = new Project(projectType);
        
        if ((projectName != null) && (projectName.length() > 0)) {
            project.name = projectName;
        }

        addProject(project);

        if (mainframe != null) {
            mainframe.panelParams.setParamSelect(-1, 0);
            mainframe.initTabs();
            //updateMainFrameTitle();
        }

    }

    public static void updateMainFrameTitle() {
        if (mainframe != null) {
            //mainframe.setTitle("D3D " + project.name);
            mainframe.setTitle(project.name);
        }
    }

    public static String initProjectClassStr(String initClassAndFunction) {

        if (initClassAndFunction == null) {
            return null;
        }

        String[] splitStr = initClassAndFunction.split("\\.");

        if (splitStr.length == 2) {
            return splitStr[0];
        }

        return null;

    }

    public static String initProjectFunctionStr(String initClassAndFunction) {

        if (initClassAndFunction == null) {
            return null;
        }

        String[] splitStr = initClassAndFunction.split("\\.");

        if (splitStr.length == 2) {
            return splitStr[1];
        }

        return null;

    }

    public static boolean initProjectClassAndFunction(String initClass, String initFunction) {

        String s;
        boolean error;

        String[] possibilities = initProjectList;

        ImageIcon icon = null;
        
        if (initClass.length() == 0) {

            s = (String) JOptionPane.showInputDialog(
                    mainframe,
                    "Choose Class:",
                    "Initialize Project",
                    JOptionPane.PLAIN_MESSAGE,
                    icon,
                    possibilities,
                    initProjectClassStr(initClassAndFunction));

            if ((s == null) || (s.length() == 0)) {
                return true; // error
            }

            initClass = s;

        }

        if (initProjectClass(initClass)) {
            return true; // error
        }

        possibilities = project.initProject.initFuncList;

        if (possibilities.length < 1) {
            return true; // error
        }
        
        if (initFunction.length() == 0) {

            s = (String) JOptionPane.showInputDialog(
                    mainframe,
                    "Choose Function:",
                    "Initialize Project",
                    JOptionPane.PLAIN_MESSAGE,
                    icon,
                    project.initProject.initFuncList,
                    possibilities[0]);

            if ((s == null) || (s.length() == 0)) {
                return true; // error
            }

            initFunction = s;

        }

        project.nullArrays();

        log("executing " + initClass + "." + initFunction + "...");

        error = project.initProject.initFunction(initFunction);

        if (error) {

            log("failed to execute " + initClass + "." + initFunction);

        } else {

            initClassAndFunction = initClass + "." + initFunction;
            project.initClassAndFunction = initClass + "." + initFunction;

            resetAllParamVectors();
            project.init();

            if (mainframe != null) {
                mainframe.panel2D.resetVoxelWidth();
                mainframe.panelParams.setParamSelect(-1, -1);
            }

            log("finished " + initClassAndFunction);

        }

        return error;

    }

    public static boolean initProject(String[][] startupArgs) {
        return Master.initProject(initClassAndFunction, startupArgs);
    }

    public static boolean initProject(String classAndFunction, String[][] startupArgsInitProject) {

        boolean error = false;

        String initClass = initProjectClassStr(classAndFunction);
        String initFunction = initProjectFunctionStr(classAndFunction);

        if (initProjectClass(initClass)) {
            return true; // error
        }

        project.nullArrays();

        log("executing " + classAndFunction + "...");

        error = project.initProject.initFunction(initFunction);

        if (error) {

            log("failed to execute " + classAndFunction);
            
        } else {

            initClassAndFunction = classAndFunction;
            project.initClassAndFunction = classAndFunction;

            log("finished " + classAndFunction);

        }

        if (mainframe != null) {
            mainframe.panel2D.resetVoxelWidth();
        }

        return error;

    }

    public static boolean initProjectClass(String initClass) {

        if (initClass.equalsIgnoreCase("InitProject")) {
            project.initProject = new InitProject(project);
        } else if (initClass.equalsIgnoreCase("InitFrapJason")) {
            project.initProject = new InitFrapJason(project);
        } else if (initClass.equalsIgnoreCase("InitAZEMJason")) {
            project.initProject = new InitAZEMJason(project);
        } else {
            log("Master.initProject error: failed to find init Class " + initClass);
            return true; // error
        }

        return false;

    }

    public static String currentDate() {

        Date date = new Date(System.currentTimeMillis());

        SimpleDateFormat dateFormat = new SimpleDateFormat("EEE, d MMM yyyy HH:mm:ss Z");

        return dateFormat.format(date);

    }

    ///////////////////////////////////////////////////////////////////
    //
    //     Export Text file functions
    //
    ///////////////////////////////////////////////////////////////////

    public static boolean writeParamVectorOpen(String directory, String fileName) {

        String fullName;

        if ((fileName == null) || (fileName.length() == 0)) {
            log("bad output parameter file name");
            return false;
        }

        fullName = directory + fileName + ".dat";

        //log("attempting to create output parameter file: " + fullName);

        try {

            bw1 = new BufferedWriter(new FileWriter(fullName));

            bw1.write("#D3D Parameter File v" + D3Dversion);
            bw1.newLine();
            bw1.write("#" + currentDate());
            bw1.newLine();

            log("created output parameter file: " + fullName);

            return true;

        } catch (IOException e) {
            log("unable to write file: " + fullName);
            return false;
        }

    }

    public static boolean writeParamVectorClose() {

        if (bw1 == null) {
            return false;
        }

        try {
            bw1.close();
            bw1 = null;
            return true;
        } catch (IOException e) {
            log("unable to close file: " + bw1.toString());
            return false;
        }

    }

    public static boolean writeParamVectorArray(Object[] pvArray) {

        ParamVector pv;

        if (bw1 == null) {
            return false;
        }

        if (pvArray == null) {
            return false;
        }

        for (int i = 0; i < pvArray.length; i++) {
            pv = (ParamVector) pvArray[i];
            writeParamVector(pv);
        }

        return true;

    }

    public static boolean writeParamVector(ParamVector pv) {

        String oname;
        String otype;

        if (bw1 == null) {
            return false;
        }

        if (pv == null) {
            return false;
        }

        ParamObject o;
        ParamObject[] v = pv.getVector();

        try {

            for (int i = 0; i < v.length; i++) {

                o = v[i];
                oname = o.getName();
                otype = o.getType();

                if (otype.compareToIgnoreCase("Class") == 0) {
                    bw1.newLine();
                    bw1.write(otype + "=" + o.getValue());
                } else if (oname.length() > 0) {
                    bw1.write(oname + "=" + o.getValue());
                } else {
                    continue;
                }

                bw1.newLine();

            }

            bw1.newLine();
            bw1.write("***************************************");
            bw1.newLine();

        } catch (IOException e) {
            log("unable to write to file: " + bw1.toString());
            return false;
        }

        return true;

    }

    ///////////////////////////////////////////////////////////////////
    //
    //     XML functions
    //
    ///////////////////////////////////////////////////////////////////

    public static boolean writeParamVectorXMLopen(String directory, String fileName) {

        String fullName1 = directory + fileName + ".xml";
        String fullName2 = directory + fileName + ".xsl";

        if ((fileName == null) || (fileName.length() == 0)) {
            return false;
        }

        try {

            bw1 = new BufferedWriter(new FileWriter(fullName1));

            bw1.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
            bw1.newLine();
            bw1.write("<?xml-stylesheet type=\"text/xsl\" href=\"" + fileName + ".xsl\"?>");
            bw1.newLine();
            bw1.write("<D3D version=\"" + D3Dversion + "\">");

            log("created output parameter file: " + fullName1);

        } catch (IOException e) {
            log("unable to write file: " + fullName1);
            return false;
        }

        try {

            bw2 = new BufferedWriter(new FileWriter(fullName2));

            bw2.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
            bw2.newLine();
            bw2.write("<xsl:stylesheet version=\"1.0\" xmlns:xsl=\"http://www.w3.org/1999/XSL/Transform\">");
            bw2.newLine();
            bw2.write("<xsl:template match=\"/\">");
            bw2.newLine();
            bw2.write("<html>");
            bw2.newLine();
            bw2.write("<body>");
            bw2.newLine();
            bw2.write("<center>");
            bw2.newLine();
            bw2.write("<h2>D3D Simulation Parameters</h2>");

            log("created output parameter file: " + fullName2);

        } catch (IOException e) {
            log("unable to write file: " + fullName2);
            return false;
        }

        return true;

    }

    public static boolean writeParamVectorXMLclose() {

        if ((bw1 == null) || (bw2 == null)){
            bw1 = null;
            bw2 = null;
            return false;
        }

        try {
            bw1.newLine();
            bw1.write("</D3D>");
            bw1.close();
            bw1 = null;
        } catch (IOException e) {
            log("unable to close file: " + bw1.toString());
            return false;
        }

        try {
            bw2.newLine();
            bw2.write("</center>");
            bw2.newLine();
            bw2.write("</body>");
            bw2.newLine();
            bw2.write("</html>");
            bw2.newLine();
            bw2.write("</xsl:template>");
            bw2.newLine();
            bw2.write("</xsl:stylesheet>");
            bw2.close();
            bw2 = null;
        } catch (IOException e) {
            log("unable to close file: " + bw2.toString());
            return false;
        }

        return true;

    }

    public static boolean writeParamVectorArrayXML(Object[] pvarray) {

        ParamVector pv;

        if ((bw1 == null) || (bw2 == null)){
            return false;
        }

        if (pvarray == null) {
            return false;
        }

        for (int i = 0; i < pvarray.length; i++) {
            pv = (ParamVector) pvarray[i];
            writeParamVectorXML(pv);
        }

        return true;

    }

    public static boolean writeParamVectorXML(ParamVector pv) {

        String oname;
        String otype;
        String classIndent = "  ";
        String[] classes = new String[10];
        String classTRcolor = "<tr bgcolor=\"#99ccff\">";

        boolean tableOpen = false;

        int classCounter = 0;

        classes[0]="D3D";

        if ((bw1 == null) || (bw2 == null)){
            return false;
        }

        if (pv == null) {
            return false;
        }

        ParamObject o;
        ParamObject[] v = pv.getVector();

        try {

            bw2.newLine();
            bw2.newLine();
            bw2.write("<table border=\"1\">");

            for (int i = 0; i < v.length; i++) {

                o = v[i];
                oname = o.getName();
                otype = o.getType();

                if (otype.compareToIgnoreCase("Class") == 0) {

                    if (tableOpen) {
                        bw2.newLine();
                        bw2.write("</table>");
                    }

                    tableOpen = true;

                    bw2.newLine();
                    bw2.newLine();
                    bw2.write("<table border=\"0\">");
                    bw2.newLine();
                    bw2.write(classTRcolor);
                    bw2.newLine();
                    bw2.write("<th width=\"300\"></th>");
                    bw2.newLine();
                    bw2.write("<th width=\"400\">" + o.getValue() + "</th>");
                    bw2.newLine();
                    bw2.write("<th width=\"100\"></th>");
                    bw2.newLine();
                    bw2.write("</tr>");

                    classTRcolor = "<tr bgcolor=\"#cccc99\">";

                    bw1.newLine();

                    for (int j = 0; j < classCounter+1; j++) {
                        bw1.write(classIndent);
                    }

                    bw1.write("<" + o.getValue() + ">"); // open new class
                    classCounter += 1;

                    classes[classCounter]=o.getValue();

                } else if (otype.compareToIgnoreCase("EndClass") == 0) {

                    bw1.newLine();
                    classCounter -= 1;

                    for (int j = 0; j < classCounter+1; j++) {
                        bw1.write(classIndent);
                    }

                    bw1.write("</" + o.getValue() + ">"); // close class

                    if (tableOpen){
                        bw2.newLine();
                        bw2.write("</table>");
                        tableOpen = false;
                    }

                } else if (oname.length() > 0) {

                    bw2.newLine();
                    bw2.write("<tr>");
                    bw2.newLine();
                    bw2.write("<td>" + oname + "</td>");
                    bw2.newLine();
                    bw2.write("<td>");
                    bw2.write("<xsl:value-of select=\"");

                    for (int j = 0; j < classCounter+1; j++) {
                        bw2.write(classes[j]+"/");
                    }

                    bw2.write(oname+"/@value");
                    bw2.write("\"/>");
                    bw2.write("</td>");
                    bw2.newLine();
                    bw2.write("<td>");
                    bw2.write("<xsl:value-of select=\"");

                    for (int j = 0; j < classCounter+1; j++) {
                        bw2.write(classes[j]+"/");
                    }

                    bw2.write(oname+"/@units");
                    bw2.write("\"/>");
                    bw2.write("</td>");
                    bw2.newLine();
                    bw2.write("</tr>");

                    bw1.newLine();

                    for (int j = 0; j < classCounter+1; j++) {
                        bw1.write(classIndent);
                    }

                    //bw1.write("<" + oname + " units=\"" + o.getUnits() + "\">" + o.getValue() + "</" + oname + ">");

                    bw1.write("<" + oname + " value=\"" + o.getValue() + "\" units=\"" + o.getUnits() + "\" />");

                }

            }

            bw2.newLine();
            bw2.newLine();
            bw2.write("</table>");
            bw2.newLine();
            bw2.write("<p></p>");

        } catch (IOException e) {
            return false;
        }

        return true;

    }

    // export 3D array to a text file
    public static void export3DArray(String file, CoordinatesVoxels c, double[][][] p) {

        int kcount = 50;

        if (p == null) {
            return;
        }

        DataOutputStream dout;
        StopWatch time = new StopWatch();

        if (file.length() == 0) {
            return;
        }

        try {

            dout = new DataOutputStream(new FileOutputStream(file));

            log("writing to file: " + file);

            time.start();

            for (int k = c.zVoxel1; k <= c.zVoxel2; k++) {

                if ((k > 0) && Math.IEEEremainder(k, kcount) == 0) {
                    time.stop();
                    log("writing k: " + k + ", t: " + time.toString());
                    time.start();
                }

                for (int j = c.yVoxel1; j <= c.yVoxel2; j++) {
                    for (int i = c.xVoxel1; i <= c.xVoxel2; i++) {
                        if ((i >= 0) && (i < p.length) && (j >= 0) && (j < p[0].length) &&
                                (k >= 0) && (k < p[0][0].length)) {
                            dout.writeBytes("\n");
                            dout.writeBytes(Integer.toString(k));
                            dout.writeBytes(",");
                            dout.writeBytes(Integer.toString(j));
                            dout.writeBytes(",");
                            dout.writeBytes(Integer.toString(i));
                            dout.writeBytes(",");
                            dout.writeBytes(Double.toString(p[i][j][k]));
                        }
                    }
                }
            }

            time.stop();

            dout.close();

            log("finished exporting PSF");

        } catch (IOException e) {
            log("unable to write file: " + file);
            return;
        }

    }

    ///////////////////////////////////////////////////////////////////
    //
    //     Diffusant Functions
    //
    ///////////////////////////////////////////////////////////////////

    public static String diffusantName() {
        return "Diffusant" + Integer.toString(project.numDiffusants());
    }

    public static int addDiffusant() {
        double C = 0;
        double D = 1;
        return addDiffusant(C, D);
    }

    public static int addDiffusant(double C, double D) {
        return addDiffusant(diffusantName(), C, D);
    }

    // create simple diffusant
    public static int addDiffusant(String name, double C, double D) {
        int dnum = project.addDiffusant(new Diffusant(project, name, C, D, null));
        updatePanel2D();
        return dnum;
    }

    ///////////////////////////////////////////////////////////////////
    //
    //     Source Functions
    //
    ///////////////////////////////////////////////////////////////////
    
    public static String sourceName() {
        return "Source" + Integer.toString(project.numSources());
    }

    public static Source addSource() {

        int diffusantNum = 0;
        double conc = 1;
        double t = 0;
        double duration = 1;

        CoordinatesVoxels c = new CoordinatesVoxels(project);

        c.matchDimensions(project.geometry);

        return addSource(diffusantNum, c, conc, t, duration);

    }

    public static Source addSource(int DiffusantNum, CoordinatesVoxels c, double conc) {
        return addSource(sourceName(), DiffusantNum, c, conc, null);
    }

    public static Source addSource(int DiffusantNum, CoordinatesVoxels c, double conc, double t,
            double duration) {
        return addSource(sourceName(), DiffusantNum, c, conc, new PulseTimer(project, t, duration));
    }

    public static Source addSource(int DiffusantNum, CoordinatesVoxels c, double conc, PulseTimer pt) {
        return addSource(sourceName(), DiffusantNum, c, conc, pt);
    }

    public static Source addSource(String name, int DiffusantNum, CoordinatesVoxels c, double conc,
            PulseTimer pt) {
        Source s = new Source(project, name, DiffusantNum, c, conc, pt);
        project.addSource(s);
        updatePanel2D();
        return s;
    }

    public static Source addSource(int DiffusantNum, CoordinatesVoxels c, PulseTimer pt, double qtotal) {
        return addSource(sourceName(), DiffusantNum, c, pt, qtotal);
    }

    public static Source addSource(String name, int DiffusantNum, CoordinatesVoxels c,
            PulseTimer pt, double qtotal) {
        Source s = new Source(project, name, DiffusantNum, c, pt, qtotal);
        project.addSource(s);
        updatePanel2D();
        return s;
    }

    public static Source addSourceImpulse() {

        int diffusantNum = 0;
        double conc = 1;
        double t = 0;

        CoordinatesVoxels c = new CoordinatesVoxels(project);

        c.matchDimensions(project.geometry);

        return addSourceImpulse(diffusantNum, c, conc, t);

    }

    public static Source addSourceImpulse(int DiffusantNum, CoordinatesVoxels c, double conc, double t) {
        return addSourceImpulse(sourceName(), DiffusantNum, c, conc, new PulseTimer(project, t));
    }

    public static Source addSourceImpulse(String name, int DiffusantNum, CoordinatesVoxels c,
            double conc, double t) {
        return addSourceImpulse(name, DiffusantNum, c, conc, new PulseTimer(project, t));
    }

    public static Source addSourceImpulse(int DiffusantNum, CoordinatesVoxels c, double conc,
            PulseTimer pt) {
        return addSourceImpulse(sourceName(), DiffusantNum, c, conc, pt);
    }

    public static Source addSourceImpulse(String name, int DiffusantNum, CoordinatesVoxels c,
            double conc, PulseTimer pt) {
        Source s = new Source(project, name, DiffusantNum, c, conc, pt);
        project.addSource(s);
        updatePanel2D();
        return s;
    }

    public static Source addSourceImpulse(int DiffusantNum, CoordinatesVoxels c, double t,
            int molecules) {
        return addSourceImpulse(sourceName(), DiffusantNum, c, new PulseTimer(project, t), molecules);
    }

    public static Source addSourceImpulse(String name, int DiffusantNum, CoordinatesVoxels c, double t,
            int molecules) {
        return addSourceImpulse(name, DiffusantNum, c, new PulseTimer(project, t), molecules);
    }

    public static Source addSourceImpulse(int DiffusantNum, CoordinatesVoxels c, PulseTimer pt,
            int molecules) {
        return addSourceImpulse(sourceName(), DiffusantNum, c, pt, molecules);
    }

    public static Source addSourceImpulse(String name, int DiffusantNum, CoordinatesVoxels c,
            PulseTimer pt, int molecules) {
        Source s = new Source(project, name, DiffusantNum, c, pt, molecules);
        project.addSource(s);
        updatePanel2D();
        return s;
    }

    ///////////////////////////////////////////////////////////////////
    //
    //     Detector Functions (Average, Weighted Average, Snapshot)
    //
    ///////////////////////////////////////////////////////////////////

    public static String detectorName() {
        return "Detector" + Integer.toString(project.numDetectors());
    }

    public static Detector addDetector() {
        int diffusantNum = 0;
        return addDetector(diffusantNum);
    }

    public static Detector addDetector(int DiffusantNum) {
        return addDetector(detectorName(), DiffusantNum);
    }

    // create average detector over entire shape
    public static Detector addDetector(String name, int DiffusantNum) {
        Detector dd = new Detector(project, name, DiffusantNum, null);
        project.addDetector(dd);
        updatePanel2D();
        return dd;
    }

    public static Detector addDetector(int DiffusantNum, CoordinatesVoxels c) {
        return addDetector(detectorName(), DiffusantNum, c);
    }

    // create average detector within coordinates
    public static Detector addDetector(String name, int DiffusantNum, CoordinatesVoxels c) {
        Detector dd = new Detector(project, name, DiffusantNum, c);
        project.addDetector(dd);
        updatePanel2D();
        return dd;
    }

    public static DetectorAvg addDetectorAvg() {
        int diffusantNum = 0;
        return addDetectorAvg(diffusantNum);
    }

    public static DetectorAvg addDetectorAvg(int DiffusantNum) {
        return addDetectorAvg(detectorName(), DiffusantNum);
    }

    // create average detector over entire shape
    public static DetectorAvg addDetectorAvg(String name, int DiffusantNum) {
        DetectorAvg dd = new DetectorAvg(project, name, DiffusantNum, null);
        project.addDetector(dd);
        updatePanel2D();
        return dd;
    }

    public static DetectorAvg addDetectorAvg(int DiffusantNum, CoordinatesVoxels c) {
        return addDetectorAvg(detectorName(), DiffusantNum, c);
    }

    // create average detector within coordinates
    public static DetectorAvg addDetectorAvg(String name, int DiffusantNum, CoordinatesVoxels c) {
        DetectorAvg dd = new DetectorAvg(project, name, DiffusantNum, c);
        project.addDetector(dd);
        updatePanel2D();
        return dd;
    }

    public static DetectorPSF addDetectorPSF(int DiffusantNum, PSF psf) {
        return addDetectorPSF(detectorName(), DiffusantNum, psf);
    }

    // create a weighted detector with PSF (over entire shape)
    public static DetectorPSF addDetectorPSF(String name, int DiffusantNum, PSF psf) {
        DetectorPSF dd = new DetectorPSF(project, name, DiffusantNum, null, psf);
        project.addDetector(dd);
        updatePanel2D();
        return dd;
    }

    public static DetectorPSF addDetectorPSF(int DiffusantNum, CoordinatesVoxels c, PSF psf) {
        return addDetectorPSF(detectorName(), DiffusantNum, c, psf);
    }

    // create a weighted detector with PSF
    public static DetectorPSF addDetectorPSF(String name, int DiffusantNum, CoordinatesVoxels c, PSF psf) {
        DetectorPSF dd = new DetectorPSF(project, name, DiffusantNum, c, psf);
        project.addDetector(dd);
        updatePanel2D();
        return dd;
    }

    public static DetectorPSF addDetectorPSF_Gauss() {
        int diffusantNum = 0;
        double xSTDV = 1;
        double ySTDV = 1;
        double zSTDV = 1;
        return addDetectorPSF_Gauss(diffusantNum, xSTDV, ySTDV, zSTDV);
    }

    public static DetectorPSF addDetectorPSF_Gauss(int DiffusantNum, double xSTDV, double ySTDV, double zSTDV) {
        return addDetectorPSF_Gauss(detectorName(), DiffusantNum, xSTDV, ySTDV, zSTDV);
    }

    public static DetectorPSF addDetectorPSF_Gauss(String name, int DiffusantNum, double xSTDV, double ySTDV, double zSTDV) {

        PSF psf = new PSFgauss(project, null, xSTDV, ySTDV, zSTDV);

        return addDetectorPSF(name, DiffusantNum, psf);

    }

    public static DetectorPSF addDetectorPSF_Gauss(int DiffusantNum, double xSTDV, double ySTDV, double zSTDV,
            double i0, double j0, double k0) {
        return addDetectorPSF_Gauss(detectorName(), DiffusantNum, xSTDV, ySTDV, zSTDV, i0, j0, k0);
    }

    public static DetectorPSF addDetectorPSF_Gauss(String name, int DiffusantNum, double xSTDV, double ySTDV, double zSTDV,
            double i0, double j0, double k0) {

        PSF psf = new PSFgauss(project, null, xSTDV, ySTDV, zSTDV);

        psf.setVoxelCenter(i0, j0, k0);

        return addDetectorPSF(name, DiffusantNum, psf);

    }

    public static DetectorPSF addDetectorPSF_Gauss(int DiffusantNum, CoordinatesVoxels c,
            double xSTDV, double ySTDV, double zSTDV) {
        return addDetectorPSF_Gauss(detectorName(), DiffusantNum, c, xSTDV, ySTDV, zSTDV);
    }

    public static DetectorPSF addDetectorPSF_Gauss(String name, int DiffusantNum, CoordinatesVoxels c,
            double xSTDV, double ySTDV, double zSTDV) {

        PSF psf = new PSFgauss(project, c, xSTDV, ySTDV, zSTDV);

        return addDetectorPSF(name, DiffusantNum, c, psf);

    }

    ///////////////////////////////////////////////////////////////////
    //
    //     Photolysis Functions (i.e. MNI uncaging or FRAP, PSFs)
    //
    ///////////////////////////////////////////////////////////////////
    
    public static boolean photoExists() {
        Diffusant[] diffusant = project.diffusants;
        if (diffusant != null) {
            for (int i = 0; i < diffusant.length; i++) {
                if (diffusant[i] instanceof DiffusantPhoto) {
                    return true;
                }
            }
        }
        return false;
    }

    public static DiffusantPhoto getPhoto() {
        Diffusant[] diffusant = project.diffusants;
        if (diffusant != null) {
            for (int i = 0; i < diffusant.length; i++) {
                if (diffusant[i] instanceof DiffusantPhoto) {
                    return (DiffusantPhoto) diffusant[i];
                }
            }
        }
        return null;
    }

    ///////////////////////////////////////////////////////////////////
    //
    //     Parameter Vector Functions
    //
    ///////////////////////////////////////////////////////////////////

    public static void updateAllParamVectors() {

        if (project != null) {

            System.out.println("updateAllParamVectors");

            project.updateVectors();

            if (project.geometry != null) {
                project.geometry.updateVectors();
            }

            updateParamVectorArray(project.diffusants);
            updateParamVectorArray(project.sources);
            updateParamVectorArray(project.detectors);
            updateParamVectorArray(project.batches);
            updateParamVectorArray(project.errors);

        }

    }

    public static void updateParamVectorArray(Object[] pv) {

        ParamVector temp;

        if (pv != null) {
            for (int i = 0; i < pv.length; i++) {
                if (pv[i] instanceof ParamVector) {
                    temp = (ParamVector) pv[i];
                    temp.updateVectors();
                }
            }
        }

    }

    public static void resetAllParamVectors() {

        if (project != null) {

            project.clearVector();
            project.createVector(true);

            if (project.geometry != null) {
                project.geometry.clearVector();
                project.geometry.createVector(true);
            }

            resetParamVectorArray(project.diffusants);
            resetParamVectorArray(project.sources);
            resetParamVectorArray(project.detectors);
            resetParamVectorArray(project.batches);
            resetParamVectorArray(project.errors);
            
        }

    }

    public static void resetParamVectorArray(Object[] pv) {

        ParamVector temp;

        if (pv != null) {
            for (int i = 0; i < pv.length; i++) {
                if (pv[i] instanceof ParamVector) {
                    temp = (ParamVector) pv[i];
                    temp.clearVector();
                    temp.createVector(true);
                }
            }
        }

    }

    // get file name
    public static String getFileName(boolean DIRECTORIES_ONLY) {

        String fileName = "";

        JFileChooser fc = new JFileChooser();

        fc.setSelectedFile(new File(""));
        
        if (DIRECTORIES_ONLY) {
            fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        }

        int approve = fc.showOpenDialog(mainframe);

        if (approve == JFileChooser.APPROVE_OPTION) {
            fileName = fc.getCurrentDirectory() + "\\" + fc.getSelectedFile().getName();
        }

        fileName = fileName.substring(2);
        fileName = fileName.replace('\\', '/');

        return fileName;

    }
    
    public static String promptForInput(String[] promptList, String promptStr, String promptTitle, String promptSelect) {

        String[] possibilities = promptList;

        ImageIcon icon = null;

        String s = (String) JOptionPane.showInputDialog(
                mainframe,
                promptStr,
                promptTitle,
                JOptionPane.PLAIN_MESSAGE,
                icon,
                possibilities,
                promptSelect);

        return s;
        
    }

}


