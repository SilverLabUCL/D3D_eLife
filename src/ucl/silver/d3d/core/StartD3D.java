package ucl.silver.d3d.core;

import javax.swing.UIManager;
import ucl.silver.d3d.utils.*;
import java.util.Date;

/**
 * <p>Title: D3D</p>
 * <p>Description: 3D Diffusion-Reaction Simulator</p>
 * <p>Copyright: Copyright (c) 2016</p>
 * <p>Company: The Silver Lab at University College London</p>
 * @author Jason S. Rothman
 * @version 1.6
 */

public class StartD3D {

    public StartD3D(String[] args) {

        boolean showGUI = true;
        boolean runSimulation = false;
        boolean packFrame = false;

        Master.startDate = new Date(System.currentTimeMillis());

        Master.log("D3D started " + Master.startDate.toString());

        String projectType = "FD";

        String initClassAndFunction = "";

        String[][] startupArgs = null;
        String[][] startupArgsProject = null;
        String[][] startupArgsInitProject = null;

        startupArgsProject = getClassStartupArgs(startupArgs, "Project");
        startupArgsInitProject = getClassStartupArgs(startupArgs, "InitProject");

        Master.createProject(projectType, "", startupArgsProject);

        if (showGUI) {
            Master.initMainFrame(packFrame);
        }

        if (initClassAndFunction.length() > 0) {
            if (Master.initProject(initClassAndFunction, startupArgsInitProject)) {
                System.exit(0);
            }
        } else if (Master.initProject(startupArgsInitProject)) {
            System.exit(0);
        }
        
        Master.resetAllParamVectors();
        Master.project.init();

        if (showGUI && (Master.mainframe != null)) {

            if (Master.mainframe.panelParams != null) {
                if (Master.project.errors != null) {
                    Master.mainframe.panelParams.setParamSelect(6, -1); // error
                } else {
                    Master.mainframe.panelParams.setParamSelect(1, -1);
                }
                Master.mainframe.panelParams.updateControls();
            }

            if (Master.mainframe.panel2D != null) {
                Master.mainframe.panel2D.resetVoxelWidth();
            }

            Master.mainframe.setVisible(true);

        }

        if (runSimulation) {
            //System.out.println("starting project simulation...");
            Master.project.simulationStart(false);
        }

    }

    private static String[][] getClassStartupArgs(String[][] startupArgs, String className) {

        int j = 0;
        String parent, parameterName;

        if ((startupArgs == null) || (startupArgs.length == 0)) {
            return null;
        }

        String[][] s = new String[startupArgs.length][2];

        for (int i = 0; i < startupArgs.length; i++) {

            parent = D3DUtility.parentClass(startupArgs[i][0]);
            parameterName = D3DUtility.parameterName(startupArgs[i][0]);

            if ((parent == null) || (parameterName == null)) {
                continue;
            }

            if (parent.equalsIgnoreCase(className)) {
                s[j][0] = startupArgs[i][0];
                s[j][1] = startupArgs[i][1];
                startupArgs[i][0] = null;
                startupArgs[i][1] = null;
                j++;
                
            }
            
        }

        return s;

    }

    private static void printUsageAndQuit() {
        System.out.println("\nUsage: java  ucl.silver.d3d.core.StartD3D [-options]");
        System.out.println("where options include: \n");
        System.out.println("-? -help        print this help message");
        System.out.println("-MonteCarlo     run Monte Carlo simulator (default finite difference)");
        System.out.println("-noGUI          run without GUI");
        System.out.println("-run            start simulation");
        System.out.println("-init Class.function   init class and function (-init InitProject.initCube");
        System.out.println("-name value     parameter name value pairs (-Project.simTime 2.5 -Geometry.xLength 0.8 -MonteCarlo.MSD.save2TextFile 1)");
        System.exit(0);
    }

    // main method, everything begins here
    public static void main(String[] args) {

        try {
            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
        } catch (Exception e) {
            //e.printStackTrace();
        }

        StartD3D d3d = new StartD3D(args);

    }

}
