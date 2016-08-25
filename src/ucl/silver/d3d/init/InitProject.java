package ucl.silver.d3d.init;

import ucl.silver.d3d.core.*;
import ucl.silver.d3d.utils.*;

/**
 * <p>Title: D3D</p>
 * <p>Description: 3D Diffusion-Reaction Simulator</p>
 * <p>Copyright: Copyright (c) 2016</p>
 * <p>Company: The Silver Lab at University College London</p>
 * @author Jason S. Rothman
 * @version 1.6
 */

public class InitProject extends ParamVector {

    Geometry geometry;
    CoordinatesVoxels coordinates; // general purpose coordinates

    public String[] initFuncList = {"initCube", "initCylinderCrank", "initSimpleMovie"};

    public InitProject(Project p) {
        super(p);
        geometry = project.geometry;
        coordinates = new CoordinatesVoxels(p);
        createVector(true);
    }

    @Override
    public boolean canEdit(String name) {
        return false;
    }

    public int initFunctionNum(String initSelect) {
        return D3DUtility.whichListItem(initFuncList, initSelect, true);
    }

    public boolean initFunction(String initSelect) {

        int i = initFunctionNum(initSelect);
        //i = 0;

        switch (i) {
            case 0:
                return initCube();
            case 1:
                return initCylinderCrank();
            case 2:
                return initSimpleMovie();
            default:
                Master.log("InitProject.initFunction error: failed to find init function " + initSelect);
                return true; // error
        }

    }

    public boolean initGeometry() {
        return false;
    }

    public boolean initCube() {

        int cubeWidth = 10;
        
        project.newFiniteDifference();

        project.name = "Simple Cube";
        project.directory = "/Jason/D3D/Simulations/";
        project.folder = "Testing";

        project.simTime = 1.0;
        project.dx = 0.05;
        project.stability = 0.4;

        geometry.resize(cubeWidth, cubeWidth, cubeWidth);
        geometry.clear();

        return false;

    }

    public boolean initCylinderCrank() { // Crank, The Mathematics of Diffusion, Eq 2.6

        double radius = Math.sqrt(1 / Math.PI);
        double diam = 2 * radius;
        double len = 20.02;
        double dx = 0.02;

        double z0, xv1, yv1, zv1, xv2, yv2, zv2;

        double source_c0 = 1.0;

        double distance_between_detectors = 0.5;

        double simTime = 1.5;

        boolean quarterSymmetry = true;

        int xVoxels = (int) (diam / dx);
        int yVoxels = xVoxels;
        int zVoxels = (int) (len / dx);
        int k, dk;
        
        project.newFiniteDifference();

        project.name = "Crank Cylinder Theoretical Solution";
        project.directory = "/Jason/D3D/Simulations/";
        project.folder = "Testing";

        project.simTime = simTime;
        project.dx = dx;
        project.stability = 0.4;
        project.set("printRate", 100);

        geometry.resize(xVoxels, yVoxels, zVoxels);
        geometry.cylinder();

        z0 = geometry.zVoxelCenter;

        xv1 = geometry.xVoxel1;
        yv1 = geometry.yVoxel1;
        zv1 = geometry.zVoxel1;
        xv2 = geometry.xVoxel2;
        yv2 = geometry.yVoxel2;
        zv2 = geometry.zVoxel2;

        if (quarterSymmetry) {
            xVoxels = (xVoxels - 2) / 2 + 2;
            yVoxels = (yVoxels - 2) / 2 + 2;
            geometry.resize(xVoxels, yVoxels, zVoxels);
        }

        int dnum = Master.addDiffusant("p", 0, 1);

        coordinates.setVoxels(xv1, yv1, z0, xv2, yv2, z0);
        Master.addSourceImpulse(dnum, coordinates, source_c0, 0.01);

        dk = (int) (distance_between_detectors / dx);

        k = (int) z0;

        for (int kcnt = 0; kcnt < 10; kcnt++) {
            coordinates.setVoxels(xv1, yv1, k, xv2, yv2, k);
            Master.addDetectorAvg(dnum, coordinates);
            k -= dk;
        }

        coordinates.setVoxels(xv1, yv1, zv1, xv2, yv2, zv2);
        Master.addDetectorAvg(dnum, coordinates); // total concentration should be step function

        return false;

    }

    public boolean initSimpleMovie() {

        double x0, y0, z0, xv1, yv1, zv1, xv2, yv2, zv2;

        double source_c0 = 1.0;

        double simTime = 2.0;

        int xVoxels = 141;
        int yVoxels = 30;
        int zVoxels = 1;

        Source s;
        
        project.newFiniteDifference();

        project.name = "Crank Cylinder Theoretical Solution";
        project.directory = "/Jason/D3D/Simulations/";
        project.folder = "Testing";

        project.simTime = simTime;
        project.dx = 0.02;
        project.stability = 0.4;
        project.set("printRate", 10);

        geometry.resize(xVoxels, yVoxels, zVoxels);
        geometry.clear();

        x0 = geometry.xVoxelCenter;
        y0 = geometry.yVoxelCenter;
        z0 = geometry.zVoxelCenter;

        xv1 = geometry.xVoxel1;
        yv1 = geometry.yVoxel1;
        zv1 = geometry.zVoxel1;
        xv2 = geometry.xVoxel2;
        yv2 = geometry.yVoxel2;
        zv2 = geometry.zVoxel2;

        int dnum = Master.addDiffusant("p", 0, 1);

        coordinates.setVoxels(x0, yv1, zv1, x0, yv2, zv2);
        s = Master.addSource(dnum, coordinates, source_c0, 0.2, 0.0003);
        s = Master.addSource(dnum, coordinates, source_c0, 0.7, 0.0003);
        s = Master.addSource(dnum, coordinates, source_c0, 1.2, 0.0003);

        return false;

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        String n = o.getName();

        if (!(o.paramVector instanceof InitProject)) {
            return false;
        }

        return false;

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, String s) {

        if ((o == null) || (s == null)) {
            return false;
        }

        if (!(o.paramVector instanceof InitProject)) {
            return false;
        }

        String n = o.getName();

        return false;

    }
    
}
