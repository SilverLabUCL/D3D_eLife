package ucl.silver.d3d.init;

import ucl.silver.d3d.core.*;
import java.util.Random;

/**
 * <p>Title: D3D</p>
 * <p>Description: 3D Diffusion-Reaction Simulator</p>
 * <p>Copyright: Copyright (c) 2016</p>
 * <p>Company: The Silver Lab at University College London</p>
 * @author Jason S. Rothman
 * @version 1.6
 */

public final class InitGeometryMito extends ParamVector {
    
    public double xdim = 8;
    public double ydim = 8;
    public double zdim = 10;
    public boolean eighth_geometry = false;
    
    public double xvcenter, yvcenter, zvcenter;
    
    public boolean mitoNonSpaceVoxels = false;
    public double mitoRadius = (0.25 / 0.89) / 2.0; // um Palay
    public double mitoAxialRatio = 2.0 / 0.25; // Palay
    public double mitoVolumeFraction = 0.28; // Zoltan average
    public double mitoVolumeFractionTolerance = 0.02;
    public double mitoVolumeFractionActual = 0;
    public boolean mitoVolumeFractionExact = false;
    public boolean mitoCylinder = true;
    public int mito_ijkselect = -1;
    
    public long seed = 8682522807148012L + System.nanoTime();
    
    private final Random ran = new Random(seed);

    public InitGeometryMito(Project p) {
        super(p);
        createVector(true);
    }
    
    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("xdim")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("ydim")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("zdim")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("mitoRadius")) {
            return project.spaceUnits;
        }
        return super.units(name);
    }

    @Override
    public boolean canEdit(String name) {
        return false;
    }

    public boolean initGeometryMC(int iPSFselect) {

        boolean error = false;

        double xdim1 = xdim;
        double ydim1 = ydim;
        double zdim1 = zdim;

        Geometry geometry = project.geometry;

        geometry.forceEvenVoxels = true;

        if (eighth_geometry) {

            xdim1 = xdim / 2;
            ydim1 = ydim / 2;

            if (iPSFselect < 2) {
                zdim1 = zdim / 2;
            }

        }

        geometry.resize(xdim1, ydim1, zdim1);

        project.geometry.clear();

        if (mitoNonSpaceVoxels) {
            error = initGeometryMitochondria();
        }

        xvcenter = geometry.xVoxelCenter;
        yvcenter = geometry.yVoxelCenter;
        zvcenter = geometry.zVoxelCenter;

        if (eighth_geometry) {

            xvcenter = geometry.xVoxels - 0.5;
            yvcenter = geometry.yVoxels - 0.5;

            if (iPSFselect < 2) {
                zvcenter = geometry.zVoxels - 0.5;
            }

        }

        geometry.update();

        return error;

    }

    private boolean initGeometryMitochondria() {

        int maxMitoTrials = 500, trialCount = 0;

        boolean error = false;

        if (mitoVolumeFraction <= 0 ) {
            return false;
        }

        Master.log("adding mitochondria...");

        for (int i = 0; i < maxMitoTrials; i++) {

            error = initGeometryMitochondria2();

            trialCount++;

            if (!error) {
                break;
            }

        }

        if (error) {
            Master.log("failed to add mitochondria, trials = " + trialCount);
            return true;
        }

        Master.log("mitochondria volume fraction = " + mitoVolumeFractionActual + " percent, trials = " + trialCount);

        return false;

    }

    private boolean initGeometryMitochondria2() {

        double dran, volumeFraction = 0;

        int linkNum = 3;

        RunMonteCarloAZ mc;

        Geometry geometry = project.geometry;

        CoordinatesVoxels c = new CoordinatesVoxels(project, geometry);

        c.setVoxels(geometry.xVoxel1, geometry.yVoxel1, geometry.zVoxel1, geometry.xVoxel2, geometry.yVoxel2, geometry.zVoxel2);

        project.geometry.clear();

        if (mito_ijkselect == -2) {

            dran = ran.nextDouble();

            if (dran < 0.3333333333) {
                mito_ijkselect = 0; // xy plane
            } else if (dran < 0.6666666666) {
                mito_ijkselect = 1; // yz plane
            } else {
                mito_ijkselect = 2; // zx plane
            }

        }

        if (mitoVolumeFractionExact) {
            volumeFraction = GeometryTools.addEllipsoids(geometry, c, mitoRadius, mitoAxialRatio, mitoVolumeFraction, -1, mitoVolumeFractionExact, mito_ijkselect, mitoCylinder, linkNum);
        } else {
            volumeFraction = GeometryTools.addEllipsoids(geometry, c, mitoRadius, mitoAxialRatio, mitoVolumeFraction - mitoVolumeFractionTolerance, -1, mitoVolumeFractionExact, mito_ijkselect, mitoCylinder, linkNum);
        }

        if ((volumeFraction < mitoVolumeFraction - mitoVolumeFractionTolerance) || (volumeFraction > mitoVolumeFraction + mitoVolumeFractionTolerance)) {
            return true; // error
        }

        if (project.monteCarlo instanceof RunMonteCarloAZ) {

            mc = (RunMonteCarloAZ) project.monteCarlo;

            double azWidthVoxels = mc.azWidth / project.dx;
            double azHeightVoxels = mc.azHeight / project.dx;
            int extra = 0;//1;

            int i1 = (int) Math.round(geometry.xVoxelCenter - (azWidthVoxels/2) - extra);
            int i2 = (int) Math.round(geometry.xVoxelCenter + (azWidthVoxels/2) + extra);
            int j1 = i1;
            int j2 = i2;
            int k1 = (int) Math.round(geometry.zVoxel2 - azHeightVoxels - extra);
            int k2 = geometry.zVoxel2;

            for (int k = k1; k <= k2; k++) {
                for (int j = j1; j <= j2; j++) {
                    for (int i = i1; i <= i2; i++) {
                        if (geometry.space[i][j][k] == -1) {
                            return true; // error, non-space voxel overlaps with active zone
                        }
                    }
                }
            }

        }

        //shape.checkSpace();

        //mitoVolumeFractionActual = (100.0 - 100.0 * shape.spaceVoxels / shape.voxels); // does not work for spherical shape
        mitoVolumeFractionActual = volumeFraction;

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
