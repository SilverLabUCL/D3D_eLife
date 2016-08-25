package ucl.silver.d3d.core;

/**
 * <p>Title: D3D</p>
 * <p>Description: 3D Diffusion-Reaction Simulator</p>
 * <p>Copyright: Copyright (c) 2016</p>
 * <p>Company: The Silver Lab at University College London</p>
 * @author Jason S. Rothman
 * @version 1.6
 */

public class PSFgauss extends PSF {

    public double xSTDV, ySTDV, zSTDV; // Gaussian standard deviations
    public double xFWHM, yFWHM, zFWHM; // Gaussian full-width half-max

    private final static double FWHMconversion = 2.0 * Math.sqrt(2.0 * Math.log(2.0)); // 2.35482

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("xSTDV")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("ySTDV")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("zSTDV")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("xFWHM")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("yFWHM")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("zFWHM")) {
            return project.spaceUnits;
        }
        return super.units(name);
    }

    public PSFgauss(Project p, CoordinatesVoxels c, double STDVx, double STDVy, double STDVz) {
        super(p, c);
        xSTDV = STDVx;
        ySTDV = STDVy;
        zSTDV = STDVz;
        xFWHM = computeFWHM(xSTDV);
        yFWHM = computeFWHM(ySTDV);
        zFWHM = computeFWHM(zSTDV);
        //xySymmetric = true; // should manually set to true
        name = "Gaussian PSF";
        createVector(true);
    }

    public static double computeSTDV(double FWHM) {
        if (FWHM < 0) {
            return -1;
        } else {
            return FWHM / FWHMconversion; // http://hyperphysics.phy-astr.gsu.edu/hbase/math/gaufcn2.html
        }
    }

    public static double computeFWHM(double STDV) {
        if (STDV < 0) {
            return -1;
        } else {
            return STDV * FWHMconversion; // Mathworld
        }
    }

    @Override
    public double computeVoxel(double xVoxel, double yVoxel, double zVoxel) {

        double x, y, z;
        double v1 = 1, v2 = 1, v3 = 1;
        double e1 = 1, e2 = 1, e3 = 1;

        dx = project.dx;
        
        xVoxel -= xVoxelCenter();
        yVoxel -= yVoxelCenter();
        zVoxel -= zVoxelCenter();

        rotateAll(xVoxel, yVoxel, zVoxel); // results saved in irotated, jrotated, krotated

        x = iRotated() * dx; // um
        y = jRotated() * dx; // um
        z = kRotated() * dx; // um

        if (xSTDV > 0) {
            v1 = Math.pow(x, 2.0);
            e1 = Math.exp(-v1 / (2.0 * Math.pow(xSTDV, 2.0)));
        }

        if (ySTDV > 0) {
            v2 = Math.pow(y, 2.0);
            e2 = Math.exp(-v2 / (2.0 * Math.pow(ySTDV, 2.0)));
        }

        if (zSTDV > 0) {
            v3 = Math.pow(z, 2.0);
            e3 = Math.exp(-v3 / (2.0 * Math.pow(zSTDV, 2.0)));
        }

        return e1 * e2 * e3;

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        if (super.setMyParams(o, v)) {
            return true;
        }

        String n = o.getName();

        if (n.equalsIgnoreCase("xSTDV")) {
            xSTDV = v;
            xFWHM = computeFWHM(xSTDV);
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("ySTDV")) {
            ySTDV = v;
            yFWHM = computeFWHM(ySTDV);
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("zSTDV")) {
            zSTDV = v;
            zFWHM = computeFWHM(zSTDV);
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("xFWHM")) {
            xFWHM = v;
            xSTDV = computeSTDV(xFWHM);
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("yFWHM")) {
            yFWHM = v;
            ySTDV = computeSTDV(yFWHM);
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("zFWHM")) {
            zFWHM = v;
            zSTDV = computeSTDV(zFWHM);
            array = null;
            return true;
        }
        return false;
    }
}
