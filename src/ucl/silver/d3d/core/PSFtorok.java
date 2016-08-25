package ucl.silver.d3d.core;

import ucl.silver.d3d.utils.*;

/**
 * <p>Title: D3D</p>
 * <p>Description: 3D Diffusion-Reaction Simulator</p>
 * <p>Copyright: Copyright (c) 2016</p>
 * <p>Company: The Silver Lab at University College London</p>
 * @author Jason S. Rothman
 * @version 1.6
 */

public class PSFtorok extends PSF {

    public int integrationSteps = 1000; // numerical integration steps
    public double numericalAperture = 1.0; // numerical aperture
    public double waveLength = 0.351; // light wavelength
    public double refractiveIndex = 1.338; // refractive index of medium
    public double ab1, ab2, ab3, ab4, ab5, ab6, ab7, ab8, ab9, ab10, abZoff, zFxnOffset; // aberration coefficients

    public double rscale = 1.0;

    private double alpha, dtheta, tk0, tk1;
    private Complex c_dtheta, c_2, i0, i1, i2, I0, I1, I2;

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("waveLength")) {
            return project.spaceUnits;
        }
        return super.units(name);
    }

    public PSFtorok(Project p, CoordinatesVoxels c) {
        super(p, c);
        //xySymmetric = true; // should manually set to true
        name = "Torok PSF";
        createVector(true);
    }

    public void setFitZoltan3() { // this is best fit to Zoltan's FRAP PSF
        waveLength = 0.488;
        numericalAperture = 0.90643;
        ab1 = -0.15507;
        ab2 = 13.488;
        ab3 = -154.08;
        ab4 = 839.65;
        ab5 = -2324.9;
        ab6 = 2568.2;
        ab7 = 1049.8;
        ab8 = -3276.1;
        abZoff = 0;
    }

    public void check(boolean recompute) {
        if (recompute || (array == null)) {
            compute();
        }
    }

    @Override
    public double getValue(int i, int j, int k) {
        if (exists()) {
            return getArrayValue(i, j, k); // already computed, get from 3D array
        } else {
            return computeValue(i, j, k, normalize); // else compute the value
        }
    }

    @Override
    public void init() {

        alpha = Math.asin(numericalAperture / refractiveIndex);
        dtheta = alpha / integrationSteps;
        tk0 = 2 * Math.PI / waveLength;
        tk1 = 2 * Math.PI * refractiveIndex / waveLength;
        c_dtheta = new Complex(dtheta, 0);
        c_2 = new Complex(2.0, 0);

        super.init();

    }

    @Override
    public double computeVoxel(double i, double j, double k) {

        double ii, jj, kk;
        double r, z, phi, rxz, axz, ryz, ayz, zoff = abZoff;
        double rotatexz = 0;
        double rotateyz = 0;

        dx = project.dx;

        //rotatexz = 2 * Math.PI / 8; // 45 degrees
        //rotateyz = 2 * Math.PI / 8; // 45 degrees

        ii = i - xVoxelCenter();
        jj = j - yVoxelCenter();
        kk = k - zVoxelCenter();

        if (rotatexz > 0) {

            // convert to cylindrical coordinates
            rxz = Math.sqrt(Math.pow(ii, 2) + Math.pow(kk, 2));
            axz = Math.atan2(kk, ii);

            axz += rotatexz; // rotate 45 degrees in xz plane

            // now convert back
            ii = rxz * Math.cos(axz);
            kk = rxz * Math.sin(axz);

        }

        if (rotateyz > 0) {

            // convert to cylindrical coordinates
            ryz = Math.sqrt(Math.pow(jj, 2) + Math.pow(kk, 2));
            ayz = Math.atan2(kk, jj);

            ayz += rotateyz; // rotate 45 degrees in yz plane

            // now convert back
            jj = ryz * Math.cos(ayz);
            kk = ryz * Math.sin(ayz);

        }

        r = Math.sqrt(Math.pow(ii, 2) + Math.pow(jj, 2)) * dx;

        z = (kk * dx) - zoff;

        phi = getPhi(ii, jj);

        integrate_Torok(r * rscale, z);

        return densityE(I0, I1, I2, phi);

    }

    double getPhi(double i, double j) {
        double phi = 0;

        if (i == 0) {
            if (j == 0) {
                return 0;
            } else if (j > 0) {
                return (Math.PI / 2.0);
            } else if (j < 0) {
                return (3 * Math.PI / 2.0);
            }
        } else if (i > 0) {
            return Math.atan(j / i);
        } else if (i < 0) {
            return (Math.PI + Math.atan(j / i));
        }

        return phi;

    }

    private void integrate_Torok(double r, double z) { // trapezoid integration

        Complex c0, c1, c2;

        I0 = new Complex(0, 0);
        I1 = new Complex(0, 0);
        I2 = new Complex(0, 0);

        Torok(r, z, 0);

        c0 = i0;
        c1 = i1;
        c2 = i2;

        for (int itheta = 0; itheta < integrationSteps; itheta += 1) {

            Torok(r, z, itheta + 1);

            c0 = c0.plus(i0);
            c0 = c0.div(c_2);

            c1 = c1.plus(i1);
            c1 = c1.div(c_2);

            c2 = c2.plus(i2);
            c2 = c2.div(c_2);

            I0 = I0.plus(c0);
            I1 = I1.plus(c1);
            I2 = I2.plus(c2);

            c0 = i0;
            c1 = i1;
            c2 = i2;

        }

    }

    // basic computation of Torok PSF
    private void Torok(double r, double z, int itheta) {
        double theta, sinT, cosT, psi;

        Complex arg, arg1, arg2;

        theta = itheta * dtheta;

        sinT = Math.sin(theta);
        cosT = Math.cos(theta);

        //taus = 2 * sinT * cosT / (Math.sin(2 * theta));
        //taup = 2 * sinT * cosT / (Math.sin(2 * theta) * Math.cos(0));

        psi = ab1 * Math.pow(sinT, 2) +
                ab2 * Math.pow(sinT, 4) +
                ab3 * Math.pow(sinT, 6) +
                ab4 * Math.pow(sinT, 8) +
                ab5 * Math.pow(sinT, 10) +
                ab6 * Math.pow(sinT, 12) +
                ab7 * Math.pow(sinT, 14) +
                ab8 * Math.pow(sinT, 16) +
                ab9 * Math.pow(sinT, 18) +
                ab10 * Math.pow(sinT, 20);

        // from Sheppard & Torok, 1996

        arg = new Complex(Math.sqrt(cosT) * sinT, 0);
        arg1 = new Complex(0, tk0 * psi);
        arg1 = arg1.exp();
        arg2 = new Complex(0, tk1 * z * cosT);
        arg2 = arg2.exp();
        arg = c_dtheta.times(arg.times(arg1.times(arg2)));

        arg1 = new Complex((1 + 1 * cosT) * MoreMath.jn(0, tk1 * r * sinT), 0);
        i0 = arg.times(arg1);

        arg1 = new Complex((1 * sinT) * MoreMath.jn(1, tk1 * r * sinT), 0);
        i1 = arg.times(arg1);

        arg1 = new Complex((1 - cosT) * MoreMath.jn(2, tk1 * r * sinT), 0);
        i2 = arg.times(arg1);

    }

    double densityEB(Complex I0, Complex I1, Complex I2) {

        return (Math.pow(I0.mod(), 2) + 2 * Math.pow(I1.mod(), 2) +
                Math.pow(I2.mod(), 2));

    }

    double densityE(Complex I0, Complex I1, Complex I2, double phi) {

        double d;

        Complex c = I2.conj();

        c = c.times(I0);

        d = Math.pow(I0.mod(), 2) + Math.pow(I2.mod(), 2);

        d += 2 * (1 + Math.cos(2 * phi)) * Math.pow(I1.mod(), 2);

        d += 2 * Math.cos(2 * phi) * c.real();

        return d;

    }

    double PoyntingVector(Complex I0, Complex I1, Complex I2) {

        double pp, q;

        pp = (Math.pow(I0.mod(), 2) - Math.pow(I2.mod(), 2)) / 2;

        Complex cc = I2.conj();
        cc = cc.minus(I0.conj());
        cc = I1.times(cc);

        q = cc.imag();

        return (Math.sqrt(Math.pow(pp, 2) + Math.pow(q, 2)));

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

        if (n.equalsIgnoreCase("integrationSteps")) {
            if (v < 5) {
                return false;
            }
            integrationSteps = (int) v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("numericalAperture")) {
            if (v < 0) {
                return false;
            }
            numericalAperture = v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("waveLength")) {
            if (v < 0) {
                return false;
            }
            waveLength = v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("refractiveIndex")) {
            if (v < 0) {
                return false;
            }
            refractiveIndex = v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("ab1")) {
            ab1 = v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("ab2")) {
            ab2 = v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("ab3")) {
            ab3 = v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("ab4")) {
            ab4 = v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("ab5")) {
            ab5 = v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("ab6")) {
            ab6 = v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("ab7")) {
            ab7 = v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("ab8")) {
            ab8 = v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("ab9")) {
            ab9 = v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("ab10")) {
            ab10 = v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("abZoff")) {
            abZoff = v;
            array = null;
            return true;
        }
        return false;
    }
}
