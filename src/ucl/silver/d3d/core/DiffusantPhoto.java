package ucl.silver.d3d.core;

/**
 * <p>Title: D3D</p>
 * <p>Description: 3D Diffusion-Reaction Simulator</p>
 * <p>Copyright: Copyright (c) 2016</p>
 * <p>Company: The Silver Lab at University College London</p>
 * @author Jason S. Rothman
 * @version 1.6
 */

public class DiffusantPhoto extends Diffusant {

    // [a] ---> [b]

    // this Diffusant is [a], and gets converted to [b] at rate kPhoto

    public String productName = null; // diffusant name for product [b]
    public int productNum = -1; // diffusant number for product [b]
    // if productNum = -1, then product is not saved

    public double kPhoto; // k of photolysis reaction (1/ms)
    private double kPhotoSave; // variable for saving kPhoto to output file

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("kPhoto")) {
            return "1/" + project.timeUnits;
        }
        return super.units(name);
    }

    public DiffusantPhoto(Project p, String NAME, double InitialConcentration, double diffusionConstant, CoordinatesVoxels c,
            int ProductNum, PulseTimer pt, PSF PSF, double Kphoto) {
        
        super(p, NAME, InitialConcentration, diffusionConstant, c);
        
        productNum = ProductNum;
        pulseTimer = pt;
        psf = PSF;
        reaction = true;
        kPhoto = Kphoto;

        createVector(true);

    }

    @Override
    public void init() {

        super.init();

        productName = project.diffusantName(productNum);

        if ((productNum >= 0) && (!project.checkDiffusantNum(productNum))) {
            error("init", "productNum", "out of range");
        }

    }

    @Override
    public void react(RunFiniteDifference fd, int aDiffusantNum) {

        int ipnt = fd.thisPnt;
        double photolysis, pt;
        double dt = project.dt;

        if (!reaction || (fd.psf == null) || (fd.it >= fd.itmax)) {
            return;
        }

        if ((pulseTimer == null) || (pulseTimer.timer == null)) {
            return;
        }

        pt = pulseTimer.timer[fd.it];
        kPhotoSave = kPhoto * pt;

        if (pt <= 0) {
            return;
        }

        photolysis = fd.diffnext[aDiffusantNum][ipnt] * kPhotoSave *
                fd.psf[aDiffusantNum][ipnt] * dt;

        fd.diffnext[aDiffusantNum][ipnt] -= photolysis;

        if ((productNum >= 0) && (productNum < fd.diffnext.length)) {
            fd.diffnext[productNum][ipnt] += photolysis;
        }

    }

    @Override
    public boolean save() {
        return save.saveData(kPhotoSave);
    }

    @Override
    public void saveDimensions() {

        if ((save == null) || (!save.autoDimensions)) {
            return;
        }

        save.xdim = project.timeUnits;
        save.ydim = "kPhoto (1/" + project.timeUnits + ")";

        if ((name != null) && (name.length() > 0)) {
            save.ydim = name + " " + save.ydim;
        }

    }

    public boolean SetProductNum(int diffusantNum) {

        Diffusant d = project.getDiffusant(diffusantNum);

        if (d != null) {
            productNum = diffusantNum;
            productName = d.name;
        } else {
            productNum = -1;
            productName = null;
        }

        updateVectors();

        return true;

    }

    public boolean SetProductName(String diffusantName) {

        int i = project.getDiffusantNum(diffusantName);

        if (i >= 0) {
            productName = diffusantName;
            productNum = i;
        } else {
            productNum = -1;
            productName = null;
        }

        updateVectors();

        return true;

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        if (!(o.paramVector instanceof DiffusantPhoto)) {
            return false;
        }
        
        if (super.setMyParams(o, v)) {
            return true;
        }
        
        String n = o.getName();

        if (n.equalsIgnoreCase("productNum")) {
            return SetProductNum((int) v);
        }
        if (n.equalsIgnoreCase("kPhoto")) {
            if (v < 0) {
                return false;
            }
            kPhoto = v;
            return true;
        }
        return false;
    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, String s) {

        if ((o == null) || (s == null)) {
            return false;
        }

        if (!(o.paramVector instanceof DiffusantPhoto)) {
            return false;
        }

        if (super.setMyParams(o, s)) {
            return true;
        }

        String n = o.getName();

        if (n.equalsIgnoreCase("productName")) {
            return SetProductName(s);
        }
        return false;
    }

}
