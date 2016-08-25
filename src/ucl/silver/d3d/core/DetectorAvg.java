package ucl.silver.d3d.core;

/**
 * <p>Title: D3D</p>
 * <p>Description: 3D Diffusion-Reaction Simulator</p>
 * <p>Copyright: Copyright (c) 2016</p>
 * <p>Company: The Silver Lab at University College London</p>
 * @author Jason S. Rothman
 * @version 1.6
 */

public class DetectorAvg extends Detector {

    public DetectorAvg(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c) {
        super(p, NAME, DiffusantNum, c);
        save.dataPoints = 1;
        createVector(true);
    }

    @Override
    public boolean saveInit() {

        if (save == null) {
            return false;
        }

        int dataPoints = 1;

        return save.init(name, coordinates(), -1, dataPoints);

    }

    @Override
    public void detect(RunFiniteDifference fd, Geometry geometry) {

        double avgC = 0;
        int index, numVoxels;

        if (fd.diffus[0].length == 1) {
            save.saveData(fd.diffus[diffusantNum][0]);
            return; // single compartment
        }

        numVoxels = coordinates().index.length;

        for (int i = 0; i < numVoxels; i++) {
            index = coordinates().index[i];
            avgC += fd.diffus[diffusantNum][index];
        }

        if (numVoxels > 1) {
            avgC /= 1.0 * numVoxels;
        }

        save.saveData(avgC);

        return;

    }

}
