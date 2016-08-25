package ucl.silver.d3d.core;

import java.util.Random;

/**
 * <p>Title: D3D</p>
 * <p>Description: 3D Diffusion-Reaction Simulator</p>
 * <p>Copyright: Copyright (c) 2016</p>
 * <p>Company: The Silver Lab at University College London</p>
 * @author Jason S. Rothman
 * @version 1.6
 */

public final class GeometryTools {

    public static CoordinatesVoxels coordinates[];
    public static boolean useCoordinates = false;

    // cannot instantiate
    private GeometryTools() {
    }

    public static int spaceVoxels(int[][][] a) {

        int spaceVoxels = 0;

        for (int k = 0; k < a[0][0].length; k++) {
            for (int j = 0; j < a[0].length; j++) {
                for (int i = 0; i < a.length; i++) {

                    if (a[i][j][k] < 0) {
                    } else {
                        spaceVoxels += 1;
                    }

                }
            }
        }

        return spaceVoxels;

    }

    // Scale values are for circular variation (0) no (1) yes
    public static void ellipsoid(double[][][] a, CoordinatesVoxels c, double value) { // create elliptical space (center a0, b0, c0)

        double xd, yd, zd;

        if ((a == null) || (c == null)) {
            return;
        }

        for (int k = c.zVoxel1; k <= c.zVoxel2; k++) {
            for (int j = c.yVoxel1; j <= c.yVoxel2; j++) {
                for (int i = c.xVoxel1; i <= c.xVoxel2; i++) {

                    if ((i < 0) || (i >= a.length)) {
                        continue;
                    }

                    if ((j < 0) || (j >= a[0].length)) {
                        continue;
                    }

                    if ((k < 0) || (k >= a[0][0].length)) {
                        continue;
                    }

                    xd = (i - c.xVoxelCenter) / (c.xVoxels / 2.0);
                    yd = (j - c.yVoxelCenter) / (c.yVoxels / 2.0);
                    zd = (k - c.zVoxelCenter) / (c.zVoxels / 2.0);

                    if ((xd * xd * c.xScale + yd * yd * c.yScale + zd * zd * c.zScale) <= 1) {
                        a[i][j][k] = value;
                    }

                }
            }
        }

    }

    public static void ellipsoid(int[][][] a, CoordinatesVoxels c, int value) { // create elliptical space (center a0, b0, c0)

        double xd, yd, zd;

        if ((a == null) || (c == null)) {
            return;
        }

        for (int k = c.zVoxel1; k <= c.zVoxel2; k++) {
            for (int j = c.yVoxel1; j <= c.yVoxel2; j++) {
                for (int i = c.xVoxel1; i <= c.xVoxel2; i++) {

                    if ((i < 0) || (i >= a.length)) {
                        continue;
                    }

                    if ((j < 0) || (j >= a[0].length)) {
                        continue;
                    }

                    if ((k < 0) || (k >= a[0][0].length)) {
                        continue;
                    }

                    xd = (i - c.xVoxelCenter) / (c.xVoxels / 2.0);
                    yd = (j - c.yVoxelCenter) / (c.yVoxels / 2.0);
                    zd = (k - c.zVoxelCenter) / (c.zVoxels / 2.0);

                    if ((xd * xd * c.xScale + yd * yd * c.yScale + zd * zd * c.zScale) <= 1) {
                        a[i][j][k] = value;
                    }

                }
            }
        }

    }

    public static double addEllipsoids(Geometry geometry, CoordinatesVoxels c, double radius, double axial_ratio, double volumeFraction, long seed, boolean exact, int ijkSelect, boolean cylinder, int linkNum) {

        int i, j, k, iradius, ihalfLength, idiameter, numCylinders = 0;
        int oldVoxels, newVoxels, cylinderVoxels, spaceVoxels, iLast;
        int maxCylinders = 500, maxTrials = 1000, trialCount = 0;
        double d = 0, rani, ranj, rank;
        boolean overlap = true;
        boolean randomIJK = false;
        boolean even = true;
        boolean link = false;
        int linkCounter = 0;

        double dx = geometry.project.dx;
        
        int[][][] spaceCopy = geometry.getSpaceCopy();

        int extra = 25;
        int xv1 = c.xVoxel1 - extra;
        int yv1 = c.yVoxel1 - extra;
        int zv1 = c.zVoxel1 - extra;
        int xv2 = c.xVoxel2 + extra;
        int yv2 = c.xVoxel2 + extra;
        int zv2 = c.xVoxel2 + extra;

        int cx1, cy1, cz1, cx2, cy2, cz2;

        CoordinatesVoxels ctemp = new CoordinatesVoxels(geometry.project);

        if (seed < 0) {
            seed = 8682522807148012L + System.nanoTime();
            //Master.log("addCylinders: random number generator seed: " + seed + "L");
        }

        Random ran = new Random(seed);

        if (useCoordinates) {
            maxCylinders = coordinates.length;
        } else {
            coordinates = new CoordinatesVoxels[maxCylinders];
        }

        iradius = (int) Math.round(radius / dx);
        idiameter = (int) Math.round(2 * radius / dx);
        ihalfLength = (int) Math.round(0.5 * axial_ratio * idiameter);

        if (Math.IEEEremainder(idiameter, 2) == 0) {
            even = true;
        } else {
            even = false;
            iradius -= 1;
        }

        spaceVoxels = geometry.spaceVoxels;
        oldVoxels = geometry.spaceVoxels;
        newVoxels = (int) (oldVoxels * (1.0 - volumeFraction));

        if (ijkSelect < 0) {
            randomIJK = true;
        }

        if (linkNum > 1) {
            link = true;
        }

        for (int ii = 0; ii < maxCylinders; ii++) {

            overlap = false;

            if (useCoordinates) {

                if (coordinates[ii] == null) {
                    continue;
                }
                
            } else {

                coordinates[ii] = new CoordinatesVoxels(geometry.project);

                if (randomIJK) {

                    d = ran.nextDouble();

                    if (d < 0.3333333333) {
                        ijkSelect = 0; // cylinderx
                    } else if (d < 0.6666666666) {
                        ijkSelect = 1; // cylindery
                    } else {
                        ijkSelect = 2; // cylinderz
                    }

                }

                trialCount = 0;

                for (int jj = 0; jj < maxTrials; jj++) {

                    if (link && (linkCounter > 0)) {

                        switch (ijkSelect) {
                            case 0: // cylinderx
                                cx1 = coordinates[ii - linkCounter].xVoxel1;
                                cy1 = coordinates[ii - linkCounter].yVoxel1 - 1 * iradius;
                                cz1 = coordinates[ii - linkCounter].zVoxel1 - 1 * iradius;
                                cx2 = coordinates[ii - linkCounter].xVoxel2;
                                cy2 = coordinates[ii - linkCounter].yVoxel2 + 1 * iradius;
                                cz2 = coordinates[ii - linkCounter].zVoxel2 + 1 * iradius;
                                break;
                            case 1: // cylindery
                                cx1 = coordinates[ii - linkCounter].xVoxel1 - 1 * iradius;
                                cy1 = coordinates[ii - linkCounter].yVoxel1;
                                cz1 = coordinates[ii - linkCounter].zVoxel1 - 1 * iradius;
                                cx2 = coordinates[ii - linkCounter].xVoxel2 + 1 * iradius;
                                cy2 = coordinates[ii - linkCounter].yVoxel2;
                                cz2 = coordinates[ii - linkCounter].zVoxel2 + 1 * iradius;
                                break;
                            case 2: // cylinderz
                                cx1 = coordinates[ii - linkCounter].xVoxel1 - 1 * iradius;
                                cy1 = coordinates[ii - linkCounter].yVoxel1 - 1 * iradius;
                                cz1 = coordinates[ii - linkCounter].zVoxel1;
                                cx2 = coordinates[ii - linkCounter].xVoxel2 + 1 * iradius;
                                cy2 = coordinates[ii - linkCounter].yVoxel2 + 1 * iradius;
                                cz2 = coordinates[ii - linkCounter].zVoxel2;
                                break;
                            default:
                                return 0;
                        }

                        rani = cx1 + ran.nextDouble() * (cx2 - cx1);
                        ranj = cy1 + ran.nextDouble() * (cy2 - cy1);
                        rank = cz1 + ran.nextDouble() * (cz2 - cz1);

                    } else {

                        rani = xv1 + ran.nextDouble() * (xv2 - xv1);
                        ranj = yv1 + ran.nextDouble() * (yv2 - yv1);
                        rank = zv1 + ran.nextDouble() * (zv2 - zv1);

                    }

                    i = (int) Math.round(rani);
                    j = (int) Math.round(ranj);
                    k = (int) Math.round(rank);

                    coordinates[ii].xScale = 1;
                    coordinates[ii].yScale = 1;
                    coordinates[ii].zScale = 1;

                    switch (ijkSelect) {

                        case 0: // cylinderx

                            if (even) {
                                coordinates[ii].setVoxels(i - ihalfLength, j - iradius + 1, k - iradius + 1, i + ihalfLength, j + iradius, k + iradius);
                            } else {
                                coordinates[ii].setVoxels(i - ihalfLength, j - iradius, k - iradius, i + ihalfLength, j + iradius, k + iradius);
                            }

                            if (cylinder) {
                                coordinates[ii].xScale = 0;
                            }

                            break;

                        case 1: // cylindery

                            if (even) {
                                coordinates[ii].setVoxels(i - iradius + 1, j - ihalfLength, k - iradius + 1, i + iradius, j + ihalfLength, k + iradius);
                            } else {
                                coordinates[ii].setVoxels(i - iradius, j - ihalfLength, k - iradius, i + iradius, j + ihalfLength, k + iradius);
                            }

                            if (cylinder) {
                                coordinates[ii].yScale = 0;
                            }

                            break;

                        case 2: // cylinderz

                            if (even) {
                                coordinates[ii].setVoxels(i - iradius + 1, j - iradius + 1, k - ihalfLength, i + iradius, j + iradius, k + ihalfLength);
                            } else {
                                coordinates[ii].setVoxels(i - iradius, j - iradius, k - ihalfLength, i + iradius, j + iradius, k + ihalfLength);
                            }

                            if (cylinder) {
                                coordinates[ii].zScale = 0;
                            }

                            break;

                    }

                    overlap = false;

                    for (int kk = 0; kk < ii; kk++) {
                        if (coordinates[ii].isInside(coordinates[kk])) {
                            overlap = true;
                        }
                    }


                    trialCount++;

                    if (!overlap) {
                        break;
                    }

                }

                if (overlap) {
                    continue;
                }

            }

            if (overlap) {
                coordinates[ii] = null;
            } else {

                GeometryTools.ellipsoid(spaceCopy, coordinates[ii], -1);
                numCylinders++;

                if (link) {

                    linkCounter++;

                    if (linkCounter >= linkNum) {
                        linkCounter = 0;
                    }

                }
                
            }
            
            spaceVoxels = spaceVoxels(spaceCopy);

            if (spaceVoxels <= newVoxels) {
                break; // finished
            }

        }
        
        if ((exact && (spaceVoxels < newVoxels))) {

            iLast = -1;

            for (int ii = 0; ii < numCylinders; ii++) {

                if (coordinates[ii] == null) {
                    continue;
                }

                iLast = ii;

            }

            for (int jj = 0; jj < geometry.xVoxels; jj++) {

                if (iLast < 0) {
                    break;
                }

                spaceCopy = geometry.getSpaceCopy();

                ctemp.matchVoxels(coordinates[iLast]);

                if (ctemp.zScale == 0) {
                    ctemp.zVoxel1 = geometry.zVoxels - jj - 1;
                    ctemp.update();
                } else if (ctemp.xScale == 0) {
                    ctemp.xVoxel1 = geometry.xVoxels - jj - 1;
                    ctemp.update();
                } else if (ctemp.yScale == 0) {
                    ctemp.yVoxel1 = geometry.yVoxels - jj - 1;
                    ctemp.update();
                }

                GeometryTools.ellipsoid(spaceCopy, ctemp, -1);

                spaceVoxels = spaceVoxels(spaceCopy);

                if (spaceVoxels <= newVoxels) {
                    break; // finished
                }

            }

            geometry.setSpace(spaceCopy);
            
        } else {

            geometry.setSpace(spaceCopy);

        }

        cylinderVoxels = oldVoxels - geometry.spaceVoxels;

        d = ( 1.0 * cylinderVoxels ) / ( 1.0  * oldVoxels); // volume fraction

        Master.log("added ellipsoids = " + numCylinders + ", space voxels = " + spaceVoxels + ", density = " + d);

        return d;

    }

}
