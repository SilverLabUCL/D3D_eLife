package ucl.silver.d3d.utils;

/**
 * <p>Title: D3D</p>
 * <p>Description: 3D Diffusion-Reaction Simulator</p>
 * <p>Copyright: Copyright (c) 2016</p>
 * <p>Company: The Silver Lab at University College London</p>
 * @author Jason S. Rothman
 * @version 1.6
 */

public final class D3DUtility {

    private D3DUtility() {
        // cannot be instantiated
    }

    public static int makeOdd(int i) {
        if (Math.IEEEremainder(i, 2) == 0) {
            return i + 1;
        }
        return i;
    }

    public static int makeEven(int i) {
        if (Math.IEEEremainder(i, 2) == 0) {
            return i;
        }
        return i + 1;
    }

    public static String parentClass(String longNameWithPeriod) { // e.g. "Detector0.PSF.numericalAperture"

        if (longNameWithPeriod == null) {
            return null;
        }

        String[] splitStr = longNameWithPeriod.split("\\.");

        if (splitStr.length >= 2) {
            return splitStr[0]; // e.g. "Detector0"
        }

        return null;

    }

    public static String childClass(String longNameWithPeriod) { // e.g. "Detector0.PSF.numericalAperture"

        if (longNameWithPeriod == null) {
            return null;
        }

        String[] splitStr = longNameWithPeriod.split("\\.");

        if (splitStr.length == 3) {
            return splitStr[1]; // e.g. "PSF"
        }

        return null;

    }

    public static String parameterName(String longName) { // e.g. "Detector0.PSF.numericalAperture"

        if (longName == null) {
            return null;
        }

        String[] splitStr = longName.split("\\.");

        if (splitStr.length > 1) {
            return splitStr[splitStr.length - 1]; // e.g. "numericalAperture"
        }

        return null;

    }

    public static int itemsInList(String strList, char stop) {

        int count = 0;

        if ((strList == null) || (strList.length() == 0)) {
            return 0;
        }

        count = 1;

        for (int i = 0; i < strList.length() - 1; i += 1) {
            if (strList.charAt(i) == stop) {
                count += 1;
            }
        }

        return count;

    }

    public static String itemFromList(String strList, int item) {

        int count = 0, start = 0, end = -1;

        if ((strList == null) || (strList.length() == 0)) {
            return "";
        }

        if (strList.indexOf(",") < 0) {
            return strList;
        }

        for (int i = 0; i < strList.length(); i += 1) {
            if (strList.charAt(i) == ',') {
                count += 1;
                start = end + 1;
                end = i;
                if (item == count - 1) {
                    return strList.substring(start, end);
                }
            }
        }

        start = end + 1;
        end = strList.length();

        return strList.substring(start, end);

    }

    public static int whichListItem(String[] strList, String itemStr, boolean ignoreCase) {

        if ((strList == null) || (strList.length == 0)) {
            return -1;
        }

        if (itemStr == null) {
            return -1;
        }

        for (int i = 0; i < strList.length; i++) {

            if (strList[i] == null) {
                continue;
            }

            if (ignoreCase) {
                if (strList[i].equalsIgnoreCase(itemStr)) {
                    return i;
                }
            } else {
                if (strList[i].equals(itemStr)) {
                    return i;
                }
            }

        }

        return -1;

    }

    // limit integer value
    public static int limit(int value, int imin, int imax) {
        return Math.min(imax, Math.max(value, imin));
    }

    // find min value of double array
    public static double getMin(double d[]) {

        double min = d[0];

        for (int i = 1; i < d.length; i++) {
            if (d[i] < min) {
                min = d[i];
            }
        }

        return min;

    }

    // find min value of double array
    public static double getMin(double d[][]) {

        double min = d[0][0];

        for (int i = 0; i < d.length; i++) {
            for (int j = 0; j < d[0].length; j++) {
                if (d[i][j] < min) {
                    min = d[i][j];
                }
            }
        }

        return min;

    }

    // find max value of double array
    public static double getMax(double d[]) {

        double max = d[0];

        for (int i = 1; i < d.length; i++) {
            if (d[i] > max) {
                max = d[i];
            }
        }

        return max;

    }

    // find max value of double array
    public static double getMax(double d[][]) {

        double max = d[0][0];

        for (int i = 0; i < d.length; i++) {
            for (int j = 0; j < d[0].length; j++) {
                if (d[i][j] > max) {
                    max = d[i][j];
                }
            }
        }

        return max;

    }

    public static int[] resizeArray(int[] array, int iNew) {

        int iNum = 0;

        if (array != null) {
            iNum = array.length;
        }

        if (iNew < 0) {
            iNew = iNum;
        }

        if (iNew <= 0) {
            return null;
        }

        int[] temp = new int[iNew];

        for (int i = 0; i < iNew; i++) {
            if (i >= iNum) {
                break;
            }
            temp[i] = array[i];
        }

        return temp;

    }

    public static byte[] resizeArray(byte[] array, int iNew) {

        int iNum = 0;

        if (array != null) {
            iNum = array.length;
        }

        if (iNew < 0) {
            iNew = iNum;
        }

        if (iNew <= 0) {
            return null;
        }

        byte[] temp = new byte[iNew];

        for (int i = 0; i < iNew; i++) {
            if (i >= iNum) {
                break;
            }
            temp[i] = array[i];
        }

        return temp;

    }

    public static double[] resizeArray(double[] array, int iNew) {

        int iNum = 0;

        if (array != null) {
            iNum = array.length;
        }

        if (iNew < 0) {
            iNew = iNum;
        }

        if (iNew <= 0) {
            return null;
        }

        double[] temp = new double[iNew];

        for (int i = 0; i < iNew; i++) {
            if (i >= iNum) {
                break;
            }
            temp[i] = array[i];
        }

        return temp;

    }

    public static double[][][] resizeArray(double[][][] array, int iNew, int jNew, int kNew) {

        int i, j, k;

        int iNum = 0, jNum = 0, kNum = 0;

        if (array != null) {
            iNum = array.length;
            jNum = array[0].length;
            kNum = array[0][0].length;
        }

        if (iNew < 0) {
            iNew = iNum;
        }

        if (jNew < 0) {
            jNew = jNum;
        }

        if (kNew < 0) {
            kNew = kNum;
        }

        if ((iNew <= 0) || (jNew <= 0) || (kNew <= 0)) {
            return null;
        }

        double[][][] temp = new double[iNew][jNew][kNew];

        for (k = 0; k < kNew; k++) {
            if (k >= kNum) {
                break;
            }
            for (j = 0; j < jNew; j++) {
                if (j >= jNum) {
                    break;
                }
                for (i = 0; i < iNew; i++) {
                    if (i >= iNum) {
                        break;
                    }
                    temp[i][j][k] = array[i][j][k];
                }
            }
        }

        return temp;

    }

    // resize byte array matrix
    public static int[][][] resizeArray(int[][][] array, int iNew, int jNew, int kNew) {

        int i, j, k;

        int iNum = 0, jNum = 0, kNum = 0;

        if (array != null) {
            iNum = array.length;
            jNum = array[0].length;
            kNum = array[0][0].length;
        }

        if (iNew < 0) {
            iNew = iNum;
        }

        if (jNew < 0) {
            jNew = jNum;
        }

        if (kNew < 0) {
            kNew = kNum;
        }

        if ((iNew <= 0) || (jNew <= 0) || (kNew <= 0)) {
            return null;
        }

        int[][][] temp = new int[iNew][jNew][kNew];

        for (k = 0; k < kNew; k++) {
            if (k >= kNum) {
                break;
            }
            for (j = 0; j < jNew; j++) {
                if (j >= jNum) {
                    break;
                }
                for (i = 0; i < iNew; i++) {
                    if (i >= iNum) {
                        break;
                    }
                    temp[i][j][k] = array[i][j][k];
                }
            }
        }

        return temp;

    }
}
