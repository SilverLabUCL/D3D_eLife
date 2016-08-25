package ucl.silver.d3d.core;

import java.io.Serializable;

/**
 * <p>Title: D3D</p>
 * <p>Description: 3D Diffusion-Reaction Simulator</p>
 * <p>Copyright: Copyright (c) 2016</p>
 * <p>Company: The Silver Lab at University College London</p>
 * @author Jason S. Rothman
 * @version 1.6
 *
 * Description: class wrapper of single parameters to be added to the Vector
 * in ParamVector class.
 */

public class ParamObject implements Serializable {

    private String name = ""; // parameter name
    private String type = ""; // "double" or "int" or "boolean" or "String" or "Class" or "EndClass"
    private String strValue = ""; // parameter value, saved as a string
    private String units = ""; // parameter units
    private String help = ""; // parameter help/description

    public String ERROR = null;
    
    public ParamVector paramVector = null; // pointer to the vector where parameter object resides

    public ParamObject(ParamVector PV, String NAME, double value, String UNITS, String HELP) { // create a double-type object
        if (NAME.length() > 0) {
            paramVector = PV;
            name = NAME;
            strValue = Double.toString(value);
            type = "double";
            units = UNITS;
            help = HELP;
        }
    }

    public ParamObject(ParamVector PV, String NAME, int ivalue, String UNITS, String HELP) { // create an integer-type object
        if (NAME.length() > 0) {
            paramVector = PV;
            name = NAME;
            strValue = Integer.toString(ivalue);
            type = "int";
            units = UNITS;
            help = HELP;
        }
    }

    public ParamObject(ParamVector PV, String NAME, boolean bvalue, String HELP) { // create a boolean-type object
        if (NAME.length() > 0) {
            paramVector = PV;
            name = NAME;
            strValue = Boolean.toString(bvalue);
            type = "boolean";
            help = HELP;
        }
    }

    public ParamObject(ParamVector PV, String NAME, String StrValue, String UNITS, String HELP) { // create a string-type object

        if (NAME.length() > 0) {

            strValue = StrValue;

            if (NAME.equalsIgnoreCase("Class")) {
                type = "Class";
            } else if (NAME.equalsIgnoreCase("EndClass")) {
                type = "EndClass";
            } else {
                paramVector = PV;
                type = "String";
                name = NAME;
                units = UNITS;
                help = HELP;
            }

        }

    }

    public boolean setValue(String sValue) {

        if (testValue(sValue)) {
            strValue = sValue;
            return true;
        }

        return false;

    }

    public boolean setValue(double value) {

        if (type.equalsIgnoreCase("int")) {
            strValue = Integer.toString((int) value);
        } else if (type.equalsIgnoreCase("double")) {
            strValue = Double.toString(value);
        } else if (type.equalsIgnoreCase("boolean")) {
            if (value == 0) {
                strValue = "false";
            } else {
                strValue = "true";
            }
        }

        return true;

    }

    public boolean testValue(String value) {

        if (value == null) {
            return false;
        }

        try {
            if (type.equalsIgnoreCase("int")) {
                Integer.valueOf(value.toString()); // may throw NumberFormatException
                return true;
            } else if (type.equalsIgnoreCase("double")) {
                Double.valueOf(value.toString()); // may throw NumberFormatException
                return true;
            } else if (type.equalsIgnoreCase("boolean")) {
                if (value.equalsIgnoreCase("true") || value.equalsIgnoreCase("false") ||
                        value.equals("1") || value.equals("0")) {
                    return true;
                }
            } else if (type.equalsIgnoreCase("String")) {
                return true;
            }
        } catch (NumberFormatException e) {
            return false;
        }

        return false;

    }

    public String getName() {
        return name;
    }

    public String getType() {
        return type;
    }

    public String getValue() {
        return strValue;
    }

    public String getUnits() {
        return units;
    }

    public String getHelp() {
        return help;
    }

    public boolean error() {
        if (ERROR == null) {
            return false;
        } else {
            return true;
        }
    }

}
