package decide;

import java.lang.Math;
import java.util.Collections;
import java.util.Arrays;

class App {

    static class Parameters{
        double LENGTH1;
        double RADIUS1 ;
        double EPSILON ;
        double AREA1;
        int Q_PTS ;
        int QUADS;
        double DIST ;
        int N_PTS ;
        int K_PTS ;
        int A_PTS ;
        int B_PTS ;
        int C_PTS ;
        int D_PTS ;
        int E_PTS ;
        int F_PTS ;
        int G_PTS ;
        double LENGTH2;
        double RADIUS2 ;
        double AREA2;
    }

    static final double PI = 3.1415926535;

    enum Connectors {
        // It is written as NOTUSED=777 in given header file, but having an equals sign doesn't seem to work in java.
        NOTUSED777,
        ORR,
        ANDD
    }


    enum Comptype {
        // It is written as LT=1111 in given header file, but having an equals sign doesn't seem to work in java.
        LT1111,
        EQ,
        GT
    }

    // GIVEN
    int numPoints;
    double[] COORDINATEX;
    double[] COORDINATEY;
    Connectors[][] LCM;
    boolean[] PUV;
    //To be created
    boolean[] CMV;
    boolean[][] PUM;
    boolean[] FUV;
    String LAUNCH;
    Parameters params;

    public static void main (String[] args) {

    }

    public App(int numPoints, double [] COORDINATEX, double [] COORDINATEY, boolean[] CMV, Parameters params, boolean[] PUV, Connectors[][] LCM){
        this.numPoints = numPoints;         this.COORDINATEX = COORDINATEX;
        this.LCM = LCM;                     this.PUV = PUV;
        this.CMV = CMV;                     this.COORDINATEY = COORDINATEY;
        this.params = params;
    }
    /*There exists at least one set of two consecutive data points
      that are a distance greater than the length, LENGTH1, apart. */
    boolean lic_0 () {
        double dist;
        for(int i = 0; i < numPoints -1; i++){
            dist = DISTANCE(COORDINATEX[i], COORDINATEX[i+1], COORDINATEY[i], COORDINATEY[i+1]);
            if(DOUBLECOMPARE (dist, params.LENGTH1) == Comptype.GT) return true;
        }
        return false;
    }
    /**
     * lic_1 checks if there are three consecutive points which cannot all be contained within or on a circle with
     * a radius = RADIUS1 where (RADIUS >= 0). Returns true if the condition is met.
     * @return true or false
     */
    boolean lic_1 () {
        if(params.RADIUS1 < 0 || numPoints < 3) return false;
        double dist1;       double dist2;       double dist3;
        double epsilon = 0.000001;
        for(int i = 0; i < numPoints - 2; i++){
            dist1 = Math.sqrt(Math.pow(COORDINATEX[i] - COORDINATEX[i+1] , 2) + Math.pow(COORDINATEY[i] - COORDINATEY[i+1] , 2));
            dist2 = Math.sqrt(Math.pow(COORDINATEX[i] - COORDINATEX[i+2] , 2) + Math.pow(COORDINATEY[i] - COORDINATEY[i+2] , 2));
            dist3 = Math.sqrt(Math.pow(COORDINATEX[i+1] - COORDINATEX[i+2] , 2) + Math.pow(COORDINATEY[i+1] - COORDINATEY[i+2] , 2));
            if(epsilon < (dist1 - 2 * params.RADIUS1) || epsilon < (dist2 - 2 * params.RADIUS1) || epsilon < (dist3 - 2 * params.RADIUS1)) return true;
        }
        return false;
    }

    boolean lic_2 () {
        if (params.EPSILON < 0 || params.EPSILON >= PI || numPoints < 3) return false;
        double distAB;       double distBC;     double distAC;      double angle;
        for (int i = 0; i < numPoints - 2; i++) {
            distAB = Math.sqrt(Math.pow(COORDINATEX[i] - COORDINATEX[i+1], 2) + Math.pow(COORDINATEY[i] - COORDINATEY[i+1], 2));
            distBC = Math.sqrt(Math.pow(COORDINATEX[i+1] - COORDINATEX[i+2], 2) + Math.pow(COORDINATEY[i+1] - COORDINATEY[i+2], 2));
            distAC = Math.sqrt(Math.pow(COORDINATEX[i] - COORDINATEX[i+2], 2) + Math.pow(COORDINATEY[i] - COORDINATEY[i+2], 2));

            // If point A or point C coincides with the vertex, point B, continue with next set of points.
            if (DOUBLECOMPARE(distAB, 0) == Comptype.EQ || DOUBLECOMPARE(distBC, 0) == Comptype.EQ) continue;
            
            angle = Math.acos((Math.pow(distAB, 2) + Math.pow(distBC, 2) - Math.pow(distAC, 2)) / (2 * distAB * distBC));
            angle = Math.toDegrees(angle);
            if (DOUBLECOMPARE(angle, (PI - params.EPSILON)) == Comptype.LT1111 || DOUBLECOMPARE(angle, (PI + params.EPSILON)) == Comptype.GT) return true;
        }
        return false;
    }

    /**
     * lic_3 checks if there are three consecutive points that are vertices of a triangle with
     * area greater than AREA1. Returns true if the condition is met.
     * @return true or false
     */
    boolean lic_3 () {
        if(params.RADIUS1 < 0 || numPoints < 3) return false;
        double dist1;       double dist2;       double dist3;
        double epsilon = 0.000001;
        for(int i = 0; i < numPoints - 2; i++){
            dist1 = Math.sqrt(Math.pow(COORDINATEX[i] - COORDINATEX[i+1] , 2) + Math.pow(COORDINATEY[i] - COORDINATEY[i+1] , 2));
            dist2 = Math.sqrt(Math.pow(COORDINATEX[i] - COORDINATEX[i+2] , 2) + Math.pow(COORDINATEY[i] - COORDINATEY[i+2] , 2));
            dist3 = Math.sqrt(Math.pow(COORDINATEX[i+1] - COORDINATEX[i+2] , 2) + Math.pow(COORDINATEY[i+1] - COORDINATEY[i+2] , 2));
            if(dist1 == 0 || dist2 == 0  || dist3 == 0) continue;
            double acos = (Math.pow(dist1, 2) + Math.pow(dist2, 2) - Math.pow(dist3, 2)) / (2 * dist1 * dist2);
            if(Math.abs(-1 - acos) < epsilon) acos = -1;    if(Math.abs(1 - acos) < epsilon) acos = 1;
            double degree = Math.acos(acos);
            double areaCalculated = Math.sin(degree) * dist1 * dist2 / 2;
            if( epsilon < (areaCalculated - params.AREA1)) return true;
        }
        return false;
    }

    /*Method for LIC 4.  */
    boolean lic_4 () {
        if (params.Q_PTS < 2 || params.Q_PTS > numPoints) return false;
        if (params.QUADS < 1 || params.QUADS > 3) return false;

        boolean[] fulfilledQuad = new boolean[4];
        int nrQuads = 0;
        for (int i = 0; i < numPoints; i++) {
            for (int j = 0; j < params.Q_PTS ; j++) {
                if (COORDINATEX[(i+j) % numPoints] >= 0 && COORDINATEY[(i+j) % numPoints] >= 0) {
                    fulfilledQuad[0] = true;
                } else if (COORDINATEX[(i+j) % numPoints] < 0 && COORDINATEY[(i+j) % numPoints] >= 0) {
                    fulfilledQuad[1] = true;
                } else if (COORDINATEX[(i+j) % numPoints] <= 0 && COORDINATEY[(i+j) % numPoints] < 0) {
                    fulfilledQuad[2] = true;
                } else {
                    fulfilledQuad[3] = true;
                }
            }
            for (boolean b : fulfilledQuad) {
                if (b)
                    nrQuads++;
            }
            if (nrQuads > params.QUADS)
                return true;
            nrQuads = 0;
            fulfilledQuad = new boolean[4];
        }
        return false;
    }
    /**Method for LIC 5. This method checks if there are consecutive data points, (X[i],Y[i]) and (X[j],Y[j]),
     such that X[j] - X[i] < 0. (where i = j-1) */
    boolean lic_5 () {
        double diffX;
        for (int i = 0; i < numPoints - 1; i++) {
            diffX = COORDINATEX[i+1] - COORDINATEX[i];
            if (diffX < 0){
                return true;
            }
        }
        return false;
    }
    /** Method for solving LIC 6 using Heron's formula. */
    boolean lic_6 () {
        double a;
        double b;
        double c;
        double s;
        double Area;
        double h;
        double distance;
        if (params.N_PTS < 3) {
            return false;
        }
        /** Solving using Heron's formula */
        for (int i = 0; i < numPoints - params.N_PTS + 1; i++) {
            if (COORDINATEX[i] == COORDINATEX[i + (params.N_PTS- 1)]) {
                for (int j = i; j < i + params.N_PTS; i++) {
                    /** This part is the distance between the point and the line*/
                    distance = Math.sqrt((COORDINATEY[i] - COORDINATEY[j]) * (COORDINATEY[i] - COORDINATEY[j]) +
                            (COORDINATEX[i] - COORDINATEX[j]) * (COORDINATEX[i] - COORDINATEX[j]));
                    if (Double.compare(distance, params.DIST) >= 0  ) {
                        return true;
                    }
                }
            }
            b = euDist(COORDINATEX[i], COORDINATEX[i + params.N_PTS - 1], COORDINATEY[i], COORDINATEY[i + params.N_PTS - 1]);
            for (int j = i; j < i + params.N_PTS - 1; j++) {
                a = euDist(COORDINATEX[i], COORDINATEX[j], COORDINATEY[i], COORDINATEY[j]);
                c = euDist(COORDINATEX[j], COORDINATEX[i + params.N_PTS - 1], COORDINATEY[j], COORDINATEY[i + params.N_PTS - 1]);
                s = (a + b + c) / 2;
                Area = triArea(a, b, c, s);
                h = (2*Area) / b;
                if (Double.compare(h, params.DIST) > 0){
                    return true;
                }
            }
        }
        return false;
    }

    boolean lic_7 () {
        if(numPoints < 3 || params.K_PTS > numPoints - 2) {return false;}
        for(int i = 0; i + params.K_PTS +1 < numPoints; i++){
            double x1 = COORDINATEX[i];
            double y1 = COORDINATEY[i];
            double x2 = COORDINATEX[i + params.K_PTS + 1];
            double y2 = COORDINATEX[i + params.K_PTS + 1];
            double dist = DISTANCE(x1, x2, y1, y2);
            if (DOUBLECOMPARE (dist, params.LENGTH1) == Comptype.GT){return true;}
        }
        return false;
    }

    boolean lic_8 () {
        if(numPoints < 5 || (params.A_PTS + params.B_PTS) > numPoints - 3) {return false;}
        for(int i = 0; i + params.A_PTS + params.B_PTS +2 < numPoints; i++){
            double x1 = COORDINATEX[i];
            double y1 = COORDINATEY[i];
            double x2 = COORDINATEX[i + params.A_PTS + 1];
            double y2 = COORDINATEY[i + params.A_PTS + 1];
            double x3 = COORDINATEX[i + params.A_PTS + params.B_PTS + 2];
            double y3 = COORDINATEY[i + params.A_PTS + params.B_PTS + 2];

            double r = findRadius(x1, y1, x2, y2, x3, y3);
            if (DOUBLECOMPARE (r, params.RADIUS1) == Comptype.GT){return true;}
        }
        return false;
    }
    /** Method for implementing LIC 9 using the law of cosine */
    boolean lic_9 () {
        if (numPoints < 5){
            return false;
        }
        double a;
        double b;
        double c;
        double cos;
        double angle;
        for (int i = 0; i < numPoints - params.C_PTS - params.D_PTS - 2; i++){
            if (COORDINATEX[i + params.C_PTS + 1] == COORDINATEX[i] && COORDINATEY[i + params.C_PTS + 1] == COORDINATEY[i]){
                continue;
            }
            if (COORDINATEX[i + params.C_PTS + 1] == COORDINATEX[i + params.C_PTS + 1 + params.D_PTS + 1] &&
                COORDINATEY[i + params.C_PTS + 1] == COORDINATEY[i + params.C_PTS + 1 + params.D_PTS + 1]){
                continue;
            }
            a = euDist(COORDINATEX[i + params.C_PTS + 1 + params.D_PTS + 1], COORDINATEX[i + params.C_PTS + 1],
                COORDINATEY[i + params.C_PTS + 1 + params.D_PTS + 1], COORDINATEY[i + params.C_PTS + 1]);
            b = euDist(COORDINATEX[i], COORDINATEX[i + params.C_PTS + 1], COORDINATEY[i], COORDINATEY[i + params.C_PTS + 1]);
            c = euDist(COORDINATEX[i], COORDINATEX[i + params.C_PTS + 1 + params.D_PTS + 1], COORDINATEY[i], COORDINATEY[i + params.C_PTS + 1 + params.D_PTS + 1]);
            //COS of the angle at second point:
            cos = (Math.pow(a, 2) + Math.pow(b, 2) - Math.pow(c, 2)) / (2 * a * b);
            angle = Math.acos(cos);
            if (Double.compare(angle, (PI - params.EPSILON)) < 0 || Double.compare(angle, (PI + params.EPSILON)) > 0 ){
                return true;
            }
        }
        return false;
    }

    boolean lic_10 () {
        if (params.E_PTS < 1 || params.F_PTS < 1 || numPoints < 5) return false;
        if (!((params.E_PTS + params.F_PTS) <= (numPoints - 3))) return false;

        double dist1;       double dist2;       double dist3;
        int second;       int third;
        double semiperimeter;       double area;
        for (int first = 0; first < numPoints; first++) {
            second = (first+params.E_PTS+1) % numPoints;
            third = (first+params.E_PTS+params.F_PTS+2) % numPoints;
            
            dist1 = Math.sqrt(Math.pow(COORDINATEX[first] - COORDINATEX[second] , 2) + Math.pow(COORDINATEY[first] - COORDINATEY[second] , 2));
            dist2 = Math.sqrt(Math.pow(COORDINATEX[second] - COORDINATEX[third] , 2) + Math.pow(COORDINATEY[second] - COORDINATEY[third] , 2));
            dist3 = Math.sqrt(Math.pow(COORDINATEX[third] - COORDINATEX[first] , 2) + Math.pow(COORDINATEY[third] - COORDINATEY[first] , 2));
            // System.out.println("Dist1: " + dist1 + ", dist2: " + dist2 + ", dist3: " + dist3);
            semiperimeter = (dist1 + dist2 + dist3)/2;
            area = Math.sqrt(semiperimeter*(semiperimeter-dist1)*(semiperimeter-dist2)*(semiperimeter-dist3));
            if (DOUBLECOMPARE(area, params.AREA1) == Comptype.GT) return true;
        }
        return false;
    }

    boolean lic_11 () {
        if(numPoints < 3 || params.G_PTS > numPoints - 2) {return false;}
        for(int i = 0; i + params.G_PTS +1 < numPoints; i++){
            double x1 = COORDINATEX[i];
            double x2 = COORDINATEX[i + params.G_PTS + 1];
            if(DOUBLECOMPARE(x2, x1) == Comptype.LT1111){return true;}
        }
        return false;
    }

    boolean lic_12 () {
        if (params.LENGTH2 < 0) return false;
        boolean condition1 = false;     boolean condition2 = false;     double dist;
        for (int i = 0; i < numPoints; i++) {
            dist = Math.sqrt(Math.pow(COORDINATEX[i] - COORDINATEX[(i+params.K_PTS+1) % numPoints] , 2) + Math.pow(COORDINATEY[i] - COORDINATEY[(i+params.K_PTS+1) % numPoints] , 2));
            if (dist > params.LENGTH1){
                condition1 = true;
            }
            if (dist < params.LENGTH2) {
                condition2 = true;
            }
        }
        if (condition1 && condition2) return true;

        return false;
    }

    boolean lic_13 () {
        double radius = 0;
        boolean flag1 = false;
        boolean flag2 = false;
        double[] line = new double[3];

        if (numPoints < 5){
            return false;
        }

        for (int i = 0; i < numPoints - params.A_PTS - params.B_PTS - 2; i++){
            if (COORDINATEX[i] == COORDINATEX[i + params.A_PTS + 1] && COORDINATEX[i] == COORDINATEX[i + params.A_PTS + 1 + params.B_PTS +1]
             && COORDINATEX[i + params.A_PTS + 1] == COORDINATEX[i + params.A_PTS + 1 + params.B_PTS + 1]
                    || COORDINATEY[i] == COORDINATEY[i + params.A_PTS + 1] && COORDINATEY[i] == COORDINATEY[i + params.A_PTS + 1 + params.B_PTS +1]
                    && COORDINATEY[i + params.A_PTS + 1] == COORDINATEY[i + params.A_PTS + 1 + params.B_PTS + 1] ){
                line[0] = euDist(COORDINATEX[i], COORDINATEX[i + params.A_PTS + 1], COORDINATEY[i], COORDINATEY[i + params.A_PTS +1]);
                line[1] = euDist(COORDINATEX[i], COORDINATEX[i + params.A_PTS + 1 + params.B_PTS + 1], COORDINATEY[i], COORDINATEY[i + params.A_PTS +1 + params.B_PTS + 1]);
                line[2] = euDist(COORDINATEX[i + params.A_PTS +1], COORDINATEX[i + params.A_PTS + 1 + params.B_PTS + 1], COORDINATEY[i + params.A_PTS + 1], COORDINATEY[i + params.A_PTS +1 + params.B_PTS + 1]);
                double mini = -1;
                for (int j = 0; j > line.length; j++){
                    if (line[j] < mini){
                        mini = line[j];
                        radius = mini;
                    }
                }
            }else {
                radius = findRadius(COORDINATEX[i], COORDINATEY[i], COORDINATEX[i + params.A_PTS + 1], COORDINATEY[i + params.A_PTS + 1],
                        COORDINATEX[i + params.A_PTS + 1 + params.B_PTS + 1], COORDINATEY[params.A_PTS + 1 + params.B_PTS + 1]);
            }
            if (radius >= params.RADIUS1){
                flag1 = true;
            }
            if (radius <= params.RADIUS2){
                flag2 = true;
            }
            if (flag1 && flag2){
                return true;
            }
        }
        return false;
    }
    /**
     * lic_14 checks if there are three points, seperated by exactly E_PTS and F_PTS, respectively, that are vertices of a triangle
     * with area greater than AREA1 and less than AREA2. (which can be the same or different from the three data points just mentioned).
     * Returns true if the condition is met.
     * @return true or false
     */
    boolean lic_14 () {
        if(params.E_PTS < 1 || params.F_PTS < 1 || params.AREA1 < 0 || params.AREA1 < 0 || numPoints < 5) return false;
        double dist1;       double dist2;       double dist3;
        double epsilon = 0.000001;  int totalRange = params.E_PTS + params.F_PTS + 2;
        boolean flag1 = false, flag2 = false;
        for(int i = 0; i < numPoints - totalRange; i++){
            dist1 = Math.sqrt(Math.pow(COORDINATEX[i] - COORDINATEX[i+params.E_PTS + 1] , 2) +
            Math.pow(COORDINATEY[i] - COORDINATEY[i+params.E_PTS + 1] , 2));

            dist2 = Math.sqrt(Math.pow(COORDINATEX[i] - COORDINATEX[i+totalRange] , 2) +
            Math.pow(COORDINATEY[i] - COORDINATEY[i+totalRange] , 2));

            dist3 = Math.sqrt(Math.pow(COORDINATEX[i + params.E_PTS + 1] - COORDINATEX[i+totalRange] , 2) +
            Math.pow(COORDINATEY[i + params.E_PTS + 1] - COORDINATEY[i+totalRange] , 2));

            if(dist1 == 0 || dist2 == 0  || dist3 == 0) continue;
            double acos = (Math.pow(dist1, 2) + Math.pow(dist2, 2) - Math.pow(dist3, 2)) / (2 * dist1 * dist2);
            if(Math.abs(-1 - acos) < epsilon) acos = -1;    if(Math.abs(1 - acos) < epsilon) acos = 1;
            double degree = Math.acos(acos);
            double areaCalculated = Math.sin(degree) * dist1 * dist2 / 2;
            if(DOUBLECOMPARE(areaCalculated, params.AREA1) == Comptype.GT) flag1 = true;
            if(DOUBLECOMPARE(areaCalculated, params.AREA2) == Comptype.LT1111) flag2 = true;
            if(flag1 && flag2) return true;
        }
        return false;

    }
    // Function to calculate the distance between 2 points
    public static double euDist (double x1, double x2, double y1, double y2){
        return Math.sqrt((y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1));
    }

    public static double triArea (double a, double b, double c, double s){
        return Math.sqrt(s * (s-a) * (s-b) * (s-c));
    }
    // Function to calculate the distance between 2 points
    public static double DISTANCE (double x1, double x2, double y1, double y2){
        return Math.sqrt(Math.pow((x2 - x1), 2) + Math.pow((y2 - y1), 2));
    }
    static Comptype DOUBLECOMPARE (double A, double B) {
        if (Math.abs(A-B)<0.000001)
            return Comptype.EQ;
        if (A<B)
            return Comptype.LT1111;
        return Comptype.GT;
    }
    public static double findRadius (double x1, double y1,double x2, double y2,double x3, double y3){
        double x12 = x1 - x2;
        double x13 = x1 - x3;
        double y12 = y1 - y2;
        double y13 = y1 - y3;
        double y31 = y3 - y1;
        double y21 = y2 - y1;
        double x31 = x3 - x1;
        double x21 = x2 - x1;
        double sx13 = (Math.pow(x1, 2) - Math.pow(x3, 2));
        double sy13 = (Math.pow(y1, 2) - Math.pow(y3, 2));
        double sx21 = (Math.pow(x2, 2) - Math.pow(x1, 2));
        double sy21 = (Math.pow(y2, 2) - Math.pow(y1, 2));
        double f = ((sx13) * (x12)+ (sy13) * (x12)+ (sx21) * (x13)+ (sy21) * (x13))/ (2 * ((y31) * (x12) - (y21) * (x13)));
        double g = ((sx13) * (y12)+ (sy13) * (y12)+ (sx21) * (y13)+ (sy21) * (y13))/ (2 * ((x31) * (y12) - (x21) * (y13)));
        double c = -Math.pow(x1, 2) - Math.pow(y1, 2) - 2 * g * x1 - 2 * f * y1;
        double h = -g;
        double k = -f;
        double sqr_of_r = h * h + k * k - c;
        double r = Math.sqrt(sqr_of_r);
        return (r);
    }

    // Should call on all the LIC-functions.
    // Should check the condition for numPoints (2 ≤ NUMPOINTS ≤ 100)
    // and other similar conditions should be checked here if they do not meet return false
    void DECIDE () {
        if(numPoints < 2 || numPoints > 100){
            LAUNCH = "";
            return;         // Because we don't want to continue in that case?
        }

        PUM = new boolean[15][15];
        FUV = new boolean[15];

        CMV[0] = lic_0();       CMV[1] = lic_1();       CMV[2] = lic_2();       CMV[3] = lic_3();   CMV[4] = lic_4();
        CMV[5] = lic_5();       CMV[6] = lic_6();       CMV[7] = lic_7();       CMV[8] = lic_8();   CMV[9] = lic_9();
        CMV[10] = lic_10();     CMV[11] = lic_11();     CMV[12] = lic_12();     CMV[13] = lic_13(); CMV[14] = lic_14();

        for(int i = 0; i < 15; i++){
            for(int j = 0; j < 15; j++){
                if(LCM[i][j] == Connectors.ANDD) PUM[i][j] = CMV[i] && CMV[j];
                else if(LCM[i][j] == Connectors.ORR) PUM[i][j] = CMV[i] || CMV[j];
                else if(LCM[i][j] == Connectors.NOTUSED777) PUM[i][j] = true;
            }
        }

        for(int k = 0; k < 15; k++){
            // Can be written as "if (!PUV[k])"
            if(PUV[k] == false) FUV[k] = true;
            else{
                FUV[k] = true;
                for(int l = 0; l < 15; l++){
                    FUV[k] = FUV[k] && PUM[k][l]; 
                }
            }
        }

        //Can be made global
        boolean launch = true;
        for(int f = 0; f < 15; f++){
            launch = launch && FUV[f];
        }
        if(launch) LAUNCH = "YES";
        else LAUNCH = "NO";

        System.out.println(LAUNCH);
    
    }
}

