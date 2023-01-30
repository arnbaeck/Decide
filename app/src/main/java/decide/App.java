package decide;

import java.lang.Math;    

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

    boolean lic_0 () {
        return false;
    }

    boolean lic_1 () {
        if(params.RADIUS1 < 0 || numPoints < 3) return false;
        double dist1;       double dist2;       double dist3;
        for(int i = 0; i < numPoints - 2; i++){
            dist1 = Math.sqrt(Math.pow(COORDINATEX[i] - COORDINATEX[i+1] , 2) + Math.pow(COORDINATEY[i] - COORDINATEY[i+1] , 2));
            dist2 = Math.sqrt(Math.pow(COORDINATEX[i] - COORDINATEX[i+2] , 2) + Math.pow(COORDINATEY[i] - COORDINATEY[i+2] , 2));
            dist3 = Math.sqrt(Math.pow(COORDINATEX[i+1] - COORDINATEX[i+2] , 2) + Math.pow(COORDINATEY[i+1] - COORDINATEY[i+2] , 2));
            if(dist1 > 2 * params.RADIUS1 || dist2 > 2 * params.RADIUS1 || dist3 > 2 * params.RADIUS1) return true; 
        }
        return false;
    }

    boolean lic_2 () {
        return false;
    }

    boolean lic_3 () {
        return false;
    }

    boolean lic_4 () {
        return false;
    }
    /*Method for LIC 5. This method checks if there are consecutive data points, (X[i],Y[i]) and (X[j],Y[j]),
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

    boolean lic_6 () {
        return false;
    }

    boolean lic_7 () {
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
    
    boolean lic_9 () {
        return false;
    }

    boolean lic_10 () {
        return false;
    }

    boolean lic_11 () {
        return false;
    }

    boolean lic_12 () {
        return false;
    }

    boolean lic_13 () {
        return false;
    }
    
    boolean lic_14 () {
       return false; 
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
        
    }
}