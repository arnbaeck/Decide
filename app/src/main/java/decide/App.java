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
        double x1; double y1;       double x2; double y2;       double x3; double y3;  
        double epsilon = 0.000001;
        for(int i = 0; i < numPoints - 2; i++){
            x1 = COORDINATEX[i]; y1 = COORDINATEY[i];       x2 = COORDINATEX[i+1]; y2 = COORDINATEY[i+1];
            x3 = COORDINATEX[i+2]; y3 = COORDINATEY[i+2];  
            double radius = findRadius(x1, y1, x2, y2, x3, y3);
            if(epsilon < (radius - params.RADIUS1) ) return true;
            
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

    // Should call on all the LIC-functions.
    // Should check the condition for numPoints (2 ≤ NUMPOINTS ≤ 100)
    // and other similar conditions should be checked here if they do not meet return false  
    void DECIDE () {
        
    }
}