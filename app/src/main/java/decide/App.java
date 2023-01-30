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