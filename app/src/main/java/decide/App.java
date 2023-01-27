package decide;

import java.lang.Math;    

    class App {

    class Parameters{
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

    boolean lic_5 () {
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
    void DECIDE () {
        
    }
}