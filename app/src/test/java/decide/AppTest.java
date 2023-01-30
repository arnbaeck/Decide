/*
 * This Java source file was generated by the Gradle 'init' task.
 */
package decide;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class AppTest {
    @Test void lic_1Test(){
        int numPoints = 10;
        double[] COORDINATEX = {1.1, -2.9, 2.1, 4.1, 5.1, 8.1, -6.9, 3.1, -1.9, 1.1};
        double[] COORDINATEY = {2.1, 3.1, -4.9, 2.1, 5.1, 4.1, -6.9, 2.1, -3.9, -5.9};
        App.Parameters params = new App.Parameters();       params.RADIUS1 = 9.30;
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_1();
        assertEquals(true, actualValue);
    }
    
    @Test void lic_1Test2(){
        int numPoints = 10;
        double[] COORDINATEX = {1.1, -2.9, 2.1, 4.1, 5.1, 8.1, -6.9, 3.1, -1.9, 1.1};
        double[] COORDINATEY = {2.1, 3.1, -4.9, 2.1, 5.1, 4.1, -6.9, 2.1, -3.9, -5.9};
        App.Parameters params = new App.Parameters();       params.RADIUS1 = 9.31;
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_1();
        assertEquals(false, actualValue);
    }

    @Test void lic_5Test1(){
        int numPoints = 10;
        double[] COORDINATEX = {1.5, 3.2, 2.5, 4.5, 1.5, -2.5, -6.9, 3.1, -1.9, 1.1};
        double[] COORDINATEY = {2.1, 3.1, -4.9, 2.1, 5.1, 4.1, -6.9, 2.1, -3.9, -5.9};
        App.Parameters params = new App.Parameters();       params.RADIUS1 = 9.31;
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_5();
        assertEquals(true, actualValue);
    }
    
    @Test void lic_5Test2(){
        int numPoints = 10;
        double[] COORDINATEX = {1.5, 3.2, 4.5, 5.5, 6.5, 6.5, 6.9, 7.1, 8, 9.1};
        double[] COORDINATEY = {2.1, 3.1, -4.9, 2.1, 5.1, 4.1, -6.9, 2.1, -3.9, -5.9};
        App.Parameters params = new App.Parameters();       params.RADIUS1 = 9.31;
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_5();
        assertEquals(false, actualValue);
    }

    @Test void lic_13Test1(){
        int numPoints = 5;
        double[] COORDINATEX = {0, 0, 2, 4, 3};
        double[] COORDINATEY = {0, 0, 1, 4, 2};
        App.Parameters params = new App.Parameters();       params.RADIUS1 = 5;
        params.RADIUS2 = 6; params.A_PTS = 1;   params.B_PTS = 1;
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_13();
        assertEquals(true, actualValue);
    }

}
