/*
 * This Java source file was generated by the Gradle 'init' task.
 */
package decide;

import org.junit.jupiter.api.Test;

import decide.App.Connectors;

import static org.junit.jupiter.api.Assertions.*;

class AppTest {
    @Test void lic_0Test(){
        int numPoints = 2;
        double[] COORDINATEX = {2.1, 0.1};
        double[] COORDINATEY = {3.3, 1.1};
        App.Parameters params = new App.Parameters();       params.RADIUS1 = 9.31; params.LENGTH1 = 1;
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_0();
        assertEquals(true, actualValue);
    }

    @Test void lic_0Test2(){
        int numPoints = 2;
        double[] COORDINATEX = {1.1, 2.9};
        double[] COORDINATEY = {2.3, 2.1};
        App.Parameters params = new App.Parameters();       params.RADIUS1 = 9.31; params.LENGTH1 = 2;
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_0();
        assertEquals(false, actualValue);
        }

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


    @Test void lic_1Test3(){
        int numPoints = 3;
        double[] COORDINATEX = {1.1, 2.1, 3.1};
        double[] COORDINATEY = {1.1, 2.1, 3.1};
        App.Parameters params = new App.Parameters();       params.RADIUS1 = 1.42;
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_1();
        assertEquals(false, actualValue);
    }

    @Test void lic_1Test4(){
        int numPoints = 3;
        double[] COORDINATEX = {1.1, 2.1, 3.1};
        double[] COORDINATEY = {1.1, 2.1, 3.1};
        App.Parameters params = new App.Parameters();       params.RADIUS1 = 1.40;
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_1();
        assertEquals(true, actualValue);
    }

    // Test that an acceptable set finds three consecutive data points with an angle within the span.
    @Test void lic_2Test1() {
        int numPoints = 5;
        double[] COORDINATEX = {1.1, -2.9, 2.1, 4.1, 5.1};
        double[] COORDINATEY = {2.1, 3.1, -4.9, 2.1, 5.1};

        App.Parameters params = new App.Parameters();
        params.EPSILON = 2;

        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_2();
        assertEquals(true, actualValue);
    }

    // Test that an acceptable set does not find three consecutive data points creating an angle within the span.
    // Also tests that a set which has two coinciding points returns false.
    @Test void lic_2Test2() {
        int numPoints = 3;
        double[] COORDINATEX = {1.6, 1.66, 1.46};
        double[] COORDINATEY = {-3.7, -3.83, -3.36};

        App.Parameters params = new App.Parameters();
        params.EPSILON = 3;

        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_2();
        assertEquals(false, actualValue);

        // Test that a set of points where one set coincides with the vertex, and
        // contains no other valid sets, returns false.
        COORDINATEX[2] = 1.66;
        COORDINATEY[2] = -3.83;

        testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        actualValue = testInstance.lic_2();
        assertEquals(false, actualValue);
    }

    @Test void lic_3Test(){
        int numPoints = 4;
        double[] COORDINATEX = {1.1, -2.9, 2.1, 4.1};
        double[] COORDINATEY = {2.1, 3.1, -4.9, 2.1};
        App.Parameters params = new App.Parameters();       params.AREA1 = 25.49;
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_3();
        assertEquals(true, actualValue);
    }


    @Test void lic_3Test2(){
        int numPoints = 4;
        double[] COORDINATEX = {1.1, -2.9, 2.1, 4.1};
        double[] COORDINATEY = {2.1, 3.1, -4.9, 2.1};
        App.Parameters params = new App.Parameters();       params.AREA1 = 25.51;
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_3();
        assertEquals(false, actualValue);
    }

    @Test void lic_3Test3(){
        int numPoints = 3;
        double[] COORDINATEX = {1.1, 2.9, 4.1};
        double[] COORDINATEY = {1.1, 2.9, 4.1};
        App.Parameters params = new App.Parameters();       params.AREA1 = 0;
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_3();
        assertEquals(false, actualValue);
    }

    // A test that finds a set of three datapoints which span over three quadrants.
    @Test void lic_4Test1(){
        int numPoints = 5;
        double[] COORDINATEX = {1.5, -3.2, 4.5, 5.5, 6.5};
        double[] COORDINATEY = {2.1, -3.1, -4.9, 2.1, 5.1};
        App.Parameters params = new App.Parameters();       params.QUADS = 2;       params.Q_PTS = 3;
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_4();
        assertEquals(true, actualValue);
    }

    // A test that does not contain a set of three datapoints which span over three quadrants.
    @Test void lic_4Test2(){
        int numPoints = 5;
        double[] COORDINATEX = {1.5, 3.2, 4.5, 5.5, 6.5};
        double[] COORDINATEY = {2.1, 3.1, -4.9, 2.1, 5.1};
        App.Parameters params = new App.Parameters();       params.QUADS = 2;       params.Q_PTS = 3;
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_4();
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

     @Test void lic_6Test1(){
        int numPoints = 3;
        double[] COORDINATEX = {20,3 ,-5};
        double[] COORDINATEY = {20,-1 ,-8};
        App.Parameters params = new App.Parameters();       params.DIST = 1.3;   params.N_PTS = 3;
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_6();
        assertEquals(true, actualValue);
    }

    @Test void lic_6Test2(){
        int numPoints = 10;
        double[] COORDINATEX = {5, -3, -5, 0, 1, 1, 2, 4, 2, 2.4};
        double[] COORDINATEY = {-5, 1, 4, 2, 2, 1, 2, 4, 0, -3};
        App.Parameters params = new App.Parameters();       params.DIST = 4.5;   params.N_PTS = 5;
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_6();
        assertEquals(false, actualValue);
    }

    @Test void lic_7Test(){
        int numPoints = 2;
        double[] COORDINATEX = {2.1, 4.1};
        double[] COORDINATEY = {3.3, 2.1};
        App.Parameters params = new App.Parameters();  params.LENGTH1 = 1; params.K_PTS = 1;
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_7();
        assertEquals(false, actualValue);
    }
    // params.K_PTS > numPoints - 2
        
    @Test void lic_7Test2(){
        int numPoints = 4;
        double[] COORDINATEX = {2.1, 4.1, 5.1, 8.1};
        double[] COORDINATEY = {3.3, 2.1, 5.1, 4.1};
        App.Parameters params = new App.Parameters(); params.RADIUS1 = 9.31; params.LENGTH1 = 1; params.K_PTS = 3;
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_7();
        assertEquals(false, actualValue);
    }
    
    @Test void lic_7Test3(){
        int numPoints = 6;
        double[] COORDINATEX = {2.1, 4.1, 5.1, 8.1, -6.9, 0.1};
        double[] COORDINATEY = {3.3, 2.1, 5.1, 4.1, -6.9, 1.1};
        App.Parameters params = new App.Parameters(); params.RADIUS1 = 9.31; params.LENGTH1 = 1; params.K_PTS = 4;
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_7();
        assertEquals(true, actualValue);
    }
    
        //numPoints < 5
    @Test void lic_8Test(){
          int numPoints = 2;
          double[] COORDINATEX = {2.1, 4.1};
          double[] COORDINATEY = {3.3, 2.1};
          App.Parameters params = new App.Parameters();  params.LENGTH1 = 1; params.K_PTS = 1;
          App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
          boolean actualValue = testInstance.lic_8();
          assertEquals(false, actualValue);
    }
    
        //params.A_PTS + params.B_PTS > numPoints - 3
    @Test void lic_8Test2(){
          int numPoints = 8;
          double[] COORDINATEX = {40.0, 5.0, 7.5, 0.0, 9.1, 1.8, 6.7, -40.0};
          double[] COORDINATEY = {0.0, 3.0, 4.3, 40.0, 14.2, 4.9, 1.3, 0.0};
          App.Parameters params = new App.Parameters(); params.RADIUS1 = 50; params.A_PTS = 2; params.B_PTS = 3;
          App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
          boolean actualValue = testInstance.lic_8();
          assertEquals(false, actualValue);
    }
    @Test void lic_8Test3(){
          int numPoints = 8;
          double[] COORDINATEX = {40.0, 5.0, 7.5, 0.0, 9.1, 1.8, 6.7, -40.0};
          double[] COORDINATEY = {0.0, 3.0, 4.3, 40.0, 14.2, 4.9, 1.3, 0.0};
          App.Parameters params = new App.Parameters(); params.RADIUS1 = 20; params.A_PTS = 2; params.B_PTS = 3;
          App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
          boolean actualValue = testInstance.lic_8();
          assertEquals(true, actualValue);
    }

    @Test void lic_9Test1(){
        int numPoints = 5;
        double[] COORDINATEX = {0, 4, 0, 4, 10};
        double[] COORDINATEY = {10, 4, 0, 4, 0};
        App.Parameters params = new App.Parameters();       params.EPSILON = 1;
        params.C_PTS = 1;   params.D_PTS = 1;
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_9();
        assertEquals(true, actualValue);
    }

    @Test void lic_9Test2(){
        int numPoints = 10;
        double[] COORDINATEX = {-5, -3, -5, 0, 1, 1, 2, 4, 2, 2.4};
        double[] COORDINATEY = {-5, 1, 4, 2, 2, 1, 2, 4, 0, -3};
        App.Parameters params = new App.Parameters();       params.EPSILON = 150;
        params.C_PTS = 2;   params.D_PTS = 3;
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_9();
        assertEquals(false, actualValue);
    }

    // A test that finds a set of three data points, separated by E_pts and F_PTS respectively, whose triangle area
    // is greater than AREA1.
    @Test void lic_10Test1(){
        int numPoints = 5;
        double[] COORDINATEX = {-1.3, 3.2, 2.5, 4.5, 1.5};
        double[] COORDINATEY = {0.8, 3.1, -4.9, 2.1, 5.1};
        App.Parameters params = new App.Parameters();       params.E_PTS = 1;      params.F_PTS = 1;
        params.AREA1 = 15;
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_10();
        assertEquals(true, actualValue);
    }

    // A test that does not find a set of three data points, separated by E_pts and F_PTS respectively, 
    // whose triangle area is greater than AREA1.
    @Test void lic_10Test2(){
        int numPoints = 5;
        double[] COORDINATEX = {-1.3, 3.2, 2.5, 4.5, 1.5};
        double[] COORDINATEY = {0.8, 3.1, -4.9, 2.1, 5.1};
        App.Parameters params = new App.Parameters();       params.E_PTS = 1;      params.F_PTS = 1;
        params.AREA1 = 20;
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_10();
        assertEquals(false, actualValue);
    }

        //numPoints < 3
    @Test void lic_11Test(){
        int numPoints = 2;
        double[] COORDINATEX = {2.1, 4.1};
        double[] COORDINATEY = {3.3, 2.1};
        App.Parameters params = new App.Parameters();  params.LENGTH1 = 1; params.K_PTS = 1;
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_11();
        assertEquals(false, actualValue);
    }
    @Test void lic_11Test2(){
        int numPoints = 6;
        double[] COORDINATEX = {2.1, 4.1, 5.1, 8.1, -6.9, 0.1};
        double[] COORDINATEY = {3.3, 2.1, 5.1, 4.1, -6.9, 1.1};
        App.Parameters params = new App.Parameters(); params.G_PTS = 3;
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_11();
        assertEquals(true, actualValue);
    }

    // A test which finds at least one set of two data points, separated by K_PTS, 
    // whose length is greater than LENGTH1. Also contains a set of two data points,
    // separated by K_PTS, whose length is lesser than LENGTH2.
    @Test void lic_12Test1(){
        int numPoints = 4;
        double[] COORDINATEX = {1.5, 3.2, 4.5, 5.5};
        double[] COORDINATEY = {2.1, 3.1, -4.9, 2.1};
        App.Parameters params = new App.Parameters();       params.K_PTS = 1;
        params.LENGTH1 = 7;        params.LENGTH2 = 5;
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_12();
        assertEquals(true, actualValue);
    }

    // A test which does not contain a set of two data points, separated by K_PTS, 
    // whose length is greater than LENGTH1. Neither does it contain a set of two 
    // data points, separated by K_PTS, whose length is lesser than LENGTH2.
    @Test void lic_12Test2(){
        int numPoints = 4;
        double[] COORDINATEX = {1.5, 3.2, 4.5, 5.5};
        double[] COORDINATEY = {2.1, 3.1, -4.9, 2.1};
        App.Parameters params = new App.Parameters();       params.K_PTS = 1;
        params.LENGTH1 = 10;        params.LENGTH2 = 2;
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_12();
        assertEquals(false, actualValue);
    }

    @Test void lic_13Test1(){
        int numPoints = 5;
        double[] COORDINATEX = {0, 0, 2, 4, 3};
        double[] COORDINATEY = {0, 0, 1, 4, 2};
        App.Parameters params = new App.Parameters();       params.RADIUS1 = 1.79;
        params.RADIUS2 = 6; params.A_PTS = 1;   params.B_PTS = 1;
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_13();
        assertEquals(true, actualValue);
    }

    @Test void lic_14Test(){
        int numPoints = 9;
        double[] COORDINATEX = {-2.9, -2.9, 1.1, 2.1, 2.1, 1.1, 1.1, 1.1, 4.1, 3.2};
        double[] COORDINATEY = {3.1,  3.1, 1.1, 8.1, -4.9, 1.1, 1.1, 1.1, 2.1, 5.3};
        App.Parameters params = new App.Parameters();       params.AREA1 = 25.49;       params.AREA2 = 15.01;
        params.E_PTS = 2;       params.F_PTS = 3;
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_14();
        assertEquals(true, actualValue);
    }

    @Test void lic_14Test2(){
        int numPoints = 9;
        double[] COORDINATEX = {-2.9, -2.9, 1.1, 2.1, 2.1, 1.1, 1.1, 1.1, 4.1, 3.2};
        double[] COORDINATEY = {3.1,  3.1, 1.1, 8.1, -4.9, 1.1, 1.1, 1.1, 2.1, 5.3};
        App.Parameters params = new App.Parameters();       params.AREA1 = 25.49;       params.AREA2 = 14.99;
        params.E_PTS = 2;       params.F_PTS = 3;
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_14();
        assertEquals(false, actualValue);
    }

    @Test void lic_14Test3(){
        int numPoints = 9;
        double[] COORDINATEX = {-2.9, -2.9, 1.1, 2.1, 2.1, 1.1, 1.1, 1.1, 4.1, 3.2};
        double[] COORDINATEY = {3.1,  3.1, 1.1, 8.1, -4.9, 1.1, 1.1, 1.1, 2.1, 5.3};
        App.Parameters params = new App.Parameters();       params.AREA1 = 25.51;       params.AREA2 = 15.01;
        params.E_PTS = 2;       params.F_PTS = 3;
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_14();
        assertEquals(false, actualValue);
    }

    @Test void lic_14Test4(){
        int numPoints = 9;
        double[] COORDINATEX = {-2.9, -2.9, 1.1, 2.1, 2.1, 1.1, 1.1, 1.1, 4.1, 3.2};
        double[] COORDINATEY = {3.1,  3.1, 1.1, 8.1, -4.9, 1.1, 1.1, 1.1, 2.1, 5.3};
        App.Parameters params = new App.Parameters();       params.AREA1 = 25.49;       params.AREA2 = 14.99;
        params.E_PTS = 2;       params.F_PTS = 3;
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        boolean actualValue = testInstance.lic_14();
        assertEquals(false, actualValue);
    }

    // Checks that an input with valid tests are launched.
    @Test void decide_test1() {
        int numPoints = 28;
        double[] COORDINATEX = {2.1, 0.1, 1.1, -2.9, 2.1, 4.1, 5.1, 1.5, -3.2, 4.5, 5.5, 6.5, 1.5, 3.2, 2.5, 4.5, 1.5, -2.5, -6.9, 3.1, -1.9, 1.1, 2.1, 4.1, 5.1, 8.1, -6.9, 0.1};
        double[] COORDINATEY = {3.3, 1.1, 2.1, 3.1, -4.9, 2.1, 5.1, 2.1, 3.1, -4.9, 2.1, 5.1, 2.1, 3.1, -4.9, 2.1, 5.1, 4.1, -6.9, 2.1, -3.9, -5.9, 3.3, 2.1, 5.1, 4.1, -6.9, 1.1};
        App.Parameters params = new App.Parameters();
        params.RADIUS1 = 9.31;
        params.LENGTH1 = 1;
        params.EPSILON = 2;
        params.QUADS = 2;
        params.K_PTS = 4;
        params.Q_PTS = 3;
        params.G_PTS = 3;
        boolean[] CMV = new boolean[15];
        boolean[] PUV = {true, false, true, false, true, true, false, true, true, false, false, true, false, false, false};
        App.Connectors[][] LCM = new App.Connectors[15][15];
        for (int i = 0; i < 15; i++) {
            for (int j = 0; j < 15; j++) {
                if (PUV[i] || PUV[j]) {
                    LCM[i][j] = Connectors.ORR;
                } else {
                    LCM[i][j] = Connectors.NOTUSED777;
                }
            }
        }

        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, CMV, params, PUV, LCM);
        testInstance.DECIDE();
        String actualValue = testInstance.LAUNCH;
        assertEquals("YES", actualValue);
    }

    // Checks that an invalid input is not handled.
    @Test void decide_test2() { 
        int numPoints = 1;
        double[] COORDINATEX = {1.2};
        double[] COORDINATEY = {3.5};
        App.Parameters params = new App.Parameters();
        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, null, params, null, null);
        testInstance.DECIDE();
        String actualValue = testInstance.LAUNCH;
        assertEquals("", actualValue);
    }

    //Checks that an input with invalid tests are not launched.
    @Test void decide_test3() {
        int numPoints = 28;
        double[] COORDINATEX = {2.1, 0.1, 1.1, -2.9, 2.1, 4.1, 5.1, 1.5, -3.2, 4.5, 5.5, 6.5, 1.5, 3.2, 2.5, 4.5, 1.5, -2.5, -6.9, 3.1, -1.9, 1.1, 2.1, 4.1, 5.1, 8.1, -6.9, 0.1};
        double[] COORDINATEY = {3.3, 1.1, 2.1, 3.1, -4.9, 2.1, 5.1, 2.1, 3.1, -4.9, 2.1, 5.1, 2.1, 3.1, -4.9, 2.1, 5.1, 4.1, -6.9, 2.1, -3.9, -5.9, 3.3, 2.1, 5.1, 4.1, -6.9, 1.1};
        App.Parameters params = new App.Parameters();
        params.RADIUS1 = 9.31;
        params.LENGTH1 = 100;
        params.EPSILON = 2;
        params.QUADS = 2;
        params.K_PTS = 4;
        params.Q_PTS = 3;
        params.G_PTS = 3;
        boolean[] CMV = new boolean[15];
        boolean[] PUV = {true, false, true, false, true, true, false, true, true, false, false, true, false, false, false};
        App.Connectors[][] LCM = new App.Connectors[15][15];
        for (int i = 0; i < 15; i++) {
            for (int j = 0; j < 15; j++) {
                if (PUV[i] || PUV[j]) {
                    LCM[i][j] = Connectors.ORR;
                } else {
                    LCM[i][j] = Connectors.NOTUSED777;
                }
            }
        }

        App testInstance = new App(numPoints, COORDINATEX, COORDINATEY, CMV, params, PUV, LCM);
        testInstance.DECIDE();
        String actualValue = testInstance.LAUNCH;
        assertEquals("NO", actualValue);
    }
    
}
