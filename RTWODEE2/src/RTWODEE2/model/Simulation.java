/*
 * Copyright (C) 2014 Robin Hankin
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package RTWODEE2.model;

import static java.lang.Math.sqrt;

/**
 * This class holds the main calculations needed to run a TWODEE simulation
 * and will run a full simulation comparable with TWODEE2013
 * 
 * Where Zalesak 1979 is mentioned in comments, this refers to:
 *  Zalesak, S. (1979). "Fully multidimensional flux-corrected transport
 *  algortihms for fluids". Journal of Computational Physics, 31, 335-362.
 * 
 * @author Stephen Denekamp and Dr. Robin Hankin
 * @version 1.0: Created based on TWODEE2013 code by Dr. Robin Hankin
 */
public class Simulation
{

    //External file objects
    protected ExternalSettingsReader settingsFile;
    protected ExternalDataHandler dataIO;
    
    //Simulation objects
    private GasCloud cloud;
    private GasCloud cloudNew; //temp cloud to storage data in between time steps
    private ArealSource source;
    private Ground ground;
    
    //genSettings.in variables (only those used in multiple methods)
    protected double maxSimulationTime; //in seconds
    protected double timeStep; //(dt in TWODEE2013)
    protected final double GRID_SIZE; //the size/scale each point on the simluation represents
    protected final double RHOA; //ambient air density
    protected final double GRAVITY;  // acceleration due to gravity
    protected final double ROUGHNESS_LENGTH; //(zRough in TWODEE2013)
    protected final double SHAPE_PARAM; //shape parameter as in Ellison and Turner.
    protected final double U_STAR;
    
    //Other Simulation variables
    private double time = 0; //real world time in simulation (in seconds)
    private int simRows; //the size of the simulation plane
    private int simCols; //the size of the simulation plane
    private final int BOUNDARY_SIZE = 3;//The size of the boundary area of the simulation plane
    private double[][] dosage; //the total dosage at each point in the simulation
    private final int[][] MONITORING_POINTS; //x and y coordinates being monitored
    
    //to be removed once separated from Entrainment
    private double[][] fhu; //high order flux - horizontal
    private double[][] fhv; //high order flux - vertical
    private double[][] flu; //low order flux - horizontal
    private double[][] flv; //low order flux - vertical

    /**
     * Constructor for the Simulation class
     * @param settingsFile reference to the external settings file
     */
    public Simulation(ExternalSettingsReader settingsFile)
    {

        this.settingsFile = settingsFile;
        
        //get variable settings from the settings file
        GRID_SIZE = this.settingsFile.getDoubleSetting("gridSize");
        maxSimulationTime = this.settingsFile.getDoubleSetting("maxSimulationTime");
        RHOA = this.settingsFile.getDoubleSetting("rhoa");
        GRAVITY = this.settingsFile.getDoubleSetting("gravity");
        ROUGHNESS_LENGTH = this.settingsFile.getDoubleSetting("roughnessLength");
        SHAPE_PARAM = this.settingsFile.getDoubleSetting("shapeParam");
        U_STAR = this.settingsFile.getDoubleSetting("uStar");
        
        //Create data handler
        dataIO = new ExternalDataHandler(settingsFile.getStringSetting("scenarioName"));
        MONITORING_POINTS = dataIO.pointsReader("points.in");
    }
    
    /**
     * Initialises the required matrices for the simulation by reading
     * them in from the external files and formating where appropriate
     * 
     */
    protected void initialSimulationSetup()
    {
        //Create cloud from files
        double[][] hVelocity = dataIO.matrixReader("horizontalVelocity.in");
        simRows = hVelocity.length;
        simCols = hVelocity[0].length;
        double[][] vVelocity = dataIO.matrixReader("verticalVelocity.in");
        double[][] height = dataIO.matrixReader("height.in");
        double[][] density = dataIO.matrixReader("density.in");
        cloud = new GasCloud
            (hVelocity, vVelocity, height, density, RHOA);
        
        //Also need an empty cloud to hold temp data during each timeStep
        cloudNew = new GasCloud
            (simRows, simCols, RHOA);
        
        //create the areal source
        double[][] sourceVelocityMatrix = dataIO.matrixReader("upwardsVelocity.in");
        double sourceDensity = settingsFile.getDoubleSetting("arealSourceDensity");
        source = new ArealSource(sourceVelocityMatrix, sourceDensity);

        //Create the gound
        double[][] groundHeight = dataIO.matrixReader("groundHeight.in");
        double[][] groundRoughness = dataIO.matrixReader("groundRoughness.in");
        ground = new Ground(groundHeight, groundRoughness);
        
        //create empty flux arrays
        fhu = new double[simRows][simCols]; //high order flux - horizontal
        fhv = new double[simRows][simCols]; //high order flux - vertical
        flu = new double[simRows][simCols]; //low order flux - horizontal
        flv = new double[simRows][simCols]; //low order flux - vertical
        
        dosage = new double[simRows][simCols];

        //Setup blank time step output files
        dataIO.createAppendFiles();
        
        /*
        Now for the ground slope. If the check value of xSlope is
        set (>999), use gen.in; otherwise set a slope of 
        (xSlope,ySlope)percent. Xslope is the theta-ex-in. 
        Same for ySlope.
        */
        double xSlope = this.settingsFile.getDoubleSetting("xSlope");
        double ySlope = this.settingsFile.getDoubleSetting("ySlope");
        for (int i = 0; i < simRows; i++) {
            for (int j = 0; j < simCols; j++) 
            { 
                if(xSlope > 999)
                {
                    //continue and read from groundHeight.in
                }
                else
                {
                    double value = (0.0+i*xSlope +j*ySlope)*GRID_SIZE*(0.01);
                    ground.setHeight(i, j, value);
                }
            }//end j loop
        }//end i loop   
    }
    
    /**
     * This method runs a full TWODEE cloud simulation
     */
    public void runSimulation()
    {
        System.out.println("Starting simulation...");
        initialSimulationSetup();        
        outputAllMatrices("");  //output initial cloud matrices 
       
        //advance through time steps
        int timeStepCounter = 1; //number of time steps passed
        int lastInterval = 0; //number of timeWrite() intervals passed
        while (time < maxSimulationTime) 
        {   
            calcTimeStepSize(); 
            
            /*
            Fill cloudNew with a copy of the cloud. This will ensure that any 
            unchanged squares are not filled with zero later on:
            */
            cloudNew.setHeight(cloud.getHeight());
            cloudNew.setDensity(cloud.getDensity());
            cloudNew.setU(cloud.getU());
            cloudNew.setV(cloud.getV());
            
            /*
            Simulation Calculations
            Comment out or change order to alter simulation
            */
            cliffs(cloud);
            iterateInteriorPoints();
            cliffs(cloudNew);
            calcEntrainment(); 
            calcHydrostaticDifference(); 
            updateArealSource();
            cliffs(cloudNew);
            updateBoundary();
            
            
            //update cloud with completed cloudNew data
            //Make sure height is not negative
            for (int i = 0; i < simRows; i++) {
                for (int j = 0; j < simCols; j++) 
                { 
                    double value = Math.max(cloudNew.getHeight(i, j), 0.0);
                    cloud.setHeight(i, j, value);
                }
            }
            cloud.setDensity(cloudNew.getDensity());
            cloud.setU(cloudNew.getU());
            cloud.setV(cloudNew.getV());
            
            updateBoundary();
            calcDosage();
            
            //calculate cloud metrics for AllData file
            double sumVolume = outputAllDataFile(timeStepCounter, time);  
            outputTerminalData(sumVolume, time, timeStepCounter);
            
            //Append time step reliant files
            outputPointsData(cloud.getHeight(), "euler_height");    //height data from points
            outputPointsData(cloud.getDensity(), "euler_density");  //density data from points
        
            //Output files at set time intervals
            lastInterval = timeWrite(lastInterval);
            
            //advance the time step
            time += timeStep;
            timeStepCounter++;
        }//end loop

        outputAllMatrices(""); //Output final cloud conditions
        System.out.println("Simulation complete");
        
    }

    /**
     * Calculates the dosage 
     */
    private void calcDosage()
    {
        double toxicExponent = settingsFile.getDoubleSetting("toxicExponent");
        
        //the density of pure chlorine gas at 15 degrees Celsius
        double chlorineGasDensity = settingsFile.getDoubleSetting("chlorineGasDensity");
        
        for (int i = 0; i < simRows; i++) {
            for (int j = 0; j < simCols; j++) 
            {          
                /*
                Calculates dosage as \int (c^toxicExponent)timeStep, where c is 
                volumetric concentration (ppm) and t is in minutes 
                (that is why there is a 60 there).
                */
                dosage[i][j] = dosage[i][j] 
                        +Math.abs((Math.pow(Math.abs((cloud.getDensity(i, j)
                        -RHOA)/(chlorineGasDensity-RHOA)
                        *1000000.0),toxicExponent))*timeStep/60);  
                
            } //end i loop
        } //end j loop
    }

    /**
     * Iterates through the common-or-garden interior points of the simulation 
     * using finite difference methods following Roache and Rottman.
     */
    private void iterateInteriorPoints()
    {   
        
        /*
        The variable `w' is the variable for each of the conservation variables
        height, horizontal velocity, vertical velocity and density.
        */
        double[][] w = new double[simRows][simCols]; //holder matrix

        /*
        The first thing to do is to make variable w hold height, 
        horizontael velocity, vertical velocity and density in turn.
        Then use the routine Zalesak to work out the new ones.
        */
        ///// HEIGHT ///////////////////////////////////////////////////////////
        //get height into w matrix
        for (int i = BOUNDARY_SIZE; i < simRows - BOUNDARY_SIZE; i++) {
            for (int j = BOUNDARY_SIZE; j < simCols - BOUNDARY_SIZE; j++) {
                w[i][j] = cloud.getHeight(i, j);
            }
        }
        w = runZalesak(w, 'h');    //run runZalesak with height

        //update cloudNew height matrix
        for (int i = BOUNDARY_SIZE; i < simRows - BOUNDARY_SIZE; i++) {
            for (int j = BOUNDARY_SIZE; j < simCols - BOUNDARY_SIZE; j++) {
                double value = Math.min(w[i][j], 1000.0);
                cloudNew.setHeight(i, j, value);
            }
        }

        ///// HORIZONTAL VELOCITY //////////////////////////////////////////////
        //get  hVelocity into w variable
        for (int i = BOUNDARY_SIZE; i < simRows - BOUNDARY_SIZE; i++) {
            for (int j = BOUNDARY_SIZE; j < simCols - BOUNDARY_SIZE; j++) {
                w[i][j] = cloud.getHeight(i, j) 
                        * cloud.getU(i, j)
                        * cloud.getDensity(i, j);
            }
        }
        w = runZalesak(w, 'u'); //run runZalesak with hVelecity

        //update cloudNew horizontal velocity matrix
        for (int i = BOUNDARY_SIZE; i < simRows - BOUNDARY_SIZE; i++) {
            for (int j = BOUNDARY_SIZE; j < simCols - BOUNDARY_SIZE; j++) {
                
                double value = 0.0; //default case

                if (cloudNew.getHeight(i, j) > ROUGHNESS_LENGTH) 
                {
                    value = w[i][j] 
                            / cloudNew.getHeight(i, j)
                            / cloudNew.getDensity(i, j);
                }

                cloudNew.setU(i, j, value);

            }
        }

        ///// VERTICAL VELOCITY ////////////////////////////////////////////////
        //get  vVelocity into w variable
        for (int i = BOUNDARY_SIZE; i < simRows - BOUNDARY_SIZE; i++) {
            for (int j = BOUNDARY_SIZE; j < simCols - BOUNDARY_SIZE; j++) 
            {    
                w[i][j] = cloud.getHeight(i, j) 
                        * cloud.getV(i, j)
                        * cloud.getDensity(i, j);
            }
        }
        w = runZalesak(w, 'v'); //run runZalesak with vVelecity

        //update cloudNew vertical velocity matrix
        for (int i = BOUNDARY_SIZE; i < simRows - BOUNDARY_SIZE; i++) {
            for (int j = BOUNDARY_SIZE; j < simCols - BOUNDARY_SIZE; j++) 
            {
                double value = 0.0; //default case

                if (cloudNew.getHeight(i, j) > ROUGHNESS_LENGTH) 
                {
                    value = w[i][j] 
                            / cloudNew.getHeight(i, j)
                            / cloudNew.getDensity(i, j);
                }

                cloudNew.setV(i, j, value);
            }
        }

        ///// DENSITY //////////////////////////////////////////////////////////
        //get  density into w variable
        for (int i = BOUNDARY_SIZE; i < simRows - BOUNDARY_SIZE; i++) {
            for (int j = BOUNDARY_SIZE; j < simCols - BOUNDARY_SIZE; j++) 
            {
                w[i][j] = cloud.getHeight(i, j)
                        * (cloud.getDensity(i, j) - RHOA);
            }
        }
        w = runZalesak(w, 'd'); //run runZalesak with density

        //update cloudNew density matrix
        for (int i = BOUNDARY_SIZE; i < simRows - BOUNDARY_SIZE; i++) {
            for (int j = BOUNDARY_SIZE; j < simCols - BOUNDARY_SIZE; j++) 
            {
                double value = RHOA; //default case

                if (cloudNew.getHeight(i, j) > ROUGHNESS_LENGTH) 
                {
                    value = Math.min(RHOA + w[i][j] 
                            / cloudNew.getHeight(i, j), source.getDensity());
                    
                    value = Math.max(value, RHOA);
                  
                    //Also update the max density matrices
                    updateMaxDensities(i, j);
                }
                
                cloudNew.setDensity(i, j, value);
                
                //Don't let the density exceed the density of the source.
                if(cloudNew.getHeight(i, j) * GRAVITY 
                        * (cloudNew.getDensity(i, j) - RHOA) 
                    / (cloudNew.getU(i, j) * cloudNew.getU(i, j) 
                        + U_STAR * U_STAR + 0.01) < 0.001)
                {
                    cloudNew.setDensity(i, j, RHOA);
                }
                    
            }
        } //end of density loop
        
    }

    /**
     * Handles the boundary conditions of the cloud
     * 
     * It uses the  zero gradient boundary conditions; meaning that a ring of 
     * cells at the edge of the computational domain will be forced to have the 
     * same values (of speed,  density and depth) as at a ring of cells 
     * slightly further in from the edge. 
     * 
     * This method will be better than having zero depth at the edge 
     * (as then fluid pours over the edge of the domain like water over a 
     * cliff, dragging fluid behind it faster and faster) or zero flux through 
     * the edge (as then fluid splashes off the impenetrable barrier of the 
     * edge and then bounces back into the centre of the domain).
     * 
     * This method is (arguably) better than that recommended by Roache, ie 
     * taking some quadratic extrapolation from the interior and forcing
     * the boundary to take the extrapolated value. This method could
     * return negative heights (or, for that matter, speeds or excess densities).
     * 
     */
    private void updateBoundary()
    {
        //copy of cloud matrices to simplify look of equations
        double[][] h = cloud.getHeight();
        double[][] u = cloud.getU();
        double[][] v = cloud.getV();
        double[][] rho = cloud.getDensity();
             
        /*
        In the following, the first part of the `if' routines handle outflow 
        and the second (that is, the `else' part) handles inflow. We set the 
        inflow conditions to zero fluid depth and excess density. 
        */
        for (int i = BOUNDARY_SIZE; i < simRows - BOUNDARY_SIZE; i++) 
        {
             int num1; //value for left side of equation
             int num2; //value for right side of equation
             
             //West side of the domain
             num1 = BOUNDARY_SIZE;  //3 for a 200x200 cloud
             num2 = BOUNDARY_SIZE+1; //4 for a 200x200 cloud
            if(u[num2][i] < 0.0) 
            {
                h[num1][i] = h[num2][i];
                u[num1][i] = u[num2][i];
                v[num1][i] = v[num2][i];
                rho[num1][i] = rho[num2][i];
            }
            else
            {
                h[num1][i] = 0.0;
                u[num1][i] = 0.0;
                v[num1][i] = v[num2][i];
                rho[num1][i] = RHOA;
            }
            
            //East side of the domain
            num1 = simRows-1 - BOUNDARY_SIZE;   //196 for a 200x200 cloud
            num2 = simRows-1 - BOUNDARY_SIZE-1; //195 for a 200x200 cloud
            if(u[num2][i] > 0.0) 
            {
                h[num1][i] = h[num2][i];
                u[num1][i] = u[num2][i];
                v[num1][i] = v[num2][i];
                rho[num1][i] = rho[num2][i];
            }
            else
            {
                h[num1][i] = 0.0;
                u[num1][i] = 0.0;
                v[num1][i] = v[num2][i];
                rho[num1][i] = RHOA;
            }
            
            //North side of the domain
            num1 = simRows-1 - BOUNDARY_SIZE;   //196 for a 200x200 cloud
            num2 = simRows-1 - BOUNDARY_SIZE-1; //195 for a 200x200 cloud
            if(v[i][num2] > 0.0) 
            {
                h[i][num1] = h[i][num2];
                u[i][num1] = u[i][num2];
                v[i][num1] = v[i][num2];
                rho[i][num1] = rho[i][num2];
            }
            else
            {
                h[i][num1] = 0.0;
                u[i][num1] = u[i][num2];
                v[i][num1] = 0.0;
                rho[i][num1] = RHOA;
            }
            
            //South side of the domain
            num1 = BOUNDARY_SIZE;   //3 for a 200x200 cloud
            num2 = BOUNDARY_SIZE+1; //4 for a 200x200 cloud
            if(v[i][num2] < 0.0) 
            {
                h[i][num1] = h[i][num2];
                u[i][num1] = u[i][num2];
                v[i][num1] = v[i][num2];
                rho[i][num1] = rho[i][num2];
            }
            else
            {
                h[i][num1] = 0.0;
                u[i][num1] = u[i][num2];
                v[i][num1] = 0.0;
                rho[i][num1] = RHOA;
            }
        }//end loop
    }
    
    /**
     * Updates a particular point in the clouds maximum density variables 
     * (maxDensityAtGround, maxDensityAtAltitude)
     * @param i Row point in cloud
     * @param j Column point in cloud
     */
    private void updateMaxDensities(int i, int j)
    {
        //Check and update maxDesnityAtGround
        if (cloudNew.getDensity(i, j) > cloud.getMaxDensityAtGround(i, j)) 
            cloud.setMaxDensityAtGround(i, j, cloudNew.getDensity(i, j));

        //Check and update maxDensityAtAltitude (calculates the exponential distn)
        double recordingAltitude = settingsFile.getDoubleSetting("recordingAltitude");
        double densityAtAltitude = RHOA
                + (cloudNew.getDensity(i, j) - RHOA)
                * 2.0 / SHAPE_PARAM
                * Math.exp(0.0 - (2 * recordingAltitude 
                        / (Math.max(cloudNew.getHeight(i, j), ROUGHNESS_LENGTH)
                                * SHAPE_PARAM)));
        
        if(densityAtAltitude > cloud.getMaxDensityAtAltitude(i, j))
            cloud.setMaxDensityAtAltitude(i, j, densityAtAltitude);
    }
    /**
     * Creates an interior matrix: a matrix which is a copy of another matrix
     * but with edge points set to 0.0
     *
     * @param inputMatrix The matrix to be copied
     * @return the interior matrix created
     */
    private double[][] createInteriorMatrix(double[][] inputMatrix)
    {
        //create an empty matrix the same size as the inputMatrix
        double[][] interiorMatrix = new double[simRows][simCols];

        //set interior points of InteriorMatrix to that of the inputMatrix
        for (int i = BOUNDARY_SIZE; i < simRows - BOUNDARY_SIZE; i++) {
            for (int j = BOUNDARY_SIZE; j < simCols - BOUNDARY_SIZE; j++) {
                double value = inputMatrix[i][j];
                interiorMatrix[i][j] = value;
            }
        }

        return interiorMatrix;
    }
    

    /**
     * Fills the high order and low oder flux array variables
     * (fhu, fhv, flu, flv)
     * The fh means high order and the u or v means the u or v flux
     * @param w The matrix where the flux originates from
     * @param matrixType The matrix that w represents (h, u, v, d)
     */
    private void fillFluxArrays(double[][] w, char matrixType)
    {
        double ubar;
        double vbar;
        double wdc; 
        
        //copy of cloud matrices to simplify look of equations
        double[][] h = cloud.getHeight();
        double[][] rho = cloud.getDensity();
        double[][] u = cloud.getU();
        double[][] v = cloud.getV();
        
        for (int i = BOUNDARY_SIZE; i < simRows - BOUNDARY_SIZE; i++) {
            for (int j = BOUNDARY_SIZE; j < simCols - BOUNDARY_SIZE; j++) {
          
                //The below is from Zalesak, JCP, 40/500/1981.
                
                fhu[i][j]=(2.0/3.0*(w[i+1][j]*u[i][j]+w[i][j]*u[i+1][j])
                        -1.0/12.0*(w[i+2][j]*u[i][j]+w[i][j]*u[i+2][j]
                        +w[i+1][j]*u[i-1][j]+w[i-1][j]*u[i+1][j]))*timeStep;
                        
                fhv[i][j]=(2.0/3.0*(w[i][j+1]*v[i][j]+w[i][j]*v[i][j+1])
                        -1.0/12.0*(w[i][j+2]*v[i][j]+w[i][j]*v[i][j+2]
                        +w[i][j+1]*v[i][j-1]+w[i][j-1]*v[i][j+1]))*timeStep;
                
                if(matrixType=='d') //d = density matrix
                {
                    double densityExchangeCoeff = 
                            settingsFile.getDoubleSetting("densityExchangeCoeff");
                    
                 
                    
                    double speedi=1/sqrt(2)*sqrt(u[i][j]*u[i][j]+v[i][j]*v[i][j]
                            +u[i+1][j]*u[i+1][j]+v[i+1][j]*v[i+1][j]);
            
                    double speedj=1/sqrt(2)*sqrt(u[i][j]*u[i][j]+v[i][j]*v[i][j]
                            +u[i][j+1]*u[i][j+1]+v[i][j+1]*v[i][j+1]);
                    
                    double h_avei=1/sqrt(2)*sqrt(h[i][j]*h[i][j]+h[i+1][j]*h[i+1][j]);
                    double h_avej=1/sqrt(2)*sqrt(h[i][j]*h[i][j]+h[i][j+1]*h[i][j+1]);
                    
                    fhu[i][j]=fhu[i][j]+densityExchangeCoeff*speedi
                            *h_avei*(rho[i][j]-rho[i+1][j]);
                    
                    fhv[i][j]=fhv[i][j]+densityExchangeCoeff*speedj
                            *h_avej*(rho[i][j]-rho[i][j+1]);
                } //end if
                
                double diffusionCoeff = 
                        settingsFile.getDoubleSetting("diffusionCoeff");

                ubar = (u[i][j]+u[i+1][j])*0.5;
                if(ubar >= 0.0)
                    wdc = w[i][j];
                else
                    wdc = w[i+1][j];
                
                flu[i][j]=ubar*wdc*timeStep-GRID_SIZE
                        *(w[i+1][j]-w[i][j])*diffusionCoeff;
                
                
                vbar = (v[i][j]+v[i][j+1])*0.5;
                if(vbar >= 0.0)
                    wdc = w[i][j];
                else
                    wdc = w[i][j+1];    
     
                flv[i][j]=vbar*wdc*timeStep-GRID_SIZE
                        *(w[i][j+1]-w[i][j])*diffusionCoeff;
                
                
                if(h[i][j]+ground.getHeight(i, j) < ground.getHeight(i+1, j))
                {
                    if(w[i+1][j] >= 0.0)
                    {
                        flu[i][j]=Math.min(flu[i][j], 0.0);
                        fhu[i][j]=Math.min(fhu[i][j], 0.0);
                    } else {
                        flu[i][j]=Math.max(flu[i][j], 0.0);
                        fhu[i][j]=Math.max(fhu[i][j], 0.0);
                    }//end if  
                }//end if
   
                if(ground.getHeight(i, j) > ground.getHeight(i+1,j)+h[i+1][j])
                {
                    if(w[i+1][j] >= 0.0)
                    {
                        flu[i][j]=Math.max(flu[i][j], 0.0);
                        fhu[i][j]=Math.max(fhu[i][j], 0.0);
                    } else {
                        flu[i][j]=Math.min(flu[i][j], 0.0);
                        fhu[i][j]=Math.min(fhu[i][j], 0.0);
                    }//end if
                }//end if
            
                if(h[i][j]+ground.getHeight(i, j) < ground.getHeight(i, j+1))
                {
                    if(w[i][j+1] >= 0.0)
                    {
                        flv[i][j]=Math.min(flv[i][j], 0.0);
                        fhv[i][j]=Math.min(fhv[i][j], 0.0);
                    } else {
                        flv[i][j]=Math.max(flv[i][j], 0.0);
                        fhv[i][j]=Math.max(fhv[i][j], 0.0);
                    }//end if
                }//end if   
              
                if(ground.getHeight(i, j) > ground.getHeight(i, j+1)+h[i][j+1])
                {
                    if(w[i][j+1] >= 0.0)
                    {
                        flv[i][j]=Math.max(flv[i][j], 0.0);
                        fhv[i][j]=Math.max(fhv[i][j], 0.0);
                    } else {
                        flv[i][j]=Math.min(flv[i][j], 0.0);
                        fhv[i][j]=Math.min(fhv[i][j], 0.0);
                    }//end if
                }//end if
                
                
                /*
                Note that the second and fourth ifthen statements above refer to
                the following type of cliff:
                                        ___
                                        |cc|
                                        |cc|
                  __________________________    c = Chlorine.
                                         ^  |
                                         |  |
                                         |  |
                        (i,j)____________/  |
                                            |__
                                            |cc|
                                            |cc|
                                            |cc|           ground
                                            |______________________
                                             ^////////////////////   
                                             |
                                             |
                         (i+1,j)_____________/


                 so the flux corresponding to the (i+1/2,j) interface
                 should be less than zero. A completely similar thing happens
                 for the (i,j+1/2) interface.
                */
                   
            }//end j loop
        }//end i loop
    }
    
    /**
     * 
     * @param w The matrix that the method will be changing
     * @param matrixType The type of that matrix w is (e.g. height, density etc)
     * @return The updated w matrix
     */
    private double[][] runZalesak(double[][] w, char matrixType)
    {

        double[][] wtd = new double[simRows][simCols]; //W matrix Transported and Diffused
        double[][] pplus = new double[simRows][simCols]; //as per equation 7 of Zalesak 1979 
        double[][] pminus = new double[simRows][simCols];
        
        //copy of cloud matrices to simplify look of equations
        double[][] h = cloud.getHeight();
        double[][] rho = cloud.getDensity();
        double[][] u = cloud.getU();
        double[][] v = cloud.getV();
        double[][] gh = ground.getHeight();
        
        //fil the four arrays flu, flv, fhu and fhv
        fillFluxArrays(w, matrixType);
        
        for (int i = BOUNDARY_SIZE; i < simRows - BOUNDARY_SIZE; i++) {
            for (int j = BOUNDARY_SIZE; j < simCols - BOUNDARY_SIZE; j++) 
            {
              
                wtd[i][j]=w[i][j]+(flu[i-1][j]-flu[i][j]+flv[i][j-1]-flv[i][j])
                        /GRID_SIZE;

                pplus[i][j]=Math.max(fhu[i-1][j]-flu[i-1][j], 0.0)
                        -Math.min(fhu[i][j]-flu[i][j], 0.0)
                        +Math.max(fhv[i][j-1]-flv[i][j-1], 0.0)
                        -Math.min(fhv[i][j]-flv[i][j], 0.0);

                pminus[i][j]=0.0-Math.min(fhu[i-1][j]-flu[i-1][j], 0.0)
                        +Math.max(fhu[i][j]-flu[i][j], 0.0)
                        -Math.min(fhv[i][j-1]-flv[i][j-1], 0.0)
                        +Math.max(fhv[i][j]-flv[i][j], 0.0);
  
            }//end j loop
        }//end i loop
        
        /*
        We now have to test whether or not this routine has been called
        from the height part or not. If it has, then the wMax and wMin
        calculations will be different as the ground height gh(i,j) must
        be taken into consideration. 
        */
        double[][] wMax = new double[simRows][simCols];
        double[][] wMin = new double[simRows][simCols];
        
        /*
        This routine will use the fluid depth plus ground height for
        calculating wMax and wMin for w=h (height). 
        Note that this routine is conceptually more difficult than
        the other one. This is because we have to consider fluid 
        elements of types, sorry, {\em cases} "C" and "D"; that is,
        the cliff type elements. The wmax and wmin of type "C" and "D"
        elements should be unaffected by the height of either the 
        cliff or the fluid on top of it. 
        */
        if(matrixType == 'h') //called from height part
        {
            for (int i = BOUNDARY_SIZE; i < simRows - BOUNDARY_SIZE; i++) {
                for (int j = BOUNDARY_SIZE; j < simCols - BOUNDARY_SIZE; j++) 
                {
                    
                    //set the upper and lower limits of w
                    wMax[i][j]=Math.max(wtd[i][j],w[i][j]) +gh[i][j];
                    wMin[i][j]=Math.min(wtd[i][j],w[i][j]) +gh[i][j];
                    
                    if(h[i+1][j]+gh[i+1][j]-gh[i][j] > 0.0 &&
                            h[i][j]+gh[i][j]-gh[i+1][j] > 0.0)
                    {
                        wMax[i][j]=Math.max(w[i+1][j]+gh[i+1][j], wMax[i][j]);
                        wMin[i][j]=Math.min(w[i+1][j]+gh[i+1][j], wMin[i][j]);
                    }
                    
                    if(h[i-1][j]+gh[i-1][j]-gh[i][j] > 0.0 &&
                            h[i][j]+gh[i][j]-gh[i-1][j] > 0.0)
                    {
                        wMax[i][j]=Math.max(w[i-1][j]+gh[i-1][j], wMax[i][j]);
                        wMin[i][j]=Math.min(w[i-1][j]+gh[i-1][j], wMin[i][j]);
                    }
                    
                    if(h[i][j+1]+gh[i][j+1]-gh[i][j] > 0.0 &&
                            h[i][j]+gh[i][j]-gh[i][j+1] > 0.0)
                    {
                        wMax[i][j]=Math.max(w[i][j+1]+gh[i][j+1], wMax[i][j]);
                        wMin[i][j]=Math.min(w[i][j+1]+gh[i][j+1], wMin[i][j]);
                    }
                    
                    if(h[i][j-1]+gh[i][j-1]-gh[i][j] > 0.0 &&
                            h[i][j]+gh[i][j]-gh[i][j-1] > 0.0)
                    {
                        wMax[i][j]=Math.max(w[i][j-1]+gh[i][j-1], wMax[i][j]);
                        wMin[i][j]=Math.min(w[i][j-1]+gh[i][j-1], wMin[i][j]);
                    }
                    
                    /*
                    Now we have to subtract the gh's which have been adding
                    to the arguments of all those 'min' and 'max' functions
                    */
                    wMax[i][j]=wMax[i][j]-gh[i][j];
                    wMin[i][j]=wMin[i][j]-gh[i][j];
                    
                    /*
                    Eliminate negative values. This is acceptable as
                    arithmetic errors in subtraction may return slightly
                    negative arguments
                    */
                    wMin[i][j]=Math.max(wMin[i][j], 0.0);       
                }//end j loop
            }//end i loop            
        } 
        else //If w is not height then wMax and wMin will be calculated as below        
        {     
            /*
            Briefly, the following lines take the max and min of non-cliff
            adjacent elements; they say ``if this one is {\em not} 
            a cliff, then take the max of wMax and that element.
            */
            for (int i = BOUNDARY_SIZE; i < simRows - BOUNDARY_SIZE; i++) {
                for (int j = BOUNDARY_SIZE; j < simCols - BOUNDARY_SIZE; j++) 
                {
                    //set the upper and lower limits of w
                    wMax[i][j]=Math.max(wtd[i][j], w[i][j]);
                    wMin[i][j]=Math.min(wtd[i][j], w[i][j]);
                    
                    if(h[i+1][j]+gh[i+1][j] > gh[i][j] &&
                            h[i][j]+gh[i][j] > gh[i+1][j] &&
                            h[i+1][j] > 0.0)
                    {
                        wMax[i][j]=Math.max(w[i+1][j]
                                *(1+(gh[i+1][j]-gh[i][j])
                                /h[i+1][j]), wMax[i][j]);
                            
                        wMin[i][j]=Math.min(w[i+1][j]
                                *(1+(gh[i+1][j]-gh[i][j])
                                /h[i+1][j]), wMin[i][j]);         
                    }
                    
                    if(h[i-1][j]+gh[i-1][j] > gh[i][j] &&
                            h[i][j]+gh[i][j] > gh[i-1][j] &&
                            h[i-1][j] > 0.0)       
                    {
                        wMax[i][j]=Math.max(w[i-1][j]
                                *(1+(gh[i-1][j]-gh[i][j])/h[i-1][j]), wMax[i][j]);
                            
                        wMin[i][j]=Math.min(w[i-1][j]
                                *(1+(gh[i-1][j]-gh[i][j])/h[i-1][j]), wMin[i][j]);     
                    }
                    
                    
                    if(h[i][j+1]+gh[i][j+1] > gh[i][j] &&
                            h[i][j]+gh[i][j] > gh[i][j+1] &&
                            h[i][j+1] > 0.0)
                    {
                        wMax[i][j]=Math.max(w[i][j+1]
                                *(1+(gh[i][j+1]-gh[i][j])/h[i][j+1]), wMax[i][j]);
                        
                        wMin[i][j]=Math.min(w[i][j+1]
                                *(1+(gh[i][j+1]-gh[i][j])/h[i][j+1]), wMin[i][j]);        
                    }
                    
                    if(h[i][j-1]+gh[i][j-1] > gh[i][j] &&
                            h[i][j]+gh[i][j] > gh[i][j-1] &&
                            h[i][j-1] > 0.0)
                    {
                        wMax[i][j]=Math.max(w[i][j-1]
                                  *(1+(gh[i][j-1]-gh[i][j])/h[i][j-1]), wMax[i][j]);
                          
                        wMin[i][j]=Math.min(w[i][j-1]
                                  *(1+(gh[i][j-1]-gh[i][j])/h[i][j-1]), wMin[i][j]);     
                    }
               
                }//end j loop
            }//end i loop
        }//end if
        
        
        double[][] rplus = new double[simRows][simCols];
        double[][] rminus = new double[simRows][simCols];
        double[][] cu = new double[simRows][simCols]; //as per equation 13 in Zalesak 1979
        double[][] cv = new double[simRows][simCols];
        
        for (int i = BOUNDARY_SIZE; i < simRows - BOUNDARY_SIZE; i++) {
            for (int j = BOUNDARY_SIZE; j < simCols - BOUNDARY_SIZE; j++) 
            {

                if(pplus[i][j] > 0.0)              
                    rplus[i][j]=Math.min(1.0, (wMax[i][j]-wtd[i][j])
                            *GRID_SIZE/pplus[i][j]);                         
                 else 
                    rplus[i][j]=0.0;

                if(pminus[i][j] > 0.0)
                    rminus[i][j]=Math.min(1.0, (wtd[i][j]-wMin[i][j])
                            *GRID_SIZE/pminus[i][j]);
                else
                    rminus[i][j]=0.0;
                
            }//end j loop            
        }//end i loop
                
        for (int i = BOUNDARY_SIZE; i < simRows - BOUNDARY_SIZE; i++) {
            for (int j = BOUNDARY_SIZE; j < simCols - BOUNDARY_SIZE; j++) 
            {
                //Horizontal Flux
                if(fhu[i][j]-flu[i][j] >= 0.0)
                    cu[i][j]=Math.min(rplus[i+1][j], rminus[i][j]);
                else
                    cu[i][j]=Math.min(rplus[i][j], rminus[i+1][j]);

                //Vertical Flux
                if(fhv[i][j]-flv[i][j] >= 0.0)
                    cv[i][j]=Math.min(rplus[i][j+1], rminus[i][j]);
                else
                    cv[i][j]=Math.min(rplus[i][j], rminus[i][j+1]);
                
            }//end j loop          
        }//end i loop  
        
        for (int i = BOUNDARY_SIZE; i < simRows - BOUNDARY_SIZE; i++) {
            for (int j = BOUNDARY_SIZE; j < simCols - BOUNDARY_SIZE; j++) 
            {
                
                double a;
                
                //The three parts to equation 14 of Zalaesak 1979
                double one, two, three; 
                
                //Horizontal Flux
                a = fhu[i][j]-flu[i][j];
                one = a * (wtd[i+1][j]-wtd[i][j]);
                two = a * (wtd[i+2][j]-wtd[i+1][j]);
                three = a * (wtd[i][j]-wtd[i-1][j]);
                
                if(one < 0.0 && (two < 0.0 || three < 0.0))
                    cu[i][j] = 0.0;
                
                //Vertical Flux
                a = fhv[i][j]-flv[i][j];
                one = a * (wtd[i][j+1]-wtd[i][j]);
                two = a * (wtd[i][j+2]-wtd[i][j+1]);
                three = a * (wtd[i][j]-wtd[i][j-1]);
                
                if(one < 0.0 && (two < 0.0 || three < 0.0))
                    cv[i][j] = 0.0;
                
            }//end j loop          
        }//end i loop
                
        for (int i = BOUNDARY_SIZE; i < simRows - BOUNDARY_SIZE; i++) {
            for (int j = BOUNDARY_SIZE; j < simCols - BOUNDARY_SIZE; j++) 
            {
                w[i][j]=wtd[i][j]+((fhu[i-1][j]-flu[i-1][j])*cu[i-1][j]
                        -(fhu[i][j]-flu[i][j])*cu[i][j]
                        +(fhv[i][j-1]-flv[i][j-1])*cv[i][j-1]
                        -(fhv[i][j]-flv[i][j])*cv[i][j])              
                        /GRID_SIZE;    
            }//end j loop          
        }//end i loop
              
        return w;
    }

    /**
     * Appends key simulation/cloud data to the all.out file
     *
     * @param counter The simulation time step counter
     * @param time The actual time in the simulation
     */
    private double outputAllDataFile(int counter, double time)
    {
        dataIO.appendDataToFile("allData", String.format("%7d", counter)); 
        dataIO.appendDataToFile("allData", String.format("%12.3f", time));

        //Cloud metrics data
        double[] cloudMetrics = cloud.calcMetrics(ROUGHNESS_LENGTH, GRID_SIZE);
        for (double metric : cloudMetrics) 
            dataIO.appendDataToFile("allData", String.format("%12.3f", metric));
        
        dataIO.appendNewLineToFile("allData");
        
        return cloudMetrics[0]; //return the sumVolume metric
    }

    /**
     * Outputs to file the data recording at monitoring points for a given
     * matrix
     *
     * @param matrix The matrix to get the data from
     * @param filename The name of the file being written to
     */
    private void outputPointsData(double[][] matrix, String filename)
    {
        
        int x;
        int y;
        

        for (int i = 0; i < MONITORING_POINTS.length; i++)
        {
            x = MONITORING_POINTS[i][0]; //get x coordinate (row)
            y = MONITORING_POINTS[i][1]; //get y coordinate (column)
 
            dataIO.appendDataToFile(filename, String.format("%8.4f", matrix[y][x]));
        }
        //add new line
        dataIO.appendNewLineToFile(filename);
    }

    /**
     * Writes cloud matrices to file after a set time has passed as specified 
 by the user in the genSettings.in fin (outputIntervalTime variable).
     * Filenames are appended with the intervalCounter number so as to not
     * overwrite existing files.
     * 
     * @param lastInterval The number of the last interval that was written
     */
    private int timeWrite(int lastInterval)
    {
        int outputIntervalTime = this.settingsFile.getIntegerSetting("outputIntervalTime");
        int currentInterval = (int)time/outputIntervalTime;
  
        //output initial interval
        if (time == 0)
        {
            outputAllMatrices("00");
            System.out.println(" - Tme inteveral outputs written");
        }
        //output if new interval reaached
        else if (currentInterval > lastInterval)
        {
            //format filenames correctly
            int leftDigit= currentInterval/10;
            int rightDigit=currentInterval-leftDigit*10;

            //output files
            outputAllMatrices(leftDigit+""+rightDigit);
            System.out.println(" - Tme inteveral outputs written");
        }
        
        return currentInterval;
    }

    /**
     * Outputs all matrix files of the simulation
     *
     * @param filenameAddition Any additional String to be added to the filename
     * this is used if you do not want to overwrite the file.
     */
    protected void outputAllMatrices(String filenameAddition)
    {    
        dataIO.matrixWriter(cloud.getU(), "horizontalVelocity" + filenameAddition);
        dataIO.matrixWriter(cloud.getV(), "verticalVelocity" + filenameAddition);
        dataIO.matrixWriter(cloud.getHeight(), "height" + filenameAddition);
        dataIO.matrixWriter(cloud.getDensity(), "density" + filenameAddition);
        dataIO.matrixWriter(dosage, "dose" + filenameAddition);
    }

    /**
     * Information that is printed to screen to help evaluate how the simulation
     * is running
     * @param volume The summed volume of the cloud
     * @param time The actual time in the simulation
     * @param counter The time step counter
     */
    protected void outputTerminalData(double volume, double time, int counter)
    {
        System.out.print("Volume = " + String.format("%.4f", volume)); 
        System.out.print("  Time = " + String.format("%.4f", time)); 
        System.out.println("  Counter = " + counter); 
    }

    /**
     * Determines the most suitable time step for the simulation at the present
     * time
     */
    protected void calcTimeStepSize()
    {

        //temp variables needed
        double umax = 0.1; //speed
        double cmax = 0.0;

        //Find the largest speed and height
        for (int i = 0; i < simRows; i++) {
            for (int j = 0; j < simCols; j++) {

                double u1 = (cloud.getU(i, j) * cloud.getU(i, j)
                        + cloud.getV(i, j) * cloud.getV(i, j));

                double c1 = (GRAVITY * (cloud.getDensity(i, j) - RHOA)
                        / cloud.getDensity(i, j) * cloud.getHeight(i, j));

                if (u1 > umax) {
                    umax = u1;
                }
                if (c1 > cmax) {
                    cmax = c1;
                }
            }
        }
        
        double courant = this.settingsFile.getDoubleSetting("courant");
        timeStep = courant * GRID_SIZE / (sqrt(umax) + sqrt(cmax));
    }

    /**
     * Simulates the behaviour of an areal source.
     * 
     * The source is a rectangle as defined by the matrix in the Source class
     * Contaminant comes out of the ground at speed 'up' and has a density 
     * of 'sourceRho'. This is thus a first approximation to an evaporating 
     * liquid pool.
     */
    private void updateArealSource()
    {
     
        double[][] up = source.getUpwardsVelocity();
        double sourceRho = source.getDensity(); //rho1 in TWODEE2013
        
        //copy of cloudNew matrices to simplify look of equations
        double[][]hNew = cloudNew.getHeight();
        double[][]uNew = cloudNew.getU();
        double[][]vNew = cloudNew.getV();
        double[][]rhoNew = cloudNew.getDensity();
        
        for(int i=0; i<simRows; i++)
        { 
            for(int j=0; j<simCols; j++)
            {
                if(up[i][j] > 0.0)
                {
                    double tempHeight = hNew[i][j]+up[i][j]*timeStep;
                
                    rhoNew[i][j]=(hNew[i][j]*rhoNew[i][j]+timeStep
                            *up[i][j]*sourceRho)/tempHeight;

                    uNew[i][j]=hNew[i][j]*rhoNew[i][j]*uNew[i][j] 
                            /(hNew[i][j]*rhoNew[i][j]+up[i][j]*timeStep*sourceRho);

                    vNew[i][j]=hNew[i][j]*rhoNew[i][j]*vNew[i][j]
                            /(hNew[i][j]*rhoNew[i][j]+up[i][j]*timeStep*sourceRho);

                    hNew[i][j] = tempHeight;
                }//end if
            }//end j loop
        }//end i loop

    }
    
    /**
     * This method searches for 'cliffs'. It suppresses motion into the cliff;
     * thus there will be zero flow over the cliff until the height of the fluid
     * exceeds that of the cliff. This may be due to the fluid simply piling up
     * at the base of the cliff, or ambient calcEntrainment; either mechanism 
     * may do.
     * 
     * @param c the cloud to be altered
     */
    private void cliffs(GasCloud c)
    {

        //copy of the ground height matrix to simplify look of equations
        double[][] gh = ground.getHeight(); 
        
        for(int i=1; i<simRows-1; i++)
        { 
            for(int j=1; j<simCols-1; j++)
            {
                double value; //temp storage
                
                //First, look for i-cliffs:
                if(c.getHeight(i, j) + gh[i][j] - gh[i+1][j] <= 0.0)   
                {
                    value = Math.min(c.getU(i, j), 0.0);
                    c.setU(i, j, value);       
                }

                if(c.getHeight(i, j) + gh[i][j] - gh[i-1][j] <= 0.0)   
                {  
                    value = Math.max(c.getU(i, j), 0.0);
                    c.setU(i, j, value);
                }
       
                //Now look for j-cliffs:
                if(c.getHeight(i, j) + gh[i][j] - gh[i][j+1] <= 0.0)   
                {
                    value = Math.min(c.getV(i, j), 0.0);
                    c.setV(i, j, value); 
                }
         
                if(c.getHeight(i, j) + gh[i][j] - gh[i][j-1] <= 0.0)   
                { 
                    value = Math.max(c.getV(i, j), 0.0);
                    c.setV(i, j, value); 
                }
            } //end of j loop
        } //end of i loop   
    }

    /**
     * Adds the effect of non-uniform height to momentum flux.
     * This is the hydrostatic difference in pressure on each side of 
     * the element.
     * 
     */
    private void calcHydrostaticDifference()
    {
        
        //copy of cloudNew matrices to simplify look of equations
        double[][]hNew = cloudNew.getHeight();
        double[][]uNew = cloudNew.getU();
        double[][]vNew = cloudNew.getV();
        double[][]rhoNew = cloudNew.getDensity();
        
        double[][]gh = ground.getHeight(); //copy of ground for same reason
        
        for (int i = BOUNDARY_SIZE; i < simRows - BOUNDARY_SIZE; i++) {
            for (int j = BOUNDARY_SIZE; j < simCols - BOUNDARY_SIZE; j++)
            {
                double uMomentum = 0.0; //momentum in the u direction (horizontal)
                double vMomentum = 0.0; //momentum in the v direction (vertical)

                double x = hNew[i][j]+gh[i][j]-hNew[i+1][j]-gh[i+1][j];
                double y = Math.min(hNew[i][j]+gh[i][j],hNew[i+1][j]+gh[i+1][j])
                          -Math.max(gh[i][j],gh[i+1][j]);
                double z = gh[i][j]-gh[i+1][j];
                    
                //loop until a case is met
                boolean isCaseFound = false;
                while(isCaseFound == false)
                {
                    //CASE A
                    if(x >= 0.0 && y >= 0.0 && z <= 0.0)
                    {
                        uMomentum = uMomentum
                                -(0.5*RHOA*x*x
                                +rhoNew[i+1][j]*x*y
                                +0.5*rhoNew[i+1][j]*y*y
                                +0.5*rhoNew[i][j]*z*z
                                -z*(RHOA*x+rhoNew[i+1][j]*y));
                        
                        isCaseFound = true;
                    }
                    //CASE B
                    else if (x <= 0.0 && y >= 0.0 && z <= 0.0)
                    {
                        uMomentum = uMomentum
                                -(rhoNew[i+1][j]*x*z
                                -(rhoNew[i+1][j]-RHOA)*x*y
                                +0.5*rhoNew[i+1][j]*y*y
                                +0.5*rhoNew[i][j]*z*z
                                -rhoNew[i][j]*y*z);
                                
                        isCaseFound = true;
                    }       
                    //CASE C
                    else if (x <= 0.0 && y <= 0.0 && z <= 0.0)
                    {
                        uMomentum = uMomentum-0.5*rhoNew[i][j]*hNew[i][j]*hNew[i][j];

                        isCaseFound = true;
                    } 
                    //CASE D
                    else if (x >= 0.0 && y <= 0.0 && z >= 0.0)
                    {
                        uMomentum = uMomentum-0.5*RHOA*hNew[i][j]*hNew[i][j];

                        isCaseFound = true;
                    }  
                    //CASE E
                    else if (x >= 0.0 && y >= 0.0 && z >= 0.0)
                    {
                        uMomentum = uMomentum
                                -(0.5*RHOA*x*x
                                +rhoNew[i+1][j]*x*y
                                +0.5*rhoNew[i+1][j]*y*y);

                        isCaseFound = true;
                    }  
                    //CASE F
                    else if (x <= 0.0 && y >= 0.0 && z >= 0.0)
                    {
                        uMomentum = uMomentum
                                -(0.0-(rhoNew[i+1][j]-RHOA)*x*y
                                +0.5*rhoNew[i+1][j]*y*y);

                        isCaseFound = true;
                    }  
                    else
                    {
                        /*
                        If the control gets to this point, it means that one 
                        of the forbidden combinations of signs has occurred. 
                        This is possibly due to the fact that y has to be of 
                        the correct sign, yet is the difference of two large 
                        numbers. The next line (repeated for the other 
                        three adjacent squares) set y to zero.
                        */
                        y = 0.0;
                    }//end of cases
                }//end of case loop
                
                x = hNew[i][j]+gh[i][j]-hNew[i-1][j]-gh[i-1][j];
                y = Math.min(hNew[i][j]+gh[i][j],hNew[i-1][j]+gh[i-1][j])
                          -Math.max(gh[i][j],gh[i-1][j]);
                z = gh[i][j]-gh[i-1][j];
                
                //loop until a case is met
                isCaseFound = false;
                while(isCaseFound == false)
                {
                    //CASE A
                    if(x >= 0.0 && y >= 0.0 && z <= 0.0)
                    {
                        uMomentum = uMomentum
                                +(0.5*RHOA*x*x
                                +rhoNew[i-1][j]*x*y
                                +0.5*rhoNew[i-1][j]*y*y
                                +0.5*rhoNew[i][j]*z*z
                                -z*(RHOA*x+rhoNew[i-1][j]*y));
                        
                        isCaseFound = true;
                    }
                    //CASE B
                    else if (x <= 0.0 && y >= 0.0 && z <= 0.0)
                    {
                        uMomentum = uMomentum
                                +(rhoNew[i-1][j]*x*z
                                -(rhoNew[i-1][j]-RHOA)*x*y
                                +0.5*rhoNew[i-1][j]*y*y
                                +0.5*rhoNew[i][j]*z*z
                                -rhoNew[i][j]*y*z);
                                
                        isCaseFound = true;
                    }       
                    //CASE C
                    else if (x <= 0.0 && y <= 0.0 && z <= 0.0)
                    {
                        uMomentum = uMomentum+0.5*rhoNew[i][j]*hNew[i][j]*hNew[i][j];

                        isCaseFound = true;
                    } 
                    //CASE D
                    else if (x >= 0.0 && y <= 0.0 && z >= 0.0)
                    {
                        uMomentum = uMomentum+0.5*RHOA*hNew[i][j]*hNew[i][j];

                        isCaseFound = true;
                    }  
                    //CASE E
                    else if (x >= 0.0 && y >= 0.0 && z >= 0.0)
                    {
                        uMomentum = uMomentum
                                +(0.5*RHOA*x*x
                                +rhoNew[i-1][j]*x*y
                                +0.5*rhoNew[i-1][j]*y*y);

                        isCaseFound = true;
                    }  
                    //CASE F
                    else if (x <= 0.0 && y >= 0.0 && z >= 0.0)
                    {
                        uMomentum = uMomentum
                                +(0.0-(rhoNew[i-1][j]-RHOA)*x*y
                                +0.5*rhoNew[i-1][j]*y*y);

                        isCaseFound = true;
                    }  
                    else
                    {
                        y = 0.0;
                    }//end of cases
                }//end of case loop
                
                x = hNew[i][j]+gh[i][j]-hNew[i][j+1]-gh[i][j+1];
                y = Math.min(hNew[i][j]+gh[i][j],hNew[i][j+1]+gh[i][j+1])
                          -Math.max(gh[i][j],gh[i][j+1]);
                z = gh[i][j]-gh[i][j+1];
                
                //loop until a case is met
                isCaseFound = false;
                while(isCaseFound == false)
                {
                    //CASE A
                    if(x >= 0.0 && y >= 0.0 && z <= 0.0)
                    {
                        vMomentum = vMomentum
                                -(0.5*RHOA*x*x
                                +rhoNew[i][j+1]*x*y
                                +0.5*rhoNew[i][j+1]*y*y
                                +0.5*rhoNew[i][j]*z*z
                                -z*(RHOA*x+rhoNew[i][j+1]*y));
                        
                        isCaseFound = true;
                    }
                    //CASE B
                    else if (x <= 0.0 && y >= 0.0 && z <= 0.0)
                    {
                        vMomentum = vMomentum
                                -(rhoNew[i][j+1]*x*z
                                -(rhoNew[i][j+1]-RHOA)*x*y
                                +0.5*rhoNew[i][j+1]*y*y
                                +0.5*rhoNew[i][j]*z*z
                                -rhoNew[i][j]*y*z);
                                
                        isCaseFound = true;
                    }       
                    //CASE C
                    else if (x <= 0.0 && y <= 0.0 && z <= 0.0)
                    {
                        vMomentum = vMomentum-0.5*rhoNew[i][j]*hNew[i][j]*hNew[i][j];

                        isCaseFound = true;
                    } 
                    //CASE D
                    else if (x >= 0.0 && y <= 0.0 && z >= 0.0)
                    {
                        vMomentum = vMomentum-0.5*RHOA*hNew[i][j]*hNew[i][j];

                        isCaseFound = true;
                    }  
                    //CASE E
                    else if (x >= 0.0 && y >= 0.0 && z >= 0.0)
                    {
                        vMomentum = vMomentum
                                -(0.5*RHOA*x*x
                                +rhoNew[i][j+1]*x*y
                                +0.5*rhoNew[i][j+1]*y*y);

                        isCaseFound = true;
                    }  
                    //CASE F
                    else if (x <= 0.0 && y >= 0.0 && z >= 0.0)
                    {
                        vMomentum = vMomentum
                                -(0.0-(rhoNew[i][j+1]-RHOA)*x*y
                                +0.5*rhoNew[i][j+1]*y*y);

                        isCaseFound = true;
                    }  
                    else
                    {
                        y = 0.0;
                    }//end of if-else cases
                }//end of loop
                
                x = hNew[i][j]+gh[i][j]-hNew[i][j-1]-gh[i][j-1];
                y = Math.min(hNew[i][j]+gh[i][j],hNew[i][j-1]+gh[i][j-1])
                          -Math.max(gh[i][j],gh[i][j-1]);
                z = gh[i][j]-gh[i][j-1];
                
                //loop until a case is met
                isCaseFound = false;
                while(isCaseFound == false)
                {
                    //CASE A
                    if(x >= 0.0 && y >= 0.0 && z <= 0.0)
                    {
                        vMomentum = vMomentum
                                +(0.5*RHOA*x*x
                                +rhoNew[i][j-1]*x*y
                                +0.5*rhoNew[i][j-1]*y*y
                                +0.5*rhoNew[i][j]*z*z
                                -z*(RHOA*x+rhoNew[i][j-1]*y));
                        
                        isCaseFound = true;
                    }
                    //CASE B
                    else if (x <= 0.0 && y >= 0.0 && z <= 0.0)
                    {
                        vMomentum = vMomentum
                                +(rhoNew[i][j-1]*x*z
                                +0.0-(rhoNew[i][j-1]-RHOA)*x*y
                                +0.5*rhoNew[i][j-1]*y*y
                                +0.5*rhoNew[i][j]*z*z
                                -rhoNew[i][j]*y*z);
                                
                        isCaseFound = true;
                    }       
                    //CASE C
                    else if (x <= 0.0 && y <= 0.0 && z <= 0.0)
                    {
                        vMomentum = vMomentum+0.5*rhoNew[i][j]*hNew[i][j]*hNew[i][j];

                        isCaseFound = true;
                    } 
                    //CASE D
                    else if (x >= 0.0 && y <= 0.0 && z >= 0.0)
                    {
                        vMomentum = vMomentum+0.5*RHOA*hNew[i][j]*hNew[i][j];

                        isCaseFound = true;
                    }  
                    //CASE E
                    else if (x >= 0.0 && y >= 0.0 && z >= 0.0)
                    {
                        vMomentum = vMomentum
                                +(0.5*RHOA*x*x
                                +rhoNew[i][j-1]*x*y
                                +0.5*rhoNew[i][j-1]*y*y);

                        isCaseFound = true;
                    }  
                    //CASE F
                    else if (x <= 0.0 && y >= 0.0 && z >= 0.0)
                    {
                        vMomentum = vMomentum
                                +(0.0-(rhoNew[i][j-1]-RHOA)*x*y
                                +0.5*rhoNew[i][j-1]*y*y);

                        isCaseFound = true;
                    }  
                    else
                    {
                        y = 0.0;
                    }//end of cases
                }//end of case loop
                
                uMomentum = uMomentum*GRAVITY*SHAPE_PARAM;
                vMomentum = vMomentum*GRAVITY*SHAPE_PARAM;
                
                if(hNew[i][j] > ROUGHNESS_LENGTH)
                {
                    uNew[i][j]=uNew[i][j]+uMomentum*timeStep
                            /(hNew[i][j]*rhoNew[i][j]*GRID_SIZE);
                    
                    vNew[i][j]=vNew[i][j]+vMomentum*timeStep
                            /(hNew[i][j]*rhoNew[i][j]*GRID_SIZE);
                } else
                {
                    uNew[i][j]=0.0;
                    vNew[i][j]=0.0;
                }
                
                //speed is the magnitude of velocity
                double speed = sqrt(uNew[i][j]*uNew[i][j]+vNew[i][j]*vNew[i][j]);
                
                double xdragu;
                double xdragv;
                
                if(hNew[i][j] > ROUGHNESS_LENGTH)
                {
                    xdragu=speed*uNew[i][j]*ground.getRoughness(i, j)/hNew[i][j];
                    xdragv=speed*vNew[i][j]*ground.getRoughness(i, j)/hNew[i][j];   
                } else
                {
                    xdragu = 0.0;
                    xdragv = 0.0;
                    uNew[i][j] = 0.0;
                    vNew[i][j] = 0.0;
                }
                
                if(uNew[i][j]*(uNew[i][j]-xdragu) > 0.0)
                    uNew[i][j] = uNew[i][j]-xdragu;
                else
                    uNew[i][j]=0.0;
                
                if(vNew[i][j]*(vNew[i][j]-xdragv) > 0.0)
                    vNew[i][j] = vNew[i][j]-xdragv;
                else
                    vNew[i][j]=0.0;
   
            }//end j loop
        }//end i loop
        
    }

    /**
     * Calculates calcEntrainment of air into the cloud
     * 
     * It calculates both the mass- and momentum
     * coupling between the dense layer and the leading edge. Also
     * the top entrainment. This routine will use the arrays fhu
     * and fhv to hold the force terms.
     */
    public void calcEntrainment()
    {
        double[][] uTop = new double[simRows][simCols]; //top calcEntrainment speed
        double uEdgeTotal = 0.0; //The total amount of calcEntrainment summed over the whole cloud
        double sumVolume = 0.0;  //The total volumne of the cloud
    
        double winduh = 0.0; //Horizontal wind velocity at reference height
        double windvh = 0.0; //Vertical wind velocity at reference height
        double windReferenceHeight = settingsFile.getDoubleSetting("windReferenceHeight");
        
        //copy of cloudNew matrices to simplify look of equations
        double[][]hNew = cloudNew.getHeight();
        double[][]uNew = cloudNew.getU();
        double[][]vNew = cloudNew.getV();
        double[][]rhoNew = cloudNew.getDensity();
        
        for (int i = BOUNDARY_SIZE; i < simRows - BOUNDARY_SIZE; i++) {
            for (int j = BOUNDARY_SIZE; j < simCols - BOUNDARY_SIZE; j++)
            {
                
                /*
                Note that a positive value for windu is to the East (from the 
                West) and a positive value for windv is to the South (from the 
                North). They are different because the file groundHeight.in 
                reads from line 1 which is the Northernmost part to line 200, 
                which is the Southernmost part. This is the opposite way round 
                from the model paradigm, which uses standard Cartesian notation, 
                with y running from down to up; with groundHeight.in being read 
                as it is, y runs from up to down (or North to South) so a 
                reflection is required to translate the two.
                */
                double windu = settingsFile.getDoubleSetting("wind_u");
                double windv = settingsFile.getDoubleSetting("wind_v");
                
                int ie;
                int je;
                
                if(hNew[i][j] > ROUGHNESS_LENGTH)
                {
                    winduh = windu*(Math.log(hNew[i][j]/ROUGHNESS_LENGTH))
                            /(Math.log(windReferenceHeight/ROUGHNESS_LENGTH));
                    
                    windvh = windv*(Math.log(hNew[i][j]/ROUGHNESS_LENGTH))
                            /(Math.log(windReferenceHeight/ROUGHNESS_LENGTH));
                }
                else
                {
                    winduh = 0.0;
                    windvh = 0.0;
                }
                
                if(hNew[i+1][j] > hNew[i-1][j])
                    ie = i-1;
                else
                    ie = i+1;
                
                if(hNew[i][j+1] > hNew[i][j-1])
                    je = j-1;
                else
                    je = j+1;
                
                double termu=(winduh*(i-ie)*((hNew[ie][j]*(uNew[ie][j]-winduh)) 
                        -(hNew[i][j]*(uNew[i][j]-winduh)))
                        +windvh*(j-je)*((hNew[i][je]*(uNew[i][je]-winduh))  
                        -(hNew[i][j]*(uNew[i][j]-winduh))));
                
                double termv=(winduh*(i-ie)*((hNew[ie][j]*(vNew[ie][j]-windvh))
                        -(hNew[i][j]*(vNew[i][j]-windvh)))
                        +windvh*(j-je)*((hNew[i][je]*(vNew[i][je]-windvh))  
                        -(hNew[i][j]*(vNew[i][j]-windvh))));
                
                double dhudt=(hNew[i][j]*(uNew[i][j]-winduh)-cloud.getHeight(i, j)
                        *(cloud.getU(i, j)-winduh));
                
                double dhvdt=(hNew[i][j]*(vNew[i][j]-windvh)-cloud.getHeight(i, j)
                        *(cloud.getV(i, j)-windvh));
                
                /*
                Note that the two lines below use fhu and fhv as the forces
                on the element (i,j). This is so that the force on a particular 
                fluid element is not dependent upon adjacent elements, some 
                of which have been affected by this routine and some of which 
                have not been. This clearly causes an asymmetry: if we consider 
                element (i,j), then element (i,j-1) has been affected by the 
                force applied but element (i,j+1) is virginal.
                */
                fhu[i][j]=dhudt-termu*timeStep/GRID_SIZE;
                fhv[i][j]=dhvdt-termv*timeStep/GRID_SIZE;
                
                
                double dhdt=(hNew[i][j]-cloud.getHeight(i, j))
                        -(winduh*(i-ie))*(hNew[ie][j]-hNew[i][j])*timeStep/GRID_SIZE
                        -(windvh*(j-je))*(hNew[i][je]-hNew[i][j])*timeStep/GRID_SIZE;  
                
                if(dhdt < 0.0)
                {
                    fhu[i][j] = 0.0;
                    fhv[i][j] = 0.0;
                }
                
                /*
                The following four kappas are kappa for the squares "up",
                "down", "right" and "left" of point (i,j).  Here, up means j-wise
                and left/right means i-wise. 
                */
                
                if(hNew[i][j] > ROUGHNESS_LENGTH)
                {
                    double kappaUp = hNew[i][j+1]
                            *sqrt(uNew[i][j+1]*uNew[i][j+1]
                                    +vNew[i][j+1]*vNew[i][j+1]);

                    double kappaDown = hNew[i][j-1]
                            *sqrt(uNew[i][j-1]*uNew[i][j-1]
                                    +vNew[i][j-1]*vNew[i][j-1]);

                    double kappaLeft = hNew[i-1][j]
                            *sqrt(uNew[i-1][j]*uNew[i-1][j]
                                    +vNew[i-1][j]*vNew[i-1][j]);

                    double kappaRight = hNew[i+1][j]
                            *sqrt(uNew[i+1][j]*uNew[i+1][j]
                                    +vNew[i+1][j]*vNew[i+1][j]);
                    
                    /*
                    The x- and y- components of the viscous force are here called
                    visX and visY. Note that, because of the dyadic in the term,
                    that visx is a function of only `u' and visy is a function only
                    of `v'. 
                    */
                    double zeta = settingsFile.getDoubleSetting("zeta");
                    
                    double visX = zeta
                            * 4.0*((kappaRight*(uNew[i+2][j]-uNew[i][j]))
                            - (kappaLeft*(uNew[i][j]-uNew[i-2][j]))
                            + (kappaUp*(uNew[i][j+2]-uNew[i][j]))
                            - (kappaDown*(uNew[i][j]-uNew[i][j-2])));
                    
                    double visY = zeta
                            * 4.0*((kappaRight*(vNew[i+2][j]-vNew[i][j]))
                            - (kappaLeft*(vNew[i][j]-vNew[i-2][j]))
                            + (kappaUp*(vNew[i][j+2]-vNew[i][j]))
                            - (kappaDown*(vNew[i][j]-vNew[i][j-2])));
                    
                    flu[i][j] = visX * timeStep;
                    flv[i][j] = visY * timeStep;   
                
                    /*
                    In the above visx and visy terms, note that I have used flu(i,j)
                    and flv(i,j) to hold the force terms. I do not use the force
                    terms straight away as to do so would result in asymmetry (see the
                    comment lines above discussing the similar use of fhu and fhv).
                     
                    Note that it is not possible to incorporate flu and fhu (as used
                    in this routine) together because fhu is the force due  to edge
                    interaction and must be multiplied by some function of the front
                    Froude number.
                    */
                    
                }//end if
                
                /*
                Now for the top calcEntrainment terms. The edge stuff will be dealt
                with at the end of this routine, separately. This is because
                the edge entrained stuff is evenly distributed throughout the 
                cloud, and not concentrated at the cloud's leading edge. Top
                calcEntrainment is local, though.
                */
                if(hNew[i][j] > ROUGHNESS_LENGTH && rhoNew[i][j] > RHOA)
                {
                    double wStar = settingsFile.getDoubleSetting("wStar"); //convection velocity
                    double alpha2 = settingsFile.getDoubleSetting("alpha2");
                    double alpha3 = settingsFile.getDoubleSetting("alpha3");
                    double alpha4 = settingsFile.getDoubleSetting("alpha4");
                    double alpha6 = settingsFile.getDoubleSetting("alpha6");
                    double alpha7 = settingsFile.getDoubleSetting("alpha7");
                    
                    double gDashStar = 
                            GRAVITY * Math.max((rhoNew[i][j]-RHOA)
                                    /RHOA, 1.0);
                    
                    //vee is the vertical entrainment velocity
                    double vee = sqrt(alpha2*wStar*wStar+U_STAR*U_STAR+ground.getRoughness(i, j)
                            /2*(alpha3*alpha3*(uNew[i][j]*uNew[i][j]
                            +vNew[i][j]*vNew[i][j]))
                            +alpha7*alpha7*((uNew[i][j]-winduh)*(uNew[i][j]-winduh)
                            +(vNew[i][j]-windvh)*(vNew[i][j]-windvh)));
                    
                    if(vee > 0.0)
                    {
                        double richard = gDashStar*hNew[i][j]/vee/vee;
                        uTop[i][j]=alpha4/(alpha4/alpha6+richard)*vee;   
                    } else {
                        uTop[i][j] = 0.0;
                    }
    
                } else {
                    
                    uNew[i][j] = 0.0;
                    vNew[i][j] = 0.0;
                    uTop[i][j] = 0.0;
                }//end if
                
                //The local edge calcEntrainment speed
                double edgeEntrainmentCoeff = settingsFile.getDoubleSetting("edgeEntrainmentCoeff");
                
                /*
                Note that uEdge is actually a velocity here (distance over time).
                Also rememember that the integral of uedge over the gravity 
                current front is equal to the volume flux into it (or, to be 
                more specific, the areal flux into it, times the cross section)
                
                Just to make it plain: uedge is the local edge entrainment speed
                and uEdgeTotal is the total amount of entrainment summed over 
                the whole cloud.
                */
                double uEdge = ((hNew[i][j]-cloud.getHeight(i, j))
                        +(winduh*(i-ie)/GRID_SIZE)*(hNew[ie][j]-hNew[i][j])*timeStep 
                        +(windvh*(j-je)/GRID_SIZE)*(hNew[i][je]-hNew[i][j])*timeStep)
                        *edgeEntrainmentCoeff;
                
                uEdge = Math.max(uEdge, 0.0)/timeStep;
                if(uEdge < 0.0)
                {
                    //stop program
                    System.out.println("");
                    System.out.println("uEdge is greater than 0.0");
                    System.out.println("Simulation terminated");
                    System.exit(0);
                    
                }
                
                uEdgeTotal = uEdgeTotal+uEdge*timeStep;
                sumVolume = sumVolume + hNew[i][j];
                
            }//end j loop
        }//end i loop
        
        /*
        Note that this part has to be done here as uEdgeTotal has to be totalled
        over all i and j. Also, sum is the total volume of the cloud and 
        will be the denominator.
        */
        if(sumVolume > 0.0)
        {
            for (int i = BOUNDARY_SIZE; i < simRows - BOUNDARY_SIZE; i++) {
                for (int j = BOUNDARY_SIZE; j < simCols - BOUNDARY_SIZE; j++)
                {
                    if(hNew[i][j] > ROUGHNESS_LENGTH)
                    {
                        rhoNew[i][j]=(rhoNew[i][j]*hNew[i][j]+RHOA
                                *hNew[i][j]*uEdgeTotal/sumVolume)
                                /(hNew[i][j]+hNew[i][j]*uEdgeTotal/sumVolume);
                        
                        hNew[i][j] = hNew[i][j]*(1.0+uEdgeTotal/sumVolume);
                    }//end if
                }//end j loop
            }//end i loop
        }//end if
        
        for (int i = BOUNDARY_SIZE; i < simRows - BOUNDARY_SIZE; i++) {
            for (int j = BOUNDARY_SIZE; j < simCols - BOUNDARY_SIZE; j++)
            {
                if(hNew[i][j] > ROUGHNESS_LENGTH)
                {
                    //beta 0100 is a Froude number of unity.
                    double beta = settingsFile.getDoubleSetting("beta");
                            
                    double k = 1/(2.0*SHAPE_PARAM*beta*beta);
                   
                    uNew[i][j]=uNew[i][j]-k*RHOA*(fhu[i][j])
                            /(rhoNew[i][j]*hNew[i][j]);
                    
                    vNew[i][j]=vNew[i][j]-k*RHOA*(fhv[i][j])
                            /(rhoNew[i][j]*hNew[i][j]);
                    
                    //Now for the viscosity terms
                    uNew[i][j]=uNew[i][j]-flu[i][j];
                    vNew[i][j]=vNew[i][j]-flv[i][j];
                    /*
                    remember that the flu and flv are already in the correct 
                    units and don't need to be multiplied by dt 
                    */
                    
                    uNew[i][j]=(uNew[i][j]*hNew[i][j]*rhoNew[i][j]
                            +winduh*timeStep*uTop[i][j]*RHOA)
                            /(hNew[i][j]*rhoNew[i][j]+timeStep*uTop[i][j]*RHOA);
                    
                    vNew[i][j]=(vNew[i][j]*hNew[i][j]*rhoNew[i][j]
                            +windvh*timeStep*uTop[i][j]*RHOA)
                            /(hNew[i][j]*rhoNew[i][j]+timeStep*uTop[i][j]*RHOA);
                    
                    rhoNew[i][j]=(rhoNew[i][j]*hNew[i][j]+RHOA*timeStep*uTop[i][j])
                            /(hNew[i][j]+timeStep*uTop[i][j]);
                    
                    hNew[i][j]=hNew[i][j]+uTop[i][j]*timeStep;
                    
                }//end if
            }//end j loop
        }//end i loop
         
    }//end of Entrainment method
}
