/*
 * Copyright (C) 2014 Dr. Robin Hankin
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
 * Subclass of Simulation. 
 * This class is currently a holder class to later be developed to run an
 * advection only simulation. Before this can be done the calculations for
 * advection will need to be separated into their own methods in the 
 * Simulation class.
 * 
 * @author Stephen Denekamp
 * @version 1.0: Created based on TWODEE2013 code by Dr. Robin Hankin
 */
public class AdvectionSimulation extends Simulation
{
    //Simulation Objects
    private ScalarArray cloud;
    private final String inputFile;
    
    //Other Simulation Variables
    private double time = 0; //real world time in simulation (in seconds)
    private int simRows; //the size of the simulation plane
    private int simCols; //the size of the simulation plane
    
    /**
     * Constructor for an Advection Simulation
     * @param settingsFile the 
     */
    public AdvectionSimulation(ExternalSettingsReader settingsFile)
    {
        super(settingsFile);
        inputFile = settingsFile.getStringSetting("advectionInputFile");
    }
    
    @Override
    public void runSimulation()
    {
        initialSimulationSetup();
        //advection simulation code to go here
        
    }
    
    @Override
    protected void initialSimulationSetup()
    {
        
        double[][] hVelocity = dataIO.matrixReader("horizontalVelocity.in");
        double[][] vVelocity = dataIO.matrixReader("verticalVelocity.in");
        double[][] advectionInput = dataIO.matrixReader("advectionInput.in");
        simRows = advectionInput.length;
        simCols = advectionInput[0].length;
        
        //create cloud
        cloud = new ScalarArray(hVelocity, vVelocity, advectionInput);
    }
    
    @Override
    protected void calcTimeStepSize()
    {
        //Cloud matrices
        double[][] v = cloud.getV();
        double[][] u = cloud.getU();
            
        //temp variables needed
        double umax=0.1;
        double u1;
    
        //Find the largest speed and height
        for(int i=0;i<simRows;i++)
        {
            for(int j=0;j<simCols;j++)
            {
                u1 = (u[i][j]*u[i][j]+v[i][j]*v[i][j]);
                if(u1 > umax) umax=u1;
            }
        }  
        
        double courant = this.settingsFile.getDoubleSetting("courant");
        timeStep = courant*GRID_SIZE/sqrt(umax);
    }
    
    @Override
    protected void outputAllMatrices(String filenameAddition)
    {
        dataIO.matrixWriter(cloud.getU(),"horizontalVelocity"+filenameAddition);
        dataIO.matrixWriter(cloud.getV(),"verticalVelocity"+filenameAddition);
        dataIO.matrixWriter(cloud.getInputVariable(), "advectionMatrix"+filenameAddition);  
    }
    
}
