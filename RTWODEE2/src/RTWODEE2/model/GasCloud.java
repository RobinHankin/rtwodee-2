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

/**
 * This class represents a cloud with multiple inputs (HEIGHT and DENSITY)
 and is used for a full R-TWODEE2 simulation
 * 
 * @author Stephen Denekamp
 * @version 1.0: Created based on TWODEE2013 code by Dr. Robin Hankin
 */
public class GasCloud extends Cloud
{

    private final double[][] HEIGHT;
    private final double[][] DENSITY;
  
    private final double[][] MAX_DENSITY_AT_GROUND;     //rhomax in TWODEE2013
    private final double[][] MAX_DENSITY_AT_ALTITUDE;   //rmaxh in TWODEE2013
    
    /**
     * Constructor for a GasCloud Cloud
     * 
     * @param hVelocity The matrix of the Clouds horizontal velocity
     * @param vVelocity The matrix of the Clouds vertical velocity
     * @param height The matrix of the Clouds HEIGHT
     * @param density  The matrix of the Clouds DENSITY
     * @param gridSize The size/scale of each point of the cloud
     * @param rhoa The ambient air DENSITY
     */
    public GasCloud(double[][] hVelocity, double[][] vVelocity,
            double[][] height, double[][] density, double rhoa)
    {
        super(hVelocity, vVelocity);
        this.HEIGHT = height;
        this.DENSITY = density;

        int cols = hVelocity.length;
        int rows = hVelocity[0].length;
        MAX_DENSITY_AT_GROUND = new double[rows][cols];
        MAX_DENSITY_AT_ALTITUDE = new double[rows][cols];
        setDensityDefaults(rhoa);        
    }
    
    /**
     * Constructor to create an empty GasCloud
     * 
     * @param gridSize The size/scale of each square in the cloud
     * @param rows Number of rows long the cloud is
     * @param cols Number of columns long the cloud is
     * @param rhoa The ambient air DENSITY
     */
    public GasCloud(int rows, int cols, double rhoa)
    {
        super(rows, cols);
        HEIGHT = new double[rows][cols];
        DENSITY = new double[rows][cols]; 
        
        MAX_DENSITY_AT_GROUND = new double[rows][cols];    //maximum DENSITY over time
        MAX_DENSITY_AT_ALTITUDE = new double[rows][cols];  //maximum DENSITY over time
        setDensityDefaults(rhoa);
    }
    
    /**
     * Sets the HEIGHT matrix of the cloud
     * 
     * @param heightMatrix matrix of new values
     */
    public void setHeight(double[][] heightMatrix)
    {     
        for (int i = 1; i < heightMatrix.length; i++) 
           for (int j = 1; j < heightMatrix[0].length; j++) 
               this.HEIGHT[i][j] = heightMatrix[i][j];
    }
    
    /**
     * Sets the HEIGHT of the cloud at a particular point
     * 
     * @param row The row coordinate (i)
     * @param col The column coordinate (j)
     * @param newValue The new newValue
     */
    public void setHeight(int row, int col, double newValue)
    {
        HEIGHT[row][col] = newValue;
    }
    
    /**
     * Sets the DENSITY matrix of the cloud
     * 
     * @param densityMatrix matrix of new values
     */
    public void setDensity(double[][] densityMatrix)
    {
        for (int i = 1; i < densityMatrix.length; i++) 
           for (int j = 1; j < densityMatrix[0].length; j++) 
               this.DENSITY[i][j] = densityMatrix[i][j];
    }
    
    /**
     * Sets the DENSITY of the cloud at a particular point
     * 
     * @param row The row coordinate (i)
     * @param col The column coordinate (j)
     * @param newValue The new newValue
     */
    public void setDensity(int row, int col, double newValue)
    {
        DENSITY[row][col] = newValue;
    }

    /**
     * Returns a reference to the HEIGHT matrix of the cloud
     * 
     * @return the matrix by reference
     */
    public double[][] getHeight()
    {
        return HEIGHT;
    }
    
    /**
     * Returns the HEIGHT of a particular point in the cloud
     * 
     * @param row The row coordinate (i)
     * @param col The column coordinate (j)
     * @return the HEIGHT value
     */
    public double getHeight(int row, int col)
    {
        return HEIGHT[row][col];
    }

    /**
     * Returns a reference to the DENSITY matrix of the cloud
     * 
     * @return the matrix by reference
     */
    public double[][] getDensity()
    {
        return DENSITY;
    }
    
    /**
     * Returns the DENSITY of a particular point in the cloud
     * 
     * @param row The row coordinate (i)
     * @param col The column coordinate (j)
     * @return the DENSITY value
     */
    public double getDensity(int row, int col)
    {
        return DENSITY[row][col];
    }
    
    /**
     * Sets the MAX_DENSITY_AT_GROUND of the cloud at a particular point
     * 
     * @param row The row coordinate (i)
     * @param col The column coordinate (j)
     * @param newValue The new newValue
     */
    public void setMaxDensityAtGround(int row, int col, double newValue)
    {
        MAX_DENSITY_AT_GROUND[row][col] = newValue;
    }
    
    /**
     * Sets the MAX_DENSITY_AT_GROUND of the cloud at a particular point
     * 
     * @param row The row coordinate (i)
     * @param col The column coordinate (j)
     * @param newValue The new newValue
     */
    public void setMaxDensityAtAltitude(int row, int col, double newValue)
    {
        MAX_DENSITY_AT_ALTITUDE[row][col] = newValue;
    }
    
    /**
     * Returns a reference to the MAX_DENSITY_AT_GROUND matrix of the cloud
     * 
     * @return the matrix by reference
     */
    public double[][] getMaxDensityAtGround()
    {
        return this.MAX_DENSITY_AT_GROUND;
    }
    
    /**
     * Returns the MAX_DENSITY_AT_GROUND of a particular point in the cloud
     * 
     * @param row The row coordinate (i)
     * @param col The column coordinate (j)
     * @return the GroundMaxDosage value
     */  
    public double getMaxDensityAtGround(int row, int col)
    {
        return MAX_DENSITY_AT_GROUND[row][col];
    }
    
    /**
     * Returns a reference to the MAX_DENSITY_AT_ALTITUDE matrix of the cloud
     * 
     * @return the matrix by reference
     */
    public double[][] getMaxDensityAtAltitude()
    {
        return this.MAX_DENSITY_AT_ALTITUDE;
    }
    
    /**
     * Returns the MAX_DENSITY_AT_ALTITUDE of a particular point in the cloud
     * 
     * @param row The row coordinate (i)
     * @param col The column coordinate (j)
     * @return the MAX_DENSITY_AT_ALTITUDE value
     */  
    public double getMaxDensityAtAltitude(int row, int col)
    {
        return MAX_DENSITY_AT_ALTITUDE[row][col];
    }

    /**
     * Calculates the following cloud metrics and returns them as an array
     * Code based on TWODEE2013 code by Dr. Robin Hankin
     * - sumVolume: the total volume of the cloud
     * - sumMomentum: the total momentum of the cloud 
     * - xCentre: the x coordinate of the cloud centre (weighted by the volume)
     * - yCentre: the y coordinate of the cloud centre (weighted by the volume)
     * - area: the area of the cloud
     * 
     * @param roughnessLength
     * @param gridSize
     * @return an array of the above metrics
     */
    public double[] calcMetrics(double roughnessLength, double gridSize)
    {
        double sumVolume = 0.0;   
        double sumMomentum = 0.0; 
        double xCentre = 0.0;      
        double yCentre = 0.0;    
        double area = 0.0;           
        
        //Caculate values for every point of the point
        for(int i=0; i<HEIGHT.length; i++)
        {
            for(int j=0; j<HEIGHT[0].length; j++)
            {
                sumVolume = sumVolume + HEIGHT[i][j];

                sumMomentum = sumMomentum + HEIGHT[i][j] * DENSITY[i][j] 
                        * super.getU()[i][j];

                xCentre = xCentre + HEIGHT[i][j] * (i+1);
                
                yCentre = yCentre + HEIGHT[i][j] * (j+1);
                
                if (HEIGHT[i][j] > roughnessLength)
                    area = area + 1.0;
                
            }//end j loop
        }//end i loop
        
        if(sumVolume > 0.001)
        {
            xCentre = xCentre / sumVolume * gridSize;
            yCentre = yCentre / sumVolume * gridSize;
        }
        
        area = area * gridSize * gridSize;
        
        //put all items into an array to be returned
        double[] metrics = new double[5];
        metrics[0] = sumVolume;
        metrics[1] = sumMomentum;
        metrics[2] = xCentre;
        metrics[3] = yCentre;
        metrics[4] = area;

        return metrics;
    }
    
    /**
     * sets the DENSITY matrices to the ambient air DENSITY (rhoa) everywhere.
     * This needs special treatment as its default value is rhoa, not zero.
     * 
     * @param amount The amount of ambient air DENSITY
     */
    private void setDensityDefaults(double amount)
    {
        for (int i = 0; i < DENSITY.length; i++) {
            for (int j = 0; j < DENSITY[0].length; j++) 
            {
                DENSITY[i][j] = amount;
                MAX_DENSITY_AT_GROUND[i][j] = amount;
                MAX_DENSITY_AT_ALTITUDE[i][j] = amount;   
            }
        }
    }
    
    

}
