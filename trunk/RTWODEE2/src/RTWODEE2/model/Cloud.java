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
 * This super class represents the cloud
 * 
 * @author Stephen Denekamp
 * @version 1.0: Created based on TWODEE2013 code by Dr. Robin Hankin
 * @version 1.1: Settings and Getters renamed to u & v at Robin's request
 */
public abstract class Cloud
{

    private final double[][] U_HORIZONTAL_VELOCITY;
    private final double[][] V_VERTICAL_VELOCITY;
    
    /**
     * Constructor for a cloud
     *
     * @param hVelocity A matrix representing the initial horizontal velocity
     * @param vVelocity A matrix representing the initial vertical velocity
     */
    public Cloud(double[][] hVelocity, double[][] vVelocity)
    {
        this.U_HORIZONTAL_VELOCITY = hVelocity;
        this.V_VERTICAL_VELOCITY = vVelocity;     
    }
    
    /**
     * Constructor to create an empty cloud
     * 
     * @param rows Number of rows long the cloud is
     * @param cols Number of columns long the cloud is
     */
    public Cloud(int rows, int cols)
    {
        U_HORIZONTAL_VELOCITY = new double[rows][cols];
        V_VERTICAL_VELOCITY = new double[rows][cols];  
    }

    /**
     * Sets the horiztonalVelocity matrix of the cloud
     * 
     * @param horizontalVelocityMatrix Matrix of new values
     */
    public void setU(double[][] horizontalVelocityMatrix)
    {
        for (int i = 1; i < horizontalVelocityMatrix.length; i++) 
           for (int j = 1; j < horizontalVelocityMatrix[0].length; j++) 
               this.U_HORIZONTAL_VELOCITY[i][j] = horizontalVelocityMatrix[i][j];
    }

    /**
     * Sets the U_HORIZONTAL_VELOCITY of the cloud at a particular point
     * 
     * @param row The row coordinate (i)
     * @param col The column coordinate (j)
     * @param newValue The new value
     */
    public void setU(int row, int col, double newValue)
    {
        U_HORIZONTAL_VELOCITY[row][col] = newValue;
        
    }
    
    /**
     * Sets the V_VERTICAL_VELOCITY matrix of the cloud
     * 
     * @param verticalVelocityMatrix Matrix of new values
     */
    public void setV(double[][] verticalVelocityMatrix)
    {
        for (int i = 1; i < verticalVelocityMatrix.length; i++) 
           for (int j = 1; j < verticalVelocityMatrix[0].length; j++) 
               this.V_VERTICAL_VELOCITY[i][j] = verticalVelocityMatrix[i][j];
    }

    /**
     * Sets the V_VERTICAL_VELOCITY of the cloud at a particular point
     * 
     * @param row The row coordinate
     * @param col The column coordinate
     * @param newValue The new value
     */
    public void setV(int row, int col, double newValue)
    {
        V_VERTICAL_VELOCITY[row][col] = newValue;
    }
    
    /**
     * Returns a reference to the horiztonalVelocity matrix of the cloud
     * 
     * @return the matrix by reference
     */
    public double[][] getU()
    {
        return U_HORIZONTAL_VELOCITY;
    }
    
    /**
     * Returns the U_HORIZONTAL_VELOCITY of a particular point in the cloud
     * 
     * @param row The row coordinate (i)
     * @param col The column coordinate (j)
     * @return the horizontalVelcoity value
     */
    public double getU(int row, int col)
    {
        return U_HORIZONTAL_VELOCITY[row][col];
    }
    
    /**
     * Returns a reference to the horiztonalVelocity matrix of the cloud
     * 
     * @return the matrix by reference
     */
    public double[][] getV()
    {
        return V_VERTICAL_VELOCITY;
    }
  
    
    /**
     * Returns the V_VERTICAL_VELOCITY of a particular point in the cloud
     * 
     * @param row The row coordinate (i)
     * @param col The column coordinate (j)
     * @return the verticalVelcoity value
     */  
    public double getV(int row, int col)
    {
        return V_VERTICAL_VELOCITY[row][col];
    }   
}
