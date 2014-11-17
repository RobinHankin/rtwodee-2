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

/**
 * Subclass of a Cloud
 * This class represents an with a single variable input and is used for 
 * an advection only simulation
 * 
 * @author Stephen Denekamp
 * @version 1.0: Created based on TWODEE2013 code by Dr. Robin Hankin
 */
public class ScalarArray extends Cloud 
{

    private final double[][] INPUT_VARIABLE;
    
    /**
     * Constructor for a ScalarArray Cloud
     * 
     * @param hVelocity The matrix of the Clouds horizontal velocity
     * @param vVelocity The matrix of the Clouds vertical velocity
     * @param simInput  The matrix that that will be monitored/changed 
     * by the simulation (e.g height)
     */
    public ScalarArray(double[][] hVelocity, double[][] vVelocity,
                            double[][] simInput) 
    {
        super(hVelocity, vVelocity);
        this.INPUT_VARIABLE = simInput;

    }
    
    /**
     * Sets the INPUT_VARIABLE matrix of the cloud
     * 
     * @param inputMatrix Matrix of new values
     */
    public void setInputVariable(double[][] inputMatrix)
    {
        for (int i = 1; i < inputMatrix.length; i++) 
           for (int j = 1; j < inputMatrix[0].length; j++) 
               this.INPUT_VARIABLE[i][j] = inputMatrix[i][j];   
    }
    
    /**
     * Sets the INPUT_VARIABLE of the cloud at a particular point
     * 
     * @param row The row coordinate (i)
     * @param col The column coordinate (j)
     * @param newValue The new newValue
     */
    public void setInputVariable(int row, int col, double newValue)
    {
        INPUT_VARIABLE[row][col] = newValue; 
    }
    
    /**
     * Returns a reference to the INPUT_VARIABLE matrix of the cloud
     * 
     * @return the matrix by reference
     */
    public double[][] getInputVariable()
    {
        return INPUT_VARIABLE;
    }
    
    /**
     * Returns the INPUT_VARIABLE of a particular point in the cloud
 This is faster than copying the entire matrix
     * 
     * @param row The row coordinate (i)
     * @param col The column coordinate (j)
     * @return the INPUT_VARIABLE value
     */
    public double getInputVariable(int row, int col)
    {
        return INPUT_VARIABLE[row][col];
    }
}
