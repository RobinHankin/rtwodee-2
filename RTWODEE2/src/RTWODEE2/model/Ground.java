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
 * Class to represent the ground
 * 
 * @author Stephen Denekamp
 * @version 1.0: Created based on TWODEE2013 code by Dr. Robin Hankin
 */
public class Ground {
    
    private final double[][] HEIGHT;    //height/altitude of ground (topography)
    private final double[][] ROUGHNESS;
   
    public Ground(double[][] height, double[][] roughness) {
        this.HEIGHT = height;
        this.ROUGHNESS = roughness;
    }
    
    /**
     * Sets the HEIGHT of the ground at a particular point
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
     * Returns a reference to the HEIGHT matrix of the Areal Source
     * 
     * @return the matrix by reference
     */
    public double[][] getHeight() {
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
     * Returns a reference to the ROUGHNESS matrix of the Areal Source
     * 
     * @return the matrix by reference
     */
    public double[][] getRoughness() {
        return ROUGHNESS;
    }
    
    /**
     * Returns the ROUGHNESS of a particular point in the cloud
     * 
     * @param row The row coordinate (i)
     * @param col The column coordinate (j)
     * @return the ROUGHNESS value
     */
    public double getRoughness(int row, int col)
    {
        return ROUGHNESS[row][col];
    }
    
}
