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
 * Class to represent the areal source
 * 
 * @author Stephen Denekamp
 * @version 1.0: Created based on TWODEE2013 code by Dr. Robin Hankin
 */
public class ArealSource {
    
    //The matrix representing the upwards velocity of the areal source
    private final double[][] UPWARDS_VELOCITY;
    
    //The DENSITY of the substance the real source is made up of
    private final double DENSITY; 
    
    /**
     * Constructor for ArealSource
     * 
     * @param upwardsVelocity The initial matrix of the upwards velocity
     * @param density The density of the stuff the areal source is made up of
     */
    public ArealSource(double[][] upwardsVelocity, double density) {
        this.UPWARDS_VELOCITY = upwardsVelocity;
        this.DENSITY = density;
    }

    /**
     * Returns a reference to the UPWARDS_VELOCITY matrix of the Areal Source
     * 
     * @return the matrix by reference
     */
    public double[][] getUpwardsVelocity() {
        return UPWARDS_VELOCITY;
    }
    
    /**
     * Returns the UPWARDS_VELOCITY of a particular Areal Source point
     * 
     * @param row The row coordinate (i)
     * @param col The column coordinate (j)
     * @return the UPWARDS_VELOCITY value
     */
    public double getUpwardsVelocity(int row, int col)
    {
        return UPWARDS_VELOCITY[row][col];
    }
    
    /**
     * 
     * @return the DENSITY variable of the areal source
     */
    public double getDensity()
    {
        return DENSITY;
    }
}
