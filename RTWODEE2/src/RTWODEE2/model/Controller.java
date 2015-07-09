/*
 * Copyright (C) 2014  Robin Hankin
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
 * This class serves as the controlling class for R-TWODEE2 that initialises
 * the correct simulation
 * 
 * @author Stephen Denekamp
 * @version 1.0: Created
 */
public class Controller 
{    
    private Simulation simulation;
    private final ExternalSettingsReader settings;

    /**
     * Constructor for the controller class
     * 
     */
    public Controller() 
    {  
        settings = new ExternalSettingsReader();
    }

  
    /**
     * Starts the correct simulation 
     * 
     */
    public void start() 
    {    
        int simulationType = settings.getIntegerSetting("simulationType");
        
        if(simulationType == 0)
            simulation = new AdvectionSimulation(settings);
        else
            simulation = new Simulation(settings); //full simulation
        
        simulation.runSimulation();
    }
    
}
