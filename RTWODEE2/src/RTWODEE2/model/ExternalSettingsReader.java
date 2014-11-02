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

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Properties;

/**
 * This class handles the reading of the external settings file that contains
 * all of the variables to instruct R-TWODEE2 how to run
 * 
 * @author Stephen Denekamp
 * @version 1.0: Created
 */
public class ExternalSettingsReader
{
    private final String inputLocation = "dataFiles/";
    private Properties prop = new Properties();
    
    /**
     * Constructor for External Settings Reader.
     * Creates a new Properties object and reads in the settings file
     */
    public ExternalSettingsReader()
    {
        try 
        {
            ClassLoader classLoader = Thread.currentThread().getContextClassLoader();
            prop.load(classLoader.getResourceAsStream("RTWODEE2/dataFiles/genSettings.in"));

        } catch (FileNotFoundException | NullPointerException ex) {
            //Logger.getLogger(ExternalSettingsReader.class.getName()).log(Level.SEVERE, null, ex);
            stopProgram("genSettings.in file not found", "");
        } catch (IOException ex) {
            stopProgram("Problem with genSettings.in file", "");
        }
    }
    
    /**
     * Use to read a String variable from the settings file
     * 
     * @param settingToGet The variable required
     * @return the variable as a String
     */
    public String getStringSetting(String settingToGet)
    { 
        String setting = "";
        try
        {
        setting = prop.getProperty(settingToGet);
        }
        catch (NullPointerException ex)
        {
            settingNotFound(settingToGet);
        }
           
        return setting;
    }
    
    /**
     * Use to read a numerical double variable from the settings file
     * 
     * @param settingToGet The variable required
     * @return the variable as a double
     */
    public double getDoubleSetting(String settingToGet)
    { 
        double setting = 0.0;
        try
        {
        setting = Double.parseDouble(prop.getProperty(settingToGet));
        }
        catch (NullPointerException ex)
        {
            settingNotFound(settingToGet);
        }
        catch (NumberFormatException ex)
        {
            settingFormatError(settingToGet, "Should be a double number.");
        }
           
        return setting;
    }
    
    /**
     * Use to read a numerical integer variable from the settings file
     * 
     * @param settingToGet The variable required
     * @return the variable as an integer
     */
    public int getIntegerSetting(String settingToGet)
    { 
        int setting = 0;
        try
        {
        setting = Integer.parseInt(prop.getProperty(settingToGet));
        }
        catch (NullPointerException ex)
        {
            settingNotFound(settingToGet);
        }
        catch (NumberFormatException ex)
        {
            settingFormatError(settingToGet, "Should be an integer.");
        }
           
        return setting;
    }
    
    /**
     * Displays error for a settings that does not exist
     * @param attemptedSetting The setting the program tried to load
     */
    private void settingNotFound(String attemptedSetting)
    {
        String message = "Tried to get '" + attemptedSetting + "' input from "
                    + "genSettings.in file but setting was not found.";
        stopProgram(message, "");
    }
    
    /**
     * Displays error for a setting input that is not of the correct data type
     * @param attemptedSetting The setting the program tried to load
     * @param messageType String stating the data type that should have been used
     */
    private void settingFormatError(String attemptedSetting, String messageType)
    {
        String message = "Tried to get '" + attemptedSetting + "' input from "
                    + "genSettings.in file but setting was not of correct type.";
        stopProgram(message, messageType);
    }
    
    /**
     * Stops RTWODEE2 from running and displays an error message to the user
     * @param message1 The first line error message to display
     * @param message2 The second line error message to display
     */
    private void stopProgram(String message1, String message2)
    {
        System.out.println("");
        System.out.println(message1);
        if(message2.length()>1)
            System.out.println(message2);
        System.out.println("Simulation terminated");
        System.exit(0);
    }
}



