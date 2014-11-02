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

import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.URLDecoder;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.List;
import java.util.Locale;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Class for reading and writing data from/to external files
 * @author Stephen Denekamp
 * @version 1.0: Created
 */
public class ExternalDataHandler
{
    private String inputLocation;
    private String outputLocation;
    private final String OUTPUT_FILE_EXTENSION = ".out";
    
    /**
     * Constructor for the External Matrix Reader
     * 
     * @param scenarioFolder The name of the folder the scenario input files
     * are located in
     */
    public ExternalDataHandler(String scenarioFolder)
    {             
        ClassLoader loader = getClass().getClassLoader();

        try {
            //Set correct locations for the input and output folders
            String tempPath = loader.getResource("RTWODEE2/dataFiles/"+scenarioFolder+"/").getPath();
            tempPath = checkForWindowsOS(tempPath);
            inputLocation = URLDecoder.decode(tempPath, "utf-8");
            outputLocation = inputLocation+"output/";
            new File(outputLocation).mkdirs(); //create output folder
        } catch (UnsupportedEncodingException ex) {
            Logger.getLogger(ExternalDataHandler.class.getName()).log(Level.SEVERE, null, ex);
            stopProgram(" ", "");
        } catch (NullPointerException ex){
            String msg1 = "Looking for RTWODEE2/dataFiles/"+scenarioFolder;
            stopProgram(msg1, "Unable to find folder '"+scenarioFolder+"'");
        }

    }
    
    /**
     * Creates a matrix from an external file
     * 
     * @param filename the name of the external file
     * @return the matrix as a double[][] 
     */
    public double[][] matrixReader(String filename)
    {
        String[][] stringMatrix = getStringMatrixFromFile(filename);
        double[][] matrix = new double[stringMatrix.length][stringMatrix[0].length];
        
        for(int i = 0; i<stringMatrix.length; i++)
                for (int j = 0; j<stringMatrix[0].length; j++)
                    matrix[i][j]=Double.parseDouble(stringMatrix[i][j]);  
        
        return matrix;
    }
     
    /**
     * Outputs a matrix to a text file. If a file already exists with the same
     * name it is overridden.
     * 
     * @param matrix The matrix to be printed to file
     * @param filename The filename to be used for the created file
     */
    public void matrixWriter(double[][] matrix, String filename) 
    {    
        createFile(filename);
        
        for (double[] matrixRow : matrix) 
        {
            String row = "";
            for(double rowItem : matrixRow)
            { 
                //format number with scientific notation
                String formattedNumber = String.format("%e", rowItem);
                
                //add current number to list
                row += formattedNumber + " ";
            }
            
            //append to external file
            appendToFile(filename, row);
            appendNewLineToFile(filename);
        } 
    }
    
    /**
     * Reads an external text file into an integer 2d array. 
     * Used for the points file to create an array of x,y coordinates.
     * 
     * @param filename the name of the external file
     * @return an integer[][] array that represents the x,y coordinates of the points file 
     */
    public int[][] pointsReader(String filename)
    {
        String[][] stringMatrix = getStringMatrixFromFile(filename);
        int[][] matrix = new int[stringMatrix.length][stringMatrix[0].length];
        
        try
        {
            for(int i = 0; i<stringMatrix.length; i++)
                for (int j = 0; j<stringMatrix[0].length; j++)
                    matrix[i][j]=Integer.parseInt(stringMatrix[i][j]);  

        } 
        catch (NumberFormatException ex) 
        {
            String errorMessage ="File '"+filename+"' not formatted correctly:";
            String systemMessage =ex.getMessage();
            stopProgram(errorMessage, systemMessage);    
        }
        
        return matrix;
    }
       
    /**
     * Format and appends data to an existing output file
     *
     * @param filename Name of file to be appended to
     * @param data String of data to be append to file
     */
    public void appendDataToFile(String filename, String data)
    {
        data = " " + data;              //add horizontal space to data
        appendToFile(filename, data);   //write data
    }
        
    /**
     * Creates new blank files for 'allData.out', 'euler_height.out', and
     * 'euler_density.out'
     *
     */
    public void createAppendFiles()
    {
        //create  blank files
        createFile("allData");
        createFile("euler_height");
        createFile("euler_density");
        
        //Add column headings to 'allData' file
        appendDataToFile("allData", String.format("%7s", "Counter"));
        String[] headings = {"Time", "sumVolume", "sumMomentum", "xCentre", "yCentre", "Area"};
        for (String heading : headings) 
        {
            appendDataToFile("allData", String.format("%12s", heading));
        }
        appendNewLineToFile("allData");
    }
    
    /**
     * Adds a new link to an external file
     * 
     * @param filename The name of the external file 
     */
    public void appendNewLineToFile(String filename)
    {
        appendToFile(filename, String.format("%n", ""));
    }
    
    /**
     * Creates an empty file in the output folder. If the file already exists
     * it is overridden
     *
     * @param filename Name of file to be created
     */
    private void createFile(String filename)
    {
        String data = "";
        String filepath = outputLocation + filename + OUTPUT_FILE_EXTENSION;
        
        try {

            Files.write(Paths.get(filepath), data.getBytes());

        } catch (IOException ex) {
            Logger.getLogger(Simulation.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    /**
     * Appends data to an existing output file
     *
     * @param filepath Name of file to be appended to
     * @param data String of data to be append to file
     */
    private void appendToFile(String filename, String data)
    {
        String filepath = outputLocation + filename + OUTPUT_FILE_EXTENSION;
        try {
            Files.write(Paths.get(filepath), data.getBytes(), StandardOpenOption.APPEND);

        } catch (IOException ex) {
            Logger.getLogger(Simulation.class.getName()).log(Level.SEVERE, null, ex);
        }    
    }

    /**
     * Creates a String[][] array from an external text file
     * 
     * @param filename The name of the external file
     * @return the String[][] array
     */
    private String[][] getStringMatrixFromFile(String filename)
    {
 
        File matrixFile = new File(inputLocation + filename);
        String[][] matrix = null;
        
        try {
            //read matrix file into List
            List<String> matrixList = Files.readAllLines(matrixFile.toPath(), Charset.forName("utf-8"));
            
            //covert list into an array of each row
            String [] matrixRows = matrixList.toArray(new String[matrixList.size()]);
            
            //covert each row in an array 
            //(so we now have a double array representing the matrix)
            String[][] matrixString = new String[matrixList.size()][1];
            for (int i = 0; i<matrixRows.length;i++)
                matrixString[i] = matrixRows[i].trim().split("\\s+");
            
            matrix = matrixString;

        } catch (IOException ex) {
            String errorMessage ="Input file '"+filename+"' not found in ";
            stopProgram(errorMessage, inputLocation);
        } 
        
        return matrix;
    }
    
    /**
     * Checks if the Operating System is Windows and if so removes leading
     * slash from the filePath. 
     * (e.g. filePath should begin with 'C:/' not '/C:/')
     * @param filePath the unedited filePath
     * @return the updated file path if OS is Windows, otherwise returns path
     * unchanged
     */
    private String checkForWindowsOS(String filePath)
    {
        String osName = System.getProperty("os.name").toLowerCase(Locale.ENGLISH);
        
        if(osName.contains("win"))
            filePath = filePath.substring(1);
        
        return filePath;
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
