/* 
 * Copyright (C) 2014 David Barry <david.barry at cancer.org.uk>
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
package DataProcessing;

import UtilClasses.Utilities;
import IAClasses.DataStatistics;
import UtilClasses.GenUtils;
import ij.IJ;
import ij.ImagePlus;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;

/**
 *
 * @author barry05
 */
public class DetectionMapAnalyser implements PlugIn {

    private static File directory;
    private final String delimiter = GenUtils.getDelimiter();

//    public static void main(String args[]) {
//        (new Detection_Map_Analyser()).run(null);
//    }
    
    public DetectionMapAnalyser() {
    }

    public void run(String arg) {
        directory = Utilities.getFolder(directory, "Select Folder", true);
        ImageProcessor signalMap = (IJ.openImage(directory.getAbsolutePath() + delimiter + "SignalMap.tif")).getProcessor();
        ImageStatistics stats = signalMap.getStatistics();
        File mapDir = new File(directory.getAbsolutePath() + delimiter + "Maps");
        File plotDir = new File(directory.getAbsolutePath() + delimiter + "PlotDataFiles");
        File maps[] = mapDir.listFiles();
        File plots[] = plotDir.listFiles();
        int size = plots.length;
//        double frameRate = getFrameRate(directory, ".params") / 60.0;
        File resultsDir = new File(GenUtils.openResultsDirectory(directory.getAbsolutePath()
                + delimiter + "Plot_SD_Sigs", delimiter));
        String headings[] = {"Time_(s)", "Normalised_Signal_SD"};
        for (int i = 0; i < size; i++) {
            int index = getNumber(plots[i].getName(), '_', '.');
            if (index >= 0) {
                ImageProcessor ip = (IJ.openImage(maps[getMapIndex(maps, index, '_', '_')].getAbsolutePath())).getProcessor();
                double sd[] = getMapSDSig(ip, stats.mean);
//                generateFile(resultsDir, sd, "Sig_SD_Plot_" + index + ".txt", headings, frameRate);
                IJ.saveAs((new ImagePlus("", generateCorrMap(ip, 5, 1))), "TIF",
                        resultsDir.getAbsolutePath() + delimiter + "CorrMap_" + index + ".tif");
            }
        }
        GenUtils.showDone(this);
    }

    /**
     * Calculates the standard deviation of pixel values for each column in the
     * specified map and returns them as an array. Only values above the
     * specified threshold are considered.
     *
     * @param map
     * @param threshold
     * @return
     */
    double[] getMapSDSig(ImageProcessor map, double threshold) {
        int width = map.getWidth();
        int height = map.getHeight();
        double sd[] = new double[width];
        ImageStatistics stats = map.getStatistics();
        map.multiply(1.0 / stats.mean);
        double threshold2 = threshold / stats.mean;
        for (int x = 0; x < width; x++) {
            int count = 0;
            for (int y = 0; y < height; y++) {
                if (map.getPixelValue(x, y) > threshold2) {
                    count++;
                }
            }
            if (count > 0) {
                double col[] = new double[count];
                for (int y = 0, i = 0; y < height; y++) {
                    if (map.getPixelValue(x, y) > threshold2) {
                        col[i] = map.getPixelValue(x, y);
                        i++;
                    }
                    DataStatistics dataStats = new DataStatistics(0.05, col, count);
                    sd[x] = dataStats.getStdDev();
                }
            } else {
                sd[x] = 0.0;
            }
        }
        return sd;
    }

    /**
     * Parses a number from a filename. The number is assumed to be demarcated
     * by c1 and c2.
     *
     * @param name
     * @param c1
     * @param c2
     * @return
     */
    int getNumber(String name, char c1, char c2) {
        int u = name.indexOf(c1) + 1;
        int p = name.indexOf(c2, u);
        return Integer.parseInt(name.substring(u, p));
    }

    /**
     * Returns the array index of the file containing the number index, expected
     * to be demarcated within the filename by c1 and c2.
     *
     * @param maps
     * @param index
     * @return
     */
    int getMapIndex(File maps[], int index, char c1, char c2) {
        int size = maps.length;
        for (int i = 0; i < size; i++) {
            int mIndex = getNumber(maps[i].getName(), c1, c2);
            if (mIndex == index) {
                return i;
            }
        }
        return -1;
    }

    /**
     * Prints the specified data in a column under the specified heading in a
     * file with the specified name in the specified directory.
     *
     * @param dir
     * @param data
     * @param index
     */
    void generateFile(File dir, double[] data, String filename, String headings[], double frameRate) {
        File resultsFile;
        PrintWriter resultsStream;
        try {
            resultsFile = new File(dir.getAbsolutePath() + delimiter
                    + filename);
            resultsStream = new PrintWriter(new FileOutputStream(resultsFile));
        } catch (FileNotFoundException e) {
            System.out.println("Error: Failed to create results file.\n");
            System.out.println(e.toString());
            return;
        }
        resultsStream.println(dir.getParent());
        for (String h : headings) {
            resultsStream.print(h + "\t");
        }
        resultsStream.println();
        for (int i = 0; i < data.length; i++) {
            resultsStream.print(i / frameRate + "\t");
            resultsStream.println(data[i]);
        }
        resultsStream.close();
    }

    /**
     * Extract the frame rate from the specified parameter file in the specified
     * directory.
     *
     * @param dir
     * @param filename
     * @return
     */
//    double getFrameRate(File dir, String filename) {
//        File params = new File(dir.getAbsolutePath() + delimiter + filename);
//        return (new FileReader()).readParam(params, 1, "timeRes");
//    }

    int correlateSingleCol(ImageProcessor map1, ImageProcessor map2, int yposc1,
            int radius, int col, int colstep) {
        double mean1 = map1.getStatistics().mean;
        double mean2 = map2.getStatistics().mean;
        int height = map1.getHeight();
        int window = 2 * radius + 1;
        double max = -Double.MAX_VALUE;
        int maxindex = 0;
        for (int yposc2 = (yposc1 < radius) ? radius : yposc1 - radius; yposc2 <= yposc1 + radius && yposc2 < height - radius; yposc2++) {
            double sum = 0.0;
            for (int i = 0; i <= window; i++) {
                int yc1 = yposc1 - radius + i;
                int yc2 = yposc2 - radius + i;
                sum += (map1.getPixelValue(col, yc2) - mean1) * (map2.getPixelValue(col - colstep, yc1) - mean2);
            }
            if (sum > max) {
                max = sum;
                maxindex = yposc1 - yposc2;
            }
        }
        return maxindex;
    }

    ImageProcessor generateCorrMap(ImageProcessor map, int radius, int colstep) {
        int width = map.getWidth();
        int height = map.getHeight();
        FloatProcessor corrMap = new FloatProcessor(width, height);
        for (int col = colstep; col < width; col++) {
            for (int row = radius; row < height - radius; row++) {
                corrMap.putPixel(col, row, correlateSingleCol(map, map, row, radius, col, colstep));
            }
        }
        return corrMap;
    }
}
