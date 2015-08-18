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
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;

/**
 *
 * @author barry05
 */
public class Data_File_Averager implements PlugIn {

    private final String TITLE = "Data Averager";
    private static int headerSize = 1;
    private static int numParams;
    private static boolean normalise = true, truncate = true, colate = true;
    private static boolean normParams[], selection[];
    public static final int VEL_INDEX = 2, TIME_INDEX = 1, LENGTH_INDEX = 8;
    private static File directory;
    private final String directoryDelim = GenUtils.getDelimiter(); // delimiter in directory strings
    private final String paramDelim = ",";
    private final String FILE_LIST = "file_list.txt";
    private final String AGG_DATA = "Aggregated_Data";
    private final String COL_DATA = "Collated_Data";
    private final String MEAN_DATA = "mean_data.csv";
    private final String headings[];
    private final String normHeadings[];

//    public static void main(String args[]) {
//        new Data_File_Averager().run(null);
//        System.exit(0);
//    }
    public Data_File_Averager(String headings[], String normHeadings[]) {
        this.headings = headings;
        this.normHeadings = normHeadings;
    }

    public void run(String arg) {
        if (arg == null) {
            directory = Utilities.getFolder(directory, "Select Directory", true);
        } else {
            directory = new File(arg);
        }
        if (directory == null) {
            return;
        }
        cleanDirectory();
        File files[] = directory.listFiles();
        int numOfFiles = files.length;
//        if (!showDialog()) {
//            return;
//        }
        FileReader reader = new FileReader(numOfFiles, headerSize);
        reader.getParamList(files, paramDelim);
        numParams = reader.getNumParams();
        ArrayList<Double> data[][] = new ArrayList[numOfFiles][numParams];
        double maxima[][] = new double[numOfFiles][numParams];
        double minima[][] = new double[numOfFiles][numParams];
        if (selection == null || selection.length != numParams) {
            selection = new boolean[numParams];
            for (int i = 0; i < numParams; i++) {
                selection[i] = true;
            }
        }
        try {
            reader.readData(data, files, paramDelim);
        } catch (Exception e) {
            IJ.log(e.toString());
            return;
        }
//        String headings[] = reader.getParamsArray();
//        selection = showSelectionDialog(headings, "Specify parameters to be output", numParams, selection);
        if (truncate) {
            truncateData(data, numParams, numOfFiles, VEL_INDEX);
        }
        getExtrema(data, minima, maxima, numParams, numOfFiles);
//        performAnalysis(data, VEL_INDEX, TIME_INDEX, numOfFiles, files);
        if (normalise) {
            if (normParams == null || normParams.length != numParams) {
                normParams = new boolean[numParams];
                for (int i = 0; i < numParams; i++) {
                    normParams[i] = false;
                    for (int j = 0; j < normHeadings.length; j++) {
                        if (headings[i].compareTo(normHeadings[j]) == 0) {
                            normParams[i] = true;
                        }
                    }
                }
            }
//            normParams = showSelectionDialog(headings, "Specify parameters to be normalised", numParams, normParams);
            normaliseData(data, numOfFiles, numParams, minima, maxima);
        }
        if (colate) {
            colateData(directory, headings, numOfFiles, data);
        }
        outputMeanData(directory, headings, numOfFiles, data);
        aggregateData(directory, headings, numOfFiles, data, numParams, selection);
//        IJ.saveAs(new ImagePlus("", buildHeatMap(data, 15, 15, VEL_INDEX, LENGTH_INDEX, 32.0, 2.1, -10.0, 0.75)), "TIF", "c:\\users\\barry05\\desktop\\heatmap.tif");
        outputFileList(directory, reader.getFilenames(), numOfFiles);
        GenUtils.showDone(this);
    }

//    boolean showDialog() {
//        GenericDialog gd = new GenericDialog(TITLE);
////        gd.addNumericField("Header Size", headerSize, 0);
////        gd.addNumericField("Number of Parameters", numParams, 0);
//        gd.addCheckbox("Normalise Data?", normalise);
////        gd.addCheckbox("Truncate Data?", truncate);
//        gd.addCheckbox("Collate Data?", colate);
//        gd.showDialog();
//        if (gd.wasCanceled()) {
//            return false;
//        }
////        headerSize = (int) Math.round(gd.getNextNumber());
////        numParams = (int) Math.round(gd.getNextNumber());
//        normalise = gd.getNextBoolean();
////        truncate = gd.getNextBoolean();
//        colate = gd.getNextBoolean();
//        return true;
//    }
    boolean[] showSelectionDialog(String[] headings, String dialogHeading, int numParams, boolean[] defaults) {
        GenericDialog gd = new GenericDialog(dialogHeading);
        boolean selection[] = new boolean[numParams];
        for (int i = 0; i < numParams; i++) {
            gd.addCheckbox(headings[i], defaults[i]);
        }
        gd.showDialog();
        if (gd.wasCanceled()) {
            return null;
        }
        for (int i = 0; i < numParams; i++) {
            selection[i] = gd.getNextBoolean();
        }
        return selection;
    }

    boolean dataAnalysis(double[] data, int n, DataAnalytics analytics) {
        /*
         * Find and first and second zero crossings of acceleration data and
         * limit calculation of mean velocity within these bounds
         */
        DataStatistics stats = new DataStatistics(0.05, data, n);
        analytics.setMaxVal(stats.getMax());
        analytics.setMinVal(stats.getMin());
        analytics.setPeakIndex(stats.getMaxIndex());
        ArrayList<Integer> zc = stats.getZerocrossings();
        if (zc == null || zc.size() < 1) {
            return false;
        }
        int zcIndex = 0;
        while (zcIndex < zc.size() && zc.get(zcIndex).intValue() < analytics.getPeakIndex()) {
            zcIndex++;
        }
        if (zcIndex < zc.size()) {
            int crossIndex = zc.get(zcIndex);
            double truncdata[] = new double[n - crossIndex - 1];
            System.arraycopy(data, crossIndex, truncdata, 0, truncdata.length);
            DataStatistics truncstats = new DataStatistics(0.05, truncdata, truncdata.length);
            analytics.setZeroCrossing(crossIndex);
            analytics.setNegMean(truncstats.getMean());
            return true;
        } else {
            return false;
        }
    }

    void performAnalysis(ArrayList<Double>[][] data, int pIndex, int tIndex, int numSets, File[] files) {
        char mu = GenUtils.mu;
        String mus = "(" + mu + "m/s)";
        System.out.println("File\tZero_Crossing_(s)\tMean_Neg_Velocity_" + mus
                + "\tMin_Velocity_" + mus + "\tMax_Velocity_" + mus + "");
        for (int i = 0; i < numSets; i++) {
            DataAnalytics analytics = new DataAnalytics();
            if (data[i][pIndex] != null) {
                int size = data[i][pIndex].size();
                double thisdata[] = new double[size];
                for (int j = 0; j < size; j++) {
                    thisdata[j] = data[i][pIndex].get(j);
                }
                if (dataAnalysis(thisdata, size, analytics)) {
                    double timeRes = data[i][tIndex].get(1) - data[i][tIndex].get(0);
                    System.out.println(files[i].getName() + "\t"
                            + ((analytics.getZeroCrossing() - analytics.getPeakIndex()) * timeRes)
                            + "\t" + analytics.getNegMean() + "\t" + analytics.getMinVal()
                            + "\t" + analytics.getMaxVal());
                }
            }
        }
    }

    void truncateData(ArrayList<Double> data[][], int numParams, int numOfFiles, int VEL_INDEX) {
        for (int i = 0; i < numOfFiles; i++) {
            boolean negvel = false, posvel = false, end = false;
            int index = 0;
            int size = data[i][0].size();
            while (!end && index < size) {
                boolean nan = false;
                for (int k = 0; k < numParams; k++) {
                    if (selection[k]) {
                        if (data[i][k] != null) {
                            nan = nan || Double.isNaN(data[i][k].get(index));
                        }
                    }
                }
                double vel = data[i][VEL_INDEX].get(index);
                if (!nan && !(negvel && vel >= 0.0)) {
                    if (posvel && vel < 0.0) {
                        negvel = true;
                    } else if (vel >= 0.0) {
                        posvel = true;
                    }
                } else {
                    for (int k = 0; k < numParams; k++) {
                        if (data[i][k] != null) {
                            removeElements(data[i][k], index, size - 1);
                        }
                    }
                    end = true;
                }
                index++;
            }
        }
    }

    void getExtrema(ArrayList<Double> data[][], double[][] minima, double[][] maxima, int numParams, int numOfFiles) {
        for (int i = 0; i < numOfFiles; i++) {
            Arrays.fill(maxima[i], -Double.MAX_VALUE);
            Arrays.fill(minima[i], Double.MAX_VALUE);
            if (data[i][0] != null) {
                int size = data[i][0].size();
                for (int m = 0; m < size; m++) {
                    for (int k = 0; k < numParams; k++) {
                        if (selection[k]) {
                            if (data[i][k] != null) {
                                double val = data[i][k].get(m);
                                if (val > maxima[i][k]) {
                                    maxima[i][k] = val;
                                }
                                if (val < minima[i][k] && val >= 0.0) {
                                    minima[i][k] = val;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * Remove elements from objects between start and end inclusive.
     *
     * @param objects
     * @param start
     * @param end
     */
    void removeElements(ArrayList objects, int start, int end) {
        for (int i = 0; i <= end - start; i++) {
            objects.remove(start);
        }
    }

    void normaliseData(ArrayList<Double>[][] data, int numOfFiles, int numParams,
            double[][] minima, double[][] maxima) {
        for (int k = 0; k < numOfFiles; k++) {
            for (int j = 0; j < numParams; j++) {
                if (data[k][j] != null) {
                    int size = data[k][j].size();
                    double norm = 1.0;
                    double offset = 0.0;
                    if (normParams[j]) {
                        offset = minima[k][j];
                        norm = maxima[k][j] - offset;
                    }
                    for (int t = 0; t < size; t++) {
                        if (data[k][j] != null && data[k][j].size() > t) {
                            data[k][j].set(t, (data[k][j].get(t) - offset) / norm);
//                        System.out.println("" + data[k][j].get(t));
                        }
                    }
                }
            }
        }
    }

    void outputMeanData(File directory, String[] headings, int numOfFiles, ArrayList<Double>[][] data) {
        File thisMeanData;
        PrintWriter thisDataStream;
        try {
            thisMeanData = new File(directory + directoryDelim + MEAN_DATA);
            thisDataStream = new PrintWriter(new FileOutputStream(thisMeanData));
        } catch (FileNotFoundException e) {
            IJ.error(e.toString());
            return;
        }
        thisDataStream.println(directory.getAbsolutePath());
        for (String h : headings) {
            thisDataStream.print(h + ", , ");
        }
        thisDataStream.println();
        boolean cont = true;
        int t = 0;
        while (cont) {
            int n = 0;
            for (int i = 0; i < numOfFiles; i++) {
                if (data[i][0] != null && data[i][0].size() > t) {
                    n++;
                }
            }
            if (n == 0) {
                thisDataStream.close();
                cont = false;
            } else {
                for (int j = 0; j < numParams; j++) {
                    if (selection[j]) {
                        double thisdata[] = new double[n];
                        for (int k = 0, index = 0; k < numOfFiles; k++) {
                            if (data[k][j] != null && data[k][j].size() > t) {
                                thisdata[index] = data[k][j].get(t);
                                index++;
                            }
                        }
                        double mean = DataStatistics.calcMean(thisdata);
                        double stdDev = DataStatistics.calcStdDev(thisdata, n, mean);
                        thisDataStream.print(mean + ", " + stdDev + ", ");
                    } else {
                        thisDataStream.print(", , ");
                    }
                }
                thisDataStream.print(n);
                thisDataStream.println();
                t++;
            }
        }
        thisDataStream.close();
    }

    void colateData(File directory, String headings[], int numOfFiles, ArrayList<Double>[][] data) {
        File colateDir = GenUtils.createDirectory(directory + directoryDelim + COL_DATA, false);
        for (int p = 0; p < numParams; p++) {
            if (selection[p]) {
                File paramFile;
                PrintWriter paramStream;
                try {
                    headings[p] = GenUtils.checkFileSep(headings[p], '-');
                    paramFile = new File(colateDir + directoryDelim + headings[p] + ".txt");
                    paramStream = new PrintWriter(new FileOutputStream(paramFile));
                } catch (FileNotFoundException e) {
                    IJ.error(e.toString());
                    return;
                }
                boolean cont = true;
                int t = 0;
                while (cont) {
                    cont = false;
                    for (int k = 0; k < numOfFiles; k++) {
                        if (data[k][p] != null && data[k][p].size() > t) {
                            paramStream.print(data[k][p].get(t) + "\t");
                            cont = true;
                        } else {
                            paramStream.print(" \t");
                        }
                    }
                    paramStream.println();
                    t++;
                }
                paramStream.close();
            }
        }
    }

    void aggregateData(File directory, String[] headings, int numOfFiles, ArrayList<Double>[][] data,
            int numParams, boolean[] selection) {
        File aggDir = GenUtils.createDirectory(directory + directoryDelim + AGG_DATA, false);
        File paramFile;
        PrintWriter paramStream;
        try {
            paramFile = new File(aggDir + directoryDelim + "aggregated_data.txt");
            paramStream = new PrintWriter(new FileOutputStream(paramFile));
        } catch (FileNotFoundException e) {
            IJ.error(e.toString());
            return;
        }
        for (int p = 0; p < numParams; p++) {
            if (selection[p]) {
                paramStream.print(headings[p] + "\t");
            }
        }
        paramStream.println();
        for (int j = 0; j < numOfFiles; j++) {
            if (data[j][0] != null) {
                int length = data[j][0].size();
                for (int l = 0; l < length; l++) {
                    for (int p = 0; p < numParams; p++) {
                        if (selection[p]) {
                            if (data[j][p] != null) {
                                paramStream.print(data[j][p].get(l) + "\t");
                            } else {
                                paramStream.print("\t\t");
                            }
                        }
                    }
                    paramStream.println();
                }
            }
        }
        paramStream.close();
    }

    void outputFileList(File directory, String[] filelist, int numOfFiles) {
        File thisFile;
        PrintWriter thisStream;
        try {
            thisFile = new File(directory + directoryDelim + FILE_LIST);
            thisStream = new PrintWriter(new FileOutputStream(thisFile));
        } catch (FileNotFoundException e) {
            IJ.error(e.toString());
            return;
        }
        for (String f : filelist) {
            thisStream.println(f);
        }
        thisStream.close();
    }

    private void cleanDirectory() {
        File files[] = directory.listFiles();
        for (int i = 0; i < files.length; i++) {
            if (!FilenameUtils.getBaseName(files[i].getAbsolutePath()).contains("bleb_data")) {
                FileUtils.deleteQuietly(files[i]);
            }
        }
    }

//    ImageStack buildHeatMap(ArrayList<Double>[][] data, int width, int height, int xIndex, int yIndex, double maxX, double maxY, double minX, double minY) {
//        int stackSize = 360;
//        ImageStack output = new ImageStack(width, height);
//        int counts[] = new int[360];
//        for (int i = 0; i < stackSize; i++) {
//            ShortProcessor bp = new ShortProcessor(width, height);
//            bp.setColor(0);
//            bp.fill();
//            output.addSlice(bp);
//        }
//        for (int i = 0; i < data.length; i++) {
//            int size = data[i][0].size();
//            for (int j = 0; j < size && j < stackSize; j++) {
//                ImageProcessor ip = output.getProcessor(j + 1);
//                int x = (int) Math.round(width * (data[i][xIndex].get(j) - minX) / (maxX - minX));
//                int y = (int) Math.round(height * (data[i][yIndex].get(j) - minY) / (maxY - minY));
//                if (x >= 0 && x < width && y >= 0 && y < height) {
//                    ip.putPixel(x, y, 1 + ip.get(x, y));
//                    counts[j]++;
//                }
//            }
//        }
//        for (int i = 0; i < stackSize; i++) {
//            (output.getProcessor(i+1)).multiply(Short.MAX_VALUE/counts[i]);
//        }
//        
//        return output;
//    }
}
