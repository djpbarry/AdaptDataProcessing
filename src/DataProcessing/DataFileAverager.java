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
import ij.gui.Plot;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FilenameFilter;
import java.io.PrintWriter;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Arrays;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.filefilter.SuffixFileFilter;

/**
 *
 * @author barry05
 */
public class DataFileAverager {

//    private final String TITLE = "Data Averager";
    private final int HEAD_SIZE = 1;
    private int numParams;
    private final boolean NORM = true, TRUNC = true, COLLATE = true;
    private boolean normParams[], selection[];
    public int velIndex = -1, timeIndex = -1;
    private static File directory;
    private final String DIR_DELIM = GenUtils.getDelimiter(); // delimiter in directory strings
    private final String PARAM_DELIM = ",";
    private final String FILE_LIST = "file_list.txt";
    private final String AGG_DATA = "Aggregated_Data";
    private final String COL_DATA = "Collated_Data";
    private final String MEAN_DATA = "mean_data.csv";
    private final String VEL;
    private final String TIME;
    private String headings[];
    private final String NORM_HEADINGS[];
    private final boolean DISPLAY_PLOTS;
    private final Charset charSet;

    public DataFileAverager(String[] headings, String[] normHeadings, boolean displayPlots, String Vel, String Time, Charset charSet) {
        this.headings = headings;
        this.NORM_HEADINGS = normHeadings;
        this.DISPLAY_PLOTS = displayPlots;
        this.VEL = Vel;
        this.TIME = Time;
        this.charSet = charSet;
    }

    public void run(String dirName) {
        try {
            if (dirName == null) {
                directory = Utilities.getFolder(directory, "Select Directory", true);
            } else {
                directory = new File(dirName);
            }
        } catch (Exception e) {
            GenUtils.error("Could not open directory.");
        }
        if (directory == null) {
            return;
        }
        cleanDirectory();
        File files[] = directory.listFiles((FilenameFilter) (new SuffixFileFilter(".csv")));
        int numOfFiles = files.length;
        FileReader reader = new FileReader(numOfFiles, HEAD_SIZE, charSet);
        try {
            reader.getParamList(files, PARAM_DELIM);
        } catch (Exception e) {
            IJ.log("Failed to read parameter list.");
            IJ.log(e.toString());
            return;
        }
        numParams = reader.getNumParams();
        ArrayList<ArrayList<ArrayList<Double>>> data = new ArrayList<>();
        double maxima[][] = new double[numOfFiles][numParams];
        double minima[][] = new double[numOfFiles][numParams];
        if (selection == null || selection.length != numParams) {
            selection = new boolean[numParams];
            for (int i = 0; i < numParams; i++) {
                selection[i] = true;
            }
        }
        try {
            reader.readData(data, files, PARAM_DELIM);
        } catch (Exception e) {
            IJ.log("Failed to read data files.");
            IJ.log(e.toString());
            return;
        }
        if (headings == null) {
            headings = reader.getParamsArray();
//            selection = showSelectionDialog(headings, "Specify parameters to be output", numParams, selection);
        }
        getParamIndices(headings);
        if (TRUNC) {
            truncateData(data, numParams, numOfFiles, velIndex);
        }
        getExtrema(data, minima, maxima, numParams, numOfFiles);
        if (NORM) {
            if (normParams == null || normParams.length != numParams) {
                normParams = new boolean[numParams];
                for (int i = 0; i < numParams; i++) {
                    normParams[i] = false;
                    for (int j = 0; j < NORM_HEADINGS.length; j++) {
                        if (headings[i].compareTo(NORM_HEADINGS[j]) == 0) {
                            normParams[i] = true;
                        }
                    }
                }
            }
//            normParams = showSelectionDialog(headings, "Specify parameters to be normalised", numParams, normParams);
            normaliseData(data, numOfFiles, numParams, minima, maxima);
        }
//        if (colate) {
//            colateData(directory, headings, numOfFiles, data);
//        }
        ArrayList<ArrayList<Double>> meanData;
        try {
            meanData = calcMeanData(directory, headings, numOfFiles, data);
        } catch (Exception e) {
            IJ.log("Failed to calculate data mean.");
            IJ.log(e.toString());
            return;
        }
        if (DISPLAY_PLOTS) {
            plotData(meanData);
        }
//        aggregateData(directory, headings, numOfFiles, data, numParams, selection);
        try {
            outputFileList(directory, reader.getFilenames(), numOfFiles);
        } catch (Exception e) {
            IJ.log("Failed to generate output files.");
            IJ.log(e.toString());
            return;
        }
        GenUtils.showDone(this);
    }

    void getParamIndices(String[] headings) {
        for (int i = 0; i < headings.length; i++) {
            if (headings[i].compareTo(VEL) == 0) {
                velIndex = i;
            }
            if (headings[i].compareTo(TIME) == 0) {
                timeIndex = i;
            }
        }
    }

//    boolean[] showSelectionDialog(String[] headings, String dialogHeading, int numParams, boolean[] defaults) {
//        GenericDialog gd = new GenericDialog(dialogHeading);
//        boolean selection[] = new boolean[numParams];
//        for (int i = 0; i < numParams; i++) {
//            gd.addCheckbox(headings[i], defaults[i]);
//        }
//        gd.showDialog();
//        if (gd.wasCanceled()) {
//            return null;
//        }
//        for (int i = 0; i < numParams; i++) {
//            selection[i] = gd.getNextBoolean();
//        }
//        return selection;
//    }
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

    void truncateData(ArrayList<ArrayList<ArrayList<Double>>> data, int numParams, int numOfFiles, int VEL_INDEX) {
        for (int i = 0; i < numOfFiles; i++) {
            boolean negvel = false, posvel = false, end = false;
            int index = 0;
            int size = data.get(i).get(0).size();
            while (!end && index < size) {
                boolean nan = false;
                for (int k = 0; k < numParams; k++) {
                    if (selection[k]) {
                        if (data.get(i).get(k) != null) {
                            nan = nan || Double.isNaN(data.get(i).get(k).get(index));
                        }
                    }
                }
                double vel = data.get(i).get(VEL_INDEX).get(index);
                if (!nan && !(negvel && vel >= 0.0)) {
                    if (posvel && vel < 0.0) {
                        negvel = true;
                    } else if (vel >= 0.0) {
                        posvel = true;
                    }
                } else {
                    for (int k = 0; k < numParams; k++) {
                        if (data.get(i).get(k) != null) {
                            removeElements(data.get(i).get(k), index, size - 1);
                        }
                    }
                    end = true;
                }
                index++;
            }
        }
    }

    void getExtrema(ArrayList<ArrayList<ArrayList<Double>>> data, double[][] minima, double[][] maxima, int numParams, int numOfFiles) {
        for (int i = 0; i < numOfFiles; i++) {
            Arrays.fill(maxima[i], -Double.MAX_VALUE);
            Arrays.fill(minima[i], Double.MAX_VALUE);
            if (data.get(i).get(0) != null) {
                int size = data.get(i).get(0).size();
                for (int m = 0; m < size; m++) {
                    for (int k = 0; k < numParams; k++) {
                        if (selection[k]) {
                            if (data.get(i).get(k) != null) {
                                double val = data.get(i).get(k).get(m);
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

    void normaliseData(ArrayList<ArrayList<ArrayList<Double>>> data, int numOfFiles, int numParams,
            double[][] minima, double[][] maxima) {
        for (int k = 0; k < numOfFiles; k++) {
            for (int j = 0; j < numParams; j++) {
                if (data.get(k).get(j) != null) {
                    int size = data.get(k).get(j).size();
                    double norm = 1.0;
                    double offset = 0.0;
                    if (normParams[j]) {
                        offset = minima[k][j];
                        norm = maxima[k][j] - offset;
                    }
                    for (int t = 0; t < size; t++) {
                        if (data.get(k).get(j) != null && data.get(k).get(j).size() > t) {
                            data.get(k).get(j).set(t, (data.get(k).get(j).get(t) - offset) / norm);
//                        System.out.println("" + data[k][j].get(t));
                        }
                    }
                }
            }
        }
    }

    ArrayList<ArrayList<Double>> calcMeanData(File directory, String[] headings, int numOfFiles, ArrayList<ArrayList<ArrayList<Double>>> data) throws FileNotFoundException {
        File thisMeanData = new File(directory + DIR_DELIM + MEAN_DATA);
        PrintWriter thisDataStream = new PrintWriter(new FileOutputStream(thisMeanData));
        thisDataStream.println(directory.getAbsolutePath());
        for (String h : headings) {
            thisDataStream.print(h + ", , ");
        }
        thisDataStream.println("N");
        for (String h : headings) {
            thisDataStream.print("Mean,Standard Deviation,");
        }
        thisDataStream.println();
        boolean cont = true;
        int t = 0;
        ArrayList<ArrayList<Double>> meanData = new ArrayList<>();
        while (cont) {
            int n = 0;
            for (int i = 0; i < numOfFiles; i++) {
                if (data.get(i).get(0) != null && data.get(i).get(0).size() > t) {
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
                            if (data.get(k).get(j) != null && data.get(k).get(j).size() > t) {
                                thisdata[index] = data.get(k).get(j).get(t);
                                index++;
                            }
                        }
                        double mean = DataStatistics.calcMean(thisdata);
                        double stdDev = DataStatistics.calcStdDev(thisdata, n, mean);
                        thisDataStream.print(mean + ", " + stdDev + ", ");
                        if (meanData.get(j) == null) {
                            meanData.add(new ArrayList<>());
                        }
                        meanData.get(j).add(mean);
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
        return meanData;
    }

    void plotData(ArrayList<ArrayList<Double>> data) {
        for (int i = 0; i < data.size(); i++) {
            if (i != timeIndex && headings[i].compareTo(TIME) != 0) {
                double xvals[] = new double[data.get(timeIndex).size()];
                double yvals[] = new double[data.get(i).size()];
                for (int j = 0; j < xvals.length; j++) {
                    xvals[j] = data.get(timeIndex).get(j);
                    yvals[j] = data.get(i).get(j);
                }
                Plot plot = new Plot(headings[i], headings[timeIndex], headings[i], xvals, yvals);
                plot.show();
            }
        }
    }

    void colateData(File directory, String headings[], int numOfFiles, ArrayList<Double>[][] data) throws FileNotFoundException {
        File colateDir = GenUtils.createDirectory(directory + DIR_DELIM + COL_DATA, false);
        for (int p = 0; p < numParams; p++) {
            if (selection[p]) {
                headings[p] = GenUtils.checkFileSep(headings[p], '-');
                File paramFile = new File(colateDir + DIR_DELIM + headings[p] + ".txt");
                PrintWriter paramStream = new PrintWriter(new FileOutputStream(paramFile));
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
            int numParams, boolean[] selection) throws FileNotFoundException {
        File aggDir = GenUtils.createDirectory(directory + DIR_DELIM + AGG_DATA, false);
        File paramFile = new File(aggDir + DIR_DELIM + "aggregated_data.txt");
        PrintWriter paramStream = new PrintWriter(new FileOutputStream(paramFile));
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

    void outputFileList(File directory, String[] filelist, int numOfFiles) throws FileNotFoundException {
        File thisFile = new File(directory + DIR_DELIM + FILE_LIST);
        PrintWriter thisStream = new PrintWriter(new FileOutputStream(thisFile));
        for (String f : filelist) {
            thisStream.println(f);
        }
        thisStream.close();
    }

    private void cleanDirectory() {
        File files[] = directory.listFiles();
        for (int i = 0; i < files.length; i++) {
            String filename = FilenameUtils.getBaseName(files[i].getAbsolutePath());
            if (filename.contains(FilenameUtils.getBaseName(MEAN_DATA))
                    || filename.contains(FilenameUtils.getBaseName(FILE_LIST))) {
                FileUtils.deleteQuietly(files[i]);
            }
        }
    }
}
