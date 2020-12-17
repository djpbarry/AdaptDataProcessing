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
package net.calm.adaptdataprocessing.DataProcessing;

import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

import java.io.File;
import java.util.ArrayList;
import java.util.Scanner;

/**
 *
 * @author barry05
 */
public class DataResampler implements PlugIn {

    private final String TITLE = "Data Resampler";
    private static int headerSize = 1;
    private static int numParams = 10;
    private static double newSampleRate = 0.25;
    private boolean normLength = false;
    private double coeffs[] = {-0.0057, 1.687};
    private int window = 5;

//    public static void main(String args[]) {
//        new Data_Resampler().run(null);
//        System.exit(0);
//    }

    public void run(String arg) {
//        String delimiter = GenUtils.getDelimiter();
//        File init = null;
//        if (IJ.getInstance() != null) {
//            init = new File(IJ.getDirectory("current"));
//        }
//        File inDir = Utilities.getFolder(init, "Select Directory", true);
//        File files[] = inDir.listFiles();
//        int numOfFiles = files.length;
//        if (!showDialog()) {
//            return;
//        }
//        File outDir = new File(GenUtils.openResultsDirectory(inDir.getAbsolutePath()
//                + delimiter + "Resampled", delimiter));
//        ArrayList<Double> data[][] = new ArrayList[numOfFiles][numParams];
//        FileReader reader = new FileReader(numOfFiles, headerSize);
//        reader.getParamList(files, ",\\s");
//        reader.readData(data, files, ",\\s");
//        String headings = reader.getParamString();
//        String filenames[] = reader.getFilenames();
//
//        for (int i = 0; i < numOfFiles; i++) {
//            File reSampledData;
//            PrintWriter reSampledDataStream;
//            try {
//                reSampledData = new File(outDir + delimiter + getFileName(files[i]) + "_resampled.txt");
//                reSampledDataStream = new PrintWriter(new FileOutputStream(reSampledData));
//            } catch (FileNotFoundException e) {
//                IJ.error(e.toString());
//                return;
//            }
//            double oldSampleRate = data[i][0].get(1) - data[i][0].get(0);
//            int size = data[i][0].size();
//            int newLength = size;
//            double newData[][];
//            double dilation = 1.0;
//            if (normLength) {
//                double measY = lastNormLength(data, i, DataFileAverager.LENGTH_INDEX, window);
//                double d = data[i][DataFileAverager.timeIndex].get(size - 1);
//                dilation = ((measY - coeffs[1]) / coeffs[0]) / d;
//                newData = new double[numParams][newLength];
//            } else {
//                newLength = (int) Math.round(size * oldSampleRate / newSampleRate);
//                newData = new double[numParams][newLength];
//            }
//            for (int j = 0; j < numParams; j++) {
//                double currentData[] = new double[size];
//                for (int k = 0; k < size; k++) {
//                    if (normLength) {
//                        if (j == DataFileAverager.timeIndex) {
//                            newData[j][k] = data[i][j].get(k) * dilation;
//                        } else {
//                            newData[j][k] = data[i][j].get(k);
//                        }
//                    } else {
//                        currentData[k] = data[i][j].get(k);
//                    }
//                }
//                if (!normLength) {
//                    newData[j] = DSPProcessor.upScale(currentData, newLength, false);
//                }
//            }
//            reSampledDataStream.println(filenames[i]);
//            reSampledDataStream.println(headings);
//            for (int l = 0; l < newLength; l++) {
//                for (int m = 0; m < numParams; m++) {
//                    reSampledDataStream.print(newData[m][l] + "\t");
//                }
//                reSampledDataStream.println();
//            }
//            reSampledDataStream.close();
//        }
    }

    boolean showDialog() {
        GenericDialog gd = new GenericDialog(TITLE);
        gd.addNumericField("Header Size", headerSize, 0);
        gd.addNumericField("Number of Parameters", numParams, 0);
        gd.addNumericField("Sample Rate", newSampleRate, 3);
        gd.showDialog();
        if (gd.wasCanceled()) {
            return false;
        }
        headerSize = (int) Math.round(gd.getNextNumber());
        numParams = (int) Math.round(gd.getNextNumber());
        newSampleRate = gd.getNextNumber();
        return true;
    }

    String getFileName(File file) {
        String filename = file.getName();
        Scanner s = new Scanner(filename).useDelimiter("[.]");
        return s.next();
    }

    double lastNormLength(ArrayList<Double> data[][], int m, int n, int w) {
        ArrayList<Double> measures = data[m][n];
        double sum = 0.0;
        int length = measures.size();
        if (w >= length) {
            w = length - 1;
        }
        for (int i = 1; i <= w; i++) {
            sum += measures.get(length - i);
        }
        return sum / w;
    }
}
