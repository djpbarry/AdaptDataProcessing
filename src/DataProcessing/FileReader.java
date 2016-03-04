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

import ij.IJ;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;

/**
 *
 * @author barry05
 */
public class FileReader {

    private int numOfFiles;
    private int headerSize;
//    private String headings;
    private String filenames[];
    private ArrayList<String> paramNames;
    private Charset charSet;

    public FileReader() {
    }

    public FileReader(int numOfFiles, int headerSize, Charset charSet) {
        this.headerSize = headerSize;
        this.numOfFiles = numOfFiles;
        this.filenames = new String[numOfFiles];
        this.paramNames = new ArrayList();
        this.charSet = charSet;
    }

    public void readData(ArrayList<Double>[][] data, File[] files, String delimiter) throws FileNotFoundException, IOException {
        for (int i = 0; i < numOfFiles; i++) {
            BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(files[i]), charSet));
            filenames[i] = br.readLine();
            ArrayList<String> thisParams = getParamsArray(br.readLine(), delimiter);
            if (thisParams != null) {
                int numThisParams = thisParams.size();
                String line = br.readLine();
                while (line != null) {
                    Scanner scan = new Scanner(line).useDelimiter(delimiter + "\\s*");
                    for (int k = 0; k < numThisParams; k++) {
                        int j = getParamIndex(thisParams.get(k), paramNames);
                        if (data[i][j] == null) {
                            data[i][j] = new ArrayList();
                        }
                        double val = scan.nextDouble();
                        data[i][j].add(val);
                    }
                    line = br.readLine();
                }
            }
            br.close();
        }
        Arrays.sort(filenames);
    }

    public String[] getFilenames() {
        return filenames;
    }

    public String[] getParamsArray() {
        int numParams = paramNames.size();
        String[] paramArray = new String[numParams];
        for (int i = 0; i < numParams; i++) {
            paramArray[i] = paramNames.get(i);
        }
        return paramArray;
    }

    public ArrayList<String> getParamsArray(String headings, String delimiter) {
        if (headings != null) {
            Scanner scanner = new Scanner(headings).useDelimiter(delimiter);
            ArrayList<String> headingsArray = new ArrayList();
            while (scanner.hasNext()) {
                headingsArray.add(scanner.next().trim());
            }
            return headingsArray;
        } else {
            return null;
        }
    }

    public String getParamString() {
        int numParams = paramNames.size();
        String paramString = "";
        for (int i = 0; i < numParams; i++) {
            paramString = paramString.concat(paramNames.get(i) + "\t");
        }
        return paramString;
    }

    public double readParam(File paramFile, int headerSize, String paramName) throws FileNotFoundException, IOException {
        BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(paramFile)));
        for (int j = 0; j < headerSize; j++) {
            br.readLine();
        }
        String line = br.readLine();
        while (line != null) {
            Scanner scan = new Scanner(line);
            if (paramName.equalsIgnoreCase(scan.next())) {
                double val;
                try {
                    val = scan.nextDouble();
                } catch (Exception e) {
                    scan.next();
                    val = scan.nextDouble();
                }
                return val;
            }
            line = br.readLine();
        }
        br.close();
        return Double.NaN;
    }

    public void getParamList(File[] files, String delimiter) throws FileNotFoundException, IOException {
        for (int i = 0; i < numOfFiles; i++) {
            BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(files[i]), charSet));
            for (int j = 0; j < headerSize; j++) {
                br.readLine();
            }
            String paramLine = br.readLine();
            if (paramLine != null) {
                Scanner scan = new Scanner(paramLine).useDelimiter(delimiter);
                while (scan.hasNext()) {
                    String thisParam = scan.next().trim();
                    if (!paramNames.contains(thisParam)) {
                        paramNames.add(thisParam);
                    }
                }
            }
            br.close();
        }
    }

    public int getNumParams() {
        return paramNames.size();
    }

    int getParamIndex(String param, ArrayList<String> refList) {
        return refList.indexOf(param);
    }

}
