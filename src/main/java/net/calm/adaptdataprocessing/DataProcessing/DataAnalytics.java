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

/**
 *
 * @author barry05
 */
public class DataAnalytics {

    private int peakIndex;
    private int zeroCrossing;
    private double negMean;
    private double posMean;
    private double maxVal;
    private double minVal;

    public double getMaxVal() {
        return maxVal;
    }

    public void setMaxVal(double maxVal) {
        this.maxVal = maxVal;
    }

    public double getMinVal() {
        return minVal;
    }

    public void setMinVal(double minVal) {
        this.minVal = minVal;
    }

    public double getNegMean() {
        return negMean;
    }

    public void setNegMean(double negMean) {
        this.negMean = negMean;
    }

    public int getPeakIndex() {
        return peakIndex;
    }

    public void setPeakIndex(int peakIndex) {
        this.peakIndex = peakIndex;
    }

    public double getPosMean() {
        return posMean;
    }

    public void setPosMean(double posMean) {
        this.posMean = posMean;
    }

    public int getZeroCrossing() {
        return zeroCrossing;
    }

    public void setZeroCrossing(int zeroCrossing) {
        this.zeroCrossing = zeroCrossing;
    }
}
