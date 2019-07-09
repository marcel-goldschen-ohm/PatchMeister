# PatchMeister
Electrophysiology patch-clamp time series analysis (NOT acquisition).

PatchMeister is...
* A free and open source graphical user interface for analyis of patch-clamp electrophysiology data or similar time series measurements.
* For analysis ONLY, it does NOT have the capability to acquire data from a digitizer/amplifier. 
* Handles multiple channels (e.g. current, voltage, etc.) and arbitrary grouping of sweeps.
* Currently primarily for analysis of ensemble measurements as opposed to, for example, single-channel analysis.
* Easily extended for custom analyses (it's all just [MATLAB](http://www.mathworks.com/products/matlab) code).

Copyright 2019 by Marcel Goldschen-Ohm <goldschen-ohm@utexas.edu>

### Install
Everything is in patchmeister.m

Baseline spline fitting requires the [splinefit](https://www.mathworks.com/matlabcentral/fileexchange/71225-splinefit) package. Easiest option is to install it via MATLAB's Add-On Explorer.

Tested with MATLAB R2019a.

### Capabilites
* *Import data.*
    * HEKA PatchMaster datafile (*NOT tested for multiple channels, otherwise working*)
    * Molecular Devices pCLAMP datafile (*NOT tested*)
    * :construction: AxoGraph or AxographX datafile
* *Organize a series of current/voltage/etc vs. time sweeps by splitting them into groups.*
    * View one or more selected channels (e.g. current, voltage, etc.)
    * Group interleaved sweeps for each channel.
    * Group blocks of sweeps for each channel.
    * Manually group arbitrary sweeps for each channel.
    * Optionally show average sweep for each group.
* *Visually step through a group of sweeps.*
    * Buttons to step forward/backward through a group of sweeps.
    * Optionally show any arbitrary range of selected sweep(s).
* *Mask sweeps to remove them from subsequent analyses.*
   * Masked sweeps are flagged, but NOT deleted.
   * Optionally show or hide masked sweeps.
* *Basic manipulation of sweep data.*
    * Sweep manipulations are stored as offsets and scale factors (**original data is NOT directly altered**).
    * **You can always optionally view the original raw data.**
    * Manipulations can be applied to one or any number of sweeps simultaneously by selecting appropriate data regions in each of the desired sweeps using MATLAB's brush tool. How the selected data points are used depends on the manipulation.
    * **Baseline Flat** - Baseline Y data to mean of selected data points.
    * **Baseline Sloping Two Region** - Baseline Y data by subtracting line through (x,y) mean of each of two regions of contiguous selected data points.
    * **Baseline Spline** - Baseline Y data by subtracting a spline fit through the selected data points.
    * **Normalize Positive/Negative/AbsPeak** - Normalize Y data to max, min or absmax of selected data points.
    * **Align to Onset** - Set time zero (X offset) to first data point following selected data points that is some number of standard deviations of the selected data points away (absolute value, positive or negative) from the mean of the selected data points. Note that this may shift traces such that the data points at the beginning or end of the trace no longer align in time with data points at the beginning or end of other traces. The group average will automatically handle these offsets by trimming each trace in the group to the set of overlapping points prior to averaging.
    * **Manually specify X and Y offsets and Y scale factor** - Set time zero (X offset), baseline (Y offset) and Y scale.
    * **Mask/Zero/Interpolate** - Mask, zero or interpolate selected segments of data points.
