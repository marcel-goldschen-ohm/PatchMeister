# PatchMeister
Electrophysiology patch-clamp time series analysis (NOT acquisition).

PatchMeister is...
* a free and open source GUI for analyis of patch-clamp electrophysiology data.
* for analysis ONLY, it does NOT have the capability to aquire data from a digitizer/amplifier. 
* currently a basic but highly useful analysis tool for analyzing series of current vs. time sweeps.
* currently primarily for ensemble current measurements, NOT for single-channel analysis.
* easily extended for custom analyses (it's all just [MATLAB](http://www.mathworks.com/products/matlab) code).

Copyright 2019 by Marcel Goldschen-Ohm <goldschen-ohm@utexas.edu>

### Install
Everything is in patch_meister.m

Tested with MATLAB R2018b.

### Capabilites
* Import data.
    * HEKA PatchMaster MATLAB export (assumes all traces belong to the same patch and experiment).
    * :construction: HEKA PatchMaster datafile
    * :construction: Molecular Devices pCLAMP datafile
    * :construction: AxoGraph or AxographX datafile
* Organize a series of current vs. time sweeps by splitting them into groups.
    * Group interleaved sweeps.
    * :construction: Group blocks of sweeps.
    * :construction: Manually group arbitrary sweeps.
    * Optionally show average of sweeps in each group.
    * :construction: Optionally overlay all groups in a single plot.
* Visually step through a group of sweeps.
    * Optionally show only the selected sweep(s).
    * Optionally show all sweeps with selected sweep(s) highlighted.
* Mask selected sweeps to remove them from analysis (e.g. average).
   * Masked sweeps are flagged, but NOT deleted.
   * Optionally show or hide masked sweeps.
* Basic manipulation of sweep data.
    * Sweep manipulations are stored as offsets and scale factors (original data is NOT directly altered).
    * Manipulations can be applied to one or any number of sweeps simultaneously by selecting appropriate data regions in each of the desired sweeps using MATLAB's brush tool. How the selected data points are used depends on the manipulation.
    * **Baseline one region (flat)** - Baseline Y data to mean of selected data points.
    * **Baseline two region (sloping)** - Baseline Y data by subtracting line through (x,y) mean of each of two regions of contiguous selected data points.
    * **Normalize** - Normalize Y data (absolute value, positive or negative) to selected data points.
    * **Align to Onset** - Set time zero (X offset) to first data point following selected data points that is some number of standard deviations of the selected data points away (absolute value, positive or negative) from the mean of the selected data points. Note that this may shift traces such that the data points at the beginning or end of the trace no longer align in time with data points at the beginning or end of other traces. The group average will automatically handle these offsets by trimming each trace in the group to the set of overlapping points prior to averaging.
    * **Manually specify X and Y offsets and Y scale factor** - :construction:
