# Tutorial for GRUTinizer

This text will describe the basic process of using GRUTinizer for analyzing
data from nuclear physics experiments at the NSCL, by providing an example of analyzing
data from a nucleon knockout experiment utilizing the S800 with GRETINA.


## Installation

* Install root6, either 6.06.08 or 6.08.02. I recommend you install it with all
  optional components enabled; however, you need at least the python and
  MathMore components installed.
* Clone and install the root6 branch of GRUTinizer 

## Benchmark GRUTinizer against previous results

In case substantial changes are made to GRUTinizer in the future, it's worth
comparing the resulting spectrum. To this end (assuming you have access to the
NSCL computers), look in:

/mnt/analysis/pecan-gade/for_grutinizer_benchmark

Contained here is raw data from Run 96 in E15127, along with  histogramming
code and root output (both the tree and histograms). Move
benchmark_grutinizer.cxx to $GRUTSYS/histos and run make again. Then we can use
GRUTinizer to re-histogram data from the raw file and therefore compare the
results. Do:

  * grutinizer -qH raw_data/BenchmarkData_raw.dat $GRUTSYS/lib/libbenchmark_grutinizer.so  -o /dev/null

This will create histograms based on the definitions in
benchmark_grutinizer.cxx from the data in BenchmarkData_raw.dat. The /dev/null
call stops the tree from being created and therefore speeds up the
histogramming process. Ensure that the number of counts in each histogram is
comparable to those in root_files/BenchmarkData_hists.root. 

## Unreacted Beam Runs
Now that you have GRUTinizer installed, let's walk through a simple analysis
process. 

First, start with data from an unreacted run for the setting of interest.
Unreacted runs center the unreacted incoming beam in the focal plane, as the
spectrograph's magnetic rigidity is determined based on the incoming beam
energy with the energy loss due to traversing the target subtracted. The
purpose of looking at this run first is it will allow us to determine an
initial shape for the incoming beam cut.

To sort the data from the raw file into a ROOT tree, use:

  * grutinizer -q /path/to/UnreactedRun.dat -o unreactedrun.root

Now let's explore the root tree a bit.

Enter the grutinizer intepreter using:

  * grutinizer unreactedrun.root

There are a few things we will want to look at in the unreacted run.
Specifically, we want to look at:

1. the outgoing particle identification (PID) plot for projectiles entering the focal plane, and
2. the incoming PID plot with gates on the outgoing PID 

If you are using the MesyTec TDCs, we need to first determine which hit for
each event is valid. To read more about the time-of-flight measurements from
the S800, see the [timing page on the S800 wiki](https://wikihost.nscl.msu.edu/S800Doc/doku.php?id=timing).

To determine the "good" time-of-flight value for both the time of flight from
the OBJ Scintillator and from the XFP scintillator, let's look at the OBJ-E1
and XFP-E1 time differences in the interpreter. *Note that histogramming from
*.cxx is MUCH prefered to the interpreter, as there can be unexpected behavior
when using the interpreter*. To use the interpreter, we can first look
at the values directly:

  * gChain->Scan("GetRawOBJ_MESY()-GetRawE1_MESY()")
  * gChain->Scan("GetRawXF_MESY()-GetRawE1_MESY()")

*gChain* is a useful global variable that includes all the different root trees
if one loads more than one file. Now to look at these visually instead, use:

  * gChain->Draw("GetRawOBJ_MESY()-GetRawE1_MESY()>>h(1000,-2000,-1000)")
  * gChain->Draw("GetRawXF_MESY()-GetRawE1_MESY()>>h2(2000, 2000, 4000)")

Based on these distributions, we can choose the "good" value for the
time-of-flight difference. We choose a value roughly in the middle of the
majority of counts in the distributions.  We then can set two **gValues** to
take these into account:

  * TARGET_MTOF_OBJE1 = *val1*
  * TARGET_MTOF_XFPE1 = *val2*

We can add these to a ".val" file. By passing this file in when sorting the
trees into histograms, we can store calibration info which will be used by the
functions in GRUTinizer. To add these to a gValue file, create a new file
called myname.val (where myname can be whatever you want) and add the following
lines with *val1*, *val2* replaced with whatever number you use:

TARGET_MTOF_OBJE1 {
  value = *val1*
}
TARGET_MTOF_XFPE1 {
  value = *val2*
}

If we only want to use this gValue for this session, we can call the following
in the interpreter:

  * GValue::AddValue(new GValue("TARGET_MTOF_OBJE1", *val1*)
  * GValue::AddValue(new GValue("TARGET_MTOF_XFPE1", *val2*)

Now, with the GValues loaded either using the AddValue function or re-opening
grutinizer with:

  * grutinizer unreactedrun.root myname.val

We can now look at the particle identification plot.  Let's draw the outgoing
PID. Note that the time-of-flight has not been corrected for angle and position
in the focal plane yet, which means the outgoing PID is not ideal at the
moment, but is ok for our purposes. To create the outgoing PID, use the typical
root syntax:

gChain->Draw("GetIonChamber().GetAve():GetMTof().GetCorrelatedObjE1()>>out(1000,-2000,-1000,4000,0,4000)", "", "colz")

This will create a 2D histogram with the ion chamber average on the y-axis and
the correlated, uncorrected obj-e1 time-of-flight on the x-axis. Note that in
root syntax the format of the draw command is: 

  * gChain->Draw("yvar:xvar>>name(bins_x,xlow,xhigh,bins_y, ylow,yhigh)","cutnames", "drawoptions", numentries)  

Create a cut on this spectrum by clicking twice to create a rough version of
the cut around the strongest component in the outgoing PID, which is the Cu-71
detected in the focal plane. Click "g" to turn the cross hatched lines to a
GCutG.

After turning the lines into a cut, we can adjust the cut to better surround
the component of the PID in which we're interested. This is done by clicking on
specific points around the cut and adjusting them. You can see all points that
can be manipulated by clicking the cut elsewhere, which also allows you to move
the entire cut.

Right click the contour after it's satisfactory, and click the SaveTo option.

After clicking SaveTo, a dialog box appears where you choose the cutname, 
filename to save it to, tagname (which allows us to use it easier in our
histogramming from the files later), and the option for opening the file, which
allows us to update or recreate the cut file (don't use recreate if you don't want
to lose everything that was in the file!). The *.cuts extension is important
for grutinizer to recognize the file when histogramming. (Note that the default
filename has a *.root ending!). 

Let's reopen GRUTinizer with the cut file (which we assume to be named
unreacted_cuts.cuts in this case) open.

  * grutinizer unreactedrun.root myname.val unreacted_cuts.cuts

Now type into the interpreter the name you gave the cut in the *.cuts file to
load it into memory.

Now that we have this outgoing PID cut, we can use it to look at the incoming
PID restricted to events with the beam component of interest in the focal
plane, and therefore determine an initial version of the incoming beam cut.

gChain->Draw("GetMTof().GetCorrelatedXfpE1():GetMTof().GetCorrelatedObjE1()>>h(4000,-4000,0,4000,0,4000)", "out_particle","colz")

We create our initial cut on this incoming beam PID, assuming it's called
"in_particle". Now that we have this cut, we can start looking at our reaction
data. 

## Initial corrections for reactions data

Now, we have an initial incoming beam cut. The following steps still need to be
completed:

  1. Use the CRDC masks to determine the calibration parameters for the 
     positions in the focal plane.
  2. Using the incoming beam cut and the calibrated CRDC values, the incoming
     beam cut needs to be refined for the reaction setting.
  3. After refining the incoming beam cut, need to do a rough gate on a few
     isotopes in the isotopic chain of interest, and these will be used to
     determine the time-of-flght corrections for the angle and position in the focal plane.
  4. After correcting the time-of-flight, we now refine our outgoing beam cut
     by utilizing the corrected time of flights. We now can further correct our
     time of flight by considering the correlation between the xfp-obj time
     difference and the corrected obj time. This will require centering the
     xfp-obj time by a shift, and then adding it with some scaling factor to
     the corrected time of flight.

Following these 4 steps, we will have completed our particle identification
process.

### Determining the CRDC calibration parameters

The CRDC mask is a mask with holes at known positions that allows us to
calibrate our CRDC X and Y values. See [the page on the GRUTinizer
wiki.](https://github.com/pcbend/GRUTinizer/wiki/CRDC-Mask-Calibration).  

The CRDC X slope and offset are always the same, and would be represented in
the GValue file as: 
 
CRDC1_X_SLOPE {
  value: 2.54
}

CRDC1_X_OFFSET {
  value: -281.94
}

(and likewise for CRDC2). Note that units of "mm" are used for the CRDC X/Y
values in GRUTinizer. For the CRDC mask positions, see the below image.  Use
the known Y values and those measured with the mask to determine the Y_SLOPE.
The Y_OFFSET will be determined by choosing it to center the Y distribution 
from the first reaction run after the mask run. The gValues for these are in the
same format as the ones listed above, with X replaced with Y. 

Note that it is also important to look at the CRDC Y value as a function of run
number and to correct the SLOPE (NOT Offset!) of each run individually to align
the centroids of the CRDC Y distributions to 0 for both CRDCs.

To plot the data for the mask run in the interpreter, use (assuming you've
already entered the X slope/offset):

gChain->Draw("GetCrdc(0).GetNonDispersiveY():GetCrdc(0).GetDispersiveX()>>h(600,-300,300,1000,0,1000)", "","colz")

where GetCrdc(0) can be replaced with GetCrdc(1) for the CRDC 2 Mask run.


### Refining the incoming beam cut
Now that the CRDCs are calibrated, we will refine the incoming beam cut in the 
first reaction setting run after the mask. 

gChain->Draw("GetMTof().GetCorrelatedXfpE1():GetMTof().GetCorrelatedObjE1()>>h(200,-1600,-1400,350,2900,3250)", "","colz")

Note that to make this run much faster, you can put a limited number of entries
to observe by rewriting the previous line as: 

gChain->Draw("GetMTof().GetCorrelatedXfpE1():GetMTof().GetCorrelatedObjE1()>>h(200,-1600,-1400,350,2900,3250)", "","colz", 500000)

where 500000 is the number of entries to loop through before drawing.

You can draw the incoming beam cut, and then extend it diagonally. After
extending it, remember to call SaveTo to save it, possibly to a different file.

### Initial rough gate on outgoing PID and first set of time-of-flight corrections
Utilizing the new incoming beam cut, we draw the outgoing PID and place a rough
gate on the region containing the isotope of interest and it's neighbors.

gChain->Draw("GetIonChamber().GetAve():GetMTof().GetCorrelatedObjE1()>>h(1000,-2000,-1000,4000,0,4000)", "in_particle", "colz", 500000)

Using this gate, assuming it's called out_particles, we can look at the following spectrum:

gChain->Draw("GetCrdc(0).GetDispersiveX():GetMTof().GetCorrelatedObjE1()>>h(1000,-2000,-1000,600,-300,300)", "in_particle&out_particles", "colz", 100000)
gChain->Draw("GetAFP():GetMTof().GetCorrelatedObjE1()>>h2(1000,-2000,-1000,1000,-1,1)", "in_particle&out_particles", "colz", 100000)

You can clearly see that there is a time-of-flight dependence on CRDC 1 X and
AFP. We correct for it in the following way: 

gChain->Draw("GetCrdc(0).GetDispersiveX():GetMTofObjE1(afp_cor, xfp_cor)>>h(1000,-2000,-1000,600,-300,300)", "in_particle&out_particles", "colz", 100000)
gChain->Draw("GetAFP():GetMTofObjE1(afp_cor, xfp_cor)>>h2(1000,-2000,-1000,1000,-1,1)", "in_particle&out_particles", "colz", 100000)

In general, we can vary the coefficients of the AFP/CRDC1X to attempt to remove
the correlation. Typically,  afp_cor ~ 1500 and xfp_cor ~  0.2 respectively,
where the gValues are: 

OBJ_MTOF_CORR_AFP	{
value:	1500
}

OBJ_MTOF_CORR_XFP	{
value:	0.2
}

These gValues will automatically be used when calling GetMTofObjE1() without
arguments. Vary these parameters for each spectrum until the correlation in
each is minized.

### Correction for XFP-OBJ time difference 
Two more useful gValues are TOFXFP_OBJ_SHIFT and TOFXFP_OBJ_CORR. To observe
the correlations these correct, plot the XFP-OBJ tof difference against the corrected TOF value.

gChain->Draw("GetMTofXfpE1(0,0)-GetMTofObjE1(0,0):GetMTofObjE1(afp_cor, xfp_cor)>>h5(1000,-2000,-1000,1000,4000,5000)", "in_particle&out_particle", "colz", 100000)

Based on the XFP-OBJ distribution, choose TOFXFP_OBJ_SHIFT as the centroid
value of xfp-obj to center the distribution. Then, using this
information, we vary the correction value TOFXFP_OBJ_CORR:

gChain->Draw("GetMTofXfpE1(0,0)-GetMTofObjE1(0,0):GetMTofObjE1(afp_cor,xfp_cor, tofxfp_obj_shift, tofxfp_obj_corr)>>h5(1000,-2000,-1000,1000,4000,5000)", "in_particle&out_particle", "colz", 100000)

The value for TOFXFP_OBJ_CORR has been seen to be typically small and negative
(~ -0.6 in two experiments recently).

### Using Histos Files
The previous corrections were rather slow, because a considerable amount of
statistics can be necessary to see enough of a component of interest to be able
to do the corrections. Also, you may want to make many histograms without
sorting through the full tree everytime. To handle this, we make histogramming
libraries in $GRUTSYS/histos, with file extensions *.cxx which will
automatically be compiled when calling make for GRUTinizer. Then these can be
passed to grutinizer based on the built library in $GRUTSYS/lib/libhistos.so
(where histos.so is replaced with the name of the histogramming library).
