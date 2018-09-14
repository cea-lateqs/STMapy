# Scampy (SCAnning tunneling Microscopy analysis in PYthon)
## Presentation
**Scampy** (previously STM_Data_Analysis) is a Graphical User Interface (GUI) to analyse CITS recorded under Nanonis (.3ds) or MATRIX (.asc).
It is written in Python using PyQt for the GUI.

## Requirements
Using Scampy requires the following to be installed :
* Python 3 (tested under Python 3.6) with the following packages
    * Numpy (tested under 1.15)
    * Scipy (tested under 1.1)
    * Matplotlib, at least version 2.0 (tested under 2.2.3)
* PyQt 5

If these are not satisfied, I advise to use the [Anaconda distribution](https://www.anaconda.com/distribution/) to install them.

## Using Scampy
### Starting
To start Scampy, run the main.py file:
```python3
python3 main.py
```
It should display the following interface.

#### File selection
To load a CITS, click on the **Open CITS** button on top-left corner. A window will appear, prompting to select a CITS of supported format (either .3ds or .asc). The filenames can be filtered according to  the format by selecting _3D binary file_ (.3ds) or _Ascii files_ (.asc).

#### Topography
Once the CITS was selected, Scampy will load the spectroscopic data and will attempt to read the topography to plot it in a seperate window. 

This always succeeds for .3ds as it plots the topography contained in the file. For .asc however, it will search for a file 'Topo.txt' in the same folder of the selected file. This 'Topo.txt' can be created by using the _Export to TXT_ method of [Gwyddion](http://gwyddion.net/). If no topographic file is found, no topography will be plotted.

<span style='color: red;'>No checks are done to see if 'Topo.txt' corresponds to the loaded CITS. Always check that the topography file was taken at the same location as the CITS.</span>

#### End of loading
If the loading succeed, a 2D plot should appear on the bottom-left widget. If not, check the console for error messages and report them (see **Further information** at the end).

### List of commands
#### Clicks on 2D plot
* **Left-click**: plots the spectrum at this location
* **Right-click**: plot a spectrum averaged around this location
* **Dragging left-click**: draws a line to plot a cut of the data along it in a separate window
* **Dragging right-click**: draws a rectangle to average the data over it and plot the resulting spectrum
#### Clicks on spectrum window
* Any click: displays the 2D plot at the clicked voltage value

#### Basic operations
* **Open CITS** (_button_): allows to select a CITS to load. If several are selected, they will be averaged.
* **Draw topo** (_button_): redraws the topo without any shape on it.
#### CITS parameters
* **Displayed channel** (_dropdown_): allows to select the working channel. It will be displayed in the 2D plot and all operations will be used on this channel.
* **Normalize current channel** (_button_): normalises the current channel. The result of this normalisation will be added as a new channel named 'Normalised' followed by the channel name.
* **V/Z index** (_spinbox_): changes the index of the voltage/altitude at which the 2D plot is displayed. The corresponding value is shown in the header of the 2D plot.
#### Averaging
##### Multiple spectra
* **Select spectra with left-click** (_checkbox_): allows to select several spectra by left-clicking to average them. Unchecking the box will display the result of the averaging in the bottom-right widget.
* **Number of pixels for right-click average** (_spinbox_): Right-clicking on the 2D plot at (X, Y) triggers an average over X - N, X + N, Y - N, Y + N. This box allows to change the value of N (default : 2).
##### Whole CITS
* **Plot average spectrum** (_button_): plots the average of all spectra of the current channel.
* **Average CITS** (_button_): averages the neighbouring spectra of the CITS along the X direction and replaces the CITS with the result. The number of spectrum to average together can be changed in the spinbox **Number of spectra to average along X direction**. Useful for line spectroscopies.
* **Average with respect to value** (_checkbox_): triggers the display of the averaging widget. This widget allows to plot two spectra: one resulting from the average of all spectra that have a value at the current voltage below a certain value (**Minimum value for the below averaging**) and the other resulting from the average of all spectra that have a value at the current voltage above a certain value (**Maximum value for the above averaging**). This operation can be started by clicking **Start the averaging** (may be lengthy).
#### Cuts
* **Waterfall** or **2D plot** (_radio buttons_): allows to choose the type of representation for cuts.
* **View selected spectra** (_checkbox_): allows to underline the pixels in the 2D plot that are chosen when doing a cut (debugging purposes).
* **Whole length cut** (_button_): does a cut of the 2D plot along the great diagonal.
#### Spectra plot
* **Clear spectra** (_button_): clears the spectra window.
* **Shift plot along X** (_text_): shifts the plotted spectra along the voltage axis by the given value. If 'topo' is given, the spectra will be shifted by the Z value at this position (unstable !).
* **Shift plot along Y** (_text_): shifts the plotted spectra along the channel axis by the given value.
* **Plot log of the spectrum** (_checkbox_): plots the logarithm of the spectrum instead of the spectrum.
* **Plot derivative** (_checkbox_): plots the derivative of the spectrum in addition to the spectrum. The derivative is computed by a [Savitzy-Golay filter](https://docs.scipy.org/doc/scipy-0.16.1/reference/generated/scipy.signal.savgol_filter.html) with a window length that can be changed in the spinbox **Window length**.
* **Fit spectra** (_button_): tries to fit the last spectra plotted by a linear function. The range on which the fit must be computed can be set by checking **Use a custom range** and setting the lower and upper limits.
#### Display parameters
* **Scale in Volts** (_checkbox_): sets voltage axis in Volts for cuts and spectra. Uses indexes otherwise.
* **Scale in metric units** (_checkbox_): sets the X and Y axis in metric units for the topo and the cuts (NOT for the 2D plot). Uses pixels otherwise.
* **Voltage index guideline** (_checkbox_): shows a dashed line in the bottom-right widget at the voltage selected in **V/Z index**.
* **Deactivate legend** (_checkbox_): removes the legend when plotting spectra. Can be useful when many spectra are plotted.
* **Colorbar settings** (_checkbox_): opens the colorbar widget where the colorbar to use can be changed. Custom limits can also be forced on the colormap by checking **Use custom limits**. <span style="color: red;">In this case, you must set both the lower AND the upper limit.</span>

## Further information
The code is available on the [Git repo](https://gitlab.com/lateqs/STM_Data_Analysis). Bugs can be reported on the repository or to directly to [me](mailto:loic.huder@gmail.com).
