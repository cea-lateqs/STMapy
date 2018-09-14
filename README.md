# Scampy
## Presentation
**Scampy** (SCAnning tunneling Microscopy analysis in PYthon) is a Graphical User Interface (GUI) to analyse CITS recorded under Nanonis (.3ds) or MATRIX (.asc).
It is written in Python using PyQt for the GUI.

## Requirements
Using **Scampy ** requires the following to be installed :
* Python 3 (tested under Python 3.6) with the following packages
    * Numpy
    * Scipy
    * Matplotlib 2.0
* PyQt 5

## Using Scampy
### Starting
To start Scampy, run the main.py file:
```python3
python3 main.py
```

### Loading a CITS
#### File selection
To load a CITS, click on the _Open CITS_ button on top-left corner. A window will appear, prompting to select a CITS of supported format (either *.3ds or *.asc). The filenames can be filtered according to  the format by selecting _3D binary file_ (*.3ds) or _Ascii files_ (*.asc).

#### Topography
Once the CITS was selected, **Scampy** will load the spectroscopic data and will attempt to read the topography to plot it in a seperate window. 

This always succeeds for *.3ds as it plots the topography contained in the file. For *.asc however, it will search for a file 'Topo.txt' in the same folder of the selected file. This 'Topo.txt' can be created by using the _Export to TXT_ method of [Gwyddion](http://gwyddion.net/). If no topographic file is found, no topography will be plotted.

<span class="warning">No checks are done to see if 'Topo.txt' corresponds to the loaded CITS. Always check that the topography file was taken at the same location as the CITS.</span>

#### End of loading
If the loading succeed, a 2D plot should appear on the bottom-left widget. If not, check the console for error messages and report to Contact Information.

### List of commands
#### Clicks on 2D plot
* Left-click: plots the spectrum at this location
* Right-click: plot a spectrum averaged around this location
* Dragging left-click: draws a line to plot a cut of the data along it in a separate window
* Dragging right-click: draws a rectangle to average the data over it and plot the resulting spectrum
#### Clicks on spectrum window
* Any click: displays the 2D plot at the clicked voltage value

### List of inputs
* _Open CITS_ (button): allows to select a CITS to load. If several are selected, they will be averaged.
* _Draw topo_ (button): redraws the topo without any shape on it.
* _Whole length cut_ (button): Do a cut of the data along the great diagonal.
#### CITS parameters
* _Displayed channel_ (dropdown): allows to select the working channel. It will be displayed in the 2D plot and all operations will be used on this channel.
* _Normalize current channel_ (button): normalises the current channel. The result of this normalisation will be added as a new channel named 'Normalised' followed by the channel name.
* _V/Z index_ (spinbox): changes the index of the voltage/altitude at which the 2D plot is displayed. The corresponding value is shown in the header of the 2D plot.
#### Averaging
##### Multiple spectra
* _Select spectra with left-click_ (checkbox): allows to select several spectra by left-clicking to average them. Unchecking the box will display the result of the averaging in the bottom-right widget.
* _Number of pixels for right-click average_ (spinbox): Right-clicking on the 2D plot at (X, Y) triggers an average over X - N, X + N, Y - N, Y + N. This box allows to change the value of N (default : 2).
* _Average with respect to value_ (checkbox): triggers the display of the averaging widget. This widget allows to plot two spectra: one resulting from the average of all spectra that have a value at the current voltage below a certain value (Minimum value for the below averaging) and the other resulting from the average of all spectra that have a value at the current voltage above a certain value (Maximum value for the below averaging). This operation can be started by clicking _Start the averaging_ (may be lengthy).
##### Whole CITS
* _Plot average spectrum_ (button): Plots the average of all spectra of the current channel.
* _Average CITS_ (button): Average the neighbouring spectra of the CITS along the X direction and replaces the CITS with the result. The number of spectrum to average together can be changed in the spinbox _Number of spectra to average along X direction_. Useful for line spectroscopies.

### Explore the data
Once the CITS is loaded, the data can be explored in several ways
#### CITS
A CITS usually contains different quantities named **channels** recorded at different voltages/altitudes (V/Z) and different positions. The bottom-left widget contains a 2D plot displaying the current channel with respect to the position at the selected V/Z index. The current channel and the selected index can be changed in the group _CITS parameters_.

The _Normalize current channel_ button will normalise the current channel. The result of this normalisation will be added as a new channel named 'Normalised' followed by the channel name.

#### Spectra
Rather than to look the data at a fixed voltage with respect to the position, one may want to look the evolution with respect to the voltage at a given position. Such data is called in **Scampy** spectra and are displayed in the bottom-right widget.

The simplest way to display a spectrum at a given position is to do a left-click on the 2D plot. A dot of given color will appear on the 2D plot (and on the topography if present) and a spectrum of the same color will be plotted on the right. This spectrum has a label of the form `[X position (in px), Y position (in px)] - channel name`.