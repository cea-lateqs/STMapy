## 0.4.0

### New features

- Added the possibility to open new channels, including calculated ones : FFT , Feenstra renormalisation, derivative
- Added display features : new mplstyle, cits aspect checkbox, fft cuts calculations
- Added statistics analysis function.
- Gifs can now be exported
- Added citation file

### Fixes

- Fixed minor issues

## 0.3.0

### New features

- Added plane fitting of the topography. The type of fitting can be chosen in the config.

### Fixes

- Fixed line fitting of the topography.

## 0.2.1

### New features

- Line fitting of the topography can be enabled/disabled in the config.

### Fixes

- Fixed crash when plotting the topo of CITS of unequal aspect ratio.
- Fixed waterfall plots not displaying correctly.

## 0.2.0

### New features

- Config is now a `JSON` file located in `scampy` folder
- Topography colormap can now be set in the configuration

### Improvements

- Lots of NumPy array calculations were rewritten by removing Python loops (e.g. using vectorization).
- The CITS map and the voltage line are now updated rather than drawn every time there is a change.
- Reading of topography files is now done with `np.genfromtxt`.
- Cut plots and the handling of 3ds Z-CITS were improved.

### Fixes

- Inconsistent ordering of pixels in cut plots was fixed.

### Clean-up

- Internals: variable renamed, functions moved in modules (`plotting`, `processing` and `reads`), tests written.
- References to Matlab file opening were removed.

## 0.1.0

Initial version.
