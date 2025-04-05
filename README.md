# xDIC: Fingertip 3D Reconstruction

[![Documentation](https://github.com/xdic-dev/xdic-dev.github.io/actions/workflows/docs.yml/badge.svg)](https://xdic-dev.github.io)


## How to run this app?

First modify the configuration files `src/global_param` and `src/dic_param`.

### For Windows users

```bash
matlab -nodisplay -nojvm -nosplash -nodesktop -r "run('xdic.m'); exit;"
```

### For Linux users

```bash
matlab -nodisplay -nojvm -nosplash -nodesktop -r "run('xdic.m'); exit;"
```

### For MacOS users

```bash
/Applications/MATLAB_RXXXa.app/bin/matlab -nodisplay -nojvm -nosplash -nodesktop -r "visu=true; run('xdic.m'); exit;"
```

### Set Parameters

if you want to set a parameter (e.g., visu) on the fly here is the how

```bash
matlab -nodisplay -nojvm -nosplash -nodesktop -r "visu=true; run('xdic.m'); exit;"
```

#### the list of possible parameters

- subject : indicates the subject identifier
- phase : indicates the phase identifier
- material : indicates the material identifier (1=glass, 2=coating, 3=coating_oil)
- calib_folder : indicates the calibration folder if it's different from the default 2
- ref_trial : indicates the reference trial number for the analysis (default is 5)
- frame_start : indicates the starting frame for the analysis
- frame_end : indicates the ending frame for the analysis
- jump : indicates the jump between the frames
- visu : indicates the visualization flag (true or false)
- debug : indicates the debug mode flag (true or false)
