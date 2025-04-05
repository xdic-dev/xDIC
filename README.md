# xDIC: Fingertip 3D Reconstruction

| **Documentation and paper** | **Build Status** |
|:-----------------:|:----------------:|
| [![DOI][paper-img]][paper-url] [![][docs-latest-img]][docs-latest-url]  [![][docs-stable-img]][docs-stable-url]        | [![Build Status][build-img]][build-url]  [![Codecov branch][codecov-img]][codecov-url]      |


[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-stable-url]: https://xdic-dev.github.io/xdic-dev.github.io/stable
[docs-latest-url]: https://xdic-dev.github.io/xdic-dev.github.io/dev
[paper-img]: https://proceedings.juliacon.org/papers/10.21105/jcon.00160/status.svg
[paper-url]: https://www.biorxiv.org/content/10.1101/2025.02.10.637366v1.full.pdf

[build-img]: https://github.com/xdic-dev/xdic-dev.github.io/actions/workflows/docs.yml/badge.svg?branch=main
[build-url]: https://github.com/xdic-dev/xdic-dev.github.io/actions?query=workflow%3ADocumentation
[codecov-img]: https://codecov.io/github/xdic-dev/xdic-dev.github.io/coverage.svg
[codecov-url]: https://app.codecov.io/github/xdic-dev/xdic-dev.github.io

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
