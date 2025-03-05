
Quick Start
===========

Installation Instructions
-------------------------
Before you can run multiDIC, follow these steps to install and configure your environment:

1. **Download the Project**  
   You can download the latest version of multiDIC from the repository or extract it from a provided zip file.

2. **MATLAB Requirements**  
   Ensure that MATLAB (version 2020a or later) is installed on your machine. You will also need the following MATLAB toolboxes:
   - **Image Processing Toolbox**: For processing and analyzing the multi-camera images.
   - **Parallel Computing Toolbox** (optional): For improving performance with parallel processing during large datasets or computationally intensive tasks.

3. **Add multiDIC to MATLAB Path**  
   Open MATLAB and add the multiDIC project folder (where the `xDIC` folder resides) to the MATLAB path.

4. **Run a Test Example**  
   To confirm everything is working, run the basic example:
   ```matlab
   runDIC;
   ```

First Steps with `runDIC.m`
---------------------------
The core script to execute the DIC process is `runDIC.m`. Before running, configure parameters in `config.m` and execute the script:
```matlab
runDIC;
```

Basic Example for Finger 3D Reconstruction
------------------------------------------
1. Prepare the data (images from multiple calibrated cameras).
2. Run the 3D reconstruction:
```matlab
runDIC;
```
3. Visualize the results using `plot_multidic_results.m`.
