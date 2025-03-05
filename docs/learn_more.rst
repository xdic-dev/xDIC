
Learn More about multiDIC
=========================

Core Concepts: 2D to 3D Reconstruction
--------------------------------------
multiDIC takes images from multiple synchronized cameras and reconstructs a 3D model using Digital Image Correlation (DIC). This involves:

1. **Camera Calibration**
2. **Digital Image Correlation**
3. **3D Reconstruction**

MultiDIC Process
----------------

.. image:: process.png
   :alt: multiDIC Process
   :align: center
   :width: 100%

Important Scripts and Their Roles
---------------------------------
- `runDIC.m`: The main script for DIC and 3D reconstruction.
- `config.m`: Configuration file for setting paths and parameters.
- `stepE_3DReconstruction.m`: Computes the 3D coordinates from 2D image points.

Detailed Workflow
-----------------
1. Data Collection
2. Set Configuration
3. Run DIC Analysis
4. Post-processing
5. Troubleshooting and Adjustments

Additional Resources
--------------------
MATLAB documentation, academic papers, and online tutorials on 3D reconstruction can further help understand multiDIC’s inner workings.
