
Dev Place
=========

File Structure Overview
-----------------------
- **xDIC/**: Main folder containing core scripts.
- **lib/**: Core algorithms and functions for 3D reconstruction.
- **utils/**: Helper scripts for data loading, error logging, and plotting.

Modifying Key Components
-------------------------
1. **Adding New Input Data Types**
   Modify data loading functions in `utils/` and calibration settings in `config.m`.

2. **Customizing the 3D Reconstruction Process**
   Modify `stepE_3DReconstruction.m` for changes in how 3D points are computed.

3. **Enhancing Plotting and Visualization**
   Modify `plot_multidic_results.m` for new visualization options.

Contribution Guidelines
-----------------------
1. Fork the repository, make changes, and submit a pull request.
