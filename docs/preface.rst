
Preface
=======

How to read this tutorial
-------------------------
Welcome to the multiDIC documentation. This tutorial is designed to help users get familiar with the multiDIC software, which enables 3D reconstruction of fingers using multi-camera data. Each section provides progressively deeper insights into how the software works, from basic usage to advanced customization. It’s recommended to start with the **Quick Start** section to get hands-on experience before diving into more complex topics.

What you will NOT find in this document
---------------------------------------
This document focuses on the technical implementation and usage of multiDIC. You will not find:
- Detailed mathematical background behind 3D reconstruction (though external references will be provided).
- Deep dive into MATLAB programming (only what's necessary for multiDIC usage).
- Guidelines for advanced hardware setups (such as complex multi-camera systems).

Installation
------------
To install multiDIC, follow these steps:
1. Clone the repository or download the project from the source.
2. Ensure you have MATLAB installed on your system (version 2020a or later recommended).
3. Add the multiDIC project folder to your MATLAB path.
4. Install any required MATLAB toolboxes such as the Image Processing Toolbox and Parallel Computing Toolbox (as needed for specific features). Mainly ensure [NCorr]() is installed
5. Run a simple test to check everything is set up correctly (use the `runDIC.m` script).

Getting Help
------------
If you encounter any issues, you can:
- Refer to the **Troubleshooting** section in the **Learn More about multiDIC** chapter.
- Check out the MATLAB documentation for any issues related to general MATLAB commands.
- Contact the developers through the official multiDIC repository or the provided communication channels.

Learning
--------
To get the most out of multiDIC, it’s important to understand the basics of 3D reconstruction, camera calibration, and MATLAB programming. The **Learn More about multiDIC** section will provide additional resources, but here are a few external references that might help you get up to speed:
- MATLAB documentation for image processing and multi-camera setups.
- Academic papers on 3D reconstruction techniques.
- Online tutorials for MATLAB basics (if you're new to MATLAB).

Miscellaneous
-------------
- **Acknowledgments**: We thank all contributors who have helped in the development of multiDIC.
- **Citation Guidelines**: If you use multiDIC in your research, please cite the following papers and repository:
  - [Include citation guidelines here].
- **Useful Links**: 
  - MATLAB Documentation: [MATLAB Documentation Link]
  - Related research papers: [Insert links here]
