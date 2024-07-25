# Cryoablation Probability Module

This 3D Slicer module calculates the probability and coverage of cryoablation treatment for tumors accounting for 1-3 needle deviations. The script performs various tasks including loading sample data, calculating iceball coverage, and generating probabilistic results for cryoablation.

## Overview

This repository contains a 3D Slicer scripted loadable module called `Probability`. The module calculates the probability and coverage of cryoablation treatment for tumors, taking into account the positioning and size of iceballs generated during the procedure.

## Installation

1. **Clone the repository**:
    ```bash
    git clone https://github.com/CR7SC3/CryoProbability.git
    ```
2. **Open 3D Slicer**.
3. **Install the module**:
    - Navigate to `Edit` > `Application Settings` > `Modules`.
    - Add the path to the cloned repository to the `Additional module paths`.
    - Restart 3D Slicer.

## Usage

1. **Load the module**:
    - Open 3D Slicer.
    - Go to the `Modules` drop-down menu and select `Examples` > `Cryoablation Probability`.

2. **Use the module**:
    - Select the input volume and probe location.
    - Adjust the error parameters.
    - Click the `Apply` button to calculate the results.

## File Structure

- `Probability.py`: Main script for the module.
- `Resources/Icons`: Folder containing icons for the module.
- `UI/Probability.ui`: UI file created with Qt Designer.

## Functions and Classes

### Main Classes

- **Probability**: Initializes the module and registers sample data.
- **ProbabilityWidget**: Manages the GUI and interactions.
- **ProbabilityLogic**: Implements the core computation logic.
- **quad**: Dictates the Gaussian Quadrature in three spatial directions

### Key Functions

- `registerSampleData()`: Registers sample data for the module.
- `ProbabilityWidget.setup()`: Sets up the GUI.
- `ProbabilityLogic.runProbability()`: Calculates the probability and coverage of the cryoablation treatment.
- `ProbabilityLogic.iceball_coverage_1()`: Calculates the iceball coverage within the tumor.
- `quad.runQuadrature()`: Runs the quadrature calculations for the probabilistic model.

## Contributors

- Alvaro Cuervo
- Pedro Moreira

## Acknowledgements

This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc., Andras Lasso, PerkLab, and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.

## License

3D Slicer is a free open source software distributed under a BSD style license.

## Contact

For any questions or issues, please contact Alvaro Cuervo at alcuervosebastian@gmail.com or raise an issue in the project's GitHub repository.
