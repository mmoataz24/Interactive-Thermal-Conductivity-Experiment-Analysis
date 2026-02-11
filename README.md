# Interactive Thermal Conductivity Experiment Analysis (MATLAB)

## Overview
This repository contains an **interactive MATLAB program** for analyzing **thermal conductivity experiments** using experimental data. The program applies **Fourier’s Law of heat conduction**, performs **complete uncertainty propagation**, and generates **statistical analysis and comparative visualizations** for multiple materials.

This project was developed for a **Heat Transfer** course in the **Mechatronics and Robotics Engineering** program and is intended for **educational and experimental analysis purposes**.

---

## Features
- Interactive, step-by-step user input
- Supports multiple materials and multiple trials
- Implements Fourier’s Law of heat conduction
- Full measurement uncertainty propagation
- Statistical analysis (mean, standard deviation, standard error)
- Professional multi-plot visualizations
- Automatically generated summary tables
- Comparative performance analysis between materials

---

## Experimental Principle
The program calculates thermal conductivity based on:
- Heat absorbed by water
- Temperature gradient along a cylindrical rod
- Rod geometry and measurement uncertainties

The analysis assumes **steady-state, one-dimensional heat conduction**.

---

## Requirements
- MATLAB (R2018 or later recommended)
- No additional toolboxes required

---

## How to Use

1. **Download the program**
   - Download the `main.m` file from this repository.

2. **Run the script**
   - Open MATLAB.
   - Navigate to the folder containing `main.m`.
   - Run the script

3. **Follow the on-screen instructions**
   - Enter the number of materials to be tested.
   - For each material, input:
     - Material name
     - Rod geometry (length, diameter) and uncertainties
     - Experimental temperature readings
     - Water temperature change
     - Time duration for each trial

4. **Enter your experimental data**
   - Type your **measured experimental outputs** as prompted in the MATLAB command window.
   - Time input supports flexible formats (e.g., seconds or minutes).

5. **View results**
   - The program automatically:
     - Computes thermal conductivity and uncertainty
     - Displays formatted numerical results
     - Generates comparative plots
     - Outputs a summary table in the command window

---

## Output
The program generates:
- Thermal conductivity values for each trial
- Average thermal conductivity with uncertainty
- Comparative bar charts between materials
- Relative performance plots
- Console summary tables

---

## Educational Value
This project demonstrates:
- Application of heat transfer theory
- Proper engineering uncertainty analysis
- Experimental data processing
- MATLAB programming with structured data
- Clear technical reporting and visualization

---

## Assumptions & Limitations
- Assumes one-dimensional heat transfer
- Neglects heat losses to surroundings
- Assumes constant material properties
- Limited to cylindrical geometries

---

## Potential Improvements
- Include convection and radiation losses
- Support non-cylindrical geometries
- Add temperature-dependent material properties
- Integrate real-time data acquisition
- Export results to Excel or PDF

---

## Suggested Topics
-matlab
-heat-transfer
-thermal-conductivity
-uncertainty-analysis
-experimental-analysis
-engineering-education
-mechatronics

---

## License
This project is intended for **educational and academic use**.  
You are free to modify and extend the code with proper attribution.
