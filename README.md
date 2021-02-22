# Robot Localization and Standard Filter Comparison

This project contains implementation of 4 main filters:
- Extended Kalman Filter
- Unscented Kalman Filter
- Particle Filter
- Invariant Kalman Filter

To see the derivations of the propogation and correction equations of the filters, look at https://github.com/nalinbendapudi/Robot-Localization/blob/master/report.pdf

## Usage Instructions

In MATLAB, enter the command: `run(100,0,0,"EKF")`. This will run the environment for 100 time-steps and localize using the EKF filter. To switch between filters enter either
"EKF", "UKF", "PF" or "InEKF" into the run function. Enter a 1 into the thrid input of the run function to create a video.

## Results

The following plots show the deviation of localized robot path from ground truth using different filters

![](https://github.com/nalinbendapudi/Robot-Localization/blob/master/EKF_deviationPlot.jpg)

Deviation Plot for Extended Kalman Filter

![](https://github.com/nalinbendapudi/Robot-Localization/blob/master/UKF_deviationPlot.jpg)

Deviation Plot for Unscented Kalman Filter

![](https://github.com/nalinbendapudi/Robot-Localization/blob/master/PF_deviationPlot.jpg)

Deviation Plot for Particle Filter

![](https://github.com/nalinbendapudi/Robot-Localization/blob/master/InEKF_deviationPlot.jpg)

Deviation Plot for Invariant Kalman Filter
