# wls_general_analysis
Shared repository for general quantification, statistical analysis, and data plotting scripts.
For use in the manuscriptt: Hydrostatic pressure shapes and canalizes semicircular canal morphology to ensure vestibular function.
doi.org/10.64898/2026.07.22.740136

Written by Jiacheng Wang and Kira L. Heikes.

Written in the Munjal Lab at Duke University.

## Repository Contents
- `AngleMeasurement.m` - Matlab function to measure pillar angles from 3D coordinates chosen in Fiji Mastodon ('spots')
- `BatchAngleMeasurements.m` - Matlab Script for batch measuring pillar angles from 3D coordinates chosen in Fiji Mastodon ('spots')
- `IndividualTimelapsePlotting.m` - Matlab Script to plot pillar angle and OV volume data in time from individually timelapsed Zebrafish OVs
- `PlotTimecourseAngleVolumeDataAndBinsAcrossAllTime.R` - R script for plotting angles vs volumes measured from pooled timecourse data
- `PlotTimelapseVSPerturbationsNormalizedAngleVolumeData.R` - R script for spearman correlation calculation and plotting of normalized angles vs volumes measured from timelapse and perturbation data
- `Plot_LargeCohort_AspectRatio_54-72hpf.R` - 
- `TimeCourse_OV-Tissue-Lumen-VolumePlot_46-72hpf.R` - 

## Dependencies
- `AngleMeasurement.m`
    -none; written in Matlab version 26.1
- `BatchAngleMeasurements.m`
    -must be stored in the same folder as the script `AngleMeasurement.m`; written in Matlab version 26.1
- `IndividualTimelapsePlotting.m`
    -none; written in Matlab version 26.1
- `PlotTimecourseAngleVolumeDataAndBinsAcrossAllTime.R`
    -requires R libraries installed: tidyverse and dplyr; written in R version 4.6.0
- `PlotTimelapseVSPerturbationsNormalizedAngleVolumeData.R`
    -requires R libraries installed: tidyverse and dplyr; written in R version 4.6.0
- `Plot_LargeCohort_AspectRatio_54-72hpf.R`
    -requires R libraries installed: tidyverse
- `TimeCourse_OV-Tissue-Lumen-VolumePlot_46-72hpf.R`
    -requires R libraries installed: tidyverse, readxl, and extrafont
