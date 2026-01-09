# Caldwell_turtle_incubation-performance

[![DOI](https://zenodo.org/badge/18202869.svg)](https://doi.org/10.5281/zenodo.18202869)

Version controlled and editable source for the data and code supporting the paper
__"Incubation and overwintering conditions influence righting performance of hatchling turtles"__ by _Molly Folkerts Caldwell, Daniel A. Warner, and Matthew E. Wolak_.

## Data

### Data citation

If you use these data, please cite the publication (preference)

>M. Folkerts Caldwell, D.A. Warner, and M.E. Wolak. _accepted_. Incubation and overwintering conditions influence righting performance of hatchling turtles. _Journal of Experimental Zoology - Part A_. [https://doi.org/](https://doi.org/).

or data package 

<!-- TODO
>Matthew E. Wolak. (2023). qgevoeco/Caldwell_turtle_nest-choice-predation: Initial release v1.0.0 (v1.0.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.7630637
-->


### Data metadata

Column header descriptions for the datasets analyzed in the paper:

#### "eggSurvival_data.txt" file

A tab-separated text file consisting of:

  - `Year` integer values indicating the replicate of the experiment which took place in either 2019 or 2020.
  
  - ` Clutch` integer values indicating a clutch identity for each egg.
   
  - `EggSurvival` integer values indicating eggs that  survived=1 or did not survive=0 to hatching.
  
  - `IncTemp` integer value for the temperature (degrees Celsius) at which an egg was incubated.
  
  - `EggID` character values assigned to give a unique identifier to each egg. Structured to indicate the `Clutch` and within-clutch egg number, separated by a dash.
  
  - `EggMass` numeric value for initial mass of eggs (grams).

#### "hatchlingAndRightingResponse_data.txt" file

A tab-separated text file contains each righting response trial measurement as a separate row and consisting of:

  - `Year`, `Clutch`, `EggMass`, and `IncTemp` as described above for the "eggSurvival_data.txt" file.
  
  - `HatchlingID` integer values assigned to give a unique identifier to each hatchling.
  
  - `HatchlingNumber` character values assigned to give a unique identifier to each hatchling within each year. Structured to indicate, within each year, the `Clutch` and within-clutch number, separated by a dash.
  
  - `MaxIncLength` 
  
  - `ScaledDevMaxIncLength`
  
  - `AgeAtTrials` integer number of days since hatching when a righting response trial took place.
  
  - `WinterHous` character denoting the "Terrestrial" or "Aquatic" treatment applied for a hatchling's overwintering conditions.
  
  - `DaysWintHous` integer number of days a hatchling spent in their winter housing treatment.
  
  - `SCL2` numeric straight carapace length in millimeters as measured with digital calipers after overwintering and at the time of the righting response trial period of the experiment.
  
  - `CW2` numeric straight carapace length in millimeters as measured with digital calipers (between the junction of the sixth and seventh marginal scutes on the carapace) after overwintering and at the time of the righting response trial period of the experiment.
  
  - `PL2` numeric straight plastron length in millimeters as measured with digital calipers after overwintering and at the time of the righting response trial period of the experiment.
  
  - `PW2` numeric straight plastron width in millimeters as measured with digital calipers after overwintering and at the time of the righting response trial period of the experiment.
  
  - `BD2` numeric shell depth in millimeters as measured with digital calipers (one caliper arm placed at the junction between the second and third vertebral scutes on the carapace and the other caliper arm placed flat along the plastron) after overwintering and at the time of the righting response trial period of the experiment.
  
  - `NV2` numeric notch-to-vent length in millimeters as measured with digital calipers spanning the posterior end/notch of the plastron to the cloaca, or vent, when the tail is firmly held straight. This measurement was taken after overwintering and at the time of the righting response trial period of the experiment.
  
  - `VT2` numeric vent-to-tip in millimeters as measured with digital calipers spanning the cloaca, or vent, to the tail tip when the tail is firmly held straight. This measurement was taken after overwintering and at the time of the righting response trial period of the experiment.
  
  - `Mass2` numeric mass of the hatchling in grams as measured with a digital balance after overwintering and at the time of the righting response trial period of the experiment.
  
  - `Round` integer value to indicate the first, second, or third righting resposne trial.
  
  - `RRTimeSec` numeric righting response time in seconds.
  
  - `LatencySec` numeric time in seconds from the start of the trial until the hatchling attempted to right itself.
  
  - `ActiveRRTimeSec` numeric time in seconds measuring the amount of time a hatchling was actively attempting to right itself.
  
  - `RRTimeMin` numeric value expressing the same quantity as `RRTimeSec`, but in minutes.
  
  - `LatencyMin` numeric value expressing the same quantity as `LatencySec`, but in minutes.
  
  - `ActiveRRTimeSec` numeric value expressing the same quantity as `RRTimeSec`, but in minutes.
  

#### "hatchling_data.txt" file

A tab-separated text file contains data on each hatchling as a single row and consisting of:

  - `Year`, `Clutch`, `EggMass`, and `IncTemp` as described above for the "eggSurvival_data.txt" file.
  
  -  `HatchlingID`, `HatchlingNumber`, `MaxIncLength`, `ScaledDevMaxIncLength`, `AgeAtTrials`, `WinterHous`, `DaysWintHous` as described above for the "hatchlingAndRightingResponse_data.txt" file.
  
  - `SCL1`, `CW1`, `PL1`, `PW1`, `BD1`, `NV1`, `VT1`, `Mass1` are the morphological traits as described above for the "hatchlingAndRightingResponse_data.txt" file for `SCL2`, `CW2`, `PL2`, `PW2`, `BD2`, `NV2`, `VT2`, `Mass2`, but the names ending in `1` indicate these measures were taken at the time of hatching.
  
  - `DaysBwMeas` integer number of days between the hatchling measurements (names ending in `1`) and the measurements after overwintering and at the time of the righting response trial period of the experiment (names ending in `2`).

#### "incubationLengths_data.txt" file

A tab-separated text file consisting of:

  - `Year`, `Clutch`, `EggMass`, and `IncTemp` as described above for the "eggSurvival_data.txt" file.
  
  - `HatchlingID`, `HatchlingNumber`, `MaxIncLength`, `ScaledDevMaxIncLength` as described above for the "hatchlingAndRightingResponse_data.txt" file.
  
  - `IncLength` integer number of days of incubation for each hatchling
  

## Changes
For ease of reference, an overview of significant changes to be noted below. Tag comments with commits or issues, where appropriate.


