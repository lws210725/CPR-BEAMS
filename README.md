# CPR-BEAMS
Example MATLAB scripts for working with CPR-BEAMS data via BCO-DMO

Scripts to extract data and generate time-series

Scripts perform the following functions:

Extract and map the locations of available samples within a specified area, for a specified range of dates, from the CPR-BEAMS resource on BCO-DMO

Extract abundance data about specified taxa from the CPR-BEAMS resource on BCO-DMO

Return a monthly time-series for abundance of a specified example taxon within an example area, for an example range of dates

Grid the data and produce localised seasonal abundance profiles for representative example taxa

Grid the data and produce localised monthly abundance timeseries for representative example taxa

Grid the data and produce localised monthly abundance timeseries for representative example taxa, highlighting missing samples and backfilling the data (with the local median value for any given month)

Initial code for the 2023/2024 version of this project is available on https://github.com/lws050721/CPR-BEAMS and includes

cprbeams_datareader.m

cprbeams_datareader_arearestriction.m

which read BCO-DMO data to end of 2021.   It is preserved here to assist with backwards compatibility.   New code for the 2024/2025 version of this project is available on https://github.com/lws210725/CPR-BEAMS and includes the above plus 

cprbeams_datareader_1958_2022.m

cprbeams_datareader_arearestriction_1958_2022.m

which reads BCO-DMO data to end of 2022.   A further update is anticipated in 2026.
