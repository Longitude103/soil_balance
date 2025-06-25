# Soil Balance Package

This package will use several of the other packages to run a soil balance on the
parcel. This soil balance will be created to return a few different results but mostly
to return what the crop water use was, the soil deep percolation and the run-off values are.

## Items Implemented

- [ ] Create a parcel with basic attributes
- [ ] Load soil parameters for parcel
- [ ] Get Crop Water Use
- [ ] Get Precipitation and calculate run-off
- [ ] Get Irrigation Amount
- [ ] Soil Evaporation (future)
- [ ] Calculate Daily Balance
- [ ] Return daily deep percolation and run-off
- [ ] Return crop water use and evaporation

## QC Process

The QC process was to compare this against a run of Hydrus-1D application version 4.17 installed on a Windows machine.
The process was to set up the same input parameters in both applications and compare the results of the bottom flux of
each
run. The result was that the solutions matched very closely (within 0.01), which is acceptable given the different
machine architectures and programming techniques used.