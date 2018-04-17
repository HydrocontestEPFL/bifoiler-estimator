# HydroContest Bifoiler state estimator

## Build example

```bash
    mkdir -p build
    cd bulid
    cmake ..
    make casadi_test
    ./casadi_test
```

## Dependencies

- CASADI 3.3 (tested with casadi-matlabR2015a-v3.3.0 under OSX High Sierra)
- yaml-cpp 0.6

Note: you need to provide CASADI_PREFIX environment variable pointing to the CASADI package.
