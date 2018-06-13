# Estimator

This is the actually deployable part of the estimator, which has no dependency on CasADi.

## Build

First go to the project root and build and run the estimator code generator.

```bash
mkdir -p estimator/codegen
mkdir -p build
cd bulid
cmake ..
make estimator_codegen
./estimator_codegen ../config/config.yaml ../estimator/codegen
```

Then go to the estimator/ directory and build the estimator test app.

```bash
cd estimator
mkdir -p build
cd bulid
cmake ..
make estimator_test
./estimator_test
```
