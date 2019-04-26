# Jet decomposition project

## Prerequisites

### root

If you type `root-config` on the command line there should be something

### gcc

We need `-std=c++17` support



## Install and run

Once cloned, type `make` to compile everything

Make a soft link from here to the sample directory in order to run the example from the makefile

```
ln -s /path/to/sample/dir Samples/
```

On lxplus one would do

```
ln -s /afs/cern.ch/user/c/chenyi/EOSBox/Share/17364_DiscretizedJetSample// Samples/
```



## Executables

### BasicJetImageAnalysis

Runs the basic jet image analysis from Yang-ting some time ago.  Type

```
make RunBasicJetImageAnalysis
```

to run the example





### TruncationTest

Runs the truncation test with Fourier transformation.  Type

```
make RunTruncationTest
```

to run the example




### FourierCatalog

Runs a simple program to demonstrate the transformed result for test jets and background shapes to build intuition

```
make RunFourierCatalog
```

to run the example

