# Simulate multi-electrode array traces

Code simulates the recorded activity of a population of neurons on a multi-electrode array. The design of these simulations incorporated realistic variations and noise characteristics seen in actual neural data. Designed for testing spike sorting algorithms.

## NOTE

Code and documention in process of being made more user-friendly for publication and use.

## Directories

### cpp

Simulation code directory

* parameters are set by command line flags at the moment, will be updated to read/write json parameters files

#### Compile

```
make simtraces
```

### SortEval
Matlab code to analyze the results of spike sorting performed on simulated traces
