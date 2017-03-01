# Modbox

## What is it?
Modbox is an algorithm for computing generic Fourier integral operators (FIOs) -- specifically, FIOs associated with canonical graphs ("graph FIOs" for short).
Some examples of these operators are:
* acoustic/elastic wave propagation (and other solution operators for hyperbolic PDEs)
* integral transforms: the Radon transform, X-ray transform (used in CT imaging), and generalizations
* differential and pseudodifferential operators

This code implements a refined version of the modified box algorithm described by de Hoop, Uhlmann, Vasy, and Wendt in ["Multiscale discrete approximations of Fourier integral operators associated with canonical transformations and caustics."](http://epubs.siam.org/doi/abs/10.1137/120889642)

You can find the mathematical details in the paper ["Computing Fourier integral operators with caustics."](http://iopscience.iop.org/article/10.1088/0266-5611/32/12/125001/meta)

## What do I need to run it?
Modbox requires a not-too-ancient version of Matlab (2012b or newer should be fine). A C++ compiler is also required for the MEX subroutines.

If you'd like to use pseudopolar FFTs, you will also need PolarLab, which can be downloaded [here](http://www.cs.technion.ac.il/~elad/software/).

## How do I install and try it out?
Compile each of the MEX C++ routines in the `mex/` folder. You will also need to add all Modbox's subfolders to your Matlab search path (most easily done with `pathtool`).

To try out Modbox, start with the driver program `driver/ApplyFIO.m`. It takes a number of arguments (describing the input function, the FIO, resolution, and so on); some convenient default settings are given by `driver/BaseArgs.m`.


## How do I compute a new FIO?
A graph FIO is specified very naturally in Modbox by its canonical transformation and amplitude. 
Some common FIOs (e.g., the Radon transform) are not graph FIOs themselves, but are sums of graph FIOs.
For this reason, Modbox's input is a sum of graph FIOs, represented by a subclass of the abstract class `FIO`.
Each graph FIO is represented as an object of a subclass of the `FIOPatch` abstract class, so an `FIO` is just a collection of `FIOPatch`es.

To write a new FIO, then, create a subclass of `FIOPatch` with implementations of the canonical transformation and amplitude functions.
For better accuracy, you can include analytical derivatives of the canonical transformation; if not included, Modbox will substitute a
finite difference approximation.
If your FIO has multiple graph FIO components, you can create a subclass of `FIO` to keep them together; otherwise, you can use the standard `SinglePatchFIO` subclass.

Note that the algorithm has two phases: a preprocessing phase (which only needs to be done once for a given FIO and input size) and the application phase.
In preprocessing, the input `FIO` is converted to a `ProcessedFIO` object containing the data needed to apply the operator to data of a given size.
