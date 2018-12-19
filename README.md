# DyConPro
DYnamic CONnectivity PROcessing (DCP) toolbox.  

MATLAB tools for computation, manipulation and processing of single-subject and group dynamic functional connectivity (dfc) tensors and matrices.  It includes a function for computing a time-varying, instantaneous measure of functional connectivity (function dcp_dcs; see [Tobia et al., 2017](http://onlinelibrary.wiley.com/doi/10.1002/hbm.23821/full)) and numerous other functions to scale, binarize, transform to matrices, and analyze dfc tensors.

Not just for dFC tensors anymore!

Denoise your 4D fMRI data.
Compute seed-based phase synchrony in 4D volumes.
Extract roi time courses from 4D volumes.
Generate correlated time series for simulations.
Dyadic and stopband Butterworth filters, or use PCA filtering on your time series matrix.
Plus so much more!

Recent update and first official release of DCP v1.1 (12/18/2018) now includes many new functions and two example files to demo how to run some of the more intricate DCP functions, such as dcp_stepmoreg for denoising, or dcp_fast_ifc_s2b for whole-brain seed-based phase synchrony. Many new functions read and write AFNI-style BRIK/HEAD files, and some will accept nifti-2 files too.

ALSO - now includes FC_SNAPSHOT, which is a crude matlab workspace gui for easy viewing of FC matrices and dFC tensors in 3 viewing modes, (1) time-dependent matrices, (2) time-dependent circle graph, and (3) time-dependent brain rendering. This is useful if you are exploring dFC tensors in matlab and want to view and scroll through time. It can scan the workspace for variables which you can then select to load/view in the gui. You can load coordinates, roi labels, and a template brain from within the gui. You can also pop out a figure to customize, save or print.

Some of these newer functionalities, especially the gui, are, of course, a ongoing WIP!

## Requirements
MATLAB version 2013b or higher.

## Usage
Command line tools. Some functions are stand-alone, some are dependencies, and some are both. Read the Input/Output information and Notes for each function file to learn how to use it. 

## References
Citing the toolbox:
- Tobia, M. J., Hayashi, K., Ballard, G., Gotlib, I. H., & Waugh, C. E. (2017). Dynamic functional connectivity and individual differences in emotions during social stress. Human brain mapping, 38(12), 6185-6205.

## Support and Communication
If you notice any bugs in the package, or would like to request features, please open an [issue](https://github.com/NBCLab/DyConPro/issues). We make no guarantee that we will incorporate requested changes, but we welcome any input and are happy to discuss.
