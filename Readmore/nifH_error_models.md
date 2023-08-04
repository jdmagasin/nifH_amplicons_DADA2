# Effects of using error models trained only on _nifH_ reads

To evaluate the impact of using error models that were trained only on _nifH_-like reads, we ran
eight _nifH_ amplicon data sets through the pipeline twice, with _nifH_ error models and with
standard error models that used all reads.

The first plot shows the percentage of reads in each sequencing run that were *gained* by using
_nifH_ error models.  Gains are shown at each stage of the pipeline from ASV inference ("Denoise")
through bimera removal, as well as after subsequent filtering for ASVs with lengths expected to
capture _nifH_ (280 - 260 nt).  Gains were usually small and occurred during ASV inference.  For
NEMO and Harding_2021 the median gains were > 3%.

![ASVs if use nifH error models](asv_abund_comparison_if_use_nifH_error_models.png)


The following plot compares ASVs found when the pipeline was run with vs. without _nifH_ error
models.  Axes indicate ASV rank abundances (#1 the most abundant).  Mostly the same ASVs were found
using either error model, indicated by points which are not on the axes.  For the top ~100 ASVs,
ranks were usualy identical (ASVs is on _y_ = _x_) and the total reads assigned to the ASV was
usually very similar (log abundance ratio ~1).  However, for rarer ASVs the ranks and abundances
differed between the two pipeline runs.  Additionally, using the standard error model found more
unique ASVs (on _y_ axis) compared to using the _nifH_ model (on _x_ axis), even though _nifH_
models resulted in using more total reads (first plot).

![Read gain if use nifH error models](read_gain_if_use_nifH_errror_models.png)