# Effects of using error models trained only on _nifH_ reads

To evaluate the impact of using error models that were trained only on _nifH_-like reads, we ran
eight _nifH_ amplicon data sets through the pipeline twice, with _nifH_ error models and with
standard error models that used all reads.

The first plot shows the percentage of reads in each sequencing run that were *gained* by using
_nifH_ error models.  Gains are shown at each stage of the pipeline from ASV inference ("Denoise")
through bimera removal, as well as after subsequent filtering for ASVs with lengths expected to
capture _nifH_ (280 - 360 nt).  Gains were usually small and occurred during ASV inference.  For
NEMO and Harding_2021 the median gains were > 3%.

![ASVs if use nifH error models](read_gain_if_use_nifH_errror_models.png)


The next plot shows that for most studies the use of _nifH_ error models resulted in more total
reads being used but in fewer total ASVs.  In other words, using the standard error models resulted
in more ASVs with very low counts.

![Read gain if use nifH error models](asv_abund_comparison2_if_use_nifH_error_models.png)


The following plot compares each ASV's abundance when run with versus without _nifH_ error models --
for ASVs that are plotted off of the axes.  The axes order the ASVs by abundance rank (#1 is the
most abundant).  Mostly the same ASVs were found using either error model (most points are off the
axes).  Using either error model, the top ~100 ASVs usually had identical ranks (ASVs on _y_ = _x_)
and very similar abundances (log ratio ~1).  However, for rarer ASVs the ranks and abundances
differed between the two error models.  Moreover, using standard models resulted in more total ASVs,
usually a few hundred to a few thousand per study. This is evident by the greater number of ASVs
found that were unique to the standard model (on the _y_ axis) compared to ASVs unique to the _nifH_
model (on the _x_ axis).

![Read gain if use nifH error models](asv_abund_comparison_if_use_nifH_error_models.png)