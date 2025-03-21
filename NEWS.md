# spatialLIBD 1.19.11

BUG FIXES

* @Nick-Eagles fixed a bug for specifying colors in `vis_clus()` and related
functions when `NA` values are present in the cluster variable. For more
details check <https://github.com/LieberInstitute/spatialLIBD/pull/102>.

# spatialLIBD 1.19.10

NEW FEATURES

* @lahuuki added the `gene_name` argument to `sig_genes_extract()` and
`sig_genes_extract_all()` to make these functions more flexible.
See <https://github.com/LieberInstitute/spatialLIBD/pull/101> for details.

# spatialLIBD 1.19.9

BUG FIXES

* @lahuuki resolved the bug we discovered on a dataset being analyzed by
@Nick-Eagles as we documented at
<https://github.com/LieberInstitute/spatialLIBD/issues/98>. Basically, an
internal function used by `layer_stat_cor_plot()` was having issues with cluster
names that were just numbers such as "1", "2", ..., "13" where
`annotate_registered_clusters()` could result in a "9/13" annotation but then
that would match with "1", "3", "9", and "13" instead of just "9" and "13", thus
introducing incorrect "X" annotations on `layer_stat_cor_plot()`. @lahuuki
resolved this at <https://github.com/LieberInstitute/spatialLIBD/pull/99>.

# spatialLIBD 1.19.8

NEW FEATURES

* `vis_gene()` now has the `cap_percentile` argument as implemented by
@Nick-Eagles. This allows you to cap the expression values at a certain
percentiles, which can be useful to exclude super high outlier values which
are common with Visium HD data. See
<https://github.com/LieberInstitute/spatialLIBD/pull/97> for details.

BUG FIXES

* Fixed colors on `layer_stat_cor_plot()` thanks to @lahuuki. Details at
<https://github.com/LieberInstitute/spatialLIBD/pull/94>.

# spatialLIBD 1.19.7

BUG FIXES

* Removed code that made the app very slow to load using `HDF5Array` objects.
See
<https://github.com/LieberInstitute/spatialLIBD/commit/b80d92c3271a6ad92859f79a3bc343f77bad9bf2>
for details.
* Fixed a bug on the Z-score calculation in `multi_gene_z_score()`. See
<https://github.com/LieberInstitute/spatialLIBD/commit/2d17ea2c3d1b73b38fbd50503765052c8487b9b1>
for details. Implemented by @Nick-Eagles.
* Fix a bug in the documentation of `run_app()` for the example stitched data.
* No longer point to Twitter, instead point to Bluesky. This is for the package
main README file as well as the default app documentation files.

# spatialLIBD 1.19.5

BUG FIXES

* Fixed internal errors in `add_qc_metrics()` on the `scuttle::isOutlier()`
function calls.

# spatialLIBD 1.19.4

NEW FEATURES

* @lahuuki fully re-implemented `gene_set_enrichment_plot()` using
`ComplexHeatmap::Heatmap()`. This new version has several new arguments that 
allow adding more annotation to the resulting heatmap. See
<https://github.com/LieberInstitute/spatialLIBD/pull/93> for more details. This
also means that `layer_matrix_plot()` has been removed from the package since
it previously served as a helper function for `gene_set_enrichment_plot()`.

# spatialLIBD 1.19.3

BUG FIXES

* Resolved <https://github.com/LieberInstitute/spatialLIBD/issues/90> which made
`add_key()` too strict and would create issues with `export_cluster()`. Reported
by @lahuuki and @manishabarse.

# spatialLIBD 1.19.2

BUG FIXES

* Merged <https://github.com/LieberInstitute/spatialLIBD/pull/92> by @lahuuki.
This fixes https://github.com/LieberInstitute/spatialLIBD/issues/72 and
https://github.com/LieberInstitute/spatialLIBD/issues/48 by making 
`registration_pseudobulk()` more robust. The original issues were reported by
@boyiguo1 and @berniejmulvey.

# spatialLIBD 1.19.1

NEW FEATURES

* Merged <https://github.com/LieberInstitute/spatialLIBD/pull/91> by @lahuuki.
This pull request fully re-implemented `layer_stat_cor_plot()` with a version
that uses `ComplexHeatmap::Heatmap()` internally. It also adds support for
incorporating the automatic annotation results from
`annotate_registered_clusters()`. NOTE that the `max` argument was renamed to
`color_max`, as well as `min` to `color_min`. Also, the default for `min` used
to be `-max` and now for `color_min` the default is the `min()` correlation
observed. The default for `max` was 0.81 and the default for `color_max()` is
the `max()` observed correlation.
* `run_app()` was also updated to match the updated in `layer_stat_cor_plot()`
and now has 2 new inputs for controlling the annotation process with
`annotate_registered_clusters()`. It also allows downloading a CSV file with
the annotation results.


# spatialLIBD 1.17.10

BUG FIXES

* `registration_wrapper()` now automatically handles the scenario where `k = 2`
by not using `registration_stats_anova()` and providing an apporpriate warning.
Implemented by @lahuuki at 
<https://github.com/LieberInstitute/spatialLIBD/pull/86>.

# spatialLIBD 1.17.9

BUG FIXES

* `read10xVisiumWrapper()` is now able to detect the GTF file used by
`SpaceRanger` for version 3.0.0+. Implemented by @nick-eagles at
<https://github.com/LieberInstitute/spatialLIBD/pull/88>.


# spatialLIBD 1.17.6

BUG FIXES

* Fixed the bug reported by @lahuuki about `vis_grid_clus()` not handling
`logical()` cluster variables.
See <https://github.com/LieberInstitute/spatialLIBD/issues/80>. To resolve this,
`sort_clusters()` and `get_colors()` had to change internally. Examples and
documentation for both functions have now been updated to showcase what happens
when you provide a `logical()` vector as an input.

# spatialLIBD 1.17.5

NEW FEATURES

* Added `add_qc_metrics()` inspired by 
<https://github.com/LieberInstitute/Visium_SPG_AD/blob/master/code/07_spot_qc/01_qc_metrics_and_segmentation.R>
which adds seven new columns to the `colData(spe)` that can be useful when
performing quality control of the data. Developed by @lahuuki.

# spatialLIBD 1.17.3

NEW FEATURES

* Added support for `SpatialExperiment` objects created with
`visiumStitched::build_spe()` 
<https://research.libd.org/visiumStitched/reference/build_spe.html> that 
stitch together multiple Visium capture areas. Developed by @Nick-Eagles.

# spatialLIBD 1.15.2

SIGNIFICANT USER-VISIBLE CHANGES

* `vis_gene()` now has a `multi_gene_method` argument which provides 3 methods
for combining multiple continuous variables: `z_score`, `pca`, and `sparsity`.
These options can now be used with `run_app()` (the interactive websites). These
methods are further illustrated and documented in a new vignette available at
<https://research.libd.org/spatialLIBD/articles/multi_gene_plots.html>. This
work was contributed by @Nick-Eagles.

# spatialLIBD 1.13.6

NEW FEATURES

* `vis_clus_p()`, `vis_clus()`, and `vis_grid_clus()` now all use implement the
`na_color` argument that was present in the `vis_gene()` functions. This
resolves https://github.com/LieberInstitute/spatialLIBD/issues/43 by @boyiguo1.

# spatialLIBD 1.13.5

NEW FEATURES

* `run_app()` now has a `auto_crop_default` argument set to `TRUE` by default.
It can be turned off in cases where you are displaying images that do not
follow the expected Visium grid dimensions, such as manually stitched images
that you don't want to automatically crop.

# spatialLIBD 1.13.4

NEW FEATURES

* Added `fetch_data("spatialDLPFC_Visium_example_subset")` which is a subset
of 3 samples with only the `lowres` images that can be used for example /
tutorial purposes.

# spatialLIBD 1.13.2

NEW FEATURES

* Louise A. Huuki-Myers @lahuuki added a vignette explaining the spatial
registration process and all related functions. See
<https://github.com/LieberInstitute/spatialLIBD/pull/46> for the full pull
request.

# spatialLIBD 1.11.13

SIGNIFICANT USER-VISIBLE CHANGES

* The vignette now has a section describing the data from the `spatialDLFPC`,
`Visium_SPG_AD`, and `locus-c` projects that were done by members of the
Keri Martinowich, Kristen R. Maynard, and Leonardo Collado-Torres LIBD teams as
well as our collaborators.

# spatialLIBD 1.11.12

SIGNIFICANT USER-VISIBLE CHANGES

* `fetch_data("Visium_SPG_AD_Visium_wholegenome_spe"")`,
`fetch_data("Visium_SPG_AD_Visium_targeted_spe")`,
`fetch_data("Visium_SPG_AD_Visium_wholegenome_pseudobulk_spe")`, and
 `fetch_data("Visium_SPG_AD_Visium_wholegenome_modeling_results")` have been
 added. Use this to access data from the
 <https://github.com/LieberInstitute/Visium_SPG_AD> project.

# spatialLIBD 1.11.11

SIGNIFICANT USER-VISIBLE CHANGES

* `fetch_data("spatialDLPFC_snRNAseq")` now works if you want to download
the snRNA-seq data used in <http://research.libd.org/spatialDLPFC/>.

# spatialLIBD 1.11.10

BUG FIXES

* `read10xVisiumAnalysis()` now supports `spaceranger` version 2023.0208.0
(internal 10x Genomics version) output files that store analysis CSVs under the
`outs/analysis_csv` directory instead of `outs/analysis` and also use the
`gene_expression_` prefix for each of the analysis directories. This was
tested with @heenadivecha on files from
<https://github.com/LieberInstitute/spatial_DG_lifespan/blob/main/code/02_build_spe/01_build_spe.R>.

# spatialLIBD 1.11.9

SIGNIFICANT USER-VISIBLE CHANGES

* `gene_set_enrichment()` now internally uses
`fisher.test(alternative = "greater")` to test for odds ratios greater than 1.
Otherwise odds ratios of 0 could be significant.

# spatialLIBD 1.11.4

SIGNIFICANT USER-VISIBLE CHANGES

* Several changes were made to the default plotting aspect of `vis_gene()`, 
`vis_clus()` and related plotting functions. This was done with input from
@lahuuki and @nick-eagles and is described in more detail at 
https://github.com/LieberInstitute/spatialLIBD/commit/8fa8459d8fa881d254824d43e52193bf2c3021c0.
Most noticeably, the aspect ratio is no longer stretched to fill the plotting
area, the NA values will be shown with a light grey that has alpha blending,
and the position of the legends has been made consistent between the plots.

# spatialLIBD 1.11.3

NEW FEATURES

* Added the function `frame_limits()` and introduced the `auto_crop` argument
to `vis_clus()`, `vis_gene()` and all related functions. This new function
enables automatically cropping the image and thus adjusting the plotting area
which is useful in cases where the image is not centered and is not a square.
This was based on work by @lahuuki at 
https://github.com/LieberInstitute/spatialDLPFC/blob/2dfb58db728c86875a86cc7b4999680ba1f34c38/code/analysis/99_spatial_plotting/01_get_frame_limits.R
and 
https://github.com/LieberInstitute/spatialDLPFC/blob/ef2952a5a0098a36b09488ebd5e36a902bb11b48/code/analysis/99_spatial_plotting/vis_gene_crop.R.


# spatialLIBD 1.9.18

BUG FIXES

* Fixed a bug related to `edgeR::filterByExpr()` inside of
`registration_pseudobulk()`.
* Moved the `min_ncells` filtering step to `registration_pseudobulk()` rather
than `registration_wrapper()` since you should drop low `ncells` before
using `edgeR::filterByExpr()`.

# spatialLIBD 1.9.15

BUG FIXES

* Fixed some bugs in `registration_stats_anova()` in cases where we only had
two different unique values to compute F-statistics with, when we need at least
3.
* Made some parts of `registration_stats_anova()` and
`registration_stats_pairwise()` more flexible.
* `registration_model()` now provides a more informative error message when you
have an empty factor level, thus leading to a non-full rank model matrix.

# spatialLIBD 1.9.12

NEW FEATURES

* Added functions for computing the modeling statistics used by the spatial
registration process. See `registration_wrapper()` and related functions.
* Added a function for using the output of `layer_stat_cor()` and for
labeling the clusters. This can help interpret the spatial registration results.
See `annotate_registered_clusters()` for more details.

# spatialLIBD 1.9.11

BUG FIXES

* Fixed bugs in `gene_set_enrichment()` for `reverse = TRUE` reported by
@sparthib.
* Added a `reverse` option on the shiny app under the gene set enrichment tab,
that we tested with the example `spe` data.

# spatialLIBD 1.9.10

SIGNIFICANT USER-VISIBLE CHANGES

* Improved the automatic color palette selector when you switch discrete
variables. It also now supports the ManualAnnotation option.
* Discrete variable (cluster) legend is no longer duplicated under the clusters
interactive tab.
* You can now search the model test, which helps if you have lots of tests
to choose from (this most likely occurs when you are looking at the pairwise
results).

# spatialLIBD 1.9.9

SIGNIFICANT USER-VISIBLE CHANGES

* Made the shiny application more memory efficient in different areas.
* Changed the default `point_size` from 1.25 to 2.
* Added the option to show or hide the spatial images on the grid panels in
the shiny web application. Turn off by default since it is more efficient.

# spatialLIBD 1.9.5

BUG FIXES

* Fix https://github.com/LieberInstitute/spatialLIBD/issues/41. Reported by
@abspangler13. Now the gene selector changes automatically when you change the
'model results' (model type) or 'model test' inputs. The gene selector is now
only shown inside the 'model boxplots' panel since it only affects that one.

# spatialLIBD 1.9.4

BUG FIXES

* Fix https://github.com/LieberInstitute/spatialLIBD/issues/40. Reported by
@Erik-D-Nelson.

# spatialLIBD 1.9.3

BUG FIXES

* Added a more informative error message when 'stats' does not have ENSEMBL
gene IDs as the `rownames()`. Reported by @abspangler13 and @sparthib at 
https://github.com/LieberInstitute/spatialLIBD/issues/33#issuecomment-1137544893

# spatialLIBD 1.7.19

SIGNIFICANT USER-VISIBLE CHANGES

* Documentation of the `layer-level data` panel at `run_app()` has been
significantly increased. You can now also visualize more than 2 reduced
dimensions computed on the pseudo-bulk level data (layer-level for the Maynard
et al, Nature Neurosci, 2021 data). 
* Users can now control the font and point size on the reduced dimension plots,
as well as the overall font size on the model boxplots.
* Image edit scenarios you might be interested
in for having a uniform color background image are now documented; for example
if you want a white or black background, or actually any valid R color name or
color HEX value.

# spatialLIBD 1.7.18

SIGNIFICANT USER-VISIBLE CHANGES

* `run_app()` now offers the option to chose any of the `paletteer::paletteer_d`
color palettes for discrete variables.
* `Polychrome` has been replaced as a dependency by `paletteer`. Note that
`Polychrome::palette36` is still the default.
* `run_app()` now looks for columns that end with '_colors' in their name
which can be used to pre-specify colors for any companion variables. For example
if you have `spe$my_groups` and `spe$my_groups_colors` then the second one
can specify the colors that will be used for visualizing `spe$my_groups`. This
makes specifying default colors more flexible than before, and the user is still
free to change them if necessary.

# spatialLIBD 1.7.17

BUG FIXES

* Fix bugs in `layer_boxplot()` where it was too specific to the Maynard et al
2021 data. We have made it more flexible now.
* Made the y-axis space more dynamic in `gene_set_enrichment_plot()` and
`layer_matrix_plot()`.

# spatialLIBD 1.7.16

BUG FIXES

* Fixed a bug in `sig_genes_extract()` when there's only one set of t-statistics
or F statistics to extract.

# spatialLIBD 1.7.12

SIGNIFICANT USER-VISIBLE CHANGES

* The visualization functions `vis_*()` of `SpatialLIBD` in this version
match the Bioconductor 3.15 version of `SpatialExperiment`.
Note that if you used `SpatialExperiment::read10xVisium()`, the names of
the spatial coordinates changed at 
https://github.com/drighelli/SpatialExperiment/commit/6710fe8b0a7919191ecce989bb6831647385ef5f
and thus you might need to switch them back if you created your
`SpatialExperiment` object before this change. You can do so with
`spatialCoordsNames(spe) <- rev(spatialCoordsNames(spe))`.
`read10xVisiumWrapper()` uses `SpatialExperiment::read10xVisium()` internally,
so this change on `SpatialExperiment` would then also affect you.

# spatialLIBD 1.7.11

NEW FEATURES

* Now `layer_stat_cor()` has the `top_n` argument which can be used for
subsetting the marker genes prior to computing the correlation as part of the
spatial registration process.

# spatialLIBD 1.7.10

NEW FEATURES

* Added the `add_key()` function to reduce code duplication and resolve
https://github.com/LieberInstitute/spatialLIBD/issues/31.

# spatialLIBD 1.7.9

NEW FEATURES

* This version is now compatible with the bioc-devel version of
SpatialExperiment where spatialData() was deprecated. Details at
https://github.com/LieberInstitute/spatialLIBD/pull/29/files.

# spatialLIBD 1.7.7

BUG FIXES

* Fixed a bug where the using the left-mouse click was not working for
annotating individual spots under the "gene (interactive)" tab.

# spatialLIBD 1.7.6

NEW FEATURES

* `vis_gene_p()`, `vis_clus_p()` and all related functions now have an argument
`point_size` which lets you control how big the points are plotted. This can
be useful for visualization purposes.
* The shiny app now has an input controlling the point size. If you increase it
to say `5`, then if you zoom in the `clusters (interactive)` panel, you can
see larger spots when zooming in.
* These features are related to
https://github.com/LieberInstitute/spatialLIBD/issues/28 although the spot
diameter is still not the true spot diameter. However, now you have more
flexibility for visualizing the spots.

# spatialLIBD 1.7.5

NEW FEATURES

* Expanded the _Using spatialLIBD with 10x Genomics public datasets_ vignette
to show how you can deploy your web application. See
https://libd.shinyapps.io/spatialLIBD_Human_Lymph_Node_10x/ for the live
example.

# spatialLIBD 1.7.4

BUG FIXES

* `vis_gene()` and `vis_grid_gene()` now support `geneid`s that are found in
the `rownames(spe)`. This makese these functions more flexible.
* `vis_grid_gene()` and `vis_grid_clus()` now have the `sample_order` argument
which gives you more control in case you want to plot a subset of samples. This
should also reduced the memory required as discovered at
https://github.com/LieberInstitute/spatialDLPFC/issues/45.

# spatialLIBD 1.7.3

NEW FEATURES

* Added support for more than one background picture per sample. This was done
through the new argument `image_id`. Resolves
https://github.com/LieberInstitute/spatialLIBD/issues/25.
* Added options for side by side visualization of the background image and the
clusters or gene expression values in the static versions. Resolves
https://github.com/LieberInstitute/spatialLIBD/issues/19.
* Allow changing the transparency level of the spots with the `alpha` argument.
Resolves https://github.com/LieberInstitute/spatialLIBD/issues/20.
* Add support for image manipulation with the `magick` package. Adds functions
`img_edit()`, `img_update()` and `img_update_all()` as well as new features
on the web application. Resolves
https://github.com/LieberInstitute/spatialLIBD/issues/21.
* Added support for more control over the gene color scale and in the web
application also added support for reversing the order of the scale.
Resolves https://github.com/LieberInstitute/spatialLIBD/issues/22 and
https://github.com/LieberInstitute/spatialLIBD/issues/23.
* Added `export_cluster()` and `import_cluster()` to help export/import
clustering results instead of having to save large `spe` objects when
exploring different clustering methods.
* Added `locate_images()` and `add_images()` for adding non-standard images
to a `spe` object.

# spatialLIBD 1.7.2

BUG FIXES

* Fixed an issue introduced by newer versions of `shiny`. This version of
`spatialLIBD` works with `shiny` version 1.7.1, though it's likely backwards
compatible. Resolves https://github.com/LieberInstitute/spatialLIBD/issues/24.
* Fix an issue where `as.data.frame(colData(spe))` uses `check.names = TRUE` by
default and then changes the column names unintentionally.

# spatialLIBD 1.7.1

NEW FEATURES

* Added `read10xVisiumWrapper()` and related functions that make it easier
to read in the SpaceRanger output files and launch a shiny web application
using `run_app()`. These new functions read in the analysis output from
SpaceRanger by 10x Genomics, in particular, the clustering and dimension
reduction (projection) results.

# spatialLIBD 1.3.19

SIGNIFICANT USER-VISIBLE CHANGES

* `spatialLIBD` has been updated to work with `SpatialExperiment` version
1.1.701 which will be released as part of Bioconductor 3.13. This changes
internal code of `spatialLIBD` which will work with any objects created with
`SpatialExperiment` version 1.1.700.

# spatialLIBD 1.3.16

SIGNIFICANT USER-VISIBLE CHANGES

* The citation information has changed now that `spatialLIBD` has a bioRxiv
pre-print at https://www.biorxiv.org/content/10.1101/2021.04.29.440149v1.

# spatialLIBD 1.3.15

SIGNIFICANT USER-VISIBLE CHANGES

* We now use `plotly::toWebGL()` to make the web application more responsive.

# spatialLIBD 1.3.14

SIGNIFICANT USER-VISIBLE CHANGES

* The documentation and help messages shown in the web application have been
revamped and improved.

# spatialLIBD 1.3.12

NEW FEATURES

* We added a new vignette that shows how you can use `spatialLIBD` with any
10x Genomics Visium dataset processed with `spaceranger`. The vignette uses
the publicly available human lymph node example from the 10x Genomics website.

# spatialLIBD 1.3.3

NEW FEATURES

* Overall the package has been updated to use `SpatialExperiment` version
1.1.427 available on Bioconductor 3.13 (bioc-devel). Several functions
were re-named such as `sce_image_gene_p()` now has a shorter name
`vis_gene_p()`. This update also changes these visualization functions to ONLY
support `SpatialExperiment` objects instead of the original modified
`SingleCellExperiment` objects.
* Updated citation information to reflect that
https://doi.org/10.1038/s41593-020-00787-0 is now public. Also added a link on
the README to https://doi.org/10.6084/m9.figshare.13623902.v1 for the
manuscript high resolution images.

# spatialLIBD 1.1.7

NEW FEATURES

* The functions `sce_image_gene_p()`, `sce_image_gene()`, `sce_image_grid()`, 
`sce_image_grid_gene()`, `sce_image_clus()`, `sce_image_clus_p()`, 
`geom_spatial()` now work with VisiumExperiment objects thanks to the new 
function `read_image()` and `ve_image_colData()`. This work was done by Brenda Pardo and Leonardo. 


# spatialLIBD 1.1.5

NEW FEATURES

* `fetch_data()` takes the data from sce object and creates a VisiumExperiment
object containing these data thanks to the function `sce_to_ve()`. VisiumExperiment object can be obtained with 
`fetch_data("ve")`. This work was done by Brenda Pardo and Leonardo. 


# spatialLIBD 1.1.4

NEW FEATURES

* `fetch_data()` now uses `BiocFileCache()` when downloading the data from
Dropbox.

# spatialLIBD 0.99.14

SIGNIFICANT USER-VISIBLE CHANGES

* Added the function `enough_ram()` which is used to control the execution of
examples. If it fails when using `fetch_data("sce")` then `fetch_data()` will
show a warning.
* `fetch_data(type = "sce_example")` is now supported and used visibly in the
vignette, eliminating the need for `eval = FALSE` chunks. This should enable
testing the vignette code on the Bioconductor Single Package Builder on
Windows (max 2.5 GB of RAM available).

BUG FIXES

* Fixed the example in `get_colors()`.
* Fixed `layer_stat_cor_plot()` for when `min` and/or `max` are specified.

# spatialLIBD 0.99.13

SIGNIFICANT USER-VISIBLE CHANGES

* Documentation website is now available at
http://LieberInstitute.github.io/spatialLIBD/. It gets updated with every
commit on the master branch (bioc-devel) using GitHub Actions and pkgdown.


# spatialLIBD 0.99.12

BUG FIXES

* Remove the spatialLIBD.Rproj file =( since BioC's SBP is asking me to do so
http://bioconductor.org/spb_reports/spatialLIBD_buildreport_20200303135350.html
* Use `system2()` instead of `system()`.
* Move the `set.seed()` call outside of `layer_boxplot()` as noted by
Martin Morgan
https://github.com/Bioconductor/Contributions/issues/1389#issuecomment-594099852
.
* Use `\linkS4class` as I see being done at
https://github.com/drisso/SingleCellExperiment/search?q=linkS4class&unscoped_q=linkS4class.
* Use `vapply()` instead of `sapply()`.
* Fix (or attempt to) some doc links.



# spatialLIBD 0.99.11

BUG FIXES

* Check if removing the `RcppAnnoy` line in the DESCRIPTION actually works now
based on Aaron Lun's comment at
https://github.com/eddelbuettel/rcppannoy/issues/57#issuecomment-594097241.


# spatialLIBD 0.99.10

SIGNIFICANT USER-VISIBLE CHANGES

* Include AWS links to the image TIFF files (~500mb each) as requested by
Qian Zhu <zqian@jimmy.harvard.edu> for visualizing the data on the Giotto
Viewer https://www.biorxiv.org/content/10.1101/701680v1.

BUG FIXES

* Fix `fetch_data()` and the vignette by specifying the `mode = "wb"` for
`utils::download.file()` in order to resolve an issue with Windows OS reported
here http://bioconductor.org/spb_reports/spatialLIBD_buildreport_20200302120158.html#tokay2_buildsrc_anchor.


# spatialLIBD 0.99.9

SIGNIFICANT USER-VISIBLE CHANGES

* Link to https://doi.org/10.1101/2020.02.28.969931 now that its
public.


# spatialLIBD 0.99.8

SIGNIFICANT USER-VISIBLE CHANGES

* https://spatial.libd.org/spatialLIBD is not supported since we
are using Shiny Server and not Shiny Server Pro. So all links have
now been updated to http://spatial.libd.org/spatialLIBD.


# spatialLIBD 0.99.7

BUG FIXES

* Run a test that might help with https://github.com/r-lib/pkgdown/issues/1230.


# spatialLIBD 0.99.6

SIGNIFICANT USER-VISIBLE CHANGES

* Add mirrors for the shiny app and change the main location.


# spatialLIBD 0.99.5

SIGNIFICANT USER-VISIBLE CHANGES

* Make `fetch_data()` more flexible. Should now work when the data is absent.


# spatialLIBD 0.99.4

BUG FIXES

* Fix Travis badges
* Fix Kristen's name on the vignette
* Add the same welcome information to the top of the vignette, since this will
be what Bioconductor users see first. Basically, we have made sure that users
will see the same information first regardless if they find the package README,
open the shiny app, or find the package vignette.


# spatialLIBD 0.99.3

SIGNIFICANT USER-VISIBLE CHANGES

* Further refine the READMEs (pkg and shiny). They now include the list of
links to the raw 10x Genomics files as well as a short description of the
project at the top. This was in response to feedback by Andrew Jaffe.


# spatialLIBD 0.99.2

SIGNIFICANT USER-VISIBLE CHANGES

* Update main package READMEs to reflect the changes to the shiny web app
README.md.


# spatialLIBD 0.99.1

NEW FEATURES

* Added Kristen R Maynard to the DESCRIPTION file.
* Improved the shiny app page footer.
* Moved around the documentation and added a new main tab with an overview in
response to the feedback by Stephanie Hicks.


# spatialLIBD 0.99.0

NEW FEATURES

* Added a `NEWS.md` file to track changes to the package.
* First full version of the package to be submitted to Bioconductor. Note that
the `ExperimentHub::ExperimentHub()` functionality won't work until they
approve the package. However, for now `fetch_data()` has a backup mechanism
in place.
* Submitted to Bioconductor
[here](https://github.com/Bioconductor/Contributions/issues/1389).
