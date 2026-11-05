This repository contains code to conduct the first two steps of the PollenPixel workflow (pollen to land cover and land cover to source area).

The pipeline can be executed by just downloading the repository and running the script main.R. Input data and packages will be downloaded and installed automatically. 

Set NODES in main.R to control workers throughout preprocessing, regression
fitting, and pixel flow. Task-based pools use fewer workers when there are fewer
tasks. When sourcing an individual stage, set NODES first in the same R session.
Legacy nrCores/n_workers function arguments are retained for compatibility;
NODES controls the worker count. Regression fitting creates and closes its own
pool and resets foreach to sequential execution when finished.

Results will be stored in a newly created directory "landcoverchange".

Time-specific ice cover
----------------------
Data preparation downloads the GLAC ice simulation (approximately 774 MiB) to
input/data/global_maps/ice_cover. Existing input installations also receive this
download. The shared configuration and ice helpers are in scripts/ice_cover.R.

The ICEM layer for an age in years is ICEM_((26100 - age) / 100). The modern
initialization uses age 0 (ICEM_261); its existing output label remains -1000.
Both ICEM_1..ICEM_261 names and GDAL's ICEM_T122KP11 coordinate names are
supported. Coordinate names must describe the complete ordered -26..0 axis;
the reconstruction still uses the ordinal age mapping above. A valid existing
.part download is checked and reused on restart instead of downloaded again.
Each requested layer must exist exactly. The pipeline rotates the source's
0-360 degree longitude grid, classifies values greater than zero as ice, and
aligns categorical masks to each site's grid using terra and nearest neighbors.
Missing ice data anywhere in the exported grid cause an error.

At each destination age, transport counts and pixel switches use only pixels
that are ice-free at that age. Vegetation targets describe proportions within
that available modeled land. Integer target counts conserve the eligible pixel
total; rounding residuals are assigned by largest remainder, with class-ID ties.

Modern ice cells receive an internal base vegetation class from the nearest
modern RDA class centroid. This state freezes while covered and re-enters
vegetation allocation when exposed. Exported ice is class 220, prescribed by the
simulation; ice is never a vegetation target. Ice is the final overlay across
the complete grid, overriding water, static land, and missing background cells.
Water, bare areas, no-data cells, and lichens/mosses retain their underlying
values and reappear unchanged when ice retreats. They never enter vegetation
allocation, even when exposed by retreating ice.

Lichens and mosses (140) are excluded by policy when class definitions are
loaded; the workbook itself is not edited. These pixels stay 140 when ice-free
and display 220 when covered. Class 140 is excluded from calibration, prediction,
RDA centroids, and transport. Only modern ice cells receive new base vegetation
seeds, chosen from the retained modeled vegetation classes.

Run main.R to rebuild and execute the complete workflow. On the first run after
this change, incompatible class-dependent preprocessing caches are invalidated.
The pixel-flow script rejects outdated preprocessing rather than reusing it.
Compatible subsequent runs reuse preprocessing, aligned masks, and flow outputs.
The static-140 revision requires recalibration and new flow outputs. Run main.R;
there is no need to delete landcoverchange. Compatible existing binary ice masks
can be reused because their geometry and ice status have not changed.

Per-site outputs include:
  summaryResults/class_140_status.csv: explicit policy exclusion, original
    pixel counts, and confirmation that no 140 pixels were converted.
  summaryResults/ice_cover_log.csv: ice coverage, source-point/source-area
    overlap, newly covered/exposed pixels, unique exposed pixels, and switches.
    newly_exposed_pixels and unique_exposed_pixels refer to modeled vegetation
    cells. Fields prefixed grid_ count the whole exported grid, including water,
    lichens/mosses, bare areas, and missing background. Water and lichen overlap
    counts show where the prescribed ice overlay covers those static classes.
  summaryResults/vegetation_proportions.csv: requested and realized ice-free
    vegetation shares and pixel counts at each age.
  psaChange/psaPixelFlow_new.rda: internal vegetation states, pixel IDs, mask
    paths, age lookup, source checksum, and configuration provenance.
  psaChange/LcovOut_new.rda and LcovOut_new.tif: maps displaying ice as 220.
  iceMasks/: aligned, cached ice-status rasters (0 = ice-free, 1 = ice).
The project-level preprocessing_status.csv records each selected site's
preprocessing status, stage, and error message. The standard run reads all sites
in ala_sites.txt and records sites with missing metadata as skipped. Other
preprocessing failures stop the workflow before pixel flow starts. A retry with
the same preprocessing fingerprint can reuse saved intermediate results. Pixel
flow requires a completion manifest. The project-level ice_flow_status.csv
records completed, cached, or failed pixel-flow sites.

Prediction omits non-finite pollen rows at ages up to 500 years, outside the
pixel-flow target ages. PSA 16238 is capped at 26,000 years for pixel flow. Its
full-age predictions are preserved separately in summaryResults. When an RDA
density grid has zero mass for a class, the grid is extended to include the
calibration scores and evaluated again with the same bandwidths.
Ice reaching a pollen source is logged; it does not currently adjust its targets.
Source-point cell indices, source-polygon boundary fractions, and area weights
are prepared once per site. Each age reuses them with the ice-status vector
already loaded for allocation, avoiding repeated projection, cell-area
calculations, and exact polygon extraction. Source-area overlap remains weighted
by the same cell areas and boundary fractions.

The current ElasticNet workflow does not require Keras or TensorFlow at startup.
The optional legacy KLD metric checks for Keras only if that metric is requested.

Validation
----------
From the repository root, run Rscript tests/ice_cover.R. Tests cover masked
transport, frozen states, integer allocation, centroid seeds, static classes,
and raster rotation. When the downloaded simulation exists, they also check real
ICEM layers and alignment. If the existing PSA 206 preprocessing files are
available, a tiny temporary crop exercises site allocation, export, and caching;
production site results are not overwritten. This fixture reuses existing
training samples to refit a test RDA without 140, using synthetic target shares.
It does not replace full production recalibration under the revised policy.

On Windows, configure_spatial_runtime() selects terra's bundled projection data
for the current R process and its workers. This prevents inherited PostgreSQL/
PostGIS PROJ_LIB or GDAL_DATA settings from causing projection failures. It does
not change persistent Windows settings or the PostgreSQL installation. Other
platforms retain their existing spatial runtime configuration.

The implementation decisions and review history are in
notes/ICE_COVER_IMPLEMENTATION_PLAN.md.
