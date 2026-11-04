This repository contains code to conduct the first two steps of the PollenPixel workflow (pollen to land cover and land cover to source area).

The pipeline can be executed by just downloading the repository and running the script main.R. Input data and packages will be downloaded and installed automatically. 

You can adjust parallelization with the nodes variable near the top of scripts/02_PixelFlow_parallel_noOptim_LS.R.

Results will be stored in a newly created directory "landcoverchange".

Time-specific ice cover
----------------------
Data preparation downloads the GLAC ice simulation (approximately 774 MiB) to
input/data/global_maps/ice_cover. Existing input installations also receive this
download. The shared configuration and ice helpers are in scripts/ice_cover.R.

The ICEM layer for an age in years is ICEM_((26100 - age) / 100). The modern
initialization uses age 0 (ICEM_261); its existing output label remains -1000.
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
preprocessing status, stage, and error message. Preprocessing failures stop the
workflow before pixel flow starts. The project-level ice_flow_status.csv records
completed, cached, or failed pixel-flow sites.
Ice reaching a pollen source is logged; it does not currently adjust its targets.

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
