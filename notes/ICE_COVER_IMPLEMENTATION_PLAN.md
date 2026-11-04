# Time-specific ice cover in PollenPixel downscaling

Status: revision approved and implemented on 4 November 2026; targeted validation passed. Full production recalibration/flow under this revision has not been run. **Section 13 defines static class 140 and ice taking precedence across the full grid.** Active summaries below reflect that revision; original user replies and debugging history are retained. Keep this plan in `notes/`.

## 1\. Intended behavior

For each transition from the previous reconstruction step to the current step:

1. Obtain the ice mask for the **current (destination) age**.
2. Count the previous step's land-cover classes only at pixels that are ice-free at the current age.
3. Calculate class-to-class transport from those restricted frequencies to the current age's predicted land-cover proportions.
4. Apply pixel switches only within that same ice-free set, retaining the previous spatial arrangement wherever no switch is needed.
5. Mark ice-covered pixels in the current output so they are not presented as reconstructed vegetation.

The first example step has no ice, so all otherwise eligible pixels participate. At the second step, newly ice-covered pixels are removed **before** calculating transport, rather than masked only after allocation.

Vegetation allocation includes modeled vegetation classes plus modern **220 (permanent snow and ice)** cells with a vegetation seed, whenever ice-free. **140 (lichens and mosses), 210 (water), bare areas, and no-data cells are static and excluded from vegetation allocation.** The full-grid ice overlay displays 220 over every covered cell, including these static classes and missing background; their original values return on retreat. Modern cropland and urban cells retain their existing initialization transfer to modeled vegetation. This addition does not include time-varying coastlines.

**Class 140 is excluded by policy from calibration, prediction, RDA, seeding, and transport.** Its original pixels remain 140 beneath the ice overlay. **Ice itself is never reconstructed as a predicted class.** Modern 220 cells receive an internal base vegetation state from the closest retained modern vegetation-class centroid; they participate in vegetation reconstruction only when ice-free. Output code 220 remains reserved for the supplied ice mask.

## 2\. Current code and proposed change locations

|Location|Current behavior|Planned change|
|-|-|-|
|`scripts/00_data_prep.R`|Downloads input archive, CHELSA, and ETOPO data.|Download and validate the ice NetCDF independently of whether the `input` directory already exists.|
|`scripts/01_Preprocessing_LS.R`, section 4.1|Builds `initRast` and handles class merging/transfer.|Retain the initialization grid as the spatial reference; expand the modeled pixel domain, seed internal vegetation states for added cells, and prepare/cache aligned ice masks.|
|`scripts/02_PixelFlow_parallel_noOptim_LS.R`, `site_pixel_flow()`, section 4.2|Builds `init_crds`, static allocation costs, and `classFlow`.|Load ice helpers/configuration, map reconstruction ages to ice layers, and retain stable pixel identities.|
|Same script, `old_p`|Counts all rows of the previous `classFlow` column and divides by `nrow(classFlow)`.|Count only current-age eligible pixels and divide by their count.|
|Same script, transport `nrPxl`|Converts transport percentages using `nrow(classFlow)`.|Use the same current-age eligible count used for `old_p`.|
|Same script, `flip_pixels()`|Selects candidates by previous land-cover class.|Require current-age ice-free eligibility as well.|
|Same script, `test_out`|Reports class percentages over the full modeled pixel count.|Compare targets with percentages over current-age ice-free modeled pixels; report ice separately.|
|Same script, section 4.3|Starts each output layer from `initRast` and overwrites modeled pixels.|Apply the full-grid ice mask after allocation, so unchanged background cells cannot show vegetation beneath ice.|
|Same script, saved results and plots|Uses file-existence checks to reuse outputs.|Include mask/configuration provenance and invalidate incompatible cached outputs.|

Prefer a small new helper file, for example `scripts/ice_cover.R`, for metadata inspection, age mapping, and mask alignment. Avoid embedding the NetCDF interpretation throughout the flow loop.

### Implementation style agreed with the user

* Limit new production scripts to the single ice helper file. Integrate setup, preprocessing, transport, and export changes into the existing scripts; avoid adding a new framework or multiple small helper scripts.
* Favor human-readable code with descriptive variable names and short, focused functions. Keep ice eligibility, internal vegetation state, and displayed output visibly distinct.
* Use line breaks and indentation more consistently than the current code: put long function arguments on separate lines, split long pipelines into named intermediate steps, and avoid nested expressions that conceal the processing order.
* Add brief comments explaining modeling choices, such as the destination-age mask, frozen vegetation state, and centroid initialization. Avoid comments that merely repeat the code.
* Improve formatting in the sections being changed; keep unrelated rewrites outside this addition.

Both data preparation and pixel flow clear the global environment with `rm(list = ls())`. Helpers/configuration must be loaded after the relevant reset and made available to parallel workers; defining them only in `main.R` before sourcing those scripts will not suffice.

## 3\. Ice data acquisition and interpretation

Source supplied for this addition:

[TOPicemsk.GLACD26kN9894GE90227A6005GGrBgic.nc](https://sharebox.lsce.ipsl.fr/index.php/s/yfuUw91ruuJXroC/download?path=/&files=TOPicemsk.GLACD26kN9894GE90227A6005GGrBgic.nc)

The supplied file has been downloaded and inspected with `terra`: it contains three subdatasets (`HDC`, `HDCB`, `ICEM`). The selected `ICEM` subdataset has 261 layers, named `ICEM_1` through `ICEM_261`, on a 360-column by 360-row grid spanning longitude 0–360 and latitude -90–90 (1 degree by 0.5 degree resolution). Sampled modern/glacial layers contain values 0 and 1. The GDAL time metadata does not supply the requested reconstruction chronology, so the explicit user-supplied index formula is used. Download size is 811,823,560 bytes; the inspected source MD5 is `e8be8847a4dabd0a446682d51f6d77c3`.

### Download changes

* Suggested local destination: `input/data/global_maps/ice_cover/TOPicemsk.GLACD26kN9894GE90227A6005GGrBgic.nc`.
* Create the destination directory and add a dedicated existence/validity check outside the initial archive-download condition.
* Follow the existing `download.file(..., method = "libcurl", mode = "wb")` pattern, with a sufficient timeout.
* Download to a temporary sibling file; validate it as a readable NetCDF containing the expected `ICEM` data before promoting it to the final filename. Reject empty files and HTML/error responses.
* Reuse valid downloads, recover from incomplete downloads, and preserve the source URL and a computed checksum for reproducibility.
* Use `terra` for new raster handling and add it explicitly to `scripts/libraries.txt` (the existing flow script uses it but it is not listed there). If metadata inspection needs `ncdf4`, add that too. Run dependency-dependent validation after package setup. Preserve existing `stars` interfaces through explicit conversions where needed; a wholesale pipeline migration is outside this addition.

### Metadata inspection before implementation

Establish whether `ICEM...` entries are separate variables, subdatasets, or raster bands; do not assume that a raster reader's band number is the simulation's layer index. Inspect dimensions, coordinate ordering, longitude convention, CRS, fill values, scale/offset, and valid value range. The user specifies a Pacific-centered raster: normalize longitude using `terra::rotate()` as appropriate before alignment. The agreed mask rule is **valid pixel value > 0 means ice**. Decode fill values as missing before applying the threshold.

## 4\. Mapping layer indices to ages

The user supplied the age-to-layer mapping:

```r
source_index <- (26100 - age) / 100
layer_name <- paste0("ICEM_", source_index)
# Inverse: age <- 26100 - 100 * source_index
```

Examples derived from this formula: age 0 maps to `ICEM_261`, age 2000 to `ICEM_241`, and age 20000 to `ICEM_61`. Interpret the formula's ages as years (0ka = 0 years), using the existing reconstruction age reference; verify metadata consistency. The supplied `seq(0, 20000, 2000)` illustrates conversion, rather than changing the pipeline's reconstruction interval. Require integer indices and verify each exact layer exists.

Build an explicit lookup table with these fields:

|Field|Purpose|
|-|-|
|`layer_name`|Exact `ICEM...` identifier in the source.|
|`source_index`|Parsed simulation index, with its indexing convention documented.|
|`age`|Age computed using the supplied formula.|
|`age_unit` / reference|Explicit units and reference convention, compatible with PollenPixel's `Age`.|

Keep the conversion in one function or configuration entry. Validate unique indices/ages and known index-to-age examples supplied with the formula.

Create a second lookup from each reconstruction age to its exact ice layer. The user confirms every reconstruction age has a mask, so no nearest-age selection or interpolation is needed. A missing layer or noninteger derived index is a validation error. Do not silently extrapolate outside the simulation range.

The existing flow uses `Age > 500` in grouped age order and labels the initial map as `-1000`. Use the **0ka mask (`ICEM_261`) for the initial/modern map**, as requested. Preserve the initialization label for compatibility while recording `ice_age = 0` explicitly; do not query a layer at age `-1000`. Make the reconstruction order explicit; “previous” means the previously processed step, which may be younger in a reconstruction proceeding into the past.

## 5\. Spatial mask and pixel identity

For each site and selected age:

1. Read only the required source layer and geographic region, allowing enough margin for reprojection.
2. Decode missing values, normalize Pacific-centered longitude with `terra::rotate()`, and classify valid values > 0 as ice.
3. Use `terra` to align the categorical mask to the complete `initRast` grid, including CRS, extent, resolution, origin, and cell order. Use nearest-neighbor resampling for a categorical mask. Resolve latitude direction explicitly and verify conversions preserve the existing raster cell mapping.
4. Distinguish ice, ice-free, and unknown/outside source coverage. Fail when unknown values affect any cell in the complete exported grid, rather than interpret missing ice data as ice-free.
5. Extract eligibility for modeled rows using `init_crds$pxlID`. Keep the full-grid mask for final map export.

Do not delete rows from the persistent pixel table as ice advances. Allocation indices are row indices into `classFlow`/`costs`; `pxlID` addresses the original raster. Preserve and validate that mapping when subsetting eligible candidates.

Cache aligned masks by source checksum, selected layer, interpretation rule, and grid signature. Share read-only caches between workers; avoid concurrent writes to the same cache path and avoid loading a global all-age stack separately for every worker.

## 6\. Transport and pixel switching

Let `eligible_t` identify modeled pixels ice-free at the destination step, and let `N_t = sum(eligible_t)`.

For class `k`, compute:

```text
old_p[k] = 100 * count(eligible_t AND previous_class == k) / N_t
new_p[k] = predicted proportion of class k at the destination age
```

Working interpretation following feedback: target land-cover proportions describe composition **within ice-free modeled land** and sum to 100% there. Ice reduces the available area; it is not added as another class in the vegetation transport matrix. Class 220 is an output ice code, not a donor or destination in vegetation transport. Do not adjust pollen-derived targets when ice reaches the source location or source area at this stage; log the overlap instead.

Align `old_p`, `new_p`, `lcovs_num`, and both axes of the transport cost matrix by class identifier, rather than relying on column position or sorting. Include zero counts for modeled classes absent from the eligible set. Validate finite, nonnegative targets and equal total transport mass.

Exclude 140 when the scripts load `input/settings/landcover_class_definitions_mod.xlsx`: `ice_class_definitions()` disables calibration and past prediction for 140, keeps it constant, and prevents transfer/merging into another class. The workbook remains unchanged. Existing sample-size filters still apply to the other candidate vegetation classes.

Refit calibration and RDA without 140, and exclude it consistently from target predictions, centroids, overlap costs, and transport. `class_140_status.csv` reports policy exclusion, observed counts, and zero converted pixels. Old 140-enabled models must not be reused. Refitting models separately at every age remains outside this addition.

Convert transport mass to pixel counts using `N_t`, never the full row count. The current `floor()` conversion can leave allocation residuals; assess this with small-domain tests and use a documented integer allocation rule that conserves eligible pixels and respects donor availability. Report any remaining difference from requested targets. Guard zero-count transfers, empty flow tables, and single-candidate cases explicitly.

Illustrative loop (conceptual, not ready-to-run R):

```text
for destination age in reconstruction order:
    ice = aligned mask for destination age
    eligible = modeled pixels known to be ice-free in ice
    previous = previous internal land-cover state
    resolve newly exposed pixels according to the agreed policy

    if count(eligible) == 0:
        skip transport and switching; emit fully ice-covered modeled output
    else:
        old_p = class frequencies of previous[eligible]
        new_p = target composition aligned by class identifier
        flow = transport(old_p, new_p, existing class costs)
        changes = allocate switches among eligible donor pixels only
        current = previous with those changes applied

    save current internal state and destination ice mask
    export land cover with ice explicitly represented
```

### Numerical example

Suppose the previous map contains 100 modeled pixels: 60 forest and 40 grass. At the next step, ice covers 30 forest pixels and 10 grass pixels.

The eligible domain is now 60 pixels: 30 forest and 30 grass. Transport must start from **50% forest / 50% grass**, rather than 60% / 40%. If the target is 25% forest / 75% grass, switch 15 eligible forest pixels to grass. The result is 15 forest, 45 grass, and 40 ice-covered pixels. No ice-covered pixel is a transport donor or a switching candidate.

## 7\. Ice retreat and initialization — frozen state accepted

Advancing ice is covered by the requested behavior. Retreating ice requires a rule for pixels that are now ice-free but had no visible vegetation reconstruction at the previous step.

**Agreed:** retain a separate internal land-cover state beneath the ice, frozen while the pixel is ice-covered. Use that stored class when the pixel becomes ice-free again, including it in `old_p` before transport. Exported maps show class 220 at covered ages. Log newly exposed pixel counts per transition, exposure events, and unique pixels exposed at least once per site. This keeps state complete and supports retreat without adding an ice class to transport, but the stored class is a modeling assumption, not evidence for vegetation beneath ice.

**Agreed initialization for modern ice pixels:** assign a base vegetation class using whichever modeled land-cover class centroid is closest to each ice pixel under modern conditions. Implementation interpretation: use Euclidean distance in the existing RDA1/RDA2 environmental space, consistent with the pipeline's class-centroid representation (`psaOvlp@centers`), rather than geographic distance. Use the stored centroids consistently; do not substitute the density-maximum points computed separately in the flow script.

1. Include modern class 220 cells in the persistent pixel table before filtering by target-class membership; also seed any modeled terrestrial cell covered by the 0ka mask that lacks a valid vegetation state.
2. Extract modern environmental predictors at those pixels and apply the existing preprocessing transformations and RDA projection/scaling used for vegetation allocation.
3. Compute distance to each modeled vegetation-class centroid in that same RDA space. Align centroids to class identifiers explicitly. Assign the nearest class, breaking exact ties by a stable class-ID order. Use centroid distance alone for this seed, rather than the combined density/distance allocation cost.
4. Save this one-time base vegetation assignment with initialization provenance. Freeze it while ice-covered. When exposed, include the stored vegetation class in `old_p` before normal transport; it is then free to change through vegetation allocation. If ice returns, freeze the latest vegetation state rather than resetting to the original seed.
5. Keep class 220 out of vegetation calibration, centroid candidates, predicted target shares, and vegetation transport. The external ice mask alone determines where output class 220 appears, including over water and every other underlying map class.

Modern 140 cells always retain their underlying class. They receive no vegetation seed and do not re-enter allocation when exposed. Only retained modeled vegetation classes are candidates for modern ice seeds; 140 is never a seed candidate.

Apply the mask to the first reconstructed step too. The user's no-ice first step should work naturally when its mask is empty; do not hard-code an assumption that every site's first actual reconstruction age is ice-free.

## 8\. Outputs, diagnostics, and cached runs

* Keep internal allocation state distinct from displayed land cover. Assign **class 220** to ice-covered pixels in exported maps, including the initial map using the 0ka mask. Preserve ordinary missing/outside-domain cells separately. Update legends/class metadata so 220 is retained as ice, rather than merged, transferred, or included in vegetation transport. Record the aligned ice mask alongside state for reproducibility.
* Apply the ice overlay last across the complete exported grid. Ice overrides water, static land, and missing background. Underlying static values are restored when ice retreats. Validate known mask status across the complete grid, and distinguish modeled-domain exposure counts from the additional full-grid counts.
* Save reconstruction ages, selected source layers/ages, source checksum, conversion formula/version, matching policy, mask rule, and grid signature with the result.
* Write a per-site/per-age log file containing modeled pixel total, ice-free count, ice-covered count, unknown count, newly covered/exposed counts, switches, and requested versus realized class proportions over the eligible domain. Also log whether ice covers the pollen source point and the ice-covered fraction of its actual source-area polygon at every age, avoiding an arbitrary threshold for “most of the source area”. Use consistent equal-area geometry for source-area fractions. Show ice separately in maps.
* Replace file-existence-only reuse with a compatibility check for ice-dependent results. Invalidate the flow state, RDA/TIFF output files, comparison plot, and map figures when the mask/configuration changes. Reuse unaffected calibration/preprocessing products where their inputs remain valid.
* If `mergeMaps()` in `scripts/LandCoverChange_ElasticNet.R` is used downstream, ensure merging cannot fill masked ice cells from a neighboring site's vegetation result; carry the ice status through merging or reapply the authoritative mask afterward.
* An explicit ice-disabled mode can support comparison with the old pipeline. Enabled mode must not silently fall back to no ice when data or age mappings are missing.

## 9\. Implementation sequence

* [x] Review this plan and settle the decisions below.
* [x] Record the supplied age/index formula and derived examples.
* [x] Inspect the NetCDF metadata and verify the exact mapped layers.
* [x] Add the setup download/validation and metadata/age lookup helpers.
* [x] Add grid-aligned masks, stable pixel extraction, and cache provenance.
* [x] Integrate eligibility into frequencies, transport counts, and switching; handle empty domains and retreat.
* [x] Update outputs, diagnostics, parallel loading, and cache compatibility checks.
* [x] Validate synthetic allocation cases and real modern/glacial mask alignment.
* [x] Exercise the actual site function on an isolated small crop of existing PSA 206 inputs, including export, diagnostics, reuse, and configuration-triggered rebuilding.
* [x] Run the same site reconstruction in a separate R worker, opening source/mask files in that worker rather than exporting raster pointers.
* [ ] Rerun full production recalibration and flow with static 140 and the full-grid overlay. Existing 140-enabled production results remain on disk until this rerun; the revised policy invalidates their cache fingerprints.
* [x] Document configuration, data provenance, and interpretation in `README.txt`.

Validation: `tests/ice_cover.R` passed with R 4.6.1, including 60 randomized allocation cases, ice advance/retreat and frozen state, fully covered/no-switch/one-pixel cases, deterministic centroid initialization, protected water/no-data, real layer mapping and alignment, source overlap, site exports and logs, cache reuse/invalidation, and worker execution. All R scripts parse successfully. The local test session unset conflicting PostgreSQL/PostGIS `PROJ_LIB`/`GDAL_DATA` environment settings before loading spatial packages; persistent user environment settings were not changed.

## 10\. Acceptance checks

1. **No ice:** reproduce existing behavior when using the original domain/class configuration and matching allocation settings. Test the expanded 140/220 domain separately, since expanding eligibility can change the reconstruction even without ice. Isolate any integer allocation changes when comparing results.
2. **Ice advance:** reproduce the 100-pixel example above, demonstrating restricted frequencies, denominator, and donor selection.
3. **Spatial continuity:** eligible pixels retain previous classes except for explicitly selected transport switches; pixel IDs remain aligned after masking.
4. **No eligible pixels / no required switches:** skip allocation cleanly and still produce valid output and diagnostics.
5. **Retreat:** newly exposed pixels re-enter according to the chosen policy and count toward both frequencies and targets.
6. **Class alignment and small domains:** handle zero-frequency classes, permuted target columns, zero rounded transfers, and one-pixel domains without changing the wrong class or overdrawing donors.
7. **Age handling:** verify numeric index parsing, the supplied formula, processing order, initial-map handling, unmatched ages, and out-of-range ages.
8. **Spatial/data handling:** verify alignment against recognizable mask features; reject unknown coverage and invalid downloads. Mask unchanged background cells in exported maps too.
9. **Caching and parallelism:** existing input directories still acquire the ice dataset; changed source/mapping/mask rules trigger recomputation; repeated and parallel runs use compatible masks.
10. **Class scope:** modern ice cells with vegetation seeds participate when ice-free. Modern 140 and water cells are static, never receive vegetation assignments, display 220 under ice, and restore their original classes on retreat. Bare and no-data cells obey the same overlay/restore rule. Class 140 must be absent from model responses, RDA, and transport regardless of sample size.
11. **Modern ice initialization:** each added modern ice pixel receives the nearest vegetation-class centroid in modern RDA space, with reproducible class alignment/tie handling. A synthetic case where the closest centroid differs from the lowest combined density/distance cost must choose the centroid. On exposure, its stored class enters transport; while covered, it remains frozen internally and displays 220. No ice target proportion or ice centroid is fitted.

## 11\. Original decisions and user replies

User replies are preserved below. Sections above incorporate these decisions; remaining clarifications are in section 12.

* **Index-to-age formula and units/reference:** pending user input. Include whether the index is zero- or one-based and a few expected examples.

  * conversion of ages to indices:
  * ages <- seq(0,20000,2000)
  * indices <- (26100-ages)/100
  * indices <-paste("ICEM",indices,sep="\_")
* **Ice value interpretation / fractional threshold:** pending metadata inspection and review.

  * every pixel value >0 is classified as ice cover
* **Unmatched reconstruction ages:** proposed exact matching with a clear failure; alternative nearest layer with an agreed tolerance.

  * all reconstruction ages will have an icemask available
* **Target proportions:** proposed interpretation is relative to ice-free modeled land. Confirm whether the pollen-derived targets require any additional treatment when ice reaches the pollen source location or most of its source area.

  * let's not treat it right now, but monitor when it happens in a log file
* **Newly exposed pixels:** proposed frozen internal class state; alternatives require an initialization rule.

  * let's use frozen internal state
  * again: monitor how often this is actually the case
* **Ice in exports:** proposed `NA` land cover plus separate ice-status raster; alternative dedicated ice code with legend/downstream updates.

  * please assign class 220 to the ice pixels
* **Domain:** proposed retention of existing modeled pixel eligibility. Confirm whether currently excluded land-cover classes, including modern permanent ice cells, should become reconstructable if the historical ice mask marks them ice-free; that would require expanding initialization and the pixel table.

  * yes, make them reconstructable 
  * Follow-up: 210 should stay constant; try making 140 and 220 reconstructable as well.
  * Clarification: ice itself is not reconstructable. Pixels exposed by retreating ice can be reconstructed with vegetation. Derive their base vegetation from the land-cover class centroid closest to each ice pixel in the modern time slice.
  * Class 140 follow-up: yes, try it as an actual predicted class; it may be rejected because there are too few pixels to reconstruct it.
* **Initial map:** confirm the physical interpretation of the `-1000` initialization label and whether the initial display should carry a separate modern ice mask.

  * there will be a 0ka step in the ice mask that can be used as the modern ice mask
* **Review notes:**
* please use terra for raster handling 
* the ice raster is pacific-centered, terra::rotate can be used to fix that 



## 12. Review decisions and remaining assumptions

Historical initial decision: try 140 as a predicted vegetation class, with water taking precedence over ice. **Superseded by the approved section 13 revision:** 140 is static, and ice overrides every exported cell. Nearest-centroid base vegetation beneath modern ice remains in effect.

No unanswered user decision currently blocks implementation. The former centroid fallback for 140 has been removed; its pixels now stay static. Source metadata and representative exact layer mappings have been verified; every requested site age is also validated at runtime.

### Follow-up: failed two-record production test

Both records (206 and 208) reported incomplete preprocessing at the pixel-flow cache check. Their preprocessing output stopped after the initial vegetation/map files, and the legacy error log contained only missing statuses. The preprocessing script was indexing the two-row status table with metadata rows 12 and 13, so errors were lost. It now matches dataset IDs, writes `landcoverchange/preprocessing_status.csv` after each record with the failing stage and message, and stops before pixel flow if preprocessing failed.

The Earth Engine land-cover image declaration in section 2.2 was commented out even though the subsequent class counts require `dataset`. That declaration has been restored. Keras and TensorFlow have been removed from startup dependencies because the current ElasticNet workflow does not use them; the optional legacy KLD metric retains its on-demand Keras check. The Python executable-discovery and glmnet build-version warnings are separate from the recorded incomplete-preprocessing errors; user interpreter settings have not been changed.

### Follow-up: ice reprojection failed after successful preprocessing

The subsequent production preprocessing completed successfully for both sites. Pixel flow then failed with `[project] Cannot do this transformation`. This was reproduced with the actual PSA 206 grid: the R process inherited `PROJ_LIB` and `GDAL_DATA` pointing at PostgreSQL/PostGIS resources, incompatible with terra's linked projection libraries. Using terra's bundled projection database made the same transformation succeed.

`configure_spatial_runtime()` now configures the current Windows R process and workers to use terra's own PROJ/GDAL resources, without changing persistent Windows settings. It is called during data setup, pixel-flow startup, and worker initialization. Calibration fingerprints now cover class-affecting code independently from ice runtime/allocation helpers; flow fingerprints still include the full helper. The two completed calibration manifests were verified against the old fingerprint and migrated to the equivalent new fingerprint, retaining their valid outputs. Modern and 20ka mask projections and source-point transformations were verified for both actual site grids.

The regression suite then passed with the inherited PostGIS environment left in place, including temporary site reconstruction with the newly calibrated class set, export, cache reuse/invalidation, and worker execution. Full production pixel flow has not been rerun; `ice_flow_status.csv` retains the previous failed-run entries until the next production run.

Working assumptions unless corrected: the supplied ages are in years on the existing reconstruction age reference; the 2000-year example does not change the reconstruction interval; target proportions remain relative to ice-free modeled land. Unknown ice data remain a validation error. Log source-point coverage, source-area ice fraction, and exposure counts at every step; do not adjust targets based on these diagnostics yet.

## 13. Approved correction: static lichens/mosses and an authoritative ice overlay

**Approved by the user on 4 November 2026. Code updated and targeted validation passed. Production models/maps have not been overwritten.**

### A. Required behavior

1. **Remove class 140 from vegetation reconstruction.** Modern lichens/mosses pixels retain their underlying class 140 throughout the reconstruction. They do not enter calibration, target prediction, RDA, vegetation transport, or pixel switching, and they are not converted to another vegetation class when 140 is excluded.
2. **Ice takes precedence over every map class.** At any destination age, every exported grid cell with a valid ice-mask value of 1 displays class 220. This includes water (210), lichens/mosses (140), bare areas, and cells whose underlying map value is 0, 9999, or missing. The exported grid extent stays the same.
3. **Preserve the underlying state beneath ice.** When ice retreats, static classes reappear unchanged: 140 returns to 140, 210 returns to 210, and bare/no-data cells return to their original values. A retreating ice mask does not make static water or lichens/mosses cells eligible for vegetation assignment.
4. **Continue reconstructing vegetation at eligible cells.** Existing modeled vegetation pixels and modern ice pixels with a seeded vegetation state remain eligible whenever ice-free. Their states freeze while ice-covered and resume allocation on exposure. Class 140 is excluded from centroid seed candidates as well as from transport destinations.

The two concepts must be separate in code: **vegetation eligibility** determines which cells can change vegetation, while the **ice overlay** determines the displayed value of every cell.

| Underlying pixel type | Ice-free output | Ice-covered output | After ice retreats | Vegetation allocation? |
| --- | --- | --- | --- | --- |
| Lichens/mosses (140) | 140 | 220 | 140 | No |
| Water (210) | 210 | 220 | 210 | No |
| Bare areas (merged 205) | 205 | 220 | 205 | No |
| No data (0, 9999, or missing) | Original no-data value | 220 | Original no-data value | No |
| Modeled vegetation | Current vegetation class | 220 | Frozen vegetation state enters allocation | Yes, when ice-free |
| Modern ice cell with a vegetation seed | Seed/current vegetation class when ice-free | 220 | Frozen vegetation state enters allocation | Yes, when ice-free |

### B. Changes to existing files

**`scripts/ice_cover.R`**

* In `ice_class_definitions()`, explicitly set 140 to excluded from calibration and past prediction, constant in the past, and not transferable. Keep the source workbook unchanged, following the current in-memory configuration approach.
* Remove 140 from `reconstructable_classes`, leaving modern ice (220) as the extra reconstructable pixel category.
* Separate the configuration for static underlying classes from output masking. A static class must never prevent ice from displaying.
* In `ice_output_values()`, first assemble the underlying map from static values and reconstructed vegetation; then assign 220 wherever the destination mask indicates ice. Remove the current `protected` exception from the ice overlay, including its exceptions for water and missing underlying values.
* Retain the current age formula, ice threshold, nearest-neighbor alignment, and Windows projection-runtime fix.

**`scripts/01_Preprocessing_LS.R`**

* Use the revised class definitions so 140 never becomes a calibration response or an RDA class. Refit the vegetation model and RDA consistently without 140; simply dropping its column from already fitted predictions is insufficient.
* Leave original 140 cells labeled 140 in `initRast`. Do not perform centroid conversion or human-class transfer on those cells.
* Keep `class_140_status.csv`, but record the explicit policy exclusion and static treatment rather than presenting it as rejection because of too few pixels. Existing counts may still be reported as context.

**`scripts/02_PixelFlow_parallel_noOptim_LS.R`**

* Build `cell_ids` from retained vegetation classes plus modern ice cells, excluding all original 140 cells. Never seed or switch static 140 pixels.
* Calculate transport source counts, target pixel counts, and vegetation comparison percentages only over ice-free pixels in this modeled subset. Static lichens/mosses and water are excluded from these denominators.
* Maintain the full original-grid values for export. Apply the full-grid ice overlay after writing reconstructed vegetation into that map, for the modern initialization layer and every subsequent layer.
* Validate that the ice mask is known across the complete exported grid, including static water and background cells. An unknown mask value must not silently behave as ice-free.
* Keep initial modern-ice seeding based on the nearest retained vegetation-class centroid, with no 140 candidate.

**Documentation and tests**

* Update `README.txt` and relevant active plan summaries after approval so they describe static 140 and ice-over-water behavior. Preserve user replies as historical decisions.
* Extend the existing test file. No new production scripts are needed; continue using the existing ice helper and readable, clearly broken-up code.

### C. Diagnostics

* Preserve `newly_exposed_pixels` and `unique_exposed_pixels` as counts for the modeled vegetation domain, and explicitly document that scope. These identify cells re-entering vegetation allocation.
* Add separate full-grid ice coverage/exposure counts so ice retreat over water, lichens/mosses, or other static cells is visible without implying that vegetation was assigned there. Suggested fields: `grid_ice_covered_pixels`, `grid_newly_exposed_pixels`, and `grid_unique_exposed_pixels`.
* Keep the water/ice overlap diagnostic; overlap now correctly displays 220. Add a lichens/mosses overlap count so the static-class override can be inspected directly.
* Report class 140 as excluded by policy in `class_140_status.csv`, with no conversion of its pixels.

### D. Cached data and rerun scope

* Class 140 was retained in both current test records. Removing it therefore requires rebuilding the class-dependent calibration, predictions, pixel sampling/environment normalization where dependent on that sampling, RDA, class costs, and initialization products. Retaining those old models would not implement this correction correctly.
* Use the existing preprocessing fingerprint/invalidation mechanism and update the flow fingerprint/version so old vegetation states, counts, maps, figures, and diagnostics cannot be reused under the new policy.
* Retain raw input downloads, the ice NetCDF, independent vegetation-source preparation, and compatible source-map crops. Reuse aligned binary ice masks where source, age mapping, threshold, and grid compatibility can be verified; rebuild masks if those checks fail.
* After implementation and validation, rerun `main.R` for the two selected test records. Deleting the whole output directory should not be necessary. Do not relabel the current models as compatible with the changed class set.

### E. Acceptance checks before asking for a production rerun

1. Run a 140 pixel through ice-free → ice-covered → ice-free: its outputs must be 140 → 220 → 140, with no vegetation seed, transport participation, or vegetation-switch event.
2. Run the same transition for water: outputs must be 210 → 220 → 210, and the pixel must never become vegetation-eligible.
3. Check bare and no-data cells: ice displays 220, and retreat restores the underlying value, including missing values.
4. Check a fully ice-covered exported grid: every cell must display 220, including cells that were water or missing in the original map. Vegetation transport is skipped when no modeled pixels are ice-free.
5. Verify that class 140 is absent from calibration responses, predicted target columns, RDA classes/centroids, overlap costs, and transport. Confirm no new class-140 pixels are generated by allocation.
6. Verify eligible vegetation counts and realized proportions on a mixed grid containing modeled vegetation, static 140, water, and modern ice cells. Static classes must not affect the vegetation denominator.
7. Confirm exposure logs distinguish modeled vegetation re-entry from full-grid retreat over static cells.
8. Confirm current 140-enabled models and old water-protected maps are rejected by cache checks, while compatible raw inputs are retained.
9. Run the existing allocation, real-mask alignment, temporary site export, and worker tests with the revised policy.

**Approval received:** sections A–E define the revised behavior. The table above records the agreed interpretation of “ice covers everything”.

### F. Implementation and validation result

* Implemented the static-140 policy, full-grid final overlay, strict full-grid mask validation, and separate modeled/full-grid exposure diagnostics in the existing scripts and ice helper.
* Updated calibration/flow fingerprints to reject old class-140 models and maps. Binary-mask fingerprints now depend only on source, layer mapping, threshold/alignment, and grid. A compatibility lookup reuses version-1 binary masks when those inputs match.
* `tests/ice_cover.R` passed: static 140/water/bare/no-data overlay and restoration, fully ice-covered output, modeled eligibility/denominators, old calibration rejection, new and legacy binary-mask reuse, class-140 policy diagnostics, allocation, temporary site exports, and worker execution. The temporary integration fixture refits an RDA without 140 from existing samples and uses synthetic target shares; it does not relabel the old production model as compatible.
* Production files are retained until the user reruns `main.R`. That run must recalibrate without 140 before generating new flow outputs. No manual deletion of the output directory is required.

## 14. Reuse source-overlap geometry (4 November 2026)

Approved performance change: compute source-point/polygon geometry once per site, then reuse it at every age. This does not change ice masking, target proportions, vegetation transport, or the area-weighted overlap definition.

* `prepare_ice_source_geometry()` projects the source point and polygon to the output grid once, records the source-point cell and polygon cell IDs, and computes cell areas times exact polygon boundary fractions. Missing underlying land cover does not remove these geometric weights.
* The pixel-flow loop passes its already loaded logical ice-status vector and this prepared geometry to `ice_source_overlap()`. Each age reads the source-point value and calculates the weighted fraction without reopening the mask, reprojecting geometries, or repeating exact polygon extraction.
* Geometry is local to each site's worker and reused in memory. No new production script or on-disk geometry cache is added. Existing grid/mask validation and missing-source-value checks remain in place.
* This change affects the flow-code fingerprint only. Calibration fingerprints depend on the class configuration and preprocessing code, which are unchanged by this optimization. Existing static-140 calibration products therefore remain compatible; the earlier class-policy revision still requires its own recalibration if not yet run.
* Regression checks compare the new weighted calculation with the original per-age raster calculation for modern, glacial, fully ice-free, and fully ice-covered masks, plus a projected grid. They also cover missing underlying land cover, invalid source masks/locations, and the complete temporary site/worker workflow.

Validation status: all R scripts parse and the regression suite passed, including comparison with the original overlap calculation to tolerance 1e-12, projected/geographic grids, missing-data cases, temporary site exports, cache reuse/invalidation, and separate worker execution. Production outputs were not overwritten by these tests. A currently running R worker retains the functions it loaded at startup; the optimization takes effect on a subsequent run.
