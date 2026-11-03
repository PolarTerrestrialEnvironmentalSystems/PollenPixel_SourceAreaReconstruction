# Time-specific ice cover in PollenPixel downscaling

Status: draft for review. This document proposes implementation; it does not change the running pipeline.

## 1. Intended behavior

For each transition from the previous reconstruction step to the current step:

1. Obtain the ice mask for the **current (destination) age**.
2. Count the previous step's land-cover classes only at pixels that are ice-free at the current age.
3. Calculate class-to-class transport from those restricted frequencies to the current age's predicted land-cover proportions.
4. Apply pixel switches only within that same ice-free set, retaining the previous spatial arrangement wherever no switch is needed.
5. Mark ice-covered pixels in the current output so they are not presented as reconstructed vegetation.

The first example step has no ice, so all otherwise eligible pixels participate. At the second step, newly ice-covered pixels are removed **before** calculating transport, rather than masked only after allocation.

“All pixels” here means the existing reconstruction domain eligible for modeled land-cover allocation. Existing exclusions, such as water or classes outside the modeled class set, remain separate from ice eligibility. Whether ice should also change the eligible domain itself is a review item below.

## 2. Current code and proposed change locations

| Location | Current behavior | Planned change |
| --- | --- | --- |
| `scripts/00_data_prep.R` | Downloads input archive, CHELSA, and ETOPO data. | Download and validate the ice NetCDF independently of whether the `input` directory already exists. |
| `scripts/01_Preprocessing_LS.R`, section 4.1 | Builds `initRast` and handles class merging/transfer. | Retain the initialization grid as the spatial reference; optionally prepare/cache aligned ice masks once the grid exists. |
| `scripts/02_PixelFlow_parallel_noOptim_LS.R`, `site_pixel_flow()`, section 4.2 | Builds `init_crds`, static allocation costs, and `classFlow`. | Load ice helpers/configuration, map reconstruction ages to ice layers, and retain stable pixel identities. |
| Same script, `old_p` | Counts all rows of the previous `classFlow` column and divides by `nrow(classFlow)`. | Count only current-age eligible pixels and divide by their count. |
| Same script, transport `nrPxl` | Converts transport percentages using `nrow(classFlow)`. | Use the same current-age eligible count used for `old_p`. |
| Same script, `flip_pixels()` | Selects candidates by previous land-cover class. | Require current-age ice-free eligibility as well. |
| Same script, `test_out` | Reports class percentages over the full modeled pixel count. | Compare targets with percentages over current-age ice-free modeled pixels; report ice separately. |
| Same script, section 4.3 | Starts each output layer from `initRast` and overwrites modeled pixels. | Apply the full-grid ice mask after allocation, so unchanged background cells cannot show vegetation beneath ice. |
| Same script, saved results and plots | Uses file-existence checks to reuse outputs. | Include mask/configuration provenance and invalidate incompatible cached outputs. |

Prefer a small new helper file, for example `scripts/ice_cover.R`, for metadata inspection, age mapping, and mask alignment. Avoid embedding the NetCDF interpretation throughout the flow loop.

Both data preparation and pixel flow clear the global environment with `rm(list = ls())`. Helpers/configuration must be loaded after the relevant reset and made available to parallel workers; defining them only in `main.R` before sourcing those scripts will not suffice.

## 3. Ice data acquisition and interpretation

Source supplied for this addition:

[TOPicemsk.GLACD26kN9894GE90227A6005GGrBgic.nc](https://sharebox.lsce.ipsl.fr/index.php/s/yfuUw91ruuJXroC/download?path=%2F&files=TOPicemsk.GLACD26kN9894GE90227A6005GGrBgic.nc)

Relevant layers have the prefix `ICEM`, as specified by the user. The remote file could not be inspected through the web tool during planning. Variable structure, coordinate system, value meanings, spatial coverage, and chronology are therefore **not yet verified**.

### Download changes

- Suggested local destination: `input/data/global_maps/ice_cover/TOPicemsk.GLACD26kN9894GE90227A6005GGrBgic.nc`.
- Create the destination directory and add a dedicated existence/validity check outside the initial archive-download condition.
- Follow the existing `download.file(..., method = "libcurl", mode = "wb")` pattern, with a sufficient timeout.
- Download to a temporary sibling file; validate it as a readable NetCDF containing the expected `ICEM` data before promoting it to the final filename. Reject empty files and HTML/error responses.
- Reuse valid downloads, recover from incomplete downloads, and preserve the source URL and a computed checksum for reproducibility.
- If metadata inspection needs `ncdf4`, add it to `scripts/libraries.txt` and run that validation after package setup. Reuse `stars`/`terra` for raster handling where practical; verify NetCDF support in the installed raster stack.

### Metadata inspection before implementation

Establish whether `ICEM...` entries are separate variables, subdatasets, or raster bands; do not assume that a raster reader's band number is the simulation's layer index. Inspect dimensions, coordinate ordering, longitude convention, CRS, fill values, scale/offset, and valid value range. Determine whether ice is binary, fractional, or encoded another way before defining the mask rule.

## 4. Mapping layer indices to ages

**Pending user input: the formula converting the ICEM index to reconstruction age.** Do not infer it from the filename, lexical ordering, or an assumed interval.

Build an explicit lookup table with these fields:

| Field | Purpose |
| --- | --- |
| `layer_name` | Exact `ICEM...` identifier in the source. |
| `source_index` | Parsed simulation index, with its indexing convention documented. |
| `age` | Age computed using the supplied formula. |
| `age_unit` / reference | Explicit units and reference convention, compatible with PollenPixel's `Age`. |

Keep the conversion in one function or configuration entry. Validate unique indices/ages and known index-to-age examples supplied with the formula.

Create a second lookup from each reconstruction age to its selected ice layer. Prefer exact matches where possible. For unmatched ages, the proposed default is to stop with a clear message until an explicit matching policy is agreed. Nearest-layer matching with a documented maximum age difference is a possible alternative; record the selected source age and mismatch. Do not silently extrapolate outside the simulation range or interpolate categorical masks.

The existing flow uses `Age > 500` in grouped age order and labels the initial map as `-1000`. Treat that initial label as an initialization convention until its physical meaning is confirmed, rather than automatically querying an ice layer at age `-1000`. Make the reconstruction order explicit; “previous” means the previously processed step, which may be younger in a reconstruction proceeding into the past.

## 5. Spatial mask and pixel identity

For each site and selected age:

1. Read only the required source layer and geographic region, allowing enough margin for reprojection.
2. Interpret source values using verified metadata. If fractional, agree an ice threshold before converting to a binary mask.
3. Align the categorical mask to the complete `initRast` grid, including CRS, extent, resolution, origin, and cell order. Use nearest-neighbor resampling for a categorical mask. Resolve longitude wrapping and latitude direction explicitly.
4. Distinguish ice, ice-free, and unknown/outside source coverage. Proposed default: fail when unknown values affect the reconstruction domain, rather than interpret missing ice data as ice-free.
5. Extract eligibility for modeled rows using `init_crds$pxlID`. Keep the full-grid mask for final map export.

Do not delete rows from the persistent pixel table as ice advances. Allocation indices are row indices into `classFlow`/`costs`; `pxlID` addresses the original raster. Preserve and validate that mapping when subsetting eligible candidates.

Cache aligned masks by source checksum, selected layer, interpretation rule, and grid signature. Share read-only caches between workers; avoid concurrent writes to the same cache path and avoid loading a global all-age stack separately for every worker.

## 6. Transport and pixel switching

Let `eligible_t` identify modeled pixels ice-free at the destination step, and let `N_t = sum(eligible_t)`.

For class `k`, compute:

```text
old_p[k] = 100 * count(eligible_t AND previous_class == k) / N_t
new_p[k] = predicted proportion of class k at the destination age
```

Proposed interpretation: target land-cover proportions describe composition **within ice-free modeled land** and sum to 100% there. Ice reduces the available area; it is not added as another class in the vegetation transport matrix. Confirm this interpretation before implementation.

Align `old_p`, `new_p`, `lcovs_num`, and both axes of the transport cost matrix by class identifier, rather than relying on column position or sorting. Include zero counts for modeled classes absent from the eligible set. Validate finite, nonnegative targets and equal total transport mass.

Keep the existing class overlap cost matrix and environmental pixel costs for the first implementation. The mask changes counts and candidate eligibility; refitting the vegetation model or RDA at each age is outside this addition.

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

## 7. Ice retreat and initialization — proposed policy for review

Advancing ice is covered by the requested behavior. Retreating ice requires a rule for pixels that are now ice-free but had no visible vegetation reconstruction at the previous step.

**Proposal:** retain a separate internal land-cover state beneath the ice, frozen while the pixel is ice-covered. Use that stored class when the pixel becomes ice-free again, including it in `old_p` before transport. Exported maps still show ice at covered ages. This keeps state complete and supports retreat without adding an ice class to transport, but the stored class is a modeling assumption, not evidence for vegetation beneath ice.

Alternatives to discuss: initialize newly exposed pixels from existing environmental suitability costs or neighboring ice-free pixels. These need a specified initialization rule and should not silently discard newly exposed pixels from the target domain.

Apply the mask to the first reconstructed step too. The user's no-ice first step should work naturally when its mask is empty; do not hard-code an assumption that every site's first actual reconstruction age is ice-free.

## 8. Outputs, diagnostics, and cached runs

- Keep internal allocation state distinct from displayed land cover. Proposed output: mask land cover to `NA` beneath ice and save a separate ice-status raster with the same grid and age labels. A dedicated ice display code is an alternative if it fits the class definitions and downstream consumers.
- Apply the full-grid mask during export, including cells copied unchanged from `initRast`. Preserve the distinction between ice and ordinary missing/outside-domain cells in the status raster and figures.
- Save reconstruction ages, selected source layers/ages, source checksum, conversion formula/version, matching policy, mask rule, and grid signature with the result.
- Per age, report modeled pixel total, ice-free count, ice-covered count, unknown count, newly covered/exposed counts, switches, and requested versus realized class proportions over the eligible domain. Show ice separately in maps.
- Replace file-existence-only reuse with a compatibility check for ice-dependent results. Invalidate the flow state, RDA/TIFF output files, comparison plot, and map figures when the mask/configuration changes. Reuse unaffected calibration/preprocessing products where their inputs remain valid.
- If `mergeMaps()` in `scripts/LandCoverChange_ElasticNet.R` is used downstream, ensure merging cannot fill masked ice cells from a neighboring site's vegetation result; carry the ice status through merging or reapply the authoritative mask afterward.
- An explicit ice-disabled mode can support comparison with the old pipeline. Enabled mode must not silently fall back to no ice when data or age mappings are missing.

## 9. Implementation sequence

- [ ] Review this plan and settle the decisions below.
- [ ] Obtain the index-to-age formula and representative expected mappings; inspect the NetCDF metadata.
- [ ] Add the setup download/validation and metadata/age lookup helpers.
- [ ] Add grid-aligned masks, stable pixel extraction, and cache provenance.
- [ ] Integrate eligibility into frequencies, transport counts, and switching; handle empty domains and retreat.
- [ ] Update outputs, diagnostics, parallel loading, and cache compatibility checks.
- [ ] Validate synthetic examples, then run a representative affected site and an unaffected site.
- [ ] Document configuration, data provenance, and interpretation in `README.txt`.

## 10. Acceptance checks

1. **No ice:** with matching inputs and allocation settings, reproduce existing behavior. If integer allocation is deliberately improved, isolate and document that change when comparing results.
2. **Ice advance:** reproduce the 100-pixel example above, demonstrating restricted frequencies, denominator, and donor selection.
3. **Spatial continuity:** eligible pixels retain previous classes except for explicitly selected transport switches; pixel IDs remain aligned after masking.
4. **No eligible pixels / no required switches:** skip allocation cleanly and still produce valid output and diagnostics.
5. **Retreat:** newly exposed pixels re-enter according to the chosen policy and count toward both frequencies and targets.
6. **Class alignment and small domains:** handle zero-frequency classes, permuted target columns, zero rounded transfers, and one-pixel domains without changing the wrong class or overdrawing donors.
7. **Age handling:** verify numeric index parsing, the supplied formula, processing order, initial-map handling, unmatched ages, and out-of-range ages.
8. **Spatial/data handling:** verify alignment against recognizable mask features; reject unknown coverage and invalid downloads. Mask unchanged background cells in exported maps too.
9. **Caching and parallelism:** existing input directories still acquire the ice dataset; changed source/mapping/mask rules trigger recomputation; repeated and parallel runs use compatible masks.

## 11. Decisions and feedback

Edit these entries directly or add comments beneath them.

- **Index-to-age formula and units/reference:** pending user input. Include whether the index is zero- or one-based and a few expected examples.
- **Ice value interpretation / fractional threshold:** pending metadata inspection and review.
- **Unmatched reconstruction ages:** proposed exact matching with a clear failure; alternative nearest layer with an agreed tolerance.
- **Target proportions:** proposed interpretation is relative to ice-free modeled land. Confirm whether the pollen-derived targets require any additional treatment when ice reaches the pollen source location or most of its source area.
- **Newly exposed pixels:** proposed frozen internal class state; alternatives require an initialization rule.
- **Ice in exports:** proposed `NA` land cover plus separate ice-status raster; alternative dedicated ice code with legend/downstream updates.
- **Domain:** proposed retention of existing modeled pixel eligibility. Confirm whether currently excluded land-cover classes, including modern permanent ice cells, should become reconstructable if the historical ice mask marks them ice-free; that would require expanding initialization and the pixel table.
- **Initial map:** confirm the physical interpretation of the `-1000` initialization label and whether the initial display should carry a separate modern ice mask.
- **Review notes:**

