#!/usr/bin/env python3

import click


@click.group(invoke_without_command=True)
@click.option("-v", "--verbosity", count=True, help="Verbosity level")
@click.option(
    "--mjd",
    default=None,
    type=int,
    help="Modified Julian date to query. Use negative values to indicate relative to current MJD",
)
@click.option("--mjd-start", default=None, type=int, help="Start of MJD range to query")
@click.option("--mjd-end", default=None, type=int, help="End of MJD range to query")
@click.option("--date", default=None, type=str, help="Date to query (e.g., 2024-01-15)")
@click.option(
    "--date-start", default=None, type=str, help="Start of date range to query"
)
@click.option("--date-end", default=None, type=str, help="End of date range to query")
@click.option("--apo", is_flag=True, help="Query Apache Point Observatory data")
@click.option("--lco", is_flag=True, help="Query Las Campanas Observatory data")
@click.option(
    "--fibers", "--fibres", is_flag=True, help="Include fibre mappings to targets"
)
@click.option(
    "--no-x-match", is_flag=True, help="Do not cross-match targets with SDSS-V database"
)
@click.option("--output", "-O", default=None, type=str, help="Output file")
@click.option(
    "--processes", "-p", default=None, type=int, help="Number of processes to use"
)
@click.pass_context
def main(
    ctx,
    verbosity,
    mjd,
    mjd_start,
    mjd_end,
    date,
    date_start,
    date_end,
    apo,
    lco,
    fibers,
    no_x_match,
    output,
    processes,
):
    """
    Almanac collects metadata from planned and actual APOGEE exposures,
    and identifies sequences of exposures that constitute epoch visits.
    """

    # This keeps the default behaviour as 'query mode' but allows for commands like 'config'.
    if ctx.invoked_subcommand is not None:
        return  # Let Click handle the subcommand

    import h5py as h5
    from itertools import product
    from rich.live import Live
    from almanac.display import ObservationsDisplay, display_exposures
    from almanac import apogee, logger, io, utils
    from contextlib import nullcontext
    from time import time, sleep

    mjds, mjd_min, mjd_max = utils.parse_mjds(
        mjd, mjd_start, mjd_end, date, date_start, date_end
    )
    observatories = utils.get_observatories(apo, lco)

    n_iterables = len(mjds) * len(observatories)
    iterable = product(mjds, observatories)
    results = []

    display = ObservationsDisplay(mjd_min, mjd_max, observatories)

    buffered_critical_logs = []
    buffered_result_rows = []

    refresh_per_second = 1
    context_manager = (
        Live(
            display.create_display(), refresh_per_second=refresh_per_second, screen=True
        )
        if verbosity >= 1
        else nullcontext()
    )
    io_kwds = dict(fibers=fibers, compression=False)
    with h5.File(output, "a") if output else nullcontext() as fp:
        with context_manager as live:
            if processes is not None:

                def initializer():
                    from sdssdb.peewee.sdss5db import database

                    if hasattr(database, "_state"):
                        database._state.closed = True
                        database._state.conn = None
                    from almanac.database import database

                # Parallel
                import os
                import signal
                import concurrent.futures

                if processes < 0:
                    processes = os.cpu_count()
                with concurrent.futures.ProcessPoolExecutor(
                    max_workers=processes, initializer=initializer
                ) as pool:

                    try:
                        futures = set()
                        for n, (mjd, observatory) in enumerate(iterable, start=1):
                            futures.add(
                                pool.submit(
                                    apogee.get_almanac_data,
                                    observatory,
                                    mjd,
                                    fibers,
                                    not no_x_match,
                                )
                            )
                            if n == processes:
                                break

                        t = time()
                        while len(futures) > 0:

                            future = next(concurrent.futures.as_completed(futures))

                            observatory, mjd, exposures, sequences = result = (
                                future.result()
                            )

                            v = mjd - mjd_min + display.offset
                            missing = [e.image_type == "missing" for e in exposures]
                            if any(missing):
                                display.missing.add(v)
                                # buffered_critical_logs.extend(missing)

                            if not exposures:
                                display.no_data[observatory].add(v)
                            else:
                                display.completed[observatory].add(v)
                                results.append(result)
                                if output:
                                    io.update(
                                        fp,
                                        observatory,
                                        mjd,
                                        exposures,
                                        sequences,
                                        **io_kwds,
                                    )

                            if (
                                live is not None
                                and (time() - t) > 1 / refresh_per_second
                            ):
                                live.update(display.create_display())
                                t = time()
                            futures.remove(future)

                            try:
                                mjd, observatory = next(iterable)
                            except StopIteration:
                                None
                            else:
                                futures.add(
                                    pool.submit(
                                        apogee.get_almanac_data,
                                        observatory,
                                        mjd,
                                        fibers,
                                        not no_x_match,
                                    )
                                )

                    except KeyboardInterrupt:
                        for pid in pool._processes:
                            os.kill(pid, signal.SIGKILL)
                        pool.shutdown(wait=False, cancel_futures=True)
                        try:
                            fp.close()
                        except:
                            None
                        raise KeyboardInterrupt
            else:
                t = time()
                for mjd, observatory in iterable:
                    *_, exposures, sequences = result = apogee.get_almanac_data(
                        observatory, mjd, fibers, not no_x_match
                    )
                    v = mjd - mjd_min + display.offset
                    if any([e.image_type == "missing" for e in exposures]):
                        display.missing.add(v)
                        # buffered_critical_logs.extend(missing)

                    if not exposures:
                        display.no_data[observatory].add(v)
                    else:
                        display.completed[observatory].add(v)
                        results.append(result)
                        if output:
                            io.update(
                                fp, observatory, mjd, exposures, sequences, **io_kwds
                            )

                    if live is not None and (time() - t) > 1 / refresh_per_second:
                        live.update(display.create_display())
                        t = time()

            if live is not None:
                live.update(display.create_display())
                if verbosity <= 1 and output is None:
                    sleep(3)

    if verbosity >= 2:
        for observatory, mjd, exposures, sequences in results:
            display_exposures(exposures, sequences)

        # Show critical logs at the end to avoid disrupting the display
        for item in buffered_critical_logs:
            logger.critical(item)


@main.command()
@click.argument("identifiers", type=int, nargs=-1)
@click.option("--output", "-O", default=None, type=str, help="Output file")
def lookup(identifiers, output, **kwargs):
    """Lookup target(s) by catalog or SDSS identifier."""

    if not identifiers:
        return

    import h5py as h5
    from contextlib import nullcontext
    from peewee import fn
    from itertools import chain, starmap
    from almanac import io
    from almanac.database import database
    from almanac.apogee import get_exposures, get_almanac_data
    from sdssdb.peewee.sdss5db.targetdb import (
        Assignment,
        AssignmentStatus,
        CartonToTarget,
        Target,
        Hole,
        Observatory,
        Design,
        DesignToField,
        Field,
    )
    from sdssdb.peewee.sdss5db.catalogdb import SDSS_ID_flat

    from rich.table import Table as RichTable
    from rich.console import Console
    from rich.live import Live

    sq = (
        SDSS_ID_flat.select(SDSS_ID_flat.sdss_id)
        .where(
            (SDSS_ID_flat.sdss_id.in_(identifiers))
            | (SDSS_ID_flat.catalogid.in_(identifiers))
        )
        .alias("sq")
    )
    q = (
        SDSS_ID_flat.select(SDSS_ID_flat.catalogid, SDSS_ID_flat.sdss_id)
        .distinct()
        .join(sq, on=(SDSS_ID_flat.sdss_id == sq.c.sdss_id))
        .order_by(SDSS_ID_flat.sdss_id.asc())
        .tuples()
    )

    catalogids = {catalogid: sdss_id for catalogid, sdss_id in q}

    if not catalogids:
        raise click.ClickException(
            f"Identifiers {identifiers} not found in SDSS-V database"
        )

    # todo hacky
    for m in (
        Target,
        Observatory,
        AssignmentStatus,
        Field,
        CartonToTarget,
        Assignment,
        Hole,
        Design,
        DesignToField,
    ):
        try:
            m.get()
        except:
            print(f"Failed {m}")
        else:
            print(f"Success {m}")

        m._meta.schema = "targetdb"

    q = (
        Target.select(
            fn.Lower(Observatory.label),
            AssignmentStatus.mjd,
            Field.field_id,
        )
        .distinct()
        .join(CartonToTarget)
        .join(Assignment)
        .join(AssignmentStatus)
        .switch(Assignment)
        .join(Hole)
        .join(Observatory)
        .switch(Assignment)
        .join(Design)
        .join(DesignToField)
        .join(Field)
        .where(
            Target.catalogid.in_(tuple(catalogids.keys()))
            & (AssignmentStatus.status == 1)
        )
        .tuples()
    )
    n, fields = (0, {})
    for obs, mjd, field_id in q:
        mjd = int(mjd)
        fields.setdefault((obs, mjd), set())
        fields[(obs, mjd)].add(field_id)
        n += 1

    console = Console()

    rich_table = RichTable(
        title=f"Exposures\n({n} obs/mjd/field combinations)",
        title_style="bold blue",
        show_header=True,
        header_style="bold cyan",
    )

    for field_name in (
        "#",
        "obs",
        "mjd",
        "exposure",
        "field",
        "fiber_id",
        "catalogid",
        "sdss_id",
    ):
        rich_table.add_column(field_name, justify="center")

    done = set()
    io_kwds = dict(fibers=True, compression=False)
    with h5.File(output, "a") if output else nullcontext() as fp:

        i = 1
        with Live(rich_table, console=console, refresh_per_second=4) as live:
            for exposure in chain(*starmap(get_exposures, fields.keys())):
                key = (exposure.observatory, exposure.mjd)
                field_ids = fields[key]

                if output and key not in done:
                    done.add(key)
                    r = get_almanac_data(*key, fibers=True, meta=True)
                    io.update(fp, *r, **io_kwds)

                if exposure.field_id in field_ids:
                    for target in exposure.targets:
                        if target.catalogid in catalogids:
                            rich_table.add_row(
                                *list(
                                    map(
                                        str,
                                        (
                                            i,
                                            exposure.observatory,
                                            exposure.mjd,
                                            exposure.exposure,
                                            exposure.field_id,
                                            target.fiber_id,
                                            target.catalogid,
                                            catalogids[target.catalogid],
                                        ),
                                    )
                                )
                            )
                            i += 1
                            break
    if output:
        console.print(f"Updated {output_path} with:")
        for obs, mjd in done:
            console.print(f"  - {obs}/{mjd}")


@main.group()
def add(**kwargs):
    """Add new information to an existing Almanac file."""
    pass


def _get_sdss_ids(fp, obs, mjd):
    group = fp.get(f"{obs}/{mjd}/fibers", [])
    sdss_ids = set()
    for config in group:
        sdss_ids.update(group[config]["sdss_id"][:])
    return sdss_ids


@add.command()
@click.argument("input_path", type=str)
@click.option(
    "--mjd",
    default=None,
    type=int,
    help="Modified Julian date to query. Use negative values to indicate relative to current MJD",
)
@click.option("--mjd-start", default=None, type=int, help="Start of MJD range to query")
@click.option("--mjd-end", default=None, type=int, help="End of MJD range to query")
@click.option("--date", default=None, type=str, help="Date to query (e.g., 2024-01-15)")
@click.option(
    "--date-start", default=None, type=str, help="Start of date range to query"
)
@click.option("--date-end", default=None, type=str, help="End of date range to query")
@click.option("--apo", is_flag=True, help="Query Apache Point Observatory data")
@click.option("--lco", is_flag=True, help="Query Las Campanas Observatory data")
@click.option("--p", default=-1, type=int, help="Number of workers to use")
def metadata(
    input_path,
    mjd,
    mjd_start,
    mjd_end,
    date,
    date_start,
    date_end,
    apo,
    lco,
    p,
    **kwargs,
):
    """Add astrometry and photometry to an existing Almanac file."""

    import os
    import h5py as h5
    import concurrent.futures
    from itertools import product
    from almanac import utils
    from almanac.catalog import query
    from tqdm import tqdm

    if p <= 0:
        p += os.cpu_count()

    observatories = utils.get_observatories(apo, lco)
    mjds, *_ = utils.parse_mjds(
        mjd, mjd_start, mjd_end, date, date_start, date_end, return_nones=True
    )
    sdss_ids = set()
    with h5.File(input_path, "r") as fp:
        if mjds is None:
            mjds = []
            for obs in observatories:
                mjds.extend(fp[obs])
            mjds = list(set(mjds))

        with concurrent.futures.ThreadPoolExecutor(max_workers=p) as executor:
            futures = [
                executor.submit(_get_sdss_ids, fp, o, m)
                for o, m in product(observatories, mjds)
            ]
            for future in tqdm(
                concurrent.futures.as_completed(futures),
                desc="Collecting SDSS identifiers",
            ):
                sdss_ids.update(future.result())

    from almanac.data_models.source import Source

    results = query(sdss_ids)
    import pickle

    with open("20251128-meta.pkl", "wb") as fp:
        pickle.dump(results, fp)


def _postprocess_chunk_worker(args):
    """
    Worker function for parallel postprocessing computations.

    This must be at module level for multiprocessing to pickle it.
    """
    import traceback
    import warnings
    import numpy as np

    task_type, chunk_indices, chunk_args = args

    # Suppress IERS auto-download warnings about polar motion data
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*IERS.*")
        warnings.filterwarnings("ignore", message=".*polar motion.*")

        try:
            if task_type == "bary":
                from almanac.postprocess import compute_barycentric_correction_chunk

                result = compute_barycentric_correction_chunk(*chunk_args)
                return (task_type, chunk_indices, result, None)
            elif task_type == "obs_meta":
                from almanac.postprocess import compute_observation_metadata_chunk

                result = compute_observation_metadata_chunk(*chunk_args)
                return (task_type, chunk_indices, result, None)
            elif task_type == "shadow":
                from almanac.postprocess import compute_shadow_heights_chunk

                result = compute_shadow_heights_chunk(*chunk_args)
                return (task_type, chunk_indices, result, None)
            else:
                return (task_type, chunk_indices, None, f"Unknown task type: {task_type}")
        except Exception as e:
            # Return NaN results on error so processing can continue
            n = len(chunk_indices)
            if task_type == "bary":
                result = np.full(n, np.nan)
            elif task_type == "obs_meta":
                result = {
                    "moon_phase": np.full(n, np.nan),
                    "moon_separation": np.full(n, np.nan),
                    "airmass": np.full(n, np.nan),
                    "alt": np.full(n, np.nan),
                    "az": np.full(n, np.nan),
                }
            elif task_type == "shadow":
                result = np.full(n, np.nan)
            else:
                result = None
            error_msg = f"{task_type}: {e}\n{traceback.format_exc()}"
            return (task_type, chunk_indices, result, error_msg)


@main.command()
@click.argument("input_path", type=str)
@click.argument("output_prefix", type=str)
@click.option(
    "--processes", "-p", default=None, type=int, help="Number of processes to use"
)
@click.option(
    "--limit", default=None, type=int, help="Limit number of spectra (for testing)"
)
def postprocess(input_path, output_prefix, processes, limit, **kwargs):
    """Post-process an existing Almanac file after reductions are complete."""

    if not input_path or not output_prefix:
        return

    import os

    if not os.path.exists(input_path):
        raise FileNotFoundError(f"The file {input_path} does not exist.")

    import h5py as h5
    import concurrent.futures
    import numpy as np
    from glob import glob
    from collections import defaultdict
    from pydantic_core import PydanticUndefined

    from almanac.display import TaskDisplay
    from almanac.data_models.source import Source
    from almanac.data_models.spectrum import Spectrum
    from almanac.io import get_hdf5_dtype, _write_models_to_hdf5_group
    from almanac.postprocess import (
        finalize_radial_velocities,
        compute_targeting_flags,
        load_almanac_file,
        load_armadgics_files,
        load_metadata_pickle,
        load_ar1d_unical_meta_batch,
        group_indices_by_keys,
        group_indices_by_array,
    )

    if processes is None or processes < 0:
        processes = os.cpu_count()

    outdir = os.path.dirname(os.path.abspath(input_path)) + "/../"
    metadata_pickle_path = "/mnt/home/acasey/almanac/20251128-meta.pkl"
    output_spectra_path = f"{output_prefix}-spectra.h5"
    output_stars_path = f"{output_prefix}-stars.h5"

    data_dict = {}

    # Find other paths
    full_list_info_path = glob(f"{outdir}/arMADGICS/raw_*/full_list_info.h5")[0]
    arMADGICS_suffix = full_list_info_path.split("/")[-2][4:]
    arMADGICS_dir = f"{outdir}/arMADGICS/wu_th_{arMADGICS_suffix}/"
    arMADGICS_paths = glob(f"{arMADGICS_dir}/arMADGICS_out_*.h5")

    with TaskDisplay() as display:

        # Add all tasks upfront
        tid_io_parallel = display.add_task(
            "Loading data in parallel (almanac, arMADGICS, metadata)", total=4
        )

        # === PARALLEL I/O PHASE ===
        # Start background I/O tasks while ar1Dunical collection runs in ProcessPool
        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as io_executor:
            # Submit I/O-bound tasks to run in parallel

            total = None
            with h5.File(full_list_info_path, "r") as full_list_info:
                keys = list(full_list_info.keys())
                tid_summary = display.add_task(
                    f"Reading {limit or len(full_list_info[keys[0]]):,} summary spectra",
                    total=len(keys),
                )
                for key in keys:
                    data_dict[key] = full_list_info[key][:]
                    display.advance(tid_summary)

            if limit is not None:
                data_dict = {k: v[:limit] for k, v in data_dict.items()}
                data_dict["mjd"] = np.array(data_dict["mjd"], dtype=int)

                almanac_future = io_executor.submit(
                    load_almanac_file,
                    input_path,
                    np.min(data_dict["mjd"]),
                    np.max(data_dict["mjd"]),
                )
            else:
                almanac_future = io_executor.submit(load_almanac_file, input_path)

            total = len(data_dict["sdss_id"])

            armadgics_future = io_executor.submit(load_armadgics_files, arMADGICS_paths)
            metadata_future = io_executor.submit(
                load_metadata_pickle, metadata_pickle_path
            )
            display.advance(tid_io_parallel)  # to trigger something to show

            star_fields = {**Source.model_fields, **Source.model_computed_fields}
            spectrum_fields = {
                **Spectrum.model_fields,
                **Spectrum.model_computed_fields,
            }

            tid_fidgeting = display.add_task("Fidgeting", total=6)

            data_dict["carton_pks"] = [None] * total

            # Group lookups by MJD for batched ar1Dunical processing (vectorized)
            data_dict["tele"] = np.array(data_dict["tele"], dtype=str)
            data_dict["mjd"] = np.array(data_dict["mjd"], dtype=int)
            data_dict["expnum"] = np.array(data_dict["expnum"], dtype=int)
            display.advance(tid_fidgeting)

            # Get unique (tele, mjd, expnum) combinations using numpy
            tele_arr = data_dict["tele"]
            mjd_arr = data_dict["mjd"]
            expnum_arr = data_dict["expnum"]

            # Stack and find unique rows
            stacked = np.column_stack(
                [
                    np.char.encode(tele_arr.astype(str), "utf-8"),
                    mjd_arr.astype("i8").view("S8"),
                    expnum_arr.astype("i8").view("S8"),
                ]
            )
            display.advance(tid_fidgeting)
            _, unique_idx = np.unique(stacked, axis=0, return_index=True)
            display.advance(tid_fidgeting)

            # Extract unique values
            unique_tele = tele_arr[unique_idx]
            unique_mjd = mjd_arr[unique_idx]
            unique_expnum = expnum_arr[unique_idx]
            display.advance(tid_fidgeting)

            # Build lookups set and lookups_by_mjd dict
            lookups = set(zip(unique_tele, unique_mjd, unique_expnum))
            lookups_by_mjd = defaultdict(list)
            for t, m, e in zip(unique_tele, unique_mjd, unique_expnum):
                lookups_by_mjd[m].append((t, e))
            display.advance(tid_fidgeting)

            n_mjds = len(lookups_by_mjd)
            n_exposures = len(lookups)
            tid_reduce_1d = display.add_task(
                f"Collecting ar1Dunical data ({n_exposures:,} exposures)",
                total=n_exposures + 1,
            )
            tid_prop_1d = display.add_task("Propagating ar1Dunical data", total=total)
            tid_fiber_lookups = display.add_task("Processing fiber lookups", total=total)
            tid_metadata = display.add_task("Assigning metadata to rows", total=5)
            tid_chunking = display.add_task("Chunking", total=total)
            tid_vrad = display.add_task("Computing radial velocities", total=1)

            tid_shadow_height = display.add_task("Computing shadow heights", total=total)
            tid_obs_meta = display.add_task(
                "Computing observation metadata (moon, airmass, alt/az)",
                total=1,
            )  # Updated later
            tid_targeting = display.add_task("Aggregating targeting cartons", total=total)
            tid_write_spectra = display.add_task(
                f"Write spectra to {os.path.basename(output_spectra_path)}",
                total=len(spectrum_fields),
            )
            tid_star_unique = display.add_task("Identifying unique stars", total=1)
            tid_write_stars = display.add_task(
                f"Write stars to {os.path.basename(output_stars_path)}",
                total=len(star_fields),
            )

            batch_args = [
                (outdir, mjd, obs_exps) for mjd, obs_exps in lookups_by_mjd.items()
            ]
            display.advance(tid_fidgeting)


            # Meanwhile, run ar1Dunical collection (CPU-bound with ProcessPool)
            display.advance(tid_reduce_1d) # to show progress while the pool maps
            ar1D_unical_meta = {}
            with concurrent.futures.ProcessPoolExecutor(
                max_workers=processes
            ) as executor:
                for batch_result in executor.map(
                    load_ar1d_unical_meta_batch, batch_args
                ):
                    ar1D_unical_meta.update(batch_result)
                    display.advance(tid_reduce_1d, len(batch_result))

            # Collect parallel I/O results
            data_dict.update(armadgics_future.result())
            display.advance(tid_io_parallel)
            source_meta = metadata_future.result()
            display.advance(tid_io_parallel)

            if limit is not None:
                data_dict = {k: v[:limit] for k, v in data_dict.items()}

            # Collect all metadata keys from source_meta (vectorized)
            sdss_ids_arr = np.asarray(data_dict["sdss_id"])
            unique_sdss_ids = np.unique(sdss_ids_arr)
            meta_keys = set()
            for sid in unique_sdss_ids:
                meta_keys.update(source_meta.get(sid, {}).keys())

            # Vectorized grouping: sdss_id -> indices (using numpy)
            sdss_id_to_indices = group_indices_by_array(sdss_ids_arr)

            # Build lookup: key (field name or alias) -> field_info for Spectrum model
            #spectrum_field_lookup = {}
            #for field_name, field_info in Spectrum.model_fields.items():
            #    spectrum_field_lookup[field_name] = field_info
            #    if field_info.alias:
            #        spectrum_field_lookup[field_info.alias] = field_info

            # Vectorized ar1D_unical_meta propagation using numpy
            tele_arr = np.array(data_dict["tele"], dtype=str)
            mjd_arr = np.array(data_dict["mjd"], dtype=int)
            expnum_arr = np.asarray(data_dict["expnum"])

            # Pre-allocate numpy arrays for ar1D metadata
            first_key = next(iter(lookups))
            ar1d_keys = list(ar1D_unical_meta[first_key].keys())
            for key in ar1d_keys:
                sample_val = ar1D_unical_meta[first_key][key]
                dtype = (
                    np.array(sample_val).dtype
                    if not isinstance(sample_val, bool)
                    else bool
                )
                data_dict[key] = np.empty(total, dtype=dtype)

            # Group rows by (tele, mjd, exp) for vectorized propagation (numpy)
            lookup_to_indices = group_indices_by_keys(tele_arr, mjd_arr, expnum_arr)

            display.tasks[tid_prop_1d].total = len(lookup_to_indices)
            for lookup_key, indices in lookup_to_indices.items():
                for k, v in ar1D_unical_meta[lookup_key].items():
                    data_dict[k][indices] = v
                display.advance(tid_prop_1d)
            # del ar1D_unical_meta, lookup_to_indices

            # Vectorized almanac lookups using numpy advanced indexing
            adj_fiber_arr = np.asarray(data_dict["adjfiberindx"])

            # Group rows by (obs, mjd) for batch processing (numpy)
            #obs_mjd_to_indices = group_indices_by_keys(tele_arr, mjd_arr)

            almanac_exposures, almanac_fibers, obs_mjd_keys = almanac_future.result()
            display.advance(tid_io_parallel)

            # Collect all field names we'll encounter (for pre-allocation)
            all_exp_keys = set()
            all_fiber_keys = set()
            for exp_data in almanac_exposures.values():
                all_exp_keys.update(exp_data.keys())
            for fiber_dict in almanac_fibers.values():
                for fiber_data in fiber_dict.values():
                    all_fiber_keys.update(fiber_data.keys())

            # Pre-allocate numpy arrays for almanac data
            for key in all_exp_keys:
                if key not in data_dict:
                    sample = next(iter(almanac_exposures.values()))[key]
                    data_dict[key] = np.empty(total, dtype=sample.dtype)

            for key in all_fiber_keys:
                if key not in data_dict:
                    for fiber_dict in almanac_fibers.values():
                        if fiber_dict:
                            sample = next(iter(fiber_dict.values())).get(key, None)
                            if sample is not None:
                                data_dict[key] = np.empty(total, dtype=sample.dtype)
                                break

            # Vectorized batch processing per (obs, mjd)
            # Update total now that we know it (couldn't know until almanac was loaded)
            display.tasks[tid_fiber_lookups].total = len(lookup_to_indices)
            for (obs, mjd, exposure), index in lookup_to_indices.items():

                exp_data = almanac_exposures[(obs, mjd)]

                for k, v in exp_data.items():
                    data_dict[k][index] = v[exposure - 1]

                reference_id = max(
                    exp_data["config_id"][exposure - 1],
                    exp_data["plate_id"][exposure - 1]
                )

                fiber_dict = almanac_fibers[(obs, mjd)][reference_id]
                fiber_index = adj_fiber_arr[index] - 1
                if fiber_index >= 300:
                    fiber_index -= 300

                for k, v in fiber_dict.items():
                    data_dict[k][index] = v[fiber_index]

                display.advance(tid_fiber_lookups)

            # Free memory
            #del almanac_exposures, almanac_fibers, obs_mjd_to_indices, tele_arr, mjd_arr

            # Vectorized metadata assignment (source_meta already loaded in parallel)
            # Build a mapping from sdss_id to row indices
            translations = [
                ("tele", "observatory"),
                ("expnum", "exposure"),
                ("starscale", "starscale0"),
                ("RV_flag", "v_rad_flags"),
                ("cartid", "cart_id"),
                ("nSkyFibers", "n_sky_fibers"),
                ("adjfiberindx", "adjusted_fiber_index"),
            ]
            for from_key, to_key in translations:
                data_dict[to_key] = data_dict.pop(from_key)
            display.advance(tid_metadata)

            for key, field_info in spectrum_fields.items():
                if (
                    key not in data_dict
                and field_info.alias is not None
                and field_info.alias in data_dict
                ):
                    data_dict[key] = data_dict.pop(field_info.alias)

            # keep this duplicate
            data_dict["adjfiberindx"] = data_dict["adjusted_fiber_index"]
            # Pre-allocate arrays for metadata fields
            for key, field_info in spectrum_fields.items():
                if key not in data_dict:
                    # Skip keys not accepted by Spectrum model (as field or alias)

                    hdf5_dtype = get_hdf5_dtype(field_info.annotation)

                    # Get default from model, with fallback for required fields
                    default = field_info.default
                    if default is PydanticUndefined:
                        if np.issubdtype(np.dtype(hdf5_dtype), np.floating):
                            default = np.nan
                        elif np.issubdtype(np.dtype(hdf5_dtype), np.integer):
                            default = -1
                        elif np.dtype(hdf5_dtype).char == "S":
                            default = ""
                        else:
                            default = 0

                    data_dict[key] = np.full(total, default, dtype=hdf5_dtype)
            display.advance(tid_metadata)

            # Vectorized metadata assignment using numpy inverse mapping
            # Instead of looping over sdss_ids and then over keys, we:
            # 1. Build an inverse mapping from row index to unique sdss_id index
            # 2. For each metadata key, build an array indexed by unique sdss_id
            # 3. Use fancy indexing to assign all rows at once

            # Get inverse mapping: row_idx -> unique_sdss_id_idx
            unique_sids, inverse_idx = np.unique(sdss_ids_arr, return_inverse=True)
            n_unique = len(unique_sids)
            display.advance(tid_metadata)

            # Build metadata arrays indexed by unique sdss_id
            # First pass: identify which keys we can vectorize
            vectorizable_keys = [
                k for k in meta_keys if k in data_dict and k != "carton_pks"
            ]

            for key in vectorizable_keys:
                # Build array of values for unique sdss_ids
                dtype = data_dict[key].dtype
                # Get default value from existing array
                default_val = data_dict[key][0] if len(data_dict[key]) > 0 else 0

                meta_values = np.full(n_unique, default_val, dtype=dtype)
                for i, sid in enumerate(unique_sids):
                    meta = source_meta.get(sid, {})
                    if meta and key in meta and meta[key] is not None:
                        meta_values[i] = meta[key]

                # Assign to all rows using inverse mapping
                data_dict[key] = meta_values[inverse_idx]
            display.advance(tid_metadata)

            # Handle carton_pks separately (list assignment)
            for sid, indices in sdss_id_to_indices.items():
                meta = source_meta.get(sid, {})
                if meta and "carton_pks" in meta and meta["carton_pks"] is not None:
                    for index in indices:
                        data_dict["carton_pks"][index] = meta["carton_pks"]

            display.advance(tid_metadata)


            # Compute targeting flags (fast, no parallelization needed)
            carton_pks_list = data_dict.pop("carton_pks", [])
            flags_array, n_unknown_carton_assignments = compute_targeting_flags(
                carton_pks_list, total
            )
            data_dict["sdss5_target_flags"] = flags_array
            display.advance(tid_targeting, total)

            # Do an intermediate write.
            def callback(name):
                def inner(*args, **kwargs):
                    try:
                        display.advance(name)
                    except:
                        pass

                return inner

            '''
            with h5.File(output_spectra_path, "w", track_order=True) as fp:
                _write_models_to_hdf5_group(
                    spectrum_fields, data_dict, fp, callback=callback(tid_write_spectra)
                )
            '''


            # === PARALLEL COMPUTATION PHASE ===
            # Use a single ProcessPoolExecutor for all parallel computations.
            # This avoids the overhead of creating/destroying process pools multiple times.

            # Pre-compute shared arrays used by multiple parallel tasks
            ra, dec, mjd_mid_exposure = map(
                np.array, (data_dict["ra"], data_dict["dec"], data_dict["mjd_mid_exposure"])
            )
            jd_mid_exposure = mjd_mid_exposure + 2400000.5
            observatory_arr = np.array(data_dict["observatory"]).astype(str)
            RV_pixoff_final = np.array(data_dict["RV_pixoff_final"])
            RV_pix_var = np.array(data_dict["RV_pix_var"])

            # Initialize output arrays
            v_bary_corr = np.full(total, np.nan)
            moon_phase = np.full(total, np.nan)
            moon_separation = np.full(total, np.nan)
            airmass = np.full(total, np.nan)
            alt_arr = np.full(total, np.nan)
            az_arr = np.full(total, np.nan)
            shadow_height = np.full(total, np.nan)

            # Create validity mask for coordinates
            valid_mask = (
                np.isfinite(ra)
                & np.isfinite(dec)
                & np.isfinite(mjd_mid_exposure)
                & (ra >= 0)
                & (ra <= 360)
                & (dec >= -90)
                & (dec <= 90)
            )

            # Build chunks for parallel processing
            chunk_size = 5000
            chunks = []  # List of (task_type, chunk_indices, args)
            n_chunks = 0
            chunk_indices = {}
            for obs in ("apo", "lco"):
                obs_mask = valid_mask & (observatory_arr == obs)
                if not np.any(obs_mask):
                    continue

                obs_indices = np.where(obs_mask)[0]
                chunk_indices[obs] = len(obs_indices)

            display.tasks[tid_chunking].total = sum(chunk_indices.values())
            for obs in ("apo", "lco"):
                obs_indices = chunk_indices.get(obs, None)
                if obs_indices is None:
                    continue

                # Split into chunks
                for chunk_start in range(0, len(obs_indices), chunk_size):
                    chunk_end = min(chunk_start + chunk_size, len(obs_indices))
                    chunk_indices = obs_indices[chunk_start:chunk_end]

                    # Barycentric correction chunk
                    chunks.append(
                        (
                            "bary",
                            chunk_indices,
                            (
                                ra[chunk_indices],
                                dec[chunk_indices],
                                mjd_mid_exposure[chunk_indices],
                                obs,
                            ),
                        )
                    )

                    # Observation metadata chunk (moon phase, separation, airmass, alt/az)
                    chunks.append(
                        (
                            "obs_meta",
                            chunk_indices,
                            (
                                ra[chunk_indices],
                                dec[chunk_indices],
                                mjd_mid_exposure[chunk_indices],
                                obs,
                            ),
                        )
                    )

                    # Shadow height chunk (uses JD, not MJD)
                    chunks.append(
                        (
                            "shadow",
                            chunk_indices,
                            (
                                ra[chunk_indices],
                                dec[chunk_indices],
                                jd_mid_exposure[chunk_indices],
                                obs,
                            ),
                        )
                    )
                    display.advance(tid_chunking, len(chunk_indices))

            # Update task totals now that we know chunk count
            n_chunks = len(chunks) // 3  # Each obs has 3 chunk types
            display.tasks[tid_vrad].total = n_chunks
            display.tasks[tid_obs_meta].total = n_chunks
            display.tasks[tid_shadow_height].total = n_chunks

            # Process all chunks in parallel using module-level worker function
            errors = []
            with concurrent.futures.ProcessPoolExecutor(max_workers=processes) as executor:
                futures = [executor.submit(_postprocess_chunk_worker, chunk) for chunk in chunks]
                for future in concurrent.futures.as_completed(futures):
                    try:
                        task_type, chunk_indices, result, error = future.result()
                        if error:
                            errors.append(error)
                        if task_type == "bary":
                            v_bary_corr[chunk_indices] = result
                            display.advance(tid_vrad)
                        elif task_type == "obs_meta":
                            moon_phase[chunk_indices] = result["moon_phase"]
                            moon_separation[chunk_indices] = result["moon_separation"]
                            airmass[chunk_indices] = result["airmass"]
                            alt_arr[chunk_indices] = result["alt"]
                            az_arr[chunk_indices] = result["az"]
                            display.advance(tid_obs_meta)
                        elif task_type == "shadow":
                            shadow_height[chunk_indices] = result
                            display.advance(tid_shadow_height)
                    except Exception as e:
                        errors.append(f"Future exception: {e}")

            if errors:
                click.echo(f"Warning: {len(errors)} errors during parallel computation")

            # Finalize radial velocities (fast, no parallelization needed)
            rv_result = finalize_radial_velocities(RV_pixoff_final, RV_pix_var, v_bary_corr)
            data_dict.update(rv_result)

            # Store observation metadata
            data_dict["moon_phase"] = moon_phase
            data_dict["moon_separation"] = moon_separation
            data_dict["airmass"] = airmass
            data_dict["alt"] = alt_arr
            data_dict["az"] = az_arr
            data_dict["shadow_height"] = shadow_height

            def callback(name):
                def inner(*args, **kwargs):
                    try:
                        display.advance(name)
                    except:
                        pass

                return inner

            with h5.File(output_spectra_path, "w", track_order=True) as fp:
                _write_models_to_hdf5_group(
                    spectrum_fields, data_dict, fp, callback=callback(tid_write_spectra)
                )

            # Now write per source
            display.tasks[tid_star_unique].total = len(data_dict)
            _, indices = np.unique(data_dict["sdss_id"], return_index=True)
            star_dict = {}
            for k, v in data_dict.items():
                star_dict[k] = v[indices]
                display.advance(tid_star_unique)

            with h5.File(output_stars_path, "w", track_order=True) as fp:
                _write_models_to_hdf5_group(
                    star_fields, star_dict, fp, callback=callback(tid_write_stars)
                )

    if n_unknown_carton_assignments > 0:
        click.echo(
            f"Warning: {n_unknown_carton_assignments:,} unknown carton assignments encountered"
        )

    # Last check point for fields
    expected = set(spectrum_fields)
    actual = set(data_dict.keys())
    for key in actual - expected:
        click.echo(f"Warning: ignored field {key} in data_dict")
    for key in expected - actual:
        click.echo(f"Warning: missing field {key}")


@main.group()
def config(**kwargs):
    """View or update configuration settings."""
    pass


@config.command()
def show(**kwargs):
    """Show all configuration settings"""

    from almanac import config, get_config_path
    from dataclasses import asdict

    click.echo(f"Configuration path: {get_config_path()}")
    click.echo(f"Configuration:")

    def _pretty_print(config_dict, indent=""):
        for k, v in config_dict.items():
            if isinstance(v, dict):
                click.echo(f"{indent}{k}:")
                _pretty_print(v, indent=indent + "  ")
            else:
                click.echo(f"{indent}{k}: {v}")

    _pretty_print(asdict(config), "  ")


@config.command
@click.argument("key", type=str)
def get(key, **kwargs):
    """Get a configuration value"""

    from almanac import config
    from dataclasses import asdict

    def traverse(config, key, provenance=None, sep="."):
        parent, *child = key.split(sep, 1)
        try:
            # TODO: Should we even allow dicts in config?
            if isinstance(config, dict):
                v = config[parent]
            else:
                v = getattr(config, parent)
        except (AttributeError, KeyError):
            context = sep.join(provenance or [])
            if context:
                context = f" within '{context}'"

            if not isinstance(config, dict):
                config = asdict(config)

            raise click.ClickException(
                f"No configuration key '{parent}'{context}. "
                f"Available{context}: {', '.join(config.keys())}"
            )

        provenance = (provenance or []) + [parent]
        return traverse(v, child[0], provenance) if child else v

    value = traverse(config, key)
    click.echo(value)


@config.command(hidden=True)
@click.argument("key")
@click.argument("value")
def update(key, value, **kwargs):
    """Update a configuration value"""
    click.echo(click.style("Deprecated: use `almanac config set`", fg="yellow"))
    return set(key, value, **kwargs)


@config.command(name="set")
@click.argument("key")
@click.argument("value")
def _set(key, value, **kwargs):
    """Set a configuration value"""

    from almanac import config, get_config_path, ConfigManager
    from dataclasses import asdict, is_dataclass

    def traverse(config, key, value, provenance=None, sep="."):
        parent, *child = key.split(sep, 1)

        try:
            scope = getattr(config, parent)
        except AttributeError:
            context = sep.join(provenance or [])
            if context:
                context = f" within '{context}'"

            if not isinstance(config, dict):
                config = asdict(config)

            raise click.ClickException(
                f"No configuration key '{parent}'{context}. "
                f"Available{context}: {', '.join(config.keys())}"
            )

        else:

            if not child:

                fields = {f.name: f.type for f in config.__dataclass_fields__.values()}
                field_type = fields[parent]
                if is_dataclass(field_type):
                    context = sep.join(provenance or [])
                    if context:
                        context = f" within '{context}'"

                    raise click.ClickException(
                        f"Key '{parent}'{context} refers to a configuration class. "
                        f"You must set the values of the configuration class individually. "
                        f"Sorry! "
                        f"Or you can directly edit the configuration file {get_config_path()}"
                    )

                setattr(config, parent, value)
            else:
                provenance = (provenance or []) + [parent]
                traverse(scope, child[0], value)

    traverse(config, key, value)
    config_path = get_config_path()
    ConfigManager.save(config, config_path)
    click.echo(f"Updated configuration {key} to {value} in {config_path}")


@main.group()
def dump(**kwargs):
    """Dump data to a summary file"""
    pass


# almanac dump star[s] almanac.h5 output.fits
def check_paths_and_format(input_path, output_path, given_format, overwrite):
    import os
    import click

    if not os.path.exists(input_path):
        raise click.ClickException(f"Input path {input_path} does not exist")

    if os.path.exists(output_path) and not overwrite:
        raise click.ClickException(
            f"Output path {output_path} already exists. Use --overwrite to overwrite."
        )

    if given_format is None:
        if output_path.lower().endswith(".fits"):
            return "fits"
        elif output_path.lower().endswith(".csv"):
            return "csv"
        elif output_path.lower().endswith(".hdf5") or output_path.lower().endswith(
            ".h5"
        ):
            return "hdf5"
        else:
            raise click.ClickException(
                "Cannot infer output format from output path. Please specify --format"
            )
    return given_format


@dump.command()
@click.argument("input_path", type=str)
@click.argument("output_path", type=str)
@click.option(
    "--format",
    "-f",
    default=None,
    type=click.Choice(["fits", "csv", "hdf5"]),
    help="Output format",
)
@click.option("--overwrite", is_flag=True, help="Overwrite existing output file")
def stars(input_path, output_path, overwrite, format, **kwargs):
    """Create a star-level summary file"""

    import h5py as h5
    from copy import deepcopy
    from collections import Counter

    stars = {}
    default = dict(
        mjds_apo=set(),
        mjds_lco=set(),
        n_visits=0,
        n_visits_apo=0,
        n_visits_lco=0,
        n_exposures=0,
        n_exposures_apo=0,
        n_exposures_lco=0,
    )

    output_format = check_paths_and_format(input_path, output_path, format, overwrite)
    assert format != "hdf5", "HDF5 output not yet supported for star summaries."
    with h5.File(input_path, "r") as fp:
        for observatory in fp:
            for mjd in fp[f"{observatory}"]:
                group = fp[f"{observatory}/{mjd}"]

                is_object = group["exposures/image_type"][:].astype(str) == "object"
                fps = is_object * (group["exposures/config_id"][:] > 0)
                plate = is_object * (group["exposures/plate_id"][:] > 0)

                if not any(fps) and not any(plate) or "fibers" not in group:
                    continue

                # fps era
                n_exposures_on_this_mjd = {}

                if any(fps):
                    config_ids = Counter(group["exposures/config_id"][:][fps])
                elif any(plate):
                    config_ids = Counter(group["exposures/plate_id"][:][plate])
                else:
                    continue

                for config_id, n_exposures in config_ids.items():
                    try:
                        config_group = group[f"fibers/{config_id}"]
                    except KeyError:
                        print(
                            f"Warning couldnt get config {config_id} for {observatory} {mjd}"
                        )
                        continue

                    ok = (
                        (config_group["catalogid"][:] > 0)
                        | (config_group["sdss_id"][:] > 0)
                        | (config_group["twomass_designation"][:].astype(str) != "")
                    ) * (
                        (config_group["category"][:].astype(str) == "science")
                        | (config_group["category"][:].astype(str) == "standard_apogee")
                        | (config_group["category"][:].astype(str) == "standard_boss")
                        | (config_group["category"][:].astype(str) == "open_fiber")
                    )
                    sdss_ids = config_group["sdss_id"][:][ok]
                    catalogids = config_group["catalogid"][:][ok]
                    for sdss_id, catalogid in zip(sdss_ids, catalogids):
                        stars.setdefault(sdss_id, deepcopy(default))
                        stars[sdss_id].setdefault(
                            "catalogid", catalogid
                        )  # this can change over time,... should we track that/
                        n_exposures_on_this_mjd.setdefault(sdss_id, 0)
                        n_exposures_on_this_mjd[sdss_id] += n_exposures

                for sdss_id, n_exposures in n_exposures_on_this_mjd.items():
                    stars[sdss_id]["n_exposures"] += n_exposures
                    stars[sdss_id][f"n_exposures_{observatory}"] += n_exposures
                    stars[sdss_id]["n_visits"] += 1
                    stars[sdss_id][f"n_visits_{observatory}"] += 1
                    stars[sdss_id][f"mjds_{observatory}"].add(int(mjd))

        rows = []
        for sdss_id, meta in stars.items():
            stars[sdss_id].update(
                mjd_min_apo=min(meta["mjds_apo"]) if meta["mjds_apo"] else -1,
                mjd_max_apo=max(meta["mjds_apo"]) if meta["mjds_apo"] else -1,
                mjd_min_lco=min(meta["mjds_lco"]) if meta["mjds_lco"] else -1,
                mjd_max_lco=max(meta["mjds_lco"]) if meta["mjds_lco"] else -1,
            )
            stars[sdss_id].pop("mjds_apo")
            stars[sdss_id].pop("mjds_lco")
            rows.append(dict(sdss_id=sdss_id, **meta))

    from astropy.table import Table

    t = Table(rows=rows)
    t.write(output_path, format=output_format, overwrite=overwrite)


@dump.command()
@click.argument("input_path", type=str)
@click.argument("output_path", type=str)
@click.option(
    "--format",
    "-f",
    default=None,
    type=click.Choice(["fits", "csv", "hdf5"]),
    help="Output format",
)
@click.option("--overwrite", is_flag=True, help="Overwrite existing output file")
def visits(input_path, output_path, format, overwrite, **kwargs):
    """Create a visit-level summary file"""

    pass


@dump.command()
@click.argument("input_path", type=str)
@click.argument("output_path", type=str)
@click.option(
    "--format",
    "-f",
    default=None,
    type=click.Choice(["fits", "csv", "hdf5"]),
    help="Output format",
)
@click.option("--overwrite", is_flag=True, help="Overwrite existing output file")
def exposures(input_path, output_path, format, overwrite, **kwargs):
    """Create an exposure-level summary file"""

    import os
    import h5py as h5
    import numpy as np

    output_format = check_paths_and_format(input_path, output_path, format, overwrite)

    from almanac.data_models import Exposure

    fields = {**Exposure.model_fields, **Exposure.model_computed_fields}
    data = dict()
    for field_name, field_spec in fields.items():
        data[field_name] = []

    with h5.File(input_path, "r") as fp:
        for observatory in ("apo", "lco"):
            for mjd in fp[observatory].keys():
                group = fp[f"{observatory}/{mjd}/exposures"]
                for key in group.keys():
                    data[key].extend(group[key][:])

    if output_format == "hdf5":
        from almanac.io import _write_models_to_hdf5_group

        fields = {**Exposure.model_fields, **Exposure.model_computed_fields}

        with h5.File(output_path, "w", track_order=True) as fp:
            _write_models_to_hdf5_group(fields, data, fp)
    else:
        from astropy.table import Table

        t = Table(data=data)
        t.write(output_path, format=output_format, overwrite=overwrite)


@dump.command()
@click.argument("input_path", type=str)
@click.argument("output_path", type=str)
@click.option(
    "--format",
    "-f",
    default=None,
    type=click.Choice(["fits", "csv", "hdf5"]),
    help="Output format",
)
@click.option("--overwrite", is_flag=True, help="Overwrite existing output file")
def fibers(input_path, output_path, format, overwrite, **kwargs):
    """Create a fiber-level summary file"""

    import os
    import h5py as h5
    import numpy as np

    output_format = check_paths_and_format(input_path, output_path, format, overwrite)

    from almanac.data_models.fps import FPSTarget
    from almanac.data_models.plate import PlateTarget

    fields = {
        **FPSTarget.model_fields,
        **FPSTarget.model_computed_fields,
        **PlateTarget.model_fields,
        **PlateTarget.model_computed_fields,
    }

    defaults = {
        name: spec.default for name, spec in fields.items() if hasattr(spec, "default")
    }
    defaults["twomass_designation"] = ""

    data = dict()
    for field_name, field_spec in fields.items():
        data[field_name] = []

    with h5.File(input_path, "r") as fp:
        for observatory in ("apo", "lco"):
            for mjd in fp[observatory].keys():
                group = fp[f"{observatory}/{mjd}/fibers"]
                for config_id in group.keys():
                    group = fp[f"{observatory}/{mjd}/fibers/{config_id}"]
                    n = len(group["sdss_id"][:])

                    for field_name in data:
                        if field_name in group.keys():
                            data[field_name].extend(group[field_name][:])
                        else:
                            data[field_name].extend([defaults[field_name]] * n)

    if output_format == "hdf5":
        from almanac.io import _write_models_to_hdf5_group

        with h5.File(output_path, "w", track_order=True) as fp:
            _write_models_to_hdf5_group(fields, data, fp)
    else:
        from astropy.table import Table

        t = Table(data=data)
        t.write(output_path, format=output_format, overwrite=overwrite)


if __name__ == "__main__":
    main()
