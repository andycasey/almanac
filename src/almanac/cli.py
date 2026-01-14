#!/usr/bin/env python3

import click


@click.group(invoke_without_command=True)
@click.option("-v", "--verbosity", count=True, help="Verbosity level")
@click.option("--mjd", default=None, type=int, help="Modified Julian date to query. Use negative values to indicate relative to current MJD")
@click.option("--mjd-start", default=None, type=int, help="Start of MJD range to query")
@click.option("--mjd-end", default=None, type=int, help="End of MJD range to query")
@click.option("--date", default=None, type=str, help="Date to query (e.g., 2024-01-15)")
@click.option("--date-start", default=None, type=str, help="Start of date range to query")
@click.option("--date-end", default=None, type=str, help="End of date range to query")
@click.option("--apo", is_flag=True, help="Query Apache Point Observatory data")
@click.option("--lco", is_flag=True, help="Query Las Campanas Observatory data")
@click.option("--fibers", "--fibres", is_flag=True, help="Include fibre mappings to targets")
@click.option("--no-x-match", is_flag=True, help="Do not cross-match targets with SDSS-V database")
@click.option("--output", "-O", default=None, type=str, help="Output file")
@click.option("--processes", "-p", default=None, type=int, help="Number of processes to use")
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

    mjds, mjd_min, mjd_max = utils.parse_mjds(mjd, mjd_start, mjd_end, date, date_start, date_end)
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
            display.create_display(),
            refresh_per_second=refresh_per_second,
            screen=True
        )
        if verbosity >= 1
        else nullcontext()
    )
    io_kwds = dict(fibers=fibers, compression=False)
    with (h5.File(output, "a") if output else nullcontext()) as fp:
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

                            observatory, mjd, exposures, sequences = result = future.result()

                            v = mjd - mjd_min + display.offset
                            missing = [e.image_type == "missing" for e in exposures]
                            if any(missing):
                                display.missing.add(v)
                                #buffered_critical_logs.extend(missing)

                            if not exposures:
                                display.no_data[observatory].add(v)
                            else:
                                display.completed[observatory].add(v)
                                results.append(result)
                                if output:
                                    io.update(fp, observatory, mjd, exposures, sequences, **io_kwds)

                            if live is not None and (time() - t) > 1 / refresh_per_second:
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
                    *_, exposures, sequences = result = apogee.get_almanac_data(observatory, mjd, fibers, not no_x_match)
                    v = mjd - mjd_min + display.offset
                    if any([e.image_type == "missing" for e in exposures]):
                        display.missing.add(v)
                        #buffered_critical_logs.extend(missing)

                    if not exposures:
                        display.no_data[observatory].add(v)
                    else:
                        display.completed[observatory].add(v)
                        results.append(result)
                        if output:
                            io.update(fp, observatory, mjd, exposures, sequences, **io_kwds)

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
        Assignment, AssignmentStatus,CartonToTarget, Target, Hole, Observatory,
        Design, DesignToField, Field
    )
    from sdssdb.peewee.sdss5db.catalogdb import SDSS_ID_flat

    from rich.table import Table as RichTable
    from rich.console import Console
    from rich.live import Live

    sq = (
        SDSS_ID_flat
        .select(SDSS_ID_flat.sdss_id)
        .where(
            (SDSS_ID_flat.sdss_id.in_(identifiers))
        |   (SDSS_ID_flat.catalogid.in_(identifiers))
        )
        .alias("sq")
    )
    q = (
        SDSS_ID_flat
        .select(SDSS_ID_flat.catalogid, SDSS_ID_flat.sdss_id)
        .distinct()
        .join(sq, on=(SDSS_ID_flat.sdss_id == sq.c.sdss_id))
        .order_by(SDSS_ID_flat.sdss_id.asc())
        .tuples()
    )

    catalogids = { catalogid: sdss_id for catalogid, sdss_id in q }

    if not catalogids:
        raise click.ClickException(f"Identifiers {identifiers} not found in SDSS-V database")

    # todo hacky
    for m in (
        Target, Observatory, AssignmentStatus, Field, CartonToTarget,
        Assignment, Hole, Design, DesignToField
    ):
        try:
            m.get()
        except:
            print(f"Failed {m}")
        else:
            print(f"Success {m}")

        m._meta.schema = "targetdb"

    q = (
        Target
        .select(
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
        &   (AssignmentStatus.status == 1)
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
        header_style="bold cyan"
    )

    for field_name in ("#", "obs", "mjd", "exposure", "field", "fiber_id", "catalogid", "sdss_id"):
        rich_table.add_column(field_name, justify="center")

    done = set()
    io_kwds = dict(fibers=True, compression=False)
    with (h5.File(output, "a") if output else nullcontext()) as fp:

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
                            rich_table.add_row(*list(map(str, (
                                i,
                                exposure.observatory,
                                exposure.mjd,
                                exposure.exposure,
                                exposure.field_id,
                                target.fiber_id,
                                target.catalogid,
                                catalogids[target.catalogid],
                            ))))
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
@click.option("--mjd", default=None, type=int, help="Modified Julian date to query. Use negative values to indicate relative to current MJD")
@click.option("--mjd-start", default=None, type=int, help="Start of MJD range to query")
@click.option("--mjd-end", default=None, type=int, help="End of MJD range to query")
@click.option("--date", default=None, type=str, help="Date to query (e.g., 2024-01-15)")
@click.option("--date-start", default=None, type=str, help="Start of date range to query")
@click.option("--date-end", default=None, type=str, help="End of date range to query")
@click.option("--apo", is_flag=True, help="Query Apache Point Observatory data")
@click.option("--lco", is_flag=True, help="Query Las Campanas Observatory data")
@click.option("--p", default=-1, type=int, help="Number of workers to use")
def metadata(input_path, mjd, mjd_start, mjd_end, date, date_start, date_end, apo, lco, p, **kwargs):
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
            for future in tqdm(concurrent.futures.as_completed(futures), desc="Collecting SDSS identifiers"):
                sdss_ids.update(future.result())

    from almanac.data_models.source import Source

    results = query(sdss_ids)
    import pickle
    with open("20251128-meta.pkl", "wb") as fp:
        pickle.dump(results, fp)

import h5py as h5
def _get_fiber_exposure(index, obs, mjd, expnum, adjusted_fiber_index, path, outdir):
    with h5.File(path, "r") as almanac_fp:

        key = f"{obs}/{mjd}"
        g = almanac_fp[key]

        d = { k: v[expnum - 1] for k, v in g["exposures"].items() }

        # Get fiber-level info.
        ref_id = max(d["config_id"], d["plate_id"])
        fiber_index = adjusted_fiber_index - 1
        if fiber_index >= 300:
            fiber_index -= 300

        gf = g[f"fibers/{ref_id}"]
        d.update({ k: gf[k][fiber_index] for k in gf.keys() })

    ar1dunical_path = f"{outdir}/apred/{mjd}/ar1Dunical_{obs}_{mjd}_{expnum:04d}_object.h5"

    ar1d_keys = (
        "bitmsk_relFluxFile",
        "cartid",
        "extraction_method",
        #"git_branch",
        #"git_clean",
        #"git_commit",
        "mjd_mid_exposure",
        "ndiff_used",
        "nread_total",
        "trace_found_match",
        "trace_orig_param_fname",
        "trace_type",
        "wavecal_type",
    )
    try:
        with h5.File(ar1dunical_path, "r") as ar1d:
            for key in ar1d_keys:
                d[key] = ar1d[f"metadata/{key}"][...].flatten()[0]
    except FileNotFoundError:
        d["flag_missing_ar1dunical"] = True
    else:
        d["flag_missing_ar1dunical"] = False
    finally:
        return (index, d)




@main.command()
@click.argument("input_path", type=str)
@click.argument("output_path", type=str)
@click.option("--processes", "-p", default=None, type=int, help="Number of processes to use")
def postprocess(input_path, output_path, processes, **kwargs):
    """Post-process an existing Almanac file after reductions are complete."""

    if not input_path or not output_path:
        return

    #print(input_path)
    #input_path = "almanac/allobs_57600_61000.h5"
    #output_path = "/mnt/ceph/users/sdssv/work/acasey/20260114-exposures.h5"
    #print(kwargs)
    import os
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"The file {input_path} does not exist.")

    if processes is None or processes < 0:
        processes = os.cpu_count()

    outdir = os.path.dirname(os.path.abspath(input_path)) + "/../"

    from glob import glob
    import h5py as h5
    import pickle
    import concurrent.futures
    import numpy as np
    import warnings
    from astropy.coordinates import SkyCoord, EarthLocation
    from astropy import units as u
    from astropy.time import Time

    from almanac.display import TaskDisplay
    from almanac.data_models.spectrum import Spectrum
    from almanac.io import _write_models_to_hdf5_group

    from sdss_semaphore.targeting import TargetingFlags

    data_dict = {}

    with TaskDisplay() as display:

        # Get full list info
        full_list_info_path = glob(f"{outdir}/arMADGICS/raw_*/full_list_info.h5")[0]
        arMADGICS_suffix = full_list_info_path.split("/")[-2][4:]

        fields = { **Spectrum.model_fields, **Spectrum.model_computed_fields }

        total = None
        with h5.File(full_list_info_path, "r") as full_list_info:
            keys = list(full_list_info.keys())
            for i, key in enumerate(keys):
                if i == 0:
                    display.add_task("summary", f"Reading {len(full_list_info[key]):,} summary spectra", total=len(keys))
                data_dict[key] = full_list_info[key][:]
                display.advance("summary")
        display.complete("summary")

        total = len(data_dict["sdss_id"])

        # Load in the arMADGICS scalar arrays
        arMADGICS_dir = f"{outdir}/arMADGICS/wu_th_{arMADGICS_suffix}/"
        arMADGICS_paths = glob(f"{arMADGICS_dir}/arMADGICS_out_*.h5")

        display.add_task("armadgics", f"Reading {len(arMADGICS_paths)} arMADGICS scalar files", total=len(arMADGICS_paths))
        display.add_task("map", "Distributing Almanac data", total=total)
        display.add_task("reduce", "Collecting Almanac data", total=total)
        display.add_task("v_rad", "Computing radial velocities", total=1)
        display.add_task("metadata", "Assigning metadata to rows", total=total)
        display.add_task("targeting", "Aggregating targeting cartons", total=total)
        display.add_task("write", f"Write to {output_path}", total=len(fields))

        for path in arMADGICS_paths:
            basename = os.path.basename(path)
            key = basename[14:-3]
            with h5.File(path, "r") as fp:
                if fp[key].ndim == 1:
                    data_dict[key] = fp[key][:]
            display.advance("armadgics")
        display.complete("armadgics")

        # Load in the almanac data for every fiber exposure
        with concurrent.futures.ProcessPoolExecutor(processes) as executor:
            futures = []
            for i in range(total):
                futures.append(
                    executor.submit(
                        _get_fiber_exposure,
                        i,
                        data_dict["tele"][i].decode(),
                        data_dict["mjd"][i].decode(),
                        data_dict["expnum"][i],
                        data_dict["adjfiberindx"][i],
                        input_path,
                        outdir
                    )
                )
                display.advance("map")
            display.complete("map")

            for future in concurrent.futures.as_completed(futures):
                i, d = future.result()
                for k, v in d.items():
                    try:
                        data_dict[k][i] = v
                    except KeyError:
                        data_dict.setdefault(k, [None] * total)
                        data_dict[k][i] = v
                display.advance("reduce")
            display.complete("reduce")

        # Load in metadata
        # TODO: This will be stored in the almanac file eventually, but for now..
        with open("/mnt/home/acasey/almanac/20251128-meta.pkl", "rb") as fp:
            source_meta = pickle.load(fp)

        for i in range(total):
            for k, v in source_meta.get(data_dict["sdss_id"][i], {}).items():
                try:
                    data_dict[k][i] = v
                except KeyError:
                    data_dict.setdefault(k, [None] * total)
                    data_dict[k][i] = v
            display.advance("metadata")
        display.complete("metadata")

        def propagate_pixels_to_z(p, δλ=6e-6):
            return δλ * np.log(10) * 10**(p * δλ)

        def propagate_z_to_v(z, c=299792.458):
            return np.abs(4 * (z + 1) / (((z + 1) ** 2 + 1) ** 2)) * c

        def pixels_to_z(x, delta=6e-6):
            return 10**(x * delta) - 1

        def z_to_v(z, c=299792.458):
            return ((z+1)**2-1)/((z+1)**2+1)*c

        def v_to_z(v, c=299792.458):
            return v / c

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            e_z = propagate_pixels_to_z(np.array(data_dict["RV_pixoff_final"])) * np.sqrt(np.array(data_dict["RV_pix_var"]))
        z = pixels_to_z(np.array(data_dict["RV_pixoff_final"]))

        ra, dec, mjd_mid_exposure = map(np.array, (
            data_dict["ra"], data_dict["dec"], data_dict["mjd_mid_exposure"]
        ))

        data_dict["v_barycentric_correction"] = np.nan * np.ones(total)
        for obs in ("apo", "lco"):
            mask = (
                (np.array(data_dict["tele"]).astype(str) == obs)
            *   (np.isfinite(ra))
            *   (np.isfinite(dec))
            *   (np.isfinite(mjd_mid_exposure))
            *   (360 >= ra) * (ra >= 0)
            *   (90 >= dec) * (dec >= -90)
            )

            data_dict["v_barycentric_correction"][mask] = (
                SkyCoord(ra=ra[mask] * u.deg, dec=dec[mask] * u.deg)
                .radial_velocity_correction(
                    kind="barycentric",
                    obstime=Time(mjd_mid_exposure[mask], format="mjd"),
                    location=EarthLocation.of_site(obs)
                )
                .to(u.km/u.s)
                .value
            )

        z_corr = v_to_z(data_dict["v_barycentric_correction"])
        data_dict.update(
            v_rel=z_to_v(z),
            v_rad=z_to_v(z + z_corr + z_corr * z),
            e_v_rad=propagate_z_to_v(z) * e_z,
        )
        display.complete("v_rad")

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


        flags = TargetingFlags(np.zeros((total, 1)))

        n_unknown_carton_assignments = 0
        for i, carton_pks in enumerate(data_dict["carton_pks"]):
            for carton_pk in (carton_pks or {}):
                try:
                    flags.set_bit_by_carton_pk(i, carton_pk)
                except KeyError:
                    n_unknown_carton_assignments += 1
            display.advance("targeting")

        data_dict["sdss5_target_flags"] = flags.array

        # Last check point for fields
        expected = set(fields)
        actual = set(data_dict.keys())
        missing_fields = expected - actual
        ignored_fields = actual - expected
        for key in missing_fields:
            data_dict[key] = [None] * total

        def callback(*args, **kwargs):
            try:
                display.advance("write")
            except:
                pass

        with h5.File(output_path, "w", track_order=True) as fp:
            _write_models_to_hdf5_group(
                fields,
                data_dict,
                fp,
                callback=callback
            )
        display.complete("write")

        if n_unknown_carton_assignments > 0:
            print(f"Warning: {n_unknown_carton_assignments:,} unknown carton assignments encountered")

        for key in ignored_fields:
            print(f"Warning: ignored field {key} in data_dict")





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
        raise click.ClickException(f"Output path {output_path} already exists. Use --overwrite to overwrite.")

    if given_format is None:
        if output_path.lower().endswith(".fits"):
            return "fits"
        elif output_path.lower().endswith(".csv"):
            return "csv"
        elif output_path.lower().endswith(".hdf5") or output_path.lower().endswith(".h5"):
            return "hdf5"
        else:
            raise click.ClickException("Cannot infer output format from output path. Please specify --format")
    return given_format


@dump.command()
@click.argument("input_path", type=str)
@click.argument("output_path", type=str)
@click.option("--format", "-f", default=None, type=click.Choice(["fits", "csv", "hdf5"]), help="Output format")
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

                is_object = (
                    (group["exposures/image_type"][:].astype(str) == "object")
                )
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
                        print(f"Warning couldnt get config {config_id} for {observatory} {mjd}")
                        continue

                    ok = (
                        (
                            (config_group["catalogid"][:] > 0)
                        |   (config_group["sdss_id"][:] > 0)
                        |   (config_group["twomass_designation"][:].astype(str) != "")
                        )
                    *   (
                            (config_group["category"][:].astype(str) == "science")
                        |   (config_group["category"][:].astype(str) == "standard_apogee")
                        |   (config_group["category"][:].astype(str) == "standard_boss")
                        |   (config_group["category"][:].astype(str) == "open_fiber")
                        )
                    )
                    sdss_ids = config_group["sdss_id"][:][ok]
                    catalogids = config_group["catalogid"][:][ok]
                    for sdss_id, catalogid in zip(sdss_ids, catalogids):
                        stars.setdefault(sdss_id, deepcopy(default))
                        stars[sdss_id].setdefault("catalogid", catalogid) # this can change over time,... should we track that/
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
@click.option("--format", "-f", default=None, type=click.Choice(["fits", "csv", "hdf5"]), help="Output format")
@click.option("--overwrite", is_flag=True, help="Overwrite existing output file")
def visits(input_path, output_path, format, overwrite, **kwargs):
    """Create a visit-level summary file"""

    pass



@dump.command()
@click.argument("input_path", type=str)
@click.argument("output_path", type=str)
@click.option("--format", "-f", default=None, type=click.Choice(["fits", "csv", "hdf5"]), help="Output format")
@click.option("--overwrite", is_flag=True, help="Overwrite existing output file")
def exposures(input_path, output_path, format, overwrite, **kwargs):
    """Create an exposure-level summary file"""

    import os
    import h5py as h5
    import numpy as np

    output_format = check_paths_and_format(input_path, output_path, format, overwrite)

    from almanac.data_models import Exposure

    fields = { **Exposure.model_fields, **Exposure.model_computed_fields }
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

        fields = { **Exposure.model_fields, **Exposure.model_computed_fields }

        with h5.File(output_path, "w", track_order=True) as fp:
            _write_models_to_hdf5_group(fields, data, fp)
    else:
        from astropy.table import Table
        t = Table(data=data)
        t.write(output_path, format=output_format, overwrite=overwrite)


@dump.command()
@click.argument("input_path", type=str)
@click.argument("output_path", type=str)
@click.option("--format", "-f", default=None, type=click.Choice(["fits", "csv", "hdf5"]), help="Output format")
@click.option("--overwrite", is_flag=True, help="Overwrite existing output file")
def fibers(input_path, output_path, format, overwrite, **kwargs):
    """Create a fiber-level summary file"""

    import os
    import h5py as h5
    import numpy as np

    output_format = check_paths_and_format(input_path, output_path, format, overwrite)

    from almanac.data_models.fps import FPSTarget
    from almanac.data_models.plate import PlateTarget

    fields = { **FPSTarget.model_fields, **FPSTarget.model_computed_fields,
              **PlateTarget.model_fields, **PlateTarget.model_computed_fields }

    defaults = { name: spec.default for name, spec in fields.items() if hasattr(spec, "default") }
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
