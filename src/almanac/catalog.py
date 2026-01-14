#from itertools import batched
import os
import tempfile
from tqdm import tqdm
from typing import List
from peewee import JOIN, BigIntegerField

from itertools import islice
import concurrent.futures

def batched(iterable, n):
    it = iter(iterable)
    while True:
        batch = tuple(islice(it, n))
        if not batch:
            return
        yield batch

def merge_dicts(*dicts):
    keys = dicts[0].keys()
    return {
        k: next((d[k] for d in dicts if d.get(k, None) is not None), None)
        for k in keys
    }


def query_targeting(sdss_ids: List[int], **kwargs):
    """
    Query the SDSS database for targeting (carton) information.

    :param sdss_ids: List[int]
        List of SDSS IDs to query
    """
    from almanac.database import catalogdb as cdb, targetdb as tdb

    cdb.database.execute_sql("DROP TABLE IF EXISTS tmp_sdss_ids")
    cdb.database.execute_sql(
        "CREATE TEMP TABLE tmp_sdss_ids (sdss_id BIGINT PRIMARY KEY)"
    )

    # Create temporary CSV file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        for sdss_id in sdss_ids:
            f.write(f"{sdss_id}\n")
        tmp_path = f.name

    with open(tmp_path, 'r') as f:
        cursor = cdb.database.connection().cursor()
        cursor.copy_from(f, 'tmp_sdss_ids', columns=('sdss_id',))
        cdb.database.connection().commit()

    os.unlink(tmp_path)

    class Source(cdb.CatalogdbModel):

        sdss_id = BigIntegerField(primary_key=True)

        class Meta:
            table_name = 'tmp_sdss_ids'
            schema = None

    q_cartons = (
        Source
        .select(
            Source.sdss_id,
            tdb.CartonToTarget.carton_pk,
        )
        .join(cdb.SDSS_ID_flat, on=(Source.sdss_id == cdb.SDSS_ID_flat.sdss_id))
        .join(tdb.Target, on=(cdb.SDSS_ID_flat.catalogid == tdb.Target.catalogid))
        .join(tdb.CartonToTarget, on=(tdb.Target.pk == tdb.CartonToTarget.target_pk))
        .tuples()
    )

    yield from q_cartons


def query(sdss_ids: List[int], batch_size: int = 10_000, tqdm_kwds=None):
    """
    Query the SDSS database for targeting, astrometry, and photometry information.

    Parameters
    ----------
    sdss_ids : List[int]
        List of SDSS IDs to query
    batch_size : int, optional
        Number of IDs per batch (default: 10,000)

    Yields
    ------
    dict
        Dictionary containing catalog data for each SDSS ID
    """

    meta = {}
    tqdm_kwds = tqdm_kwds or {}
    with tqdm(desc="Querying catalog", total=len(sdss_ids), **tqdm_kwds) as pb:
        with concurrent.futures.ProcessPoolExecutor(max_workers=32) as executor:
            futures = [
                executor.submit(_query_catalog, batch, i)
                for i, batch in enumerate(batched(sorted(sdss_ids), batch_size))
            ]
            for future in concurrent.futures.as_completed(futures):
                meta.update(future.result())
                pb.update(batch_size)

    raise a
    return meta


def _query_catalog(sdss_ids: List[int], suffix=""):
    """
    Query the SDSS database for targeting, astrometry, and photometry information.

    Uses a temporary table for efficient querying of large ID lists.

    Parameters
    ----------
    sdss_ids : List[int]
        List of SDSS IDs to query
    """

    from almanac.database import catalogdb as cdb, targetdb as tdb

    # Create temporary CSV file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        for sdss_id in sdss_ids:
            f.write(f"{sdss_id}\n")
        tmp_path = f.name

    cdb.database.execute_sql(
        f"CREATE TEMP TABLE tmp_sdss_ids{suffix} (sdss_id BIGINT PRIMARY KEY)"
    )
    with open(tmp_path, 'r') as f:
        cursor = cdb.database.connection().cursor()
        cursor.copy_from(f, f'tmp_sdss_ids{suffix}', columns=('sdss_id',))
        cdb.database.connection().commit()

    cdb.database.execute_sql("SET enable_seqscan = off")
    cdb.database.execute_sql("SET enable_hashjoin = off")
    cdb.database.execute_sql("SET enable_mergejoin = off")
    cdb.database.execute_sql(f"ANALYZE tmp_sdss_ids{suffix}")

    os.unlink(tmp_path)

    class Source(cdb.CatalogdbModel):

        sdss_id = BigIntegerField(primary_key=True)

        class Meta:
            table_name = f'tmp_sdss_ids{suffix}'
            schema = None

    # Query astrometry and photometry using temp table
    q = (
        Source
        .select(
            cdb.SDSS_ID_To_Catalog.sdss_id,
            # TODO: catalogid and version_id may not be reliable in that we
            #       will get multiple versions per sdss_id. There is no good
            #       way to handle this because some other surveys are crossmatched
            #       to an early version, and some to a later version.
            #       We need a thinko here. We need Zach Way.
            cdb.SDSS_ID_To_Catalog.catalogid,
            cdb.SDSS_ID_To_Catalog.version_id,
            cdb.SDSS_ID_To_Catalog.lead,
            cdb.SDSS_ID_To_Catalog.allstar_dr17_synspec_rev1,
            cdb.SDSS_ID_To_Catalog.allwise,
            cdb.SDSS_ID_To_Catalog.catwise,
            cdb.SDSS_ID_To_Catalog.catwise2020,
            cdb.SDSS_ID_To_Catalog.gaia_dr2_source,
            cdb.SDSS_ID_To_Catalog.gaia_dr3_source,
            cdb.SDSS_ID_To_Catalog.glimpse,
            cdb.SDSS_ID_To_Catalog.guvcat,
            cdb.SDSS_ID_To_Catalog.panstarrs1,
            cdb.SDSS_ID_To_Catalog.ps1_g18,
            cdb.SDSS_ID_To_Catalog.sdss_dr13_photoobj,
            cdb.SDSS_ID_To_Catalog.sdss_dr17_specobj,
            cdb.SDSS_ID_To_Catalog.skymapper_dr2,
            cdb.SDSS_ID_To_Catalog.supercosmos,
            cdb.SDSS_ID_To_Catalog.tic_v8,
            cdb.SDSS_ID_To_Catalog.twomass_psc,
            cdb.SDSS_ID_To_Catalog.tycho2,
            cdb.SDSS_ID_To_Catalog.unwise,
            cdb.Gaia_DR3.source_id.alias('gaia_source_id'),
            cdb.Gaia_DR3.ra.alias('gaia_ra'),
            cdb.Gaia_DR3.ra_error.alias('gaia_ra_error'),
            cdb.Gaia_DR3.dec.alias('gaia_dec'),
            cdb.Gaia_DR3.dec_error.alias('gaia_dec_error'),
            cdb.Gaia_DR3.parallax.alias('gaia_parallax'),
            cdb.Gaia_DR3.parallax_error.alias('gaia_parallax_error'),
            cdb.Gaia_DR3.pm.alias('gaia_pm'),
            cdb.Gaia_DR3.pmra.alias('gaia_pmra'),
            cdb.Gaia_DR3.pmra_error.alias('gaia_pmra_error'),
            cdb.Gaia_DR3.pmdec.alias('gaia_pmdec'),
            cdb.Gaia_DR3.pmdec_error.alias('gaia_pmdec_error'),
            cdb.Gaia_DR3.ruwe.alias('gaia_ruwe'),
            cdb.Gaia_DR3.duplicated_source.alias('gaia_duplicated_source'),
            cdb.Gaia_DR3.phot_g_mean_mag.alias('gaia_phot_g_mean_mag'),
            cdb.Gaia_DR3.phot_bp_mean_mag.alias('gaia_phot_bp_mean_mag'),
            cdb.Gaia_DR3.phot_rp_mean_mag.alias('gaia_phot_rp_mean_mag'),
            cdb.Gaia_DR3.phot_bp_rp_excess_factor.alias('gaia_phot_bp_rp_excess_factor'),
            cdb.Gaia_DR3.radial_velocity.alias('gaia_radial_velocity'),
            cdb.Gaia_DR3.radial_velocity_error.alias('gaia_radial_velocity_error'),
            cdb.Gaia_DR3.rv_nb_transits.alias('gaia_rv_nb_transits'),
            cdb.Gaia_DR3.rv_nb_deblended_transits.alias('gaia_rv_nb_deblended_transits'),
            cdb.Gaia_DR3.rv_visibility_periods_used.alias('gaia_rv_visibility_periods_used'),
            cdb.Gaia_DR3.rv_expected_sig_to_noise.alias('gaia_rv_expected_sig_to_noise'),
            cdb.Gaia_DR3.rv_renormalised_gof.alias('gaia_rv_renormalised_gof'),
            cdb.Gaia_DR3.rv_chisq_pvalue.alias('gaia_rv_chisq_pvalue'),
            cdb.Gaia_DR3.rv_time_duration.alias('gaia_rv_time_duration'),
            cdb.Gaia_DR3.rv_amplitude_robust.alias('gaia_rv_amplitude_robust'),
            cdb.Gaia_DR3.rv_template_teff.alias('gaia_rv_template_teff'),
            cdb.Gaia_DR3.rv_template_logg.alias('gaia_rv_template_logg'),
            cdb.Gaia_DR3.rv_template_fe_h.alias('gaia_rv_template_fe_h'),
            cdb.Gaia_DR3.rv_atm_param_origin.alias('gaia_rv_atm_param_origin'),
            cdb.Gaia_DR3.vbroad.alias('gaia_vbroad'),
            cdb.Gaia_DR3.vbroad_error.alias('gaia_vbroad_error'),
            cdb.Gaia_DR3.vbroad_nb_transits.alias('gaia_vbroad_nb_transits'),
            cdb.Gaia_DR3.grvs_mag.alias('gaia_grvs_mag'),
            cdb.Gaia_DR3.grvs_mag_error.alias('gaia_grvs_mag_error'),
            cdb.Gaia_DR3.grvs_mag_nb_transits.alias('gaia_grvs_mag_nb_transits'),
            cdb.Gaia_DR3.rvs_spec_sig_to_noise.alias('gaia_rvs_spec_sig_to_noise'),
            cdb.Gaia_DR3.teff_gspphot.alias('gaia_teff_gspphot'),
            cdb.Gaia_DR3.logg_gspphot.alias('gaia_logg_gspphot'),
            cdb.Gaia_DR3.mh_gspphot.alias('gaia_mh_gspphot'),
            cdb.Gaia_DR3.distance_gspphot.alias('gaia_distance_gspphot'),
            cdb.Gaia_DR3.azero_gspphot.alias('gaia_azero_gspphot'),
            cdb.Gaia_DR3.ag_gspphot.alias('gaia_ag_gspphot'),
            cdb.TwoMassPSC.designation.alias('twomass_designation'),
            cdb.TwoMassPSC.j_m.alias('twomass_j_m'),
            cdb.TwoMassPSC.j_cmsig.alias('twomass_j_cmsig'),
            cdb.TwoMassPSC.j_msigcom.alias('twomass_j_msigcom'),
            cdb.TwoMassPSC.j_snr.alias('twomass_j_snr'),
            cdb.TwoMassPSC.h_m.alias('twomass_h_m'),
            cdb.TwoMassPSC.h_cmsig.alias('twomass_h_cmsig'),
            cdb.TwoMassPSC.h_msigcom.alias('twomass_h_msigcom'),
            cdb.TwoMassPSC.h_snr.alias('twomass_h_snr'),
            cdb.TwoMassPSC.k_m.alias('twomass_k_m'),
            cdb.TwoMassPSC.k_cmsig.alias('twomass_k_cmsig'),
            cdb.TwoMassPSC.k_msigcom.alias('twomass_k_msigcom'),
            cdb.TwoMassPSC.k_snr.alias('twomass_k_snr'),
            cdb.TwoMassPSC.ph_qual.alias('twomass_ph_qual'),
            cdb.TwoMassPSC.rd_flg.alias('twomass_rd_flg'),
            cdb.TwoMassPSC.bl_flg.alias('twomass_bl_flg'),
            cdb.TwoMassPSC.cc_flg.alias('twomass_cc_flg'),
        )
        .join(cdb.SDSS_ID_To_Catalog, on=(Source.sdss_id == cdb.SDSS_ID_To_Catalog.sdss_id))
        .join(cdb.Gaia_DR3, join_type=JOIN.LEFT_OUTER, on=(cdb.Gaia_DR3.source_id == cdb.SDSS_ID_To_Catalog.gaia_dr3_source))
        .switch(cdb.SDSS_ID_To_Catalog)
        .join(cdb.TwoMassPSC, join_type=JOIN.LEFT_OUTER, on=(cdb.TwoMassPSC.pts_key == cdb.SDSS_ID_To_Catalog.twomass_psc))
        .switch(cdb.SDSS_ID_To_Catalog)
        .dicts()
    )

    meta = {}
    for item in q.iterator():
        sdss_id = item['sdss_id']
        meta[sdss_id] = merge_dicts(item, meta.get(sdss_id, {}))

    q_cartons = (
        Source
        .select(
            Source.sdss_id,
            tdb.CartonToTarget.carton_pk,
        )
        .join(cdb.SDSS_ID_flat, on=(Source.sdss_id == cdb.SDSS_ID_flat.sdss_id))
        .join(tdb.Target, on=(cdb.SDSS_ID_flat.catalogid == tdb.Target.catalogid))
        .join(tdb.CartonToTarget, on=(tdb.Target.pk == tdb.CartonToTarget.target_pk))
        .tuples()
    )
    cartons = {}
    for sdss_id, carton_pks in q_cartons:
        try:
            cartons[sdss_id].add(carton_pks)
        except KeyError:
            cartons[sdss_id] = {carton_pks}

    for sdss_id in meta.keys():
        meta[sdss_id]["carton_pks"] = cartons.pop(sdss_id, set())

    return meta
