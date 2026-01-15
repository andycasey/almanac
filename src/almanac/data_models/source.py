import numpy as np
from pydantic import BaseModel, Field
from almanac.data_models.types import *

from pydantic import BaseModel, field_validator
from pydantic_core import PydanticUseDefault

class Source(BaseModel):
    #> Catalog Cross-Matches (SDSS_ID_To_Catalog)
    sdss_id: Int64 = Field()
    catalogid: Int64 = Field()
    version_id: int = Field(description="Version ID for SDSS ID to catalog mapping")
    lead: str = Field(description="Lead catalog for this SDSS ID")
    allstar_dr17_synspec_rev1: str = Field(default="", description="APOGEE allstar DR17 synspec rev1 ID")
    allwise: Int64 = Field(default=-1, description="AllWISE catalog ID")
    catwise: str = Field(default="", description="CatWISE catalog ID")
    catwise2020: str = Field(default="", description="CatWISE 2020 catalog ID")
    gaia_dr2_source: Int64 = Field(default=-1, description="Gaia DR2 source ID")
    gaia_dr3_source: Int64 = Field(default=-1, description="Gaia DR3 source ID")
    glimpse: Int64 = Field(default=-1, description="GLIMPSE catalog ID")
    guvcat: Int64 = Field(default=-1, description="GALEX GUVcat ID")
    panstarrs1: Int64 = Field(default=-1, description="Pan-STARRS1 catalog ID")
    ps1_g18: Int64 = Field(default=-1, description="PS1 G18 catalog ID")
    sdss_dr13_photoobj: Int64 = Field(default=-1, description="SDSS DR13 photo object ID")
    sdss_dr17_specobj: str = Field(default="", description="SDSS DR17 spec object ID")
    skymapper_dr2: Int64 = Field(default=-1, description="SkyMapper DR2 object ID")
    supercosmos: Int64 = Field(default=-1, description="SuperCOSMOS object ID")
    tic_v8: Int64 = Field(default=-1, description="TESS Input Catalog v8 ID")
    twomass_id: str = Field(default="", description="Input 2MASS ID to config file (use with caution)")
    twomass_psc: Int64 = Field(default=-1, description="2MASS Point Source Catalog ID")
    tycho2: str = Field(default="", description="Tycho-2 catalog designation")
    unwise: str = Field(default="", description="unWISE object identifier")

    #> Gaia DR3 Data
    gaia_source_id: Int64 = Field(default=-1, alias="source_id", description="Gaia DR3 source ID")
    gaia_ra: float = Field(default=float('NaN'), description="Gaia DR3 right ascension [deg]")
    gaia_ra_error: float = Field(default=float('NaN'), description="Gaia DR3 RA error [mas]")
    gaia_dec: float = Field(default=float('NaN'), description="Gaia DR3 declination [deg]")
    gaia_dec_error: float = Field(default=float('NaN'), description="Gaia DR3 Dec error [mas]")
    gaia_parallax: float = Field(default=float('NaN'), alias="parallax", description="Gaia DR3 parallax [mas]")
    gaia_parallax_error: float = Field(default=float('NaN'), alias="parallax_error", description="Gaia DR3 parallax error [mas]")
    gaia_pm: float = Field(default=float('NaN'), alias="pm", description="Gaia DR3 total proper motion [mas/yr]")
    gaia_pmra: float = Field(default=float('NaN'), alias="pmra", description="Gaia DR3 proper motion in RA [mas/yr]")
    gaia_pmra_error: float = Field(default=float('NaN'), alias="pmra_error", description="Gaia DR3 PM RA error [mas/yr]")
    gaia_pmdec: float = Field(default=float('NaN'), alias="pmdec", description="Gaia DR3 proper motion in Dec [mas/yr]")
    gaia_pmdec_error: float = Field(default=float('NaN'), alias="pmdec_error", description="Gaia DR3 PM Dec error [mas/yr]")
    gaia_ruwe: float = Field(default=float('NaN'), alias="ruwe", description="Gaia DR3 RUWE")
    gaia_duplicated_source: bool = Field(default=False, alias="duplicated_source", description="Gaia DR3 duplicated source flag")
    gaia_phot_g_mean_mag: float = Field(default=float('NaN'), alias="phot_g_mean_mag", description="Gaia DR3 G-band mean magnitude [mag]")
    gaia_phot_bp_mean_mag: float = Field(default=float('NaN'), alias="phot_bp_mean_mag", description="Gaia DR3 BP-band mean magnitude [mag]")
    gaia_phot_rp_mean_mag: float = Field(default=float('NaN'), alias="phot_rp_mean_mag", description="Gaia DR3 RP-band mean magnitude [mag]")
    gaia_phot_bp_rp_excess_factor: float = Field(default=float('NaN'), alias="phot_bp_rp_excess_factor", description="Gaia DR3 BP/RP excess factor")
    gaia_radial_velocity: float = Field(default=float('NaN'), alias="radial_velocity", description="Gaia DR3 radial velocity [km/s]")
    gaia_radial_velocity_error: float = Field(default=float('NaN'), alias="radial_velocity_error", description="Gaia DR3 radial velocity error [km/s]")
    gaia_rv_nb_transits: int = Field(default=-1, alias="rv_nb_transits", description="Gaia DR3 number of RV transits")
    gaia_rv_nb_deblended_transits: int = Field(default=-1, alias="rv_nb_deblended_transits", description="Gaia DR3 number of deblended RV transits")
    gaia_rv_visibility_periods_used: int = Field(default=-1, alias="rv_visibility_periods_used", description="Gaia DR3 RV visibility periods used")
    gaia_rv_expected_sig_to_noise: float = Field(default=float('NaN'), alias="rv_expected_sig_to_noise", description="Gaia DR3 expected RV S/N")
    gaia_rv_renormalised_gof: float = Field(default=float('NaN'), alias="rv_renormalised_gof", description="Gaia DR3 RV renormalized GoF")
    gaia_rv_chisq_pvalue: float = Field(default=float('NaN'), alias="rv_chisq_pvalue", description="Gaia DR3 RV chi-squared p-value")
    gaia_rv_time_duration: float = Field(default=float('NaN'), alias="rv_time_duration", description="Gaia DR3 RV time duration [days]")
    gaia_rv_amplitude_robust: float = Field(default=float('NaN'), alias="rv_amplitude_robust", description="Gaia DR3 RV amplitude robust [km/s]")
    gaia_rv_template_teff: float = Field(default=float('NaN'), alias="rv_template_teff", description="Gaia DR3 RV template Teff [K]")
    gaia_rv_template_logg: float = Field(default=float('NaN'), alias="rv_template_logg", description="Gaia DR3 RV template log(g) [dex]")
    gaia_rv_template_fe_h: float = Field(default=float('NaN'), alias="rv_template_fe_h", description="Gaia DR3 RV template [Fe/H] [dex]")
    gaia_rv_atm_param_origin: int = Field(default=-1, alias="rv_atm_param_origin", description="Gaia DR3 RV atmospheric parameter origin")
    gaia_vbroad: float = Field(default=float('NaN'), alias="vbroad", description="Gaia DR3 spectral line broadening [km/s]")
    gaia_vbroad_error: float = Field(default=float('NaN'), alias="vbroad_error", description="Gaia DR3 spectral line broadening error [km/s]")
    gaia_vbroad_nb_transits: int = Field(default=-1, alias="vbroad_nb_transits", description="Gaia DR3 vbroad number of transits")
    gaia_grvs_mag: float = Field(default=float('NaN'), alias="grvs_mag", description="Gaia DR3 G_RVS magnitude [mag]")
    gaia_grvs_mag_error: float = Field(default=float('NaN'), alias="grvs_mag_error", description="Gaia DR3 G_RVS magnitude error [mag]")
    gaia_grvs_mag_nb_transits: int = Field(default=-1, alias="grvs_mag_nb_transits", description="Gaia DR3 G_RVS number of transits")
    gaia_rvs_spec_sig_to_noise: float = Field(default=float('NaN'), alias="rvs_spec_sig_to_noise", description="Gaia DR3 RVS spectrum S/N")
    gaia_teff_gspphot: float = Field(default=float('NaN'), alias="teff_gspphot", description="Gaia DR3 GSP-Phot Teff [K]")
    gaia_logg_gspphot: float = Field(default=float('NaN'), alias="logg_gspphot", description="Gaia DR3 GSP-Phot log(g) [dex]")
    gaia_mh_gspphot: float = Field(default=float('NaN'), alias="mh_gspphot", description="Gaia DR3 GSP-Phot [M/H] [dex]")
    gaia_distance_gspphot: float = Field(default=float('NaN'), alias="distance_gspphot", description="Gaia DR3 GSP-Phot distance [pc]")
    gaia_azero_gspphot: float = Field(default=float('NaN'), alias="azero_gspphot", description="Gaia DR3 GSP-Phot A0 [mag]")
    gaia_ag_gspphot: float = Field(default=float('NaN'), alias="ag_gspphot", description="Gaia DR3 GSP-Phot A_G [mag]")
    gaia_ebpminrp_gspphot: float = Field(default=float('NaN'), alias="ebpminrp_gspphot", description="Gaia DR3 GSP-Phot E(BP-RP) [mag]")

    #> 2MASS Point Source Catalog Data
    twomass_designation: str = Field(default="", description="2MASS PSC designation")
    twomass_j_m: float = Field(default=float('NaN'), alias="j_m", description="2MASS J-band magnitude [mag]")
    twomass_j_cmsig: float = Field(default=float('NaN'), alias="j_cmsig", description="2MASS J-band corrected photometric uncertainty [mag]")
    twomass_j_msigcom: float = Field(default=float('NaN'), alias="j_msigcom", description="2MASS J-band combined uncertainty [mag]")
    twomass_j_snr: float = Field(default=float('NaN'), alias="j_snr", description="2MASS J-band signal-to-noise ratio")
    twomass_h_m: float = Field(default=float('NaN'), alias="h_m", description="2MASS H-band magnitude [mag]")
    twomass_h_cmsig: float = Field(default=float('NaN'), alias="h_cmsig", description="2MASS H-band corrected photometric uncertainty [mag]")
    twomass_h_msigcom: float = Field(default=float('NaN'), alias="h_msigcom", description="2MASS H-band combined uncertainty [mag]")
    twomass_h_snr: float = Field(default=float('NaN'), alias="h_snr", description="2MASS H-band signal-to-noise ratio")
    twomass_k_m: float = Field(default=float('NaN'), alias="k_m", description="2MASS K-band magnitude [mag]")
    twomass_k_cmsig: float = Field(default=float('NaN'), alias="k_cmsig", description="2MASS K-band corrected photometric uncertainty [mag]")
    twomass_k_msigcom: float = Field(default=float('NaN'), alias="k_msigcom", description="2MASS K-band combined uncertainty [mag]")
    twomass_k_snr: float = Field(default=float('NaN'), alias="k_snr", description="2MASS K-band signal-to-noise ratio")
    twomass_ph_qual: str = Field(default="", alias="ph_qual", description="2MASS photometric quality flag")
    twomass_rd_flg: str = Field(default="", alias="rd_flg", description="2MASS read flag")
    twomass_bl_flg: str = Field(default="", alias="bl_flg", description="2MASS blend flag")
    twomass_cc_flg: str = Field(default="", alias="cc_flg", description="2MASS contamination and confusion flag")

    #> Target of Opportunity (FPS)
    too: bool = Field(default=False, description="Target of opportunity (FPS era)")
    too_id: int = Field(default=-1, description="Target of opportunity ID (FPS era)")
    too_program: str = Field(default="", description="Target of opportunity program (FPS era)")

    #> SDSS-V Targeting Cartons
    sdss5_target_flags: np.ndarray = Field(
        default_factory=lambda: np.zeros((1, 1), dtype=np.uint64),
        description="SDSS-V target flags bitmask array"
    )

    @field_validator("*", mode="before")
    @classmethod
    def none_to_default(cls, v):
        if v is None:
            raise PydanticUseDefault()
        return v

    class Config:
        validate_by_name = True
        validate_assignment = True
        arbitrary_types_allowed = True
