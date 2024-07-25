import sys
from pyraf import iraf

from utilities import *

iraf.set(stdimage="imtgmos")

# Completed observations:
# 702S0060
# 702S0063
# 609S0067
# 706S0081
# 706S0082
# N20240609S0068.fits
# N20240609S0069.fits

# Silly
# N20240609S0072.fits
# N20240609S0073.fits 	

# See config in utilities.py too
config = {
    # MDF and observatory name
    "mdf": "gnifu_slitr_mdf.fits",
    "observatory": "Gemini-North",
    # Biases
    "biasesRefs": [
        "N20240630S0004",
        "N20240630S0005",
        "N20240630S0006",
        "N20240630S0007",
        "N20240630S0008",
        "N20240710S0017",
        "N20240710S0018",
        "N20240710S0019",
        "N20240710S0020",
        "N20240710S0021",
        "N20240607S0003",
        "N20240607S0004",
        "N20240607S0005",
        "N20240607S0006",
        "N20240607S0007",
    ],
    # Science
    "science": {
        "refs": ["N20240702S0060"],
        "flatRefs": ["N20240702S0061"],
        "arcRefs": ["N20240702S0124"],
        "bpmRef": "../data/bpm_20230729_gmos-n_Ham_11_full_12amp.fits",
    },
    # Standard Star
    "standardStar": {
        "name": "bd284211", #BD+28 4211 
        "stdRoot": "bd284211_",
        "caldir": 'onedstds$spec50cal/',
        "extinction": 'gmos$calib/mkoextinct.dat',
        "refs": ["N20240611S0118"],
        "flatRefs": ["N20240611S0119"],
        "arcRefs": ["N20240611S0149"],
        "bpmRef": "../data/bpm_20230729_gmos-n_Ham_11_full_12amp.fits",
    },
}

housekeeping()
# create_master_bias(config)
# create_MDF(config["mdf"], config["science"]["flatRefs"])

def standard_star():
    refs = config["standardStar"]["refs"]
    flats = config["standardStar"]["flatRefs"]
    arcs = config["standardStar"]["arcRefs"]
    "Flats and sensitivity function"
    wavelength(config["standardStar"])
    flat_bundle_gaps(flats)
    remove_scatter(flats, interactive=False, xorder=[6], yorder=[5])
    qe_correct(flats, arcs)
    response_function(flats)
    # view_response(flats)
    
    "Sensitivity function"
    sci_trace_reference(refs) # rg
    remove_scatter(refs, gapSolution=flats[0], xorder=[4], yorder=[4], interactive=False) # b    
    view_scatter(refs)
    skip_step(refs, "x", "brg")
    sci_qe_correct(refs, arcs, flats) # eq
    skip_step(refs, "x", "eqxbrg")
    dw = angstroms_per_pixel(refs, arcs)
    rectify_spectra(refs, arcs, dw) # tx
    subtract_sky(refs) # s
    sensitivity_function(config["standardStar"])
    
def science_flats_arc(show=False):
    flats = config["science"]["flatRefs"]
    arcs = config["science"]["arcRefs"]
    wavelength(config["science"])
    flat_bundle_gaps(flats)
    remove_scatter(flats, interactive=False, xorder=[6], yorder=[6])
    # view_scatter(flats)
    # sys.exit()
    qe_correct(flats, arcs)
    response_function(flats)
    if show:
        view_response(flats)

def science():
    flats = config["science"]["flatRefs"]
    arcs = config["science"]["arcRefs"]
    refs = config["science"]["refs"]
    sci_trace_reference(refs)
    remove_scatter(refs, gapSolution=flats[0], xorder=[4], yorder=[3], interactive=False)
    reject_cosmic_rays(refs) # long step!
    sci_qe_correct(refs, arcs, flats)
    skip_step(refs, "x", "eqxbrg")
    dw = angstroms_per_pixel(refs, arcs)
    rectify_spectra(refs, arcs, dw) # long step :(
    subtract_sky(refs)
    spectrophotometric(refs, config)
    encubenate(refs)

standard_star()
science_flats_arc(False)
science()
print("Finished :)")