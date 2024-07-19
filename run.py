from pyraf import iraf

from utilities import *

iraf.set(stdimage="imtgmos")

# See config in utilities.py too
config = {
    # MDF and observatory name
    "mdf": "gnifu_slitr_mdf.fits",
    "observatory": "Gemini-North",
    # Biases
    "biasesRefs": [
        "N20240609S0158",
        "N20240609S0157",
        "N20240609S0159",
        "N20240609S0160",
        "N20240609S0161",
        "N20240610S0227",
        "N20240610S0228",
        "N20240610S0229",
        "N20240610S0230",
        "N20240610S0231",
        "N20240611S0150",
        "N20240611S0151",
        "N20240611S0152",
        "N20240611S0153",
        "N20240611S0154",
        "N20240612S0142",
        "N20240613S0168",
        "N20240613S0169",
        "N20240613S0170",
        "N20240613S0171",
        "N20240613S0172",
        "N20240613S0267",
        "N20240613S0268",
        "N20240613S0269",
        "N20240614S0001",
        "N20240614S0002"
    ],
    # Science
    "science": {
        "refs": ["N20240609S0067"],
        "flatRefs": ["N20240609S0066"],
        "arcRefs": ["N20240609S0151"],
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

housekeeping(config)
# create_MDF(config["mdf"], config["standardStar"]["flatRefs"])

def standard_star():
    refs = config["standardStar"]["refs"]
    flats = config["standardStar"]["flatRefs"]
    arcs = config["standardStar"]["arcRefs"]
    "Flats and sensitivity function"
    wavelength(config["standardStar"])
    flat_bundle_gaps(flats)
    remove_scatter(flats, interactive=False, xorder=[5], yorder=[6])
    qe_correct(flats, arcs)
    response_function(flats)
    view_response(flats)
    
    "Sensitivity function"
    sci_trace_reference(refs) # rg
    remove_scatter(refs, gapSolution=flats[0], xorder=[5], yorder=[6], interactive=False) # b    
    reject_cosmic_rays(refs) # x
    # TODO: We are here
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
    remove_scatter(flats, interactive=False, xorder=[5], yorder=[6])
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

# standard_star()
science_flats_arc(False)
science()
print("Finished :)")