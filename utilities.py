from time import sleep
from pyraf import iraf
from pyraf.iraf import gemini, gmos
from astropy.io import fits
import os

# Config
rawDir = "../data/"
biasDir = "../data/biases/"
calDir = "../calibrations/"
masterBiasName = "bias.fits"
bpmRef = "../data/bpm_20230729_gmos-n_Ham_11_full_12amp.fits"

biasPath = calDir + masterBiasName
mdfPath = "gnifu_slitr_mdf.fits"

def housekeeping():
    # iraf.delete("database/*")
    iraf.delete("tmp*")
    iraf.delete("*.log")
        
    

# FUNCTIONS THAT HELP OUT

def skip_step(inRefs, addPrefix, inPrefix):
    for ref in inRefs:
        iraf.copy(inPrefix+ref+'.fits', addPrefix+inPrefix+ref+'.fits')

def iraf_list(lst, prefix="", suffix=""):
    return ",".join(map(lambda x: prefix+x+suffix, lst))
  
def display_image(image):
    try:
        iraf.gdisplay(image, 1, fl_paste="no")
    except:
        print("Problem displaying image :\(")
      
def view_flats(flatRefs):
    print("Verifying flats")
    for flat in flatRefs:
        flat = flat.strip()
        for i in range(3):
            iraf.imexamine("brg"+flat+"[sci,"+str(i+1)+"]", 1)
              
def view_qe(flatRefs):
    for flat in flatRefs:
        iraf.gfdisplay("eqbrg"+flat, 1, version=1)
             
def view_response(flatRefs):
    for flat in flatRefs:
        iraf.gfdisplay(flat+"_resp", 1, version="1")

def view_scatter(refs):
    for ref in refs:
        for i in range(3):
            iraf.imexamine('brg'+ref+'[sci,'+str(i+1)+']', 1)
# FUNCTIONS THAT PROCESS THINGS
def create_master_bias(config):
    biasRefs = config["biasesRefs"]
    iraf.imdelete(biasPath, verify="no")
    iraf.gbias(iraf_list(biasRefs), masterBiasName, rawpath=biasDir, fl_vardq="yes")
    iraf.copy(masterBiasName, calDir)

def create_MDF(mdf, flatRefs):
    print("Viewing MDF")
    iraf.imdelete(iraf_list(flatRefs, "g"), verify="no")
    iraf.imdelete(iraf_list(flatRefs, "rg"), verify="no")
    if (not os.path.isfile(mdf)):
        print("Copying brand new MDF")
        iraf.copy("gmos$data/"+mdf, ".", verbose="no")
    iraf.gfreduce(iraf_list(flatRefs), rawpath=rawDir, fl_extract="no",    bias=biasPath, \
        fl_over="yes", fl_trim="yes", mdffile=mdfPath, mdfdir="./",  \
        slits="red", fl_fluxcal="no", fl_gscrrej="no", \
        fl_wavtran="no", fl_skysub="no", fl_inter="no", fl_vardq="no", Stdout=1)
    iraf.imdelete(iraf_list(flatRefs, "erg"), verify="no")
    iraf.gfextract(iraf_list(flatRefs, "rg"), fl_inter="yes")
    
def sensitivity_function(stdStar):
    """Calculates the sensitivity function for a standard star. Make sure to answer "yes" to the bandpass question.
    """
    outflux = stdStar["stdRoot"]+"std"
    sensfunc = calDir+stdStar["stdRoot"]+'sens'
    ref = stdStar["refs"][0]
    
    iraf.imdelete(iraf_list(stdStar["refs"], "astxeqxbrg"))
    iraf.delete(outflux, verify="no")
    iraf.imdelete(sensfunc, verify="no")
    
    iraf.gfapsum(iraf_list(stdStar["refs"], "stxeqxbrg"), combine="sum", fl_inter="no") 

    iraf.gsstandard("astxeqxbrg"+ref, outflux, sensfunc, \
            starname=stdStar["name"], observatory="Gemini-North", \
            caldir=stdStar["caldir"], extinction=stdStar["extinction"], fl_inter="yes", \
            function="spline3", order=7)

def wavelength(inObs):
    flatRefs = inObs["flatRefs"]
    arcRefs = inObs["arcRefs"]
    iraf.imdelete(iraf_list(flatRefs, "g"), verify="no")
    iraf.imdelete(iraf_list(flatRefs, "rg"), verify="no")
    iraf.imdelete(iraf_list(flatRefs, "erg"), verify="no")
    iraf.imdelete(iraf_list(arcRefs, "g"), verify="no")
    iraf.imdelete(iraf_list(arcRefs, "rg"), verify="no")
    iraf.imdelete(iraf_list(arcRefs, "erg"), verify="no")
    print("Finding wavelength solution...")
    print("Extracting flats")
    iraf.gfreduce(iraf_list(flatRefs), rawpath=rawDir, fl_extract="yes", bias=biasPath, \
        fl_over="yes", fl_trim="yes", mdffile=mdfPath, mdfdir="./",  \
        slits="red", fl_fluxcal="no", fl_gscrrej="no", \
        fl_wavtran="no", fl_skysub="no", fl_inter="no", fl_vardq="yes", bpmfile=bpmRef)
    print("Extracting arc")
    iraf.gfreduce(iraf_list(arcRefs), rawpath=rawDir, fl_extract="yes", recenter="no", \
        trace="no", reference=iraf_list(flatRefs, "erg"), fl_bias="no", \
        fl_over="yes", fl_trim="yes", mdffile=mdfPath, mdfdir="./", \
        slits="red", fl_fluxcal="no", fl_gscrrej="no", \
        fl_wavtran="no", fl_skysub="no", fl_inter="no")
    print("Finding solution (interactive)")
    iraf.gswavelength(iraf_list(arcRefs, "erg"), \
        nlost=10, ntarget=15, threshold=25, \
        coordlis="gmos$data/GCALcuar.dat", fl_inter="yes")    
    
def sci_trace_reference(inRefs):
    print("Making trace reference...")
    iraf.imdelete(iraf_list(inRefs, "g"), verify="no")
    iraf.imdelete(iraf_list(inRefs, "rg"), verify="no")
    
    iraf.gfreduce(iraf_list(inRefs), rawpath=rawDir, fl_extract="no",    bias=biasPath, \
        fl_over="yes", fl_trim="yes", mdffile=mdfPath, mdfdir="./",  \
        slits="red", fl_fluxcal="no", fl_gscrrej="no", \
        fl_wavtran="no", fl_skysub="no", fl_inter="no", fl_vardq="yes", Stdout=1)
            
def flat_bundle_gaps(inRefs):
    iraf.delete("blkmask_*", verify="no")
    
    for flat in inRefs:
        iraf.gffindblocks("rg"+flat, "erg"+flat, "blkmask_"+flat)
    
def remove_scatter(inRefs, gapSolution = "", xorder = (3), yorder = (3), interactive = True):
    iraf.imdelete(iraf_list(inRefs, "brg"))
    xorderStr = ",".join(map(str, xorder))
    yorderStr = ",".join(map(str, yorder))
    print("Remove scatter time!")
    for img in inRefs:
        gaps = gapSolution if gapSolution else img
        iraf.gfscatsub("rg"+img, "blkmask_"+gaps, outimage="", \
                prefix="b", xorder=xorderStr, yorder=yorderStr, \
                cross="yes", fl_inter=interactive)
      
def qe_correct(inRefs, arcRefs):
    iraf.imdelete(iraf_list(inRefs, "qbrg"))
    iraf.imdelete(iraf_list(inRefs, "eqbrg"))
    iraf.gfreduce(iraf_list(inRefs, "brg"), fl_extract=True, fl_qecorr="yes", \
         qe_refim="erg"+arcRefs[0], fl_addmdf="no", fl_bias="no", \
         fl_over="no", fl_trim="no", mdffile=mdfPath, mdfdir="./", \
         slits="red", fl_fluxcal="no", fl_gscrrej="no", \
         fl_wavtran="no", fl_skysub="no", fl_inter="no", \
         fl_vardq="yes")
    
def response_function(flatRefs):
    for flat in flatRefs:
        iraf.imdelete(flat+"_resp")
        iraf.gfresponse("eqbrg"+flat, outimage=flat+"_resp", sky="", \
                        order=45, func="spline3", sample="*", \
                        fl_fit="yes", fl_inter="yes")
    
def reject_cosmic_rays(inRefs):
    iraf.imdelete(iraf_list(inRefs, "xbrg"))
    for img in inRefs:
        iraf.gemcrspec("brg"+img, "xbrg"+img, logfile="crrej.log", \
                key_gain="GAIN", key_ron="RDNOISE", xorder=9, \
                yorder=-1, sigclip=4.5, sigfrac=0.5, objlim=1., \
                    # TODO: maybe bump up niter!
                niter=4, verbose="yes", fl_vardq="yes")

def sci_qe_correct(inRefs, refArcs, refFlats):
    iraf.imdelete(iraf_list(inRefs, "qxbrg"))
    iraf.imdelete(iraf_list(inRefs, "eqxbrg"))
    arc = refArcs[0]
    response = refFlats[0] + "_resp"
    refflat = "eqbrg" + refFlats[0]
    for img in inRefs:
        iraf.gqecorr("xbrg"+img, refimage="erg"+arc, fl_correct="yes", \
                fl_vardq="yes", verbose="yes")
        iraf.gfextract("qxbrg"+img, response=response, recenter="no", \
                trace="no", reference=refflat, weights="none", \
                fl_vardq="yes", bpmfile=bpmRef)
       
# fxhead eqxbrgS20060327S0043 to find width/height
def bad_column_fix(inRefs, badCols, width, height):
    # write badCols to maskbadcol.txt first
    with open("maskbadcol.txt", "w") as f:
        for col in badCols:
            f.write(str(col) + "\n")
    iraf.text2mask("maskbadcol.txt", "maskbadcol.pl", width, height)
    
    for img in inRefs:
        iraf.copy("eqxbrg"+img+".fits", "tmp"+img+".fits")
        iraf.proto.fixpix("tmp"+img+".fits[sci,1]", "maskbadcol.pl", \
                        linterp="1,2,3,4")
        iraf.copy("tmp"+img+".fits", "xeqxbrg"+img+".fits")
        iraf.imarith("maskbadcol.pl", "+", "xeqxbrg"+img+".fits[dq,1]", \
                    "tmpdq"+img)
        iraf.imcopy("tmpdq"+img+"[*,*]", "xeqxbrg"+img+".fits[dq,1][*,*]")

    iraf.imdelete("tmp*.fits")
    iraf.delete("mask*.pl") 
    
def angstroms_per_pixel(inRefs, arcRefs):
    test = inRefs[0]
    arc = arcRefs[0]
    iraf.imdelete("txeqxbrg"+test, verify="no")
    iraf.gftransform("xeqxbrg"+test, wavtraname="erg"+arc, fl_vardq="no")
    return iraf.hselect("txeqxbrg"+test+"[sci,1]", "CD1_1", "yes", Stdout=1 )[0]

def rectify_spectra(inRefs, arcRefs, dw):
    arc = arcRefs[0]
    for sci in inRefs:
        print("Rectifying " + sci)
        iraf.imdelete("txeqxbrg"+sci, verify="no")
        iraf.gftransform("xeqxbrg"+sci, wavtraname="erg"+arc, dw=dw, \
                        fl_vardq="yes")

def subtract_sky(inRefs):
    print("Subtracting sky")
    iraf.imdelete(iraf_list(inRefs, "stxeqxbrg"), verify="no")

    for sci in inRefs:
        iraf.gfskysub("txeqxbrg"+sci, fl_inter="no")
        
def spectrophotometric(inRefs, config):
    extinction = config["standardStar"]["extinction"]
    observatory = config["observatory"]
    sensfunc = calDir+config["standardStar"]["stdRoot"]+'sens'
    
    iraf.imdelete(iraf_list(inRefs, "cstxeqxbrg"), verify="no")

    for sci in inRefs:
        iraf.gscalibrate("stxeqxbrg"+sci, sfunction=sensfunc, \
                    obs=observatory, extinction=extinction, \
                    fl_ext="yes", fl_vardq="yes")

def encubenate(inRefs):
    for img in inRefs:
        iraf.imdelete("cstxeqxbrg"+img+"_3D", verify="no")
        iraf.gfcube("cstxeqxbrg"+img,
            outimage="cstxeqxbrg"+img+"_3D", fl_atmdisp="yes", \
            fl_var="yes", fl_dq="yes")
        
def seeing(inRefs):
    iraf.imdelete(iraf_list(inRefs, suffix="_collapsed"), verify="no")
    for sci in inRefs:
        iraf.imdelete(sci+'_collapsed')
        iraf.imcombine('cstxeqxbrg'+sci+'_3D.fits[sci,1][*,*,400:4747]', 
                       sci+'_collapsed', project='yes')
        iraf.imexamine(sci+'_collapsed', 1)