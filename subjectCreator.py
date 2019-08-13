import os
import shutil
import json
from numpy import sum, max, array, copy
import astropy.units as u
from astropy.io import fits
import sdssCutoutGrab as scg
import createSubjectsFunctions as csf
import montage_wrapper as montage
import requests
from astropy.nddata.utils import NoOverlapError


def make_dir(d):
    if not os.path.exists(d):
        os.mkdir(d)


# stolen and edited from SDSS_cutoutGrab
def queryFromRaDec(ra, dec, radius=0.5, limit=10, verbose=False):
    queryUrl = 'http://skyserver.sdss.org/dr13/en/tools/search/x_results.aspx'
    res = requests.get(queryUrl, params={
        'searchtool': 'Radial',
        'TaskName': 'Skyserver.Search.Radial',
        'whichphotometry': 'optical',
        'coordtype': 'equatorial',
        'ra': ra,
        'dec': dec,
        'radius': radius,
        'limit': limit,
        'format': 'json',
    })
    if res.status_code == 200:
        try:
            result = res.json()[0]['Rows']
            return sorted(
                result,
                key=lambda i: (
                    (i['ra'] - 160.65883)**2 + (i['dec'] - 23.95189)**2
                )
            )
        except json.decoder.JSONDecodeError:
            if verbose:
                print(res.url)
                print('Could not parse returned JSON: ' + res.url)
            return []
    else:
        if verbose:
            print('Could not connect to SDSS skyserver: ' + res.url)
        return []


def reduceFrameList(frameList):
    s = set()
    out = []

    def makeHashable(f):
        return '{run},{camcol},{field},{rerun}'.format(**f)

    for frame in frameList:
        foo = makeHashable(frame)
        if foo not in s:
            s |= {foo}
            out += [{
                a: b for a, b in zip(
                    ['run', 'camcol', 'field', 'rerun'],
                    foo.split(',')
                )
            }]
    return s, out


def getFramesForMontage(frameList, ra, dec):
    montageFolder = 'montageGroups/{}+{}'.format(ra, dec)
    make_dir('montageGroups')
    make_dir(montageFolder)
    for frame in frameList:
        downloadedFile = os.path.abspath(scg.getBandFits(frame))
        try:
            os.symlink(
                downloadedFile,
                '{}/{}'.format(montageFolder, downloadedFile.split('/')[-1])
            )
        except FileExistsError:
            pass
    return montageFolder


def doMontage(folder, outFile, header=None):
    if (os.path.exists(outFile)):
        shutil.rmtree(outFile)
    folderContents = os.listdir(folder)
    if len(folderContents) > 1:
        montage.wrappers.mosaic(
            folder,
            outFile,
            header=header,
            background_match=True
        )
        return os.path.abspath(outFile + '/mosaic.fits')
    else:
        return os.path.abspath(folder + '/' + os.listdir(folder)[0])


def makeModel(imageData, psfCut, resizeTo=(512, 512)):
    return {
        'psf': psfCut.tolist(),
        'psfWidth': psfCut.shape[1],
        'psfHeight': psfCut.shape[0],
        'width': imageData.shape[1],
        'height': imageData.shape[0],
        'imageWidth': resizeTo[0],
        'imageHeight': int(
            imageData.shape[0] / imageData.shape[1] * resizeTo[0]
        )
    }


def makeDifference(imageData, psfCut, mask, resizeTo=(512, 512)):
    return {
        'psf': psfCut.tolist(),
        'psfWidth': psfCut.shape[1],
        'psfHeight': psfCut.shape[0],
        'mask': array(mask, dtype=float).tolist(),
        'imageData': (imageData / max(imageData)).tolist(),
        'multiplier': float(max(imageData)),
        'width': imageData.shape[1],
        'height': imageData.shape[0],
        'imageWidth': resizeTo[0],
        'imageHeight': int(
            imageData.shape[0] / imageData.shape[1] * resizeTo[0]
        ),
    }


def makeMetadata(galObj, extra_metadata={}):
    return {
        'ra': '{:05f}'.format(galObj['ra']),
        'dec': '{:05f}'.format(galObj['dec']),
        'redshift': '{:05f}'.format(galObj['z']),
        'SDSS dr7 id': str(galObj['dr7objid']),
        'Common name': (
            str(galObj['IAUNAME'])
            if galObj['IAUNAME'] else 'Unknown'
        ),
        'NSA id': str(galObj['NSAID']),
        'Estimated distance': '{} * c / H_0'.format(galObj['ZDIST']),
        'Petrosian radius (degrees)': '{:05f}'.format(galObj['petrotheta']),
        'Run': str(galObj['RUN']),
        'Camcol': str(galObj['CAMCOL']),
        'Field': str(galObj['FIELD']),
        'Rerun': str(galObj['RERUN']),
        'Sersic axis ratio': '{:05f}'.format(galObj['SERSIC_BA']),
        'Url to view': "[View on SkyServer](+tab+http://skyserver.sdss.org/dr14/en/tools/chart/navi.aspx?ra={ra}&dec={dec}&opt=F)".format(
            ra=float(galObj['ra']),
            dec=float(galObj['dec'])
        ),
        '#isModelling': True,
        '#models': [
            {'frame': 0, 'model': 'GALAXY_BUILDER_DIFFERENCE'},
            {'frame': 2, 'model': 'GALAXY_BUILDER_MODEL'},
        ],
        **extra_metadata
    }


# galObj required keys:
# ra, dec, petrotheta for cutout
# dr7objid, IAUNAME, NSAID, ZDIST, RUN, CAMCOL, FIELD, RERUN, SERSIC_BA for
# metadata
def pipeline(galObj, outputFolder, subjectName, extra_metadata={},
             source_extraction_kwargs={}, verbose=False):
    make_dir(outputFolder)
    # get all frames out to 10 * petrosian radius
    # (with conversion to asec to amin)
    n = 10
    frame_list = queryFromRaDec(
        galObj['ra'],
        galObj['dec'],
        radius=galObj['petrotheta'] * (n / 60),
        limit=300
    )
    # Combine frames to improve signal to noise
    # reducedFrameList = reduceFrameList(frameList)

    dirToMontage = getFramesForMontage(frame_list, galObj['ra'], galObj['dec'])
    make_dir('montageOutputs')

    output_dir = 'montageOutputs/{}+{}'.format(galObj['ra'], galObj['dec'])
    file_loc = doMontage(dirToMontage, output_dir)

    fitsFile = fits.open(file_loc)

    # read it in and crop out around the galaxy
    cutout_result = scg.cutFits(
        fitsFile,
        galObj['ra'], galObj['dec'],
        size=(
            4 * galObj['petrotheta'] * u.arcsec,
            4 * galObj['petrotheta'] * u.arcsec
        )
    )
    if cutout_result is False and verbose:
        print('\t💀  \033[31mReturned False from image Data\033[0m')
        print('\tRa: {} Dec: {}'.format(galObj['ra'], galObj['dec']))
        raise(NoOverlapError(
            'Ra: {} Dec: {}'.format(galObj['ra'], galObj['dec'])
        ))
    image_data = cutout_result.data

    # Use source extractor to identify objects TODO proper deblending
    try:
        background = fitsFile[2].data[0][0]
    except IndexError:
        # mosaiced images don't have the background, pass a blank array
        background = None
    objects, segmentation_map = csf.sourceExtractImage(
        image_data,
        background,
        **source_extraction_kwargs,
    )
    # create a true/false masking array
    mask = csf.maskArr(image_data, segmentation_map, objects[-1][0] + 1)

    maskedImageData = copy(image_data)
    maskedImageData[mask] = 0

    # Now we find the PSF
    psf = scg.getPSF(
        (galObj['ra'], galObj['dec']),
        frame_list[0],
        fitsFile,
        fname="{}/psf_{}.png".format(outputFolder, subjectName)
    )
    c = 20
    # crop out most of the 0-ish stuff
    psfCut = psf[c:-c, c:-c]
    # normalise so we don't lose flux
    psfCut = psfCut / sum(psfCut)
    # psfs += [[str(i) for i in psfCut.reshape(psfCut.size)]]

    # generate the model json
    model = makeModel(maskedImageData, psfCut)

    # and the difference json
    difference = makeDifference(maskedImageData, psfCut, mask)

    # and the metadata
    metadata = makeMetadata(galObj, extra_metadata)

    # apply an asinh stretch and save the image to the outfolder
    resizeTo = (512, 512)
    csf.saveImage(
        csf.stretchArray(maskedImageData[:, ::-1]),
        fname="{}/image_{}.png".format(outputFolder, subjectName),
        resize=True,
        size=resizeTo
    )

    # now save the model json
    modelFileName = '{}/model_{}.json'.format(outputFolder, subjectName)
    with open(modelFileName, 'w') as f:
        json.dump(model, f)

    # write out the difference
    diffFileName = '{}/difference_{}.json'.format(outputFolder, subjectName)
    with open(diffFileName, 'w') as f:
        json.dump(difference, f)

    # and the metadata!
    metaFileName = '{}/metadata_{}.json'.format(outputFolder, subjectName)
    with open(metaFileName, 'w') as f:
        json.dump(metadata, f)

    return [
        "{}/image_{}.png".format(outputFolder, subjectName),
        modelFileName,
        diffFileName,
        metaFileName
    ]
