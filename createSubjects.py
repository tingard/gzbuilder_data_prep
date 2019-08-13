import os
import json
import gzip
import shutil
from numpy import sum, max, array
import astropy.units as u
from astropy.io import fits
import sdssCutoutGrab as scg
import createSubjectsFunctions as csf


def main(objList, outFolder='subjects', metadataList=None):
    # objList is list of (ra, dec, petrotheta)s for target galaxies
    # make sure the outfolder exists
    if not os.path.exists(outFolder):
        os.mkdir(outFolder)
    # cycle through the input objects
    for i, (ra, dec, petrotheta) in enumerate(objList):
        # search by ra, dec
        print('🛰  Looking for galaxy at {}, {}'.format(ra, dec))
        frame = scg.queryFromRaDec(ra, dec)
        if not len(frame):
            print("💩  Couldn\'t find any galaxies")
            continue
        fileLoc = scg.getBandFits(frame)
        fitsFile = fits.open(fileLoc)
        # read it in and crop out around the galaxy
        imageData = scg.cutFits(
            fitsFile,
            ra, dec,
            size=(4 * petrotheta * u.arcsec, 4 * petrotheta * u.arcsec)
        )
        if imageData is False:
            print('\t💀  \033[31mReturned False from image Data\033[0m')
            print('\tRa: {} Dec: {}'.format(ra, dec))
            continue

        # Use source extractor to identify objects TODO proper deblending
        objects, segmentation_map = csf.sourceExtractImage(
            imageData,
            fitsFile[2].data[0][0]
        )
        # create a true/false masking array
        mask = csf.maskArr(imageData, segmentation_map, objects[-1][0] + 1)

        # create the masked image
        maskedImageData = imageData[:]
        maskedImageData[mask] = 0

        # apply an asinh stretch and save the image to the outfolder
        resizeTo = (512, 512)
        csf.saveImage(
            csf.stretchArray(maskedImageData[:, ::-1]),
            fname="{}/image_{}.png".format(outFolder, i),
            resize=True,
            size=resizeTo
        )
        # Now we find the PSF
        psf = scg.getPSF((ra, dec), frame, fitsFile)
        c = 20
        # crop out most of the 0-ish stuff
        psfCut = psf[c:-c, c:-c]
        # normalise so we don't lose flux
        psfCut = psfCut / sum(psfCut)

        # generate the model json
        model = {
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
        # and the difference json
        difference = {
            'psf': psfCut.tolist(),
            'psfWidth': psfCut.shape[1],
            'psfHeight': psfCut.shape[0],
            'mask': array(mask, dtype=float).tolist(),
            'imageData': (imageData / max(imageData)).tolist(),
            'width': imageData.shape[1],
            'height': imageData.shape[0],
            'imageWidth': resizeTo[0],
            'imageHeight': int(
                imageData.shape[0] / imageData.shape[1] * resizeTo[0]
            )
        }
        # and the subject metadata
        metadata = {
            '#originalBrightness': float(max(imageData)),
            'Ra': ra,
            'Dec': dec,
            'petrotheta': petrotheta,
            'SDSS_ID': int(frame.get('objID', 0)),
            '#isModelling': True,
            '#models': [
                {'frame': 1, 'model': 'GALAXY_BUILDER_DIFFERENCE'},
                {'frame': 2, 'model': 'GALAXY_BUILDER_MODEL'},
            ]
        }
        print('metadata:', metadata)
        # write out the model (saving a gzipped and non-gzipped version)
        modelFileName = '{}/model_{}.json'.format(outFolder, i)
        with open(modelFileName, 'w') as f:
            json.dump(model, f)
        with open(modelFileName, 'rb') as f_in, \
                gzip.open(modelFileName + '.gz', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

        # write out the difference
        diffFileName = '{}/difference_{}.json'.format(outFolder, i)
        with open(diffFileName, 'w') as f:
            json.dump(difference, f)
        with open(diffFileName, 'rb') as f_in, \
                gzip.open(diffFileName + '.gz', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

        metaFileName = '{}/metadata_{}.json'.format(outFolder, i)
        with open(metaFileName, 'w') as f:
            json.dump(metadata, f)


lucyGals = (
    (160.65883, 23.95189, 22.3782),  # In beta
    (119.06931414139731, 11.662177891345076, 25.510479),
    (236.14108, 10.29315, 27.98376),
    (216.53496423791432, 5.237946739507734, 27.81166),  # In beta
    (248.7370584891894, 25.69259397115592, 26.418024),
    (239.5076772567969, 14.963535027843466, 26.63623),
    (178.44322153341767, 10.40313979646676, 18.855566),  # In beta
)

betaGals = [
    [41.641447536432004, -0.24713745138021195, 60.431831],
    [141.84790761271455, 30.44079266056109, 29.399462],
    [160.91220560347105, 14.871805138233954, 75.733269],
    [168.16590912088591, 9.0558683286936947, 40.021893],
    [189.15865415391502, 24.42881785060105, 12.993234],
    [196.30898397322505, 31.999692657501768, 31.009161],
    [205.16589100978419, 54.332963367499858, 59.585342],
    [213.10179825210173, 18.411757808899239, 27.811661],
    [216.16885767118126, 26.625148863396365, 25.35359],
    [230.75661842627977, -1.3471502025661983, 53.713493],
    [259.27702903505087, 40.844937775866171, 63.939758],
] + [lucyGals[0], lucyGals[3], lucyGals[6]]

if __name__ == '__main__':
    main(betaGals)
