import sys
import numpy as np
import argparse
import subjectCreator as subC
from astropy import log
from tqdm import tqdm

log.setLevel('ERROR')

parser = argparse.ArgumentParser(
    description=(
        'Creat Zooniverse Subject data for a given Galaxy Builder subject set'
    )
)
parser.add_argument('n', metavar='N', type=int, help='Subject Set number')

args = parser.parse_args()

try:
    set_indices = np.load('original_subjectSetIndices/subjectSet{}.npy'.format(args.n))
except FileNotFoundError:
    print('Could not find indices file')
    parser.help()
    sys.exit(0)

finalOrderedByRedshift = np.load('original_finalOrderedByRedshift.npy')

files, success, err, montageFailures = [], [], [], []
with tqdm(range(0, set_indices.shape[0]), desc='Making Subjects') as bar:
    for i in bar:
        gal = finalOrderedByRedshift[set_indices[i]]
        try:
            res = subC.pipeline(
                gal,
                outputFolder='subject_sets/set_{}'.format(args.n),
                subjectName='subject{}'.format(i),
                verbose=True,
            )
            if res:
                success += [i]
                files += [res]
            else:
                montageFailures += [i]
        except Exception as e:
            print(e, i, set_indices[i])
            err += [i]

print('Subjects made: {}'.format(len(success)))
print('Subjects failed: {}'.format(len(err)))
