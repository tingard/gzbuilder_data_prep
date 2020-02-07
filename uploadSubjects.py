from panoptes_client import SubjectSet, Subject, Project, Panoptes
import getpass
import os
import sys
import time
import json
import re
from tqdm import tqdm
# import numpy as np


def addLocation(subjectObj, locationDict):
    fPath = next(iter(locationDict.values()))
    if os.path.isfile(fPath):
        with open(fPath, 'rb') as f:
            media_data = f.read()
        media_type = next(iter(locationDict.keys()))
        subjectObj.locations.append(media_type)
        subjectObj._media_files.append(media_data)
        return True
    return False


def uploadSubjectToSet(project, subjectSet, locationsList, metadataList):
    print('Uploading {} subjects to {}'.format(len(locationsList), subjectSet))
    # imagePath can be string or list, metadata must be same dimension
    if not len(locationsList) == len(metadataList):
        print(
            '\t\033[31mInvalid arguments, locationsList and metadataList',
            'must have same length\033[0m'
        )
        return
    subjects = []
    for locations, meta in tqdm(zip(locationsList, metadataList)):
        # the json subjects need to be added in a more manual way so we can
        # specify a MIME type
        subjects.append(Subject())
        subjects[-1].links.project = project
        # comparison between model and image
        addLocation(subjects[-1], {'application/json': locations[1]})
        # actual galaxy image
        subjects[-1].add_location(locations[0])
        # and now just the model
        addLocation(subjects[-1], {'application/json': locations[2]})
        for k, v in meta.items():
            subjects[-1].metadata[k] = v
        try:
            subjects[-1].save()
        except RuntimeError:
            pass
    subjectSet.add(subjects)
    return subjectSet


def main(production=False):
    uname = input('Enter your username: ')
    pwd = getpass.getpass()
    Panoptes.connect(
        username=uname,
        password=pwd,
        # endpoint='https://panoptes-staging.zooniverse.org',
        admin=True
    )
    pId = 5733  # if production else 1820
    project = Project.find(pId)
    subject_set = SubjectSet()
    subject_set.links.project = project
    subject_set.display_name = 'Test_subject_set_' + str(int(time.time()))
    subject_set.save()

    loc = os.path.abspath(os.path.dirname(__file__))
    subjects = os.listdir(loc + '/subjects')
    images, differences, model, metadata = [
        sorted((
            int(re.match(r'{}_([0-9]+)\.(?:json|png)$'.format(s), i).group(1))
            for i in subjects
            if re.match(r'{}_([0-9]+)\.(?:json|png)$'.format(s), i)
        ))
        for s in ('difference', 'image', 'model', 'metadata')
    ]
    if not images == differences == model == metadata:
        print(
            'Images, differences, model and metadata '
            + 'must all have same length'
        )

    # TODO: change subject directory structure to be more efficient
    #       (not having 12,000+ files in a folder...)
    for i in images:
        try:
            with open('{}/subjects/metadata_{}.json'.format(loc, i)) as f:
                metadata = json.load(f)
        except IOError:
            metadata = {}
        subject_set = uploadSubjectToSet(
            project, subject_set,
            [[j.format(loc, i) for j in (
                '{}/subjects/image_{}.png',
                '{}/subjects/difference_{}.json',
                '{}/subjects/model_{}.json'
            )]],  # locations
            [metadata],
        )


def create_subject_set(folder_name, set_name='test_subject_set'):
    nSubjects = len([
        f for f in os.listdir(folder_name)
        if '.DS' not in f
    ])//4
    subject_names = [
        i.group(1)
        for i in (
            re.match(r'image_(.*?).png', f)
            for f in os.listdir('affirmation_subjects')
        )
        if i is not None
    ]
    files = [
        (
            '{}/image_{}.png'.format(folder_name, name),
            '{}/difference_{}.json'.format(folder_name, name),
            '{}/model_{}.json'.format(folder_name, name),
            '{}/metadata_{}.json'.format(folder_name, name),
        )
        for name in subject_names
    ]
    assert all(os.path.exists(j) for i in files for j in i), 'Missing files!'
    uname = input('Enter your username: ')
    pwd = getpass.getpass()
    Panoptes.connect(
        username=uname,
        password=pwd,
        admin=True
    )
    pId = 5590
    project = Project.find(pId)
    subject_set = SubjectSet()
    subject_set.links.project = project
    subject_set.display_name = set_name
    subject_set.save()
    metadata_list = []
    for fs in files:
        try:
            with open(fs[3]) as metaF:
                metadata = json.load(metaF)
        except IOError:
            metadata = {}
        metadata_list.append(metadata)
    subject_set = uploadSubjectToSet(
        project, subject_set,
        [i[:3] for i in files],
        metadata_list,
    )


if __name__ == '__main__':
    try:
        folder_name = sys.argv[1]
        set_name = sys.argv[2]
    except IndexError:
        print('Please provide the path to input files and a subject set name')
        sys.exit(0)
    create_subject_set(folder_name, set_name)
