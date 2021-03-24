# coding: utf-8
# Distributed under the terms of the MIT License.
'''
Currently the default collection is: db.vasp.tasks
'''

import os
import shutil
import time
import logging
from pymongo import MongoClient


def get_mongodb_entries(dir_prefix):
    '''
    Find documents using the criterion
    
    Args:
        dir_prefix: part of the dir_name as the search criterion
    Return:
        retrieved documents
    '''
    db = MongoClient()
    criteria = {'dir_name': {'$regex': dir_prefix}}
    return list(db.vasp.tasks.find(criteria))


def field_modify(docs, mgd_field, rep_list):
    '''
    Change the field value of the documents.  

    Args:
        docs: documents to be modified
        field: field to be modified, e.g 'dir_name'
        values: two item list for the replace function, e.g. ['-', '_']
    Return:
        retrieved documents
    '''
    db = MongoClient()
    for doc in docs:
        mongo_id = doc['_id']
        new = doc[mgd_field].replace(rep_list[0], rep_list[1])
        db.vasp.tasks.update_one({'_id':mongo_id}, {"$set": {mgd_field: new}}, upsert=False)


def rerun_killed(docs):
    '''
    Rerun the unfinished jobs.  
    Should be done when the xml files were saved and previous folders were deleted.  
    Do the following:  
        1. Rename the saved xml, add suffix 'not_finish'
        2. Current contcar?

    Args:
        docs: documents to be modified

    Return:
        retrieved documents
    '''
    db = MongoClient()
    for doc in docs:
        if doc['state'] == 'killed':
            old_dir = doc['dir_name'].replace('pbs-master:', '')
            try:
                shutil.move(old_dir+'/CONTCAR', old_dir+'/POSCAR')
            except:
                logging.info('Cannot move {} CONTCAR to POSCAR'.format(old_dir))
            try:
                os.remove(old_dir+'/POTCAR')
            except:
                logging.info('Cannot remove {}/POTCAR'.format(old_dir))
            try:
                os.remove(old_dir+'/OUTCAR')
            except:
                logging.info('Cannot remove {}/OUTCAR'.format(old_dir))

            old_dir = old_dir.replace('/workspace/RZ', '/backup_1/xml_files')
            old_file = old_dir + '.xml'
            new_file = old_dir + '_not_finish.xml'
            try:
                shutil.move(old_file, new_file)
            except:
                logging.info('Cannot move {} to {}'.format(old_file, new_file))

            mongo_id = doc['_id']
            new = doc['dir_name'] + '_not_finish'
            db.vasp.tasks.update_one({'_id':mongo_id}, {"$set": {'dir_name': new}}, upsert=False)
        else:
            old_dir = doc['dir_name'].replace('pbs-master:', '')
            try:
                shutil.rmtree(old_dir)
            except:
                logging.info('Cannot remove {}'.format(old_dir))


