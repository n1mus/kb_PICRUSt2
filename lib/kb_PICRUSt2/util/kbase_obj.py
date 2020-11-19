import shutil
import logging
import pandas as pd
import numpy as np
import os
import sys
import gzip
import json
import functools

from .dprint import dprint
from .config import var
from .error import * # custom Exceptions


pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1000)
pd.set_option('display.max_colwidth', 20)


####################################################################################################
####################################################################################################
MISSING_VALS = [None, '', 'nan', 'None', 'null', np.nan] # Below functions should implicitly handle these missing values

####################################################################################################
####################################################################################################
def replace_missing(a: np.ndarray, rep=np.nan):
    '''
    Use np.ndarray because
    (1) Easier to detect `None`
    (2) Common denominator
    '''    
    if np.issubdtype(a.dtype, int) and rep in [np.nan, None]:
        raise Exception('Cannot assign special missing values to numpy integer array')

    for missing in MISSING_VALS:
        """
        dprint( 
                'missing',
                'type(missing)',
                'a', 
                'a.dtype', 
                'a == np.array(missing, dtype=object)', 
                run={**globals(),**locals()}
                )
        """
        try:
            a[a == np.array(missing, dtype=object)] = rep
        except:
            raise Exception(missing)

    return a

####################################################################################################
####################################################################################################
def get_num_missing(a: np.ndarray):
    num = 0
    for missing in MISSING_VALS:
        
        """
        dprint( 
                'missing',
                'type(missing)',
                'a', 
                'a.dtype', 
                'a == np.array(missing, dtype=object)', 
                run={**globals(),**locals()}
                )
        """
        
        num = num + (a == np.array(missing, dtype=object)).sum() # cast missing as obj np.array to get element-wise comparison
    
    return num


####################################################################################################
####################################################################################################
def as_numeric(a: np.ndarray, rep=np.nan, dtype=float):
    '''
    Warning: NumPy integer arrays do not support missing values (np.nan etc.)
    '''

    if np.issubdtype(dtype, int) and rep in [np.nan, None]:
        raise Exception('Cannot assign special missing values to numpy integer array')

    if np.issubdtype(a.dtype, int) and rep in [np.nan, None]:
        a = a.astype(float) # cast to float first to allow assigning special values

    a = replace_missing(a)

    try:
        a = a.astype(dtype)
        return a
    except ValueError:
        return None

    raise


####################################################################################################
####################################################################################################
####################################################################################################

class AmpliconMatrix:

    def __init__(self, upa):
        self.upa = upa
        self._get_obj()


    def _get_obj(self):
        logging.info('Loading AmpliconMatrix object')

        obj = var.dfu.get_objects({
            'object_refs': [self.upa]
            })

        self.name = obj['data'][0]['info'][1]
        self.row_attrmap_upa = obj['data'][0]['data'].get('row_attributemapping_ref')
        self.obj = obj['data'][0]['data']
        
        if 'run_dir' in var: # comment directory with AmpMat name. optional since unit tests may not have run_dir
            dprint('touch %s' % os.path.join(var.run_dir, '#' + self.name), run='cli')


    def to_seq_abundance_table(self, flpth):
        '''
        Prerequisite: validate first with `validate_seq_abundance_data`
        '''
        logging.info(f"Writing sequence abundance table to %s" % flpth)

        data = as_numeric(np.array(self.obj['data']['values']), rep=np.nan, dtype=float) # Rep as floats so you can accept any numeric/missing then write with %g 
                                                                # np integer arrays don't take missing values
        row_ids = self.obj['data']['row_ids']
        col_ids = self.obj['data']['col_ids']

        data = pd.DataFrame(
            data, 
            index=row_ids, # ASV Ids 
            columns=col_ids # sample names
            )
        data.index.name = var.amplicon_header_name
        data.to_csv(flpth, sep='\t', float_format='%g')

    def to_fasta(self, flpth):
        fetched_flpth = var.gapi.fetch_sequence(self.upa)
        shutil.copyfile(fetched_flpth, flpth)


    def validate_seq_abundance_data(self):
        '''
        Can't be all missing
        Should be count data, missing allowed

        Because of KBase float types, which this is composed of,
        don't have to worry about complex, inf, etc.
        '''
        a = np.array(self.obj['data']['values'])

        # Can't be all missing
        if get_num_missing(a) == a.size:
            raise ValidationException(
                'Input AmpliconMatrix cannot be all missing values'
            )

        a = as_numeric(a, rep=np.nan, dtype=float) # cast as float because that can pick up float-like strings like '1.0'

        base_msg = (
            'Input AmpliconMatrix must have count data (missing values allowed). '
        )

        # Numeric, missing allowed
        if a is None:
            raise ValidationException(
                base_msg + 'Non-numeric detected'
            )

        #dprint('a', 'a.dtype', 'as_numeric(a, dtype=float)', run={**globals(), **locals()})

        # Integer
        allclose_ = functools.partial(
            np.allclose, 
            equal_nan=True, # NaNs equal 
            atol=1e-8, # same
            rtol=1e-8, # decrease otherwise 4.00001 == 4
        )
        a_round = np.round(a)
        if not allclose_(a, a_round) or not allclose_(a_round, a):
            #dprint('a', 'a_round', 'allclose_(a, a_round)', 'allclose_(a_round, a)', run={**globals(), **locals()})
            raise ValidationException(
                base_msg + 'Non-integer detected after attemping to round/truncate'
            )

        # Gte 0
        if np.any(a_round < 0):
            raise ValidationException(
                base_msg + 'Negative value detected'
            )


    def _map_id2attr_ids(self, id2attr, axis='row'):
        '''
        Parameters
        ----------
        id2attr - AmpliconMatrix row_ids to attribute you want to give AttributeMapping


        Behavior
        --------
        Swap out ids in id2attr so they end up mapping AttributeMapping ids to attributes
        '''
        if f'{axis}_attributemapping_ref' not in self.obj:
            raise Exception(
                'Trying to map AmpliconMatrix %s_ids to %s AttributeMapping ids '
                "when AmpliconMatrix doesn't have %s AttributeMapping"
                % (axis, axis, axis)
            )
        elif f'{axis}_mapping' not in self.obj:
            msg = (
                'Dude this object has a %s_attributemapping_ref '
                'and needs a %s_mapping. Letting it slide for now.'
                % (axis, axis)
            )
            logging.warning(msg)
            var.warnings.append(msg)
            return id2attr

        id2attr = {
            self.obj[f'{axis}_mapping'][id]: attr
            for id, attr in id2attr.items()
        }

        return id2attr


    def save(self, name=None):
        logging.info('Saving AmpliconMatrix')

        upa_new = var.gapi.save_object({
            'obj_type': 'KBaseMatrices.AmpliconMatrix', # TODO version
            'obj_name': name if name is not None else self.name,
            'data': self.obj,
            'workspace_id': var.params['workspace_id'],
        })['obj_ref']

        return upa_new



####################################################################################################
####################################################################################################
####################################################################################################

class AttributeMapping:

    def __init__(self, upa, amp_mat):
        '''
        Needs amp_mat for root ref and id mapping
        '''
        self.upa = upa
        self.amp_mat = amp_mat
        self._get_obj()

    def _get_obj(self):
        logging.info('Loading AttributeMapping object')

        obj = var.dfu.get_objects({
            'object_refs': ['%s;%s' %(self.amp_mat.upa, self.upa)]
            })

        self.name = obj['data'][0]['info'][1]
        self.obj = obj['data'][0]['data']


    def map_update_attribute(self, ind: int, id2attr: dict, map_ids_first=True):
        '''
        Update attribute at index `ind` using mapping `id2attr`
        '''
        if map_ids_first is True:
            id2attr = self.amp_mat._map_id2attr_ids(id2attr)

        for id, attr in id2attr.items():
            self.obj['instances'][id][ind] = attr


    def get_attribute_slot_warn(self, attribute, source) -> int:
        '''
        Get attribute slot matching both `attribute` and `source`
        
        Return the index of that slot, 
        which corresponds to both obj['attributes'] and each list in obj['instances']

        If slot matching both `attribute` and `source` does not exist
        * if `create=True` add slot for it
        * if `create=False` return -1
        '''
        
        # check if already exists
        for ind, attr_d in enumerate(self.obj['attributes']):
            if attr_d['attribute'] == attribute and attr_d['source'] == source:
                msg = (
                    'Overwriting attribute `%s` with source `%s` '
                    'in row AttributeMapping with name `%s`'
                    % (attribute, source, self.name)
                )
                logging.warning(msg)
                var.warnings.append(msg)
                return ind

        # append slot to `attributes`
        self.obj['attributes'].append({
            'attribute': attribute,
            'source': source,
            })

        # append slots to `instances` 
        for attr_l in self.obj['instances'].values():
            attr_l.append('')

        return len(attr_l) - 1


    def save(self):
        logging.info('Saving AttributeMapping')
        
        info = var.dfu.save_objects(
            {'id': var.params['workspace_id'],
             "objects": [{
                 "type": "KBaseExperiments.AttributeMapping", # TODO version
                 "data": self.obj,
                 "name": self.name,
                 "extra_provenance_input_refs": [self.upa]
             }]})[0]

        upa_new = "%s/%s/%s" % (info[6], info[0], info[4])

        return upa_new



####################################################################################################
####################################################################################################
####################################################################################################
class Report:
    '''For facilitating testing'''

    def __init__(self, upa):
        self.upa = upa
        self._get_obj()

    def _get_obj(self):

        logging.info('Loading AttributeMapping object')

        obj = var.dfu.get_objects({
            'object_refs': [self.upa]
            })

        self.obj = obj['data'][0]['data']




