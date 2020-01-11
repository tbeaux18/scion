#!/usr/bin/env python3
"""
@author: Timothy Baker
@version: 1.0.0

zumi_config_builder.py

Must include yaml representer to preserve order
zUMI requires all paths to be absolute; will need to take this into account
for creating the barcode text file and trimmed fastq samples
"""


from collections import OrderedDict
import yaml
from yaml.resolver import Resolver

def represent_dictionary_order(self, dict_data):
    """ instantiates yaml dict mapping """
    return self.represent_mapping('tag:yaml.org,2002:map', dict_data.items())

def setup_yaml():
    """ adds the representer to the yaml instance """
    yaml.add_representer(OrderedDict, represent_dictionary_order)

    # remove resolver entries for On/Off/Yes/No
    for char_bool in "OoYyNn":
        if len(Resolver.yaml_implicit_resolvers[char_bool]) == 1:
            del Resolver.yaml_implicit_resolvers[char_bool]
        else:
            Resolver.yaml_implicit_resolvers[char_bool] = [
                x for x in Resolver.yaml_implicit_resolvers[char_bool]
                if x[0] != 'tag:yaml.org,2002:bool']


class ZumiConfigBuilder:
    """
    ZumiConfigBuilder object instantiates the appropriate dictionaries for
    easy value updates in the format provided by the zUMIs pipeline.

    Order matters for how the yaml file is produced. This is a hacked way of
    getting the proper format for the correct structure so zUMIs will work.

    Attributes:
        top_yaml_dict
        sequence_file_dict
        reference_dict
        filter_cutoff_dict
        barcode_dict
        counting_opts_dict

    Instance Methods:
        set_nested_dict
        return_top_yaml
        return_sequence_dict
        return_reference_dict
        return_filter_dict
        return_barcode_dict
        return_count_dict
        write_2_yaml

    """

    def __init__(self):

        self.top_yaml_dict = OrderedDict([
            ('project', None),
            ('sequence_files', None),
            ('reference', None),
            ('out_dir', None),
            ('num_threads', None),
            ('mem_limit', None),
            ('filter_cutoffs', None),
            ('barcodes', None),
            ('counting_opts', None),
            ('make_stats', "yes"),
            ('which_Stage', None),
            ('samtools_exec', 'samtools'),
            ('Rscript_exec', 'Rscript'),
            ('STAR_exec', 'STAR'),
            ('pigz_exec', 'pigz')
        ])

        self.sequence_file_dict = OrderedDict([
            ('file1', OrderedDict([
                ('name', None), # update from SampleSheetParser
                ('base_definition', ['BC(7-12)', 'UMI(1-6)'])
            ])
            ),
            ('file2', OrderedDict([
                ('name', None), # update from SampleSheetParser
                ('base_definition', ['cDNA(1-50)'])
            ])
            )
        ])

        self.reference_dict = OrderedDict([
            ('STAR_index', None), # update from SampleSheetParser
            ('GTF_file', None), # update from SampleSheetParser
            ('additional_files', None),
            ('additional_STAR_params', None)
        ])

        self.filter_cutoff_dict = OrderedDict([
            ('BC_filter', OrderedDict([
                ('num_bases', None), # update from SampleSheetParser
                ('phred', None) # update from SampleSheetParser
            ])
            ),
            ('UMI_filter', OrderedDict([
                ('num_bases', None), # update from SampleSheetParser
                ('phred', None) # update from SampleSheetParser
            ])
            )
        ])

        self.barcode_dict = OrderedDict([
            ('barcode_num', None),
            ('barcode_file', None), # file path from SampleSheetParser
            ('automatic', "yes"),
            ('BarcodeBinning', None), # int from SampleSheetParser
            ('nReadsperCell', 100)
        ])

        self.counting_opts_dict = OrderedDict([
            ('introns', "yes"),
            ('downsampling', 0),
            ('strand', 0),
            ('Ham_Dist', None), # must stay at 0 for now until zUMI fixes multi-thread issue
            ('velocyto', "no"),
            ('primaryHit', "yes"),
            ('twoPass', "yes")
        ])

    def update_top_level(self, run_name, out_path, threads, stage):
        """ updates top level values """
        self.top_yaml_dict['project'] = run_name
        self.top_yaml_dict['out_dir'] = out_path
        self.top_yaml_dict['num_threads'] = threads
        self.top_yaml_dict['which_Stage'] = stage

    def update_file_names(self, fastq_read1, fastq_read2):
        """ updates fastq file name values """
        self.sequence_file_dict['file1']['name'] = fastq_read1
        self.sequence_file_dict['file2']['name'] = fastq_read2

    def update_reference_files(self, star_index, gtf_file, add_params):
        """ updates reference dict values """
        self.reference_dict['STAR_index'] = star_index
        self.reference_dict['GTF_file'] = gtf_file
        self.reference_dict['additional_STAR_params'] = add_params

    def update_filter_cutoffs(self, bc_num, bc_phred, umi_num, umi_phred):
        """ updates filter cutoffs values """
        self.filter_cutoff_dict['BC_filter']['num_bases'] = bc_num
        self.filter_cutoff_dict['BC_filter']['phred'] = bc_phred
        self.filter_cutoff_dict['UMI_filter']['num_bases'] = umi_num
        self.filter_cutoff_dict['UMI_filter']['phred'] = umi_phred

    def update_barcodes(self, barcode_path, barcode_ham):
        """ updates barcode options """
        self.barcode_dict['barcode_file'] = barcode_path
        self.barcode_dict['BarcodeBinning'] = barcode_ham

    def update_count_opts(self, umi_ham):
        """ updates umi counting opts """
        self.counting_opts_dict['Ham_Dist'] = umi_ham

    def set_nested_dict(self):
        """ setting nested keys with their respective values """
        # update methods should be ran before this method
        self.top_yaml_dict['sequence_files'] = self.sequence_file_dict
        self.top_yaml_dict['reference'] = self.reference_dict
        self.top_yaml_dict['filter_cutoffs'] = self.filter_cutoff_dict
        self.top_yaml_dict['barcodes'] = self.barcode_dict
        self.top_yaml_dict['counting_opts'] = self.counting_opts_dict

    def write_yaml(self, zumi_config_path):
        """ writes the top level yaml dict to yaml file """

        # zumi_yaml_full_path = current_dir + '/' + basename + '_zumi_config.yaml'
        with open(zumi_config_path, 'w') as outfile:
            yaml.dump(
                self.top_yaml_dict,
                outfile,
                default_flow_style=False,
                default_style=None
            )
        return True

def main():
    """ run main for testing """
    pass

if __name__ == '__main__':
    main()
