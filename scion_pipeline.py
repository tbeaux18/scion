#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 01/11/2020

Need to add a lot more error handling if anything other than Filtering is chosen
"""

import re
import csv
import argparse
import pathlib
import collections
import shlex
import subprocess
import logging
import multiprocessing as mp
import pandas as pd
import yaml
from yaml.resolver import Resolver
from zumi_config_builder import ZumiConfigBuilder


logging.basicConfig(
    filename="main_pipeline.log",
    level=logging.DEBUG,
    format="%(asctime)s:%(levelname)s:%(message)s"
    )

# path of directory where main_pipeline.py runs
PARENT_PATH = pathlib.Path(__file__).resolve().parent
logging.info('parent path to file: %s', PARENT_PATH)

def arg_parser():
    """ Argument input from command line """

    parser = argparse.ArgumentParser(
        description='Runs the single-cell RNA-seq pipeline.'
    )
    parser.add_argument(
        '-e',
        '--exp-design',
        type=str,
        help='path/to/exp_design.csv'
    )
    parser.add_argument(
        '-s',
        '--sample-sheet',
        type=str,
        help='path/to/sample_info.csv'
    )

    return parser.parse_args()


def represent_dictionary_order(self, dict_data):
    """ instantiates yaml dict mapping """
    return self.represent_mapping('tag:yaml.org,2002:map', dict_data.items())


def setup_yaml():
    """ adds the representer to the yaml instance """
    yaml.add_representer(collections.OrderedDict, represent_dictionary_order)

    # remove resolver entries for On/Off/Yes/No
    for char_bool in "OoYyNn":
        if len(Resolver.yaml_implicit_resolvers[char_bool]) == 1:
            del Resolver.yaml_implicit_resolvers[char_bool]
        else:
            Resolver.yaml_implicit_resolvers[char_bool] = [
                x for x in Resolver.yaml_implicit_resolvers[char_bool]
                if x[0] != 'tag:yaml.org,2002:bool']



def run_fastqc(**kwargs):
    """ runs fastqc with the provided args from build_fastqc_args()
        params:
            kwargs : dict specified in build_fastqc_args
            --extract : qc log files will be unzipped
            --format : fastq default
        returns:
            None
    """

    fastqc_command = shlex.split(
        """fastqc
            --extract
            --threads {threads}
            --format fastq
            --outdir {fastqc_dir}
            {fastq_read_1}
            {fastq_read_2}""".format(**kwargs)
        )

    logging.info('fastqc command entered: %s', fastqc_command)

    with open(kwargs['fastqc_log'], 'ab+') as fastqc_stdout:
        try:
            subprocess.run(
                fastqc_command,
                bufsize=20,
                stdout=fastqc_stdout,
                stderr=fastqc_stdout
            )

            logging.info('fastqc command completed.')

        except subprocess.CalledProcessError as call_error:
            logging.exception('command failed: %s', call_error)



def run_cutadapt(**kwargs):
    """ runs cutadapt with provided arguments from build_cutadapt_args()
        params:
            kwargs : dict specified in build_cutadapt_args
        returns:
            None
    """

    cutadapt_command = shlex.split(
        """cutadapt
            -j {threads}
            -a {adapter_3}
            -A {adapter_5}
            -m {min_length_r1}:{min_length_r2}
            -o {trimmed_r1}
            -p {trimmed_r2}
            {fastq_read_1}
            {fastq_read_2}""".format(**kwargs)
        )

    logging.info("cutadapt args: %s", cutadapt_command)

    with open(kwargs['cutadapt_log'], 'ab') as cutadapt_stdout:
        try:
            subprocess.run(
                cutadapt_command,
                bufsize=20,
                stdout=cutadapt_stdout,
                stderr=cutadapt_stdout
            )

            logging.info('fastqc command completed')

        except subprocess.CalledProcessError as call_error:
            logging.exception('command failed: %s', call_error)


def build_star_index(**kwargs):
    """ takes kwargs and builds the STAR index

        params:
            threads : int for threading
            dmel_index_dir : output dir to put index files, needs to be created
                                prior to running STAR
            ref_fasta_file : absolute path to fasta file, grabs from SampleSheetParser
        returns:
            None
    """
    star_build_index_cmd = shlex.split(
        """STAR
            --runMode genomeGenerate
            --runThreadN {threads}
            --genomeDir {star_idx_dir}
            --genomeFastaFiles {ref_genome_fasta}""".format(**kwargs)
        )

    logging.info('STAR index build args: %s', star_build_index_cmd)

    with open(kwargs['star_idx_log'], 'ab+') as staridx_stdout:

        try:
            subprocess.run(
                star_build_index_cmd,
                bufsize=50,
                stdout=staridx_stdout,
                stderr=staridx_stdout
            )

            logging.info('star idx built successfully')

        except subprocess.CalledProcessError as call_error:
            logging.exception('star idx failed: %s', call_error)


def check_star_idx_exists(star_idx_dir):
    """ checks to see if all files in star idx are present """

    staridx_files = [
        'chrLength.txt',
        'chrName.txt',
        'Genome',
        'SA',
        'chrNameLength.txt',
        'chrStart.txt',
        'genomeParameters.txt',
        'SAindex'
    ]

    # boolean array
    file_result = []

    for file in staridx_files:

        if pathlib.Path(star_idx_dir).joinpath(file).exists():

            file_result.append(True)

        else:
            file_result.append(False)

    if all(file_result):
        return True

    return False



def run_zumi_pipeline(**kwargs):
    """ runs the zumi pipeline. Steps include filtering, alignment, and counting

        params:
            zumi_yaml : path/to/zumi_config.yaml file constructed from the
                        ZumiConfigBuilder
        returns:
            None
    """

    zumi_cmd = shlex.split(
        """bash {zumi_pipeline_path}/zUMIs-master.sh
            -y {zumi_yaml}""".format(**kwargs)
        )

    logging.info("zUMI args: %s", zumi_cmd)

    if pathlib.Path(kwargs['zumi_yaml']).is_file():

        with open(kwargs['zumi_stdout'], 'ab+') as zumi_output:
            try:
                subprocess.run(
                    zumi_cmd,
                    bufsize=50,
                    stdout=zumi_output,
                    stderr=zumi_output
                )

            except subprocess.CalledProcessError as call_error:
                logging.exception('zumi failed: %s', call_error)

    else:
        logging.info("zumi yaml file is not made.")



def parse_sample_exp_args(sample_sheet_path):
    """ parses sample info into dict
        param:
            sample_sheet_path (str) : sample info path from args
        return:
            sample_info_dict (dict) : dict of info to pass to pipeline
    """

    sample_info_fields_int = [
        'bc_filter_num_bases',
        'bc_filter_phred',
        'bc_ham_dist',
        'umi_filter_num_bases',
        'umi_filter_phred',
        'umi_ham_dist',
        'min_length_r1',
        'min_length_r2'
    ]

    pipeline_path = [
        'fastq_read_1',
        'fastq_read_2',
        'ref_genome_fasta',
        'ref_annotation_gtf',
        'zumi_pipeline_path'
    ]

    zumi_stage_options = ['Filtering', 'Mapping', 'Counting', 'Summarising']

    pipeline_regex = re.compile(r'^/pipeline/')

    # dict will hold all info and paths
    sample_info_dict = {}

    # need to add more checks to make sure no unsuspecting inputs
    with open(sample_sheet_path, 'r') as input_file:

        reader = csv.reader(input_file)

        for row in reader:

            # setting numbers in sampleinfo to integers
            if row[0] in sample_info_fields_int:
                sample_info_dict[row[0]] = int(row[1])

            # ensuring that /pipeline/ is first in path due to docker paths
            elif row[0] in pipeline_path:
                if pipeline_regex.match(row[1]):
                    sample_info_dict[row[0]] = row[1]
                    logging.info('/pipeline/ detected in path: %s', row[1])

                else:
                    logging.error('/pipeline/ is missing from path')

            # ensuring that a valid zumi stage was entered
            elif row[0] == 'zumi_start_stage':
                if row[1] in zumi_stage_options:
                    sample_info_dict[row[0]] = row[1]
                    logging.info('correct stage detected: %s', row[1])
                else:
                    print('Error: Filtering, Mapping, Counting, Summarising are the only options.')
                    logging.error('incorrect params entered.')

            else:
                sample_info_dict[row[0]] = row[1]

    return sample_info_dict


def add_paths_to_sample_dict(sample_info_dict, exp_design_path):
    """ adding paths to the main sample dict """

    path_array = [
        ('fastqc_dir', 'fastqc_output'),
        ('barcode_file', 'barcode_whitelist.txt'),
        ('star_idx_dir', 'dmel_star_idx_NOGTF')
    ]

    log_paths = [
        ('fastqc_log', 'fastqc_log.txt'),
        ('cutadapt_log', 'cutadapt_log.txt'),
        ('star_idx_log', 'star_idx_log.txt'),
        ('zumi_stdout', 'zumi_stdout.txt')
    ]

    host_threads = mp.cpu_count()
    logging.info('host threads avail: %s', host_threads)

    sample_info_dict['threads'] = int(host_threads - 5)
    logging.info('threads to use: %s', sample_info_dict['threads'])

    # adding necessary paths to pass to zumi yaml builder
    sample_info_dict['current_wd'] = pathlib.Path.cwd()


    for path_key, path_val in path_array:
        sample_info_dict[path_key] = str(PARENT_PATH.joinpath(path_val))


    for log_key, log_file in log_paths:
        sample_info_dict[log_key] = str(PARENT_PATH.joinpath('logs', log_file))


    sample_info_dict['trimmed_r1'] = str(
        sample_info_dict['current_wd'].joinpath(
            '.'.join([sample_info_dict['base_name'], 'trimmed', 'R1', 'fastq.gz'])
        )
    )
    sample_info_dict['trimmed_r2'] = str(
        sample_info_dict['current_wd'].joinpath(
            '.'.join([sample_info_dict['base_name'], 'trimmed', 'R2', 'fastq.gz'])
        )
    )

    sample_info_dict['zumi_output_dir'] = '_'.join(
        [sample_info_dict['base_name'], 'zumi_output']
    )

    sample_info_dict['zumi_yaml'] = sample_info_dict['current_wd'].joinpath(
        '_'.join(
            [sample_info_dict['base_name'], 'zumi_config.yaml']
        )
    )

    # load exp design in pandas path
    expdesign_df = pd.read_csv(exp_design_path)
    # need to ensure barcode_sequence is present
    with open(sample_info_dict['barcode_file'], 'w+') as barcode_output:
        barcode_output.write(
            '\n'.join(expdesign_df['barcode_sequence'].tolist())
        )
        logging.info('created barcode whitelist text file: %s', sample_info_dict['barcode_file'])


    # need to creating star index directory
    pathlib.Path(sample_info_dict['star_idx_dir']).mkdir(exist_ok=True)
    pathlib.Path(sample_info_dict['zumi_output_dir']).mkdir(exist_ok=True)


    logging.info('sample info before passing to zumi yaml build: %s', sample_info_dict)
    return sample_info_dict, expdesign_df




def build_zumi_yaml(sample_info_dict):
    """ building the zumi yaml config sheet """

    # create representer
    setup_yaml()
    zumi_config_obj = ZumiConfigBuilder()

    additional_star_params = "--limitOutSJcollapsed 3000000 --limitSjdbInsertNsj 3000000"

    zumi_config_obj.update_top_level(
        sample_info_dict['run_name'],
        str(sample_info_dict['current_wd'].joinpath(sample_info_dict['zumi_output_dir'])),
        sample_info_dict['threads'],
        sample_info_dict['zumi_start_stage']
    )

    zumi_config_obj.update_file_names(
        str(sample_info_dict['current_wd'].joinpath(sample_info_dict['trimmed_r1'])),
        str(sample_info_dict['current_wd'].joinpath(sample_info_dict['trimmed_r2']))
    )

    # need to better handle star index build path
    # setting reference paths for zumi output
    zumi_config_obj.update_reference_files(
        str(
            sample_info_dict['current_wd'].joinpath(
                sample_info_dict['star_idx_dir']
            )
        ),
        str(
            sample_info_dict['current_wd'].joinpath(
                sample_info_dict['ref_annotation_gtf']
            )
        ), # gtf path
        additional_star_params
    )

    # passing filter cutoffs
    zumi_config_obj.update_filter_cutoffs(
        int(sample_info_dict['bc_filter_num_bases']),
        int(sample_info_dict['bc_filter_phred']),
        int(sample_info_dict['umi_filter_num_bases']),
        int(sample_info_dict['umi_filter_phred']),
    )

    # need to handle barcode whitelist path better
    zumi_config_obj.update_barcodes(
        str(sample_info_dict['current_wd'].joinpath(sample_info_dict['barcode_file'])),
        int(sample_info_dict['bc_ham_dist']),
    )

    # updating collapsing ham distance
    zumi_config_obj.update_count_opts(
        int(sample_info_dict['umi_ham_dist'])
    )

    # creates the nested dicts in the proper order to be dumped to YAML file
    # must run first before writing to yaml file
    zumi_config_obj.set_nested_dict()

    # contained in root directory of the repo
    zumi_written_bool = zumi_config_obj.write_yaml(sample_info_dict['zumi_yaml'])

    if zumi_written_bool:
        logging.info('zumi config file successfully written')
    else:
        logging.info('no zumi config file was written')

    return sample_info_dict


def run_main_pipeline(sample_info_dict):
    """ function that runs all pipeline """

    print('running fastqc')
    run_fastqc(**sample_info_dict)


    if (pathlib.Path(sample_info_dict['trimmed_r1']).exists() and
            pathlib.Path(sample_info_dict['trimmed_r2']).exists()):
        logging.info('trimmed fastqs detected, skipping cutadapt')
    else:
        print('trimmed files not detected. running cutadapt')
        run_cutadapt(**sample_info_dict)


    if check_star_idx_exists(sample_info_dict['star_idx_dir']):
        logging.info('star idx files detected')

    else:
        print('building star index')
        build_star_index(**sample_info_dict)

    print('running zumi')
    run_zumi_pipeline(**sample_info_dict)

def main():
    """ main function """

    args = arg_parser()

    # setting up logging directory
    pathlib.Path('logs').mkdir(exist_ok=True)

    logging.info('sample sheet path: %s', args.sample_sheet)
    logging.info('exp design path: %s', args.exp_design)

    sample_info = parse_sample_exp_args(args.sample_sheet)

    sample_info_w_paths, exp_df = add_paths_to_sample_dict(sample_info, args.exp_design)

    new_sample_info = build_zumi_yaml(sample_info_w_paths)

    run_main_pipeline(new_sample_info)

if __name__ == '__main__':
    main()
