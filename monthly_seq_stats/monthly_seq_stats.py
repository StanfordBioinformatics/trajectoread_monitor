#!/usr/bin/env python

import os
import re
import sys
import pdb
import dxpy
import argparse
import datetime
import subprocess

from collections import defaultdict

def parse_lane_html(html_dxfile, lane_index, lane_name):
    html = html_dxfile.read()
    lines = html.split('\n')

    for i in range(0, len(lines)):
        lane_match = re.search('<td>%d</td>' % lane_index, lines[i])

        if lane_match:
            pf_clusters_line = lines[i + 1]
            perc_perfect_barcode_line = lines[i + 3]
            perc_one_mismatch_barcode_line = lines[i + 4]
            yield_line = lines[i + 5]
            perc_pf_line = lines[i + 6]
            perc_q30_bases_line = lines[i + 7]
            mean_quality_line = lines[i + 8]
            
            perc_perfect_barcode = parse_html_value(perc_perfect_barcode_line, lane_name)
            perc_one_mismatch_barcode = parse_html_value(perc_one_mismatch_barcode_line, lane_name)
            perc_pf = parse_html_value(perc_pf_line, lane_name)
            perc_q30_bases = parse_html_value(perc_q30_bases_line, lane_name)
            mean_quality = parse_html_value(mean_quality_line, lane_name)
            pf_clusters = parse_html_value(pf_clusters_line, lane_name)

            yield_mbases_raw = parse_html_value(yield_line, lane_name)
            yield_elements = yield_mbases_raw.split(',')
            yield_mbases = ''.join(yield_elements) # Remove commas from mbase value

            lane_output = {
                           "yield_mbases": yield_mbases,
                           "perc_pf": perc_pf,
                           "pf_clusters": pf_clusters,
                           "perc_perfect_barcode": perc_perfect_barcode,
                           "perc_one_mismatch_barcode": perc_one_mismatch_barcode,
                           "perc_q30_bases": perc_q30_bases,
                           "mean_quality": mean_quality
                          }

            return (lane_output)

def parse_html_value(html_line, lane_name):
    html_match = re.search('<td>(.+)</td>', html_line)
    if html_match:
        value = html_match.group(1)
        return value
    else:
        print 'Error: could not find value in html line: %s, Lane %s' % (html_line, lane_name)
        return 'NA'

def classify_instrument(name):
    '''Get sequencer type based on current instrument names.
    Used in the case of old sequencing records without sequencer type.
    
    Args:
        name (str): Name of the sequencing machine in the LIMS
    '''

    hiseq_4000s = ['Gadget', 'Cooper']
    # HiSeq 2500 : Software Version: HCS 2.0.5
    hiseq_2500s = ['Briscoe', 'Marple']
    miseqs = ['Holmes', 'Spenser']

    if name in hiseq_4000s:
        sequencer_type = 'HiSeq_4000'
    elif name in miseqs:
        sequencer_type = 'MiSeq'
    elif name in hiseq_2500s:
        sequencer_type = 'HiSeq_2500'
    else:
        sequencer_type = 'HiSeq_2000'
    return sequencer_type

def parse_args(args):

    parser = argparse.ArgumentParser()
    parser.add_argument(
                        '-m', 
                        '--month', 
                        dest = 'month', 
                        type = int,
                        required = False,
                        help = 'Get metrics for all sequencing lanes processed on DNAnexus during this month.')
    parser.add_argument(
                        '-y', 
                        '--year', 
                        dest = 'year', 
                        type = int,
                        required = False,
                        help = 'Get metrics for sequencing lanes processed during this year.')
    parser.add_argument(
                        '-c',
                        '--cron',
                        dest = 'cron',
                        action = 'store_true',
                        default = False,
                        help = 'Use CRON logic to determine month to process.')
    parser.add_argument(
                        '-o',
                        '--outfile',
                        dest = 'outfile',
                        type = str,
                        default = 'seq_stats.txt',
                        help = 'Specify outfile to append results & use to ' + 
                               'generate figures.')
    args = parser.parse_args(args)

    if args.cron:
        print 'Processing monthly seq stats in CRON mode'
    elif args.year and args.month:
        print 'Processing monthly seq stats for {}-{}'.format(args.year, args.month)
    else:
        print 'ERROR: You did not specify year/month nor CRON option'
        parser.print_help()
        sys.exit()
    return(args)

def main():
    '''Get lane, base, and read stats for the given month.
    '''

    args = parse_args(sys.argv[1:])
    outfile = args.outfile

    if args.cron:
        now = datetime.datetime.now()
        year = now.year
        month = now.month - 1
        print 'Info: Collecting metrics for %d-%d' % (year, month)
    else:   
        year = args.year
        month = args.month
    monthly_outfile = '%d-%d_seq-stats.txt' % (year, month)

    monthly_metrics = defaultdict()
    '''
    monthly_metrics = {
                       'lane_count' : 0,
                       'base_count' : 0,
                       'read_count' : 0
                      }
    '''

    # Dev: NEED TO QC THIS
    after_date = '%d-%d-01' % (year, month)
    before_date = '%d-%d-01' % (year, int(month + 1))
    print 'After: %s' % after_date
    print 'Before: %s ' % before_date

    monthly_records = dxpy.find_data_objects(
                                             classname = 'record',
                                             project = 'project-BY82j6Q0jJxgg986V16FQzjx',
                                             folder = '/',
                                             typename = 'SCGPMRun',
                                             created_after = after_date,
                                             created_before = before_date)

    MOUT = open(monthly_outfile, 'w')
    for record in monthly_records:
        sequencer_type = None
        
        print record['id']
        dxrecord = dxpy.DXRecord(record['id'], record['project']) 
        lane_details = dxrecord.get_details()
        lane_properties = dxrecord.get_properties()

        try:
            production = lane_properties['production']
            if production == 'false':
                continue
        except:
            print 'Skipping: not production'
            continue

        dxproject = lane_details['laneProject']
        lane_index = int(lane_details['lane'])
        run_name = str(lane_details['run'])
        paired_end = bool(lane_properties['paired_end'])

        #try:
        #    sequencer_type = str(lane_properties['sequencer_type'])
        #except:
        seq_instrument = str(lane_properties['seq_instrument'])
        sequencer_type = classify_instrument(seq_instrument)
        #sequencer_type = sequencer_type.replace(' ', '_')

        # Create monthly_records defaultdict to store cumulative stats
        if not sequencer_type in monthly_metrics.keys():
            monthly_metrics[sequencer_type] = {
                                               'lane_count' : 0,
                                               'base_count' : 0,
                                               'read_count' : 0
                                              }

        lane_name = '%s_L%d' % (run_name, lane_index)
        print 'Processing %s' % lane_name

        try:
            html_file = dxpy.find_one_data_object(
                                                  classname = 'file',
                                                  name = '*.lane.html',
                                                  name_mode = 'glob',
                                                  project = dxproject,
                                                  folder = '/stage0_bcl2fastq/miscellany',
                                                  more_ok = True,
                                                  zero_ok = False)
        except:
            print 'Warning: Could not get lane.html file. Skipping'
            continue

        html_dxfile = dxpy.DXFile(html_file['id'], html_file['project'])
        lane_metrics = parse_lane_html(html_dxfile, lane_index, lane_name)

        pf_clusters = int(lane_metrics['pf_clusters'].replace(',',''))
        if paired_end:
            read_count = pf_clusters * 2
        else:
            read_count = pf_clusters

        mbase_count = int(lane_metrics['yield_mbases'].replace(',',''))
        base_count = mbase_count * 1000000

        monthly_metrics[sequencer_type]['lane_count'] += 1
        monthly_metrics[sequencer_type]['base_count'] += base_count
        monthly_metrics[sequencer_type]['read_count'] += read_count

        # Write individual record data to the monthly out file.
        # One record/lane per line.
        #pdb.set_trace()
        mout_str = (
                    '{}\t'.format(year) +
                    '{}\t'.format(month) +
                    '{}\t'.format(run_name) +
                    '{}\t'.format(lane_index) +
                    '{}\t'.format(read_count) +
                    '{}\t'.format(base_count) +
                    '{}\n'.format(sequencer_type))
        MOUT.write(mout_str)
    MOUT.close()

    # Add header to new outfile
    if not os.path.isfile(outfile):
        with open(outfile, 'w') as OUT:
            OUT.write('Year\tMonth\tLane_Count\tRead_Count\tBase_Count\tSeq_Type\n')
    # Write monthly metrics to outfile
    with open(outfile, 'a') as OUT:
        '''Old string formatting method
        out_str = '%d\t%d\t%d\t%d\t%d\n' % (
                                            year,
                                            month,
                                            monthly_metrics['lane_count'],
                                            monthly_metrics['read_count'],
                                            monthly_metrics['base_count'],
                                            sequencer_type)
        '''
        #pdb.set_trace()
        for sequencer_type in monthly_metrics.keys():
            out_str = (
                       '{}\t'.format(year) +
                       '{}\t'.format(month) +
                       '{}\t'.format(monthly_metrics[sequencer_type]['lane_count']) +
                       '{}\t'.format(monthly_metrics[sequencer_type]['read_count']) +
                       '{}\t'.format(monthly_metrics[sequencer_type]['base_count']) +
                       '{}\n'.format(sequencer_type))
            OUT.write(out_str)
    out_prefix = outfile.split('.')[0]
    #subprocess.call(['Rscript','monthly_seq_stats.R','-f', outfile, '-o', out_prefix])


if __name__ == "__main__":
    main()