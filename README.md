# trajectoread monitor 

## Overview
The monitor scripts are designed to track and analyze Stanford Sequencing Center operations on DNAnexus.

- monthly_seq_stats.py: Collect metrics on monthly sequencing output based on DNAnexus records and lane.html files generated by bcl2fastq2. On completion it will call "monthly_seq_stats.R" to generate summary plots. Designed to be called monthly by a Cron job.
     
- monthly_seq_stats.R: Create R plots describing monthly sequencing output. Plots are generated from files created by monthly_seq_stats.py.

## Output

![Sequencing Center Output 2013-2016](https://cloud.githubusercontent.com/assets/14796101/21828654/8415c95a-d746-11e6-93c5-1b6abbb5d384.png)

## Usage
To run the python script to gather data and then automatically generate plots:

    $ python monthly_seq_stats.py -y 2016 -m 12 -o seq_stats.txt
Or to run the Rscript indepedently:
    
    $ Rscript monthly_seq_stats.R -f stats_files/seq_stats.txt -o seq_stats

### Crontab

    0       12      1       *       *       bash -l -c "python /usr/local/autocopy/cron_scripts/monthly_seq_stats.py -c

