# trajectoread monitor 

## Overview
The monitor scripts are designed to track and analyze Stanford Sequencing Center operations on DNAnexus.

- monthly_seq_stats.py: Collect metrics on monthly sequencing output based on DNAnexus records and lane.html files generated by bcl2fastq2. On completion it will call "monthly_seq_stats.R" to generate summary plots. Designed to be called monthly by a Cron job.

    ```0       12      1       *       *       bash -l -c "python /usr/local/autocopy/cron_scripts/monthly_seq_stats.py -c```

- monthly_seq_stats.R: Create R plots describing monthly sequencing output. Plots are generated from files created by monthly_seq_stats.py.


