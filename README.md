# plot-nanopore-run
Generates plots reporting on nanopore runs

## General usage

usage: plot_nanopore_run.py [-h] -s/--summary SEQ_SUMMARIES [SEQ_SUMMARIES ...] [-o/--output OUTPUT_PREFIX]
                -r/--run_id RUN_ID

Calculate nanopore run report

optional arguments:
  -h, --help            show this help message and exit
  -s SEQ_SUMMARIES [SEQ_SUMMARIES ...], --summary SEQ_SUMMARIES [SEQ_SUMMARIES ...]
                        sequencing summary tab-separated file(s)
  -o OUTPUT_PREFIX, --output OUTPUT_PREFIX
                        prefix of output directory
  -r RUN_ID, --run_id RUN_ID
                        run id will be displayed as part of titles of plots.
                        NOTE: wrap in quotes if spaces included.

### Input

Example sequence summary file:

filename    read_id    run_id    channel    start_time    duration    num_events    passes_filtering    template_start    num_events_template    template_duration    num_called_template    sequence_length_template    mean_qscore_template    strand_score_template    calibration_strand_genome_template    calibration_strand_identity_template    calibration_strand_accuracy_template    aligned_speed_bps_template
p_102_20180123_0004A30B001B5D4B_PH_p_102_00_sequencing_run_180123_Minden9693_run3_21050_read_10208_ch_1531_strand.fast5    075fb561-9f3c-4248-8382-1d79df8977ba    75052db6272cde941f37a391e84c912fd8f459c9    1531    6548.61575    1.11875    895    False    0.0    895    1.11875    895    443    5.091    -0.0013    filtered_out    -1.0    -1.0    0.0
p_102_20180123_0004A30B001B5D4B_PH_p_102_00_sequencing_run_180123_Minden9693_run3_21050_read_10194_ch_1531_strand.fast5    12e7dd12-ab0b-4bab-abab-e03c595b96da    75052db6272cde941f37a391e84c912fd8f459c9    1531    6545.75575    1.14    912    False    0.0    912    1.14    912    511    4.7    -0.0012    filtered_out    -1.0    -1.0    0.0

### Output 

1. Total yield as a function of sequencing time
2. Maximum read length as a function of cumulative sequencing time
3. Average read length as a function of cumulative sequencing time
4. Maximum read length per hour
5. Average read length per hour
