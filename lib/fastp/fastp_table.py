import os,sys,json
fastp_json = sys.argv[1]
outdir = sys.argv[2]
samplename = sys.argv[3]
summary = json.load(open(fastp_json))['summary']
header = 'Sample_name,Raw_reads,Raw_bases(G),Clean_reads,Clean_bases(G),Clean_ratio(%),Q20(%),Q30(%),GC_content(%)'
data = ['{}'.format(samplename),
        '{}'.format(int(summary['before_filtering']['total_reads']/2)),
        '{:.2f}'.format(summary['before_filtering']['total_bases']/10**9),
        '{}'.format(int(summary['after_filtering']['total_reads']/2)),
        '{:.2f}'.format(summary['after_filtering']['total_bases']/10**9),
        '{:.2f}'.format(summary['after_filtering']['total_bases']/summary['before_filtering']['total_bases']*100),
        '{:.2f}'.format(summary['after_filtering']['q20_rate']*100),
        '{:.2f}'.format(summary['after_filtering']['q30_rate']*100),
        '{:.2f}'.format(summary['after_filtering']['gc_content']*100)
       ]
table_csv = os.path.join(outdir, '{}_table.csv'.format(samplename))
with open(table_csv, 'w') as fh:
    fh.write(header + '\n')
    fh.write(','.join(data) + '\n')
