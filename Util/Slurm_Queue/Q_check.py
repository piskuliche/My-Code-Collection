#!/usr/bin/env python
import os

lines=os.popen('squeue -u "e924p726"').read()
job_printed=[]
for line in lines.split('\n'):
    if "JOBID" in line:
        print(line)
    elif "R" in line:
        jid=line.split()[0]
        if "_" in jid:
            job=jid.split('_')[0]
            if job not in job_printed:
                print(line)
                job_printed.append(job)


