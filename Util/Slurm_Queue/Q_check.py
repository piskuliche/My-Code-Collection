#!/usr/bin/env python
import os
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-u', default='e924p726', type=str, help='Username')
parser.add_argument('-n', default=0, type=int, help='Print only unique names')
args = parser.parse_args()

user=args.u
uniquename = False
if args.n == 1:
    uniquename=True


print(uniquename,user)

lines=os.popen('squeue -u "%s"' % user).read()
job_printed=[]
names_printed=[]
for line in lines.split('\n'):
    if "JOBID" in line:
        print(line)
    elif "R" in line:
        jid=line.split()[0]
        if "_" in jid:
            job=jid.split('_')[0]
            if job not in job_printed:
                job_printed.append(job)
                if uniquename is True:
                    name=line.split()[2]
                    if name not in names_printed:
                        print(line)
                        names_printed.append(name)
                else:
                    print(line)


