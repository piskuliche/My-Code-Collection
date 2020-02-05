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

lines=os.popen('squeue -u "%s" -o  "%%.18i %%.9P %%.12j %%.8u %%.2t %%.10M %%.6D %%R"' % user).read().split('\n')
lines.pop()
print(len(lines))

job_printed=[]
names_printed=[]
running=[]
pending=[]
start=[]
rcount=0
pdcount=0
for line in lines:
    if "JOBID" in line:
        start=line
    elif "R" in line.split()[4]:
        jid=line.split()[0]
        rcount+=1
        if "_" in jid:
            job=jid.split('_')[0]
            if job not in job_printed:
                job_printed.append(job)
                if uniquename is True:
                    name=line.split()[2]
                    if name not in names_printed:
                        running.append(line)
                        names_printed.append(name)
                else:
                    running.append(line)
        else:
            running.append(line)
            job_printed.append(jid)
    elif "PD" in line.split()[4]:
        jid=line.split()[0]
        pdcount+=1
        if "_" in jid:
            job=jid.split('_')[0]
            if job not in job_printed:
                job_printed.append(job)
                if uniquename is True:
                    name=line.split()[2]
                    if name not in names_printed:
                        pending.append(line)
                        names_printed.append(name)
                else:
                    pending.append(line)
        else:
            pending.append(line)
            job_printed.append(jid)

print("----------------")
print("%d RUNNING      " % rcount)
print("----------------")
print(start)
for job in running:
    print(job)
print("----------------")
print("%d PENDING      " % pdcount)
print("----------------")
print(start)
for job in pending:
    print(job)
    

