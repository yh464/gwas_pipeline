#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-06-17

A quick tool to check the SLURM queue
'''

import os
from subprocess import check_output
from collections import Counter

def parse_entry(line):
  line = line.split()
  
  # number of jobs
  jobid = line[0]
  if jobid.find('[') == -1: njobs = 1
  else:
    array = jobid[jobid.find('[')+1:].replace(']','')
    array = array.split(',')
    njobs = 0
    for x in array:
      if x.find('-') == -1: njobs += 1
      else: 
        last_job = int(x.split('-')[1].split('%')[0])
        first_job = int(x.split('-')[0].split('%')[0])
        if last_job < first_job: last_job *= 10 # if last digit is not displayed
        njobs += last_job - first_job + 1
  
  # crsid
  userid = line[3]
  # number of nodes
  nnodes = int(line[6])
  # running or not
  # running = (line[5] != '0:00')
  return userid, njobs, nnodes*njobs

def partition_info(partition):
  # get queue info
  squeue = check_output(['squeue','-p', partition]).decode().splitlines()
  jobs = Counter(); nodes = Counter(); arrays = Counter()
  for line in squeue[1:]:
    uid, j, n = parse_entry(line)
    jobs.update({uid: j}); nodes.update({uid: n}); arrays.update({uid: 1})

  sinfo = check_output(['sinfo','-p', partition]).decode().splitlines()
  up = 0; down = 0; idle = 0
  for line in sinfo[1:]:
    line = line.split()
    if line[4] in ['mix','alloc','idle']: up += int(line[3])
    else: down += int(line[3])
    if line[4] == 'idle': idle += int(line[3])
  
  print('#' * 100)
  print(f'Partition: {partition:15}Idle: {idle!s:6}Available: {up!s:6}Down: {down}')
  print(f'Total jobs: {jobs.total()} in {arrays.total()} arrays    Total nodes: {nodes.total()}')
  print('Top users:')
  print('User          Jobs /arrays   Nodes')
  for uid, n in nodes.most_common(5):
    print(f'{uid:15}{jobs[uid]!s:7}/{arrays[uid]!s:7}{n}')
  print('#'*100)
  return up, down, idle, jobs.total(), nodes.total(), arrays.total()

info_table = ['Partition      Idle  Avail Down  Jobs / arrays   Nodes']
for partition in ['cclake','cclake-himem','sapphire','icelake','icelake-himem','desktop']:
  up, down, idle, jobs, nodes, arrays = partition_info(partition)
  info_table.append(f'{partition:15}{idle!s:6}{up!s:6}{down!s:6}{jobs!s:7}/{arrays!s:7}{nodes}')
for x in info_table: print(x)