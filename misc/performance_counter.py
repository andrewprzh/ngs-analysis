#!/usr/bin/env python3
#

import os
import sys
import psutil
import subprocess
from traceback import print_exc
from time import time
from time import sleep
from collections import defaultdict


def get_current_stats(pids):
    total_vms = 0
    total_rss = 0
    total_cpu = 0.0
    for pid in pids.keys():
        proc = pids[pid]
        try:
            mem_info = proc.memory_info()
            total_vms += mem_info.vms
            total_rss += mem_info.rss
            total_cpu += proc.cpu_percent(interval=None)
        except psutil.NoSuchProcess:
            pass

    return total_rss, total_vms, total_cpu


def get_children(pid):
    pids = {pid}
    for p in psutil.Process(pid).children(recursive=True):
        pids.add(p.pid)
    return pids


def main():
    if sys.argv[0].startswith("python"):
        cmd = sys.argv[2:]
    else:
        cmd = sys.argv[1:]

    outf = open("performance_stats.tsv", "w")
    outf.write("Time\tRSS\tVMS\tCPU %\n")
    launched = subprocess.Popen(cmd, shell=False)
    main_pid = launched.pid
    main_proc = psutil.Process(main_pid)
    all_pids = {main_pid : main_proc}
    start_time = time()
    max_rss = 0
    cpu_times = defaultdict(float)
    while main_proc.status() != 'zombie':
        children_pids = get_children(main_pid)
        for child in children_pids:
            if child not in all_pids:
                all_pids[child] = psutil.Process(child)
        total_rss, total_vms, total_cpu = get_current_stats(all_pids)
        max_rss = max(max_rss, total_rss)
        current_time = time() - start_time
        outf.write("%d\t%d\t%d\t%.1f\n" % (current_time, total_rss, total_vms, total_cpu))
        for p in all_pids.keys():
            try:
                cpu_times[p] = max(cpu_times[p], all_pids[p].cpu_times().user)
            except psutil.NoSuchProcess:
                pass
        sleep(1)

    outf.close()
    cpu_time = sum(cpu_times.values())
    print("Max RSS: %d, CPU time: %d, wall clock time: %d" % (max_rss, cpu_time, time() - start_time))


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)


