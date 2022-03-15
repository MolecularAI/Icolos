#!/bin/env python2
#
# jobinfo - collect job information from slurm in nicely readable format
#
# Copyright 2015 Anders Halager <aeh@birc.au.dk>
#
# LICENSE: MIT

import math
import os
import re
import subprocess
import sys
import time
from collections import namedtuple as NT

def parse_time(t):
    # Format: [DD-[hh:]]mm:ss
    time_parts = re.compile(r'(((?P<days>\d+)-)?(?P<hours>\d\d):)?' +
                            r'(?P<minutes>\d\d):(?P<seconds>\d\d(\.\d+)?)')
    m = time_parts.match(t)
    if m is None:
        return 0.0, 0, 0, 0
    ss = float(m.group('seconds'))
    mm = int(m.group('minutes'))
    hh = int(m.group('hours') or '0')
    dd = int(m.group('days') or '0')
    return ss, mm, hh, dd

def elapsed_to_seconds(elapsed):
    ss, mm, hh, dd = parse_time(elapsed)
    return dd*24*60*60 + hh*60*60 + mm*60 + ss

def format_bs(x):
    postfix = ' KMGTPE'
    e = int(math.log(x+1, 2)/10)
    return "%.2f%s" % (x / 2**(10*e), postfix[e])

def whoami():
    import pwd
    return pwd.getpwuid(os.getuid()).pw_name

# Constructors / parsers:
def str_set(x=None):
    if x in [None, '']:
        return set()
    return set([x])

def gpu_str(s=None):
    if s is None or "gres/gpu=" not in s.lower():
        return 0
    m = re.match(r".*gres/gpu=(?P<num_gpus>\d+)", s.lower())
    return int(m.group('num_gpus'))

def byte_size(s=None):
    if s in [None, "", "16?"]:
        return -1.0
    m = {'K': 10, 'M': 20, 'G': 30, 'T': 40, 'P': 50, 'E': 60}
    scale = 2 ** m.get(s[-1], 0)
    if scale != 1:
        s = s[:-1]
    return scale * float(s)

def date_str(s=None):
    if s is None or s.strip() == "":
        return "9999-01-01T00:00:00"
    return s

# Combinators:
def keep_first(a, b):
    return a == '' and b or a

def time_max(a, b):
    if 'UNLIMITED' in [a, b]:
        return 'UNLIMITED'
    if a in ['', 'INVALID']:
        return b
    if b in ['', 'INVALID']:
        return a
    return max(a, b)

# Formatters:
def f_rss(x, meta):
    if x < 0:
        return "--"
    return "%s (%s)" % (format_bs(x), ",".join(meta.MaxRSSNode))

def f_dw(x, meta):
    if x < 0:
        return "--"
    return "%s (%s)" % (format_bs(x), ",".join(meta.MaxDiskWriteNode))

def f_dr(x, meta):
    if x < 0:
        return "--"
    return "%s (%s)" % (format_bs(x), ",".join(meta.MaxDiskReadNode))

def f_cpu(x, meta):
    total = elapsed_to_seconds(meta.TotalCPU)
    if total == 0:
        return "--"
    xp = elapsed_to_seconds(x)
    return "%5.2f%%" % (xp/total*100)

def f_mem(x, meta):
    if x.endswith('c'):
        return "%s/core" % (x[:-1])
    elif x.endswith('n'):
        return "%s/node" % (x[:-1])
    else:
        return x

def f_time(x, meta):
    all_times = [meta.timelimit, meta.elapsed, meta.TotalCPU, '-']
    days_len = max(len(y.split('-')[0]) for y in all_times if '-' in y)
    ss, mm, hh, dd = parse_time(x)
    if days_len == 0:
        dd = ""
    else:
        if dd > 0:
            dd = ("%i-" % dd).rjust(days_len)
        else:
            dd = " "*(days_len + 1)
    res = "%s%02i:%02i:%02i" % (dd, hh, mm, ss)
    if res.strip() == "00:00:00":
        return "--"
    return res

def f_str(x, meta):
    return str(x)

def f_date(x, meta):
    if str(x).lower() == "unknown":
        return "--"
    return str(x)

def f_state(states, meta):
    if len(states) > 1:
        states = states - set(["COMPLETED", ""])
    reason = meta.reason
    if reason != '':
        reason = ' ' + reason
    deps = meta.dependencies
    if deps != '':
        deps = " (%s)" % deps
    return ','.join(states) + reason + deps

hide, show = False, True
Field = NT('Field', 'name ctor combinator shown prefer_sstat formatter desc')
FIELDS = [
        Field("JobName",             str,        keep_first,   show,  False, f_str,     "Name"),
        Field("User",                str,        keep_first,   show,  False, f_str,     "User"),
        Field("Partition",           str,        keep_first,   show,  False, f_str,     "Partition"),
        Field("NodeList",            str,        keep_first,   show,  False, f_str,     "Nodes"),
        Field("ncpus",               int,        max,          show,  False, f_str,     "Cores"),
        Field("ReqTRES",             gpu_str,    max,          show,  False, f_str,     "GPUs"),
        Field("State",               str_set,    set.union,    show,  False, f_state,   "State"),
        Field("Submit",              str,        keep_first,   show,  False, f_str,     "Submit"),
        Field("start",               date_str,   min,          show,  False, f_date,    "Start"),
        Field("end",                 str,        time_max,     show,  False, f_date,    "End"),
        Field("timelimit",           str,        time_max,     show,  False, f_time,    "Reserved walltime"),
        Field("elapsed",             str,        time_max,     show,  False, f_time,    "Used walltime"),
        Field("TotalCPU",            str,        max,          show,  False, f_time,    "Used CPU time"),
        Field("UserCPU",             str,        max,          show,  False, f_cpu,     "% User (Computation)"),
        Field("SystemCPU",           str,        max,          show,  False, f_cpu,     "% System (I/O)"),
        Field("ReqMem",              str,        keep_first,   show,  False, f_mem,     "Mem reserved"),
        Field("MaxRSS",              byte_size,  max,          show,  True,  f_rss,     "Max Mem used"),
        Field("MaxDiskWrite",        byte_size,  max,          show,  True,  f_dw,      "Max Disk Write"),
        Field("MaxDiskRead",         byte_size,  max,          show,  True,  f_dr,      "Max Disk Read"),

        Field("MaxRSSNode",          str_set,    set.union,    hide,  True,  None,      ""),
        Field("MaxDiskWriteNode",    str_set,    set.union,    hide,  True,  None,      ""),
        Field("MaxDiskReadNode",     str_set,    set.union,    hide,  True,  None,      ""),
        ]

FIELD_NAMES = [f.name for f in FIELDS]
FIELD_NAMES_SSTAT = [f.name for f in FIELDS if f.prefer_sstat]
FIELD_CTORS = [f.ctor for f in FIELDS]
FIELD_COMB = [f.combinator for f in FIELDS]
FORMAT_STR = "--format=%s" % (",".join(FIELD_NAMES))
FORMAT_SSTAT_STR = "--format=%s" % (",".join(FIELD_NAMES_SSTAT))
Meta = NT('Meta', FIELD_NAMES + ['dependencies', 'reason'])

def combine(xs):
    r = xs[0]
    for x in xs[1:]:
        for i, comb in enumerate(FIELD_COMB):
            r[i] = comb(r[i], x[i])
    return r

def get_values(jobid):
    info = subprocess.Popen(['sacct', FORMAT_STR, '--parsable', '--noheader', '-j', jobid], stdout=subprocess.PIPE)
    xs = []
    for line in info.stdout:
        fields = line.strip().split('|')
        xs.append([ctor(s) for ctor, s in zip(FIELD_CTORS, fields)])
    if len(xs) == 0:
        print >>sys.stderr, "No such job"
        sys.exit(1)
    return xs

def get_sstat_values(jobid):
    info = subprocess.Popen(['sstat', FORMAT_SSTAT_STR, '--parsable', '--noheader', '-a', '-j', jobid], stdout=subprocess.PIPE)
    xs = []
    for line in info.stdout:
        j = 0
        vals = line.strip().split('|')
        x = []
        for f in FIELDS:
            if f.prefer_sstat:
                x.append(f.ctor(vals[j]))
                j += 1
            else:
                x.append(f.ctor())
            xs.append(x)
    return xs

def main(jobid):
    y = combine(get_values(jobid))
    meta = Meta._make(y + ['', ''])
    ys = [y]
    if "RUNNING" in meta.State and (os.getuid() == 0 or meta.User == whoami()):
        # get more info from sstat
        tmp = get_sstat_values("%s,%s.batch" % (jobid, jobid))
        if len(tmp) != 0:
            ys.append(combine(tmp))
    if "PENDING" in meta.State:
        info = subprocess.Popen(['squeue', '--format=%E;%R', '--noheader', '-a', '-j', jobid], stdout=subprocess.PIPE)
        deps, reason = info.stdout.readline().strip().split(";")
        dependencies = deps
    else:
        dependencies = ""
        reason = ""
    y = combine(ys)
    meta = Meta._make(y + [dependencies, reason])

    for i,(name,parse,comb,show,prefer_sstat,format,desc) in enumerate(FIELDS):
        val = y[i]
        if show:
            print "%-20s: %s" % (desc, format(val, meta))

def usage(pipe):
    print >>pipe, \
"""jobinfo - collates job information from the 'sstat', 'sacct' and
'squeue' SLURM commands to give a uniform interface for both current
and historical jobs.

Usage:
    jobinfo <job id>

Report bugs to Anders Halager <aeh@birc.au.dk>
  or on GitHub (https://github.com/birc-aeh/slurm-utils)"""

if __name__ == "__main__":
    if "-h" in sys.argv or "--help" in sys.argv:
        usage(sys.stdout)
        sys.exit(0)
    if len(sys.argv) != 2:
        usage(sys.stderr)
        sys.exit(1)
    jobid = sys.argv[1]
    if len(set(jobid) - set("0123456789_.")) > 0:
        print >>sys.stderr, "The argument does not look like a valid job id"
        usage(sys.stderr)
        sys.exit(1)
    main(jobid)