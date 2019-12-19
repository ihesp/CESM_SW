#!/usr/bin/env python3
import argparse
import sys
import os
#import exceptions
parser = argparse.ArgumentParser(description='Easier case creator')
parser.add_argument('-l', '--phys-loadbalance', help='Set phys_loadbalance option for new case.', action='store', type=int, default=None)
parser.add_argument('-r', '--resolution', help='Resolution of new case.', choices=['ne120', 'ne30', 't12', 'g16', 's11'], action='store', default='ne30')
parser.add_argument('-o', '--output', help='Enable data output for new case.', action='store_true', default=False)
parser.add_argument('-P', '--pes-file', help='env_mach_pes file for this case.', action='store', default=None)
parser.add_argument('-p', '--pecount', help='Number of pe to be used.', action='store', default=None)
parser.add_argument('-s', '--src-branch', help='Branch for SourceMods.', action='store', default=None, required=True)
parser.add_argument('-c', '--compset', help='Compset to be used for case', choices=['B', 'D', 'G', 'F'], action='store', default='B')
parser.add_argument('-C', '--compiler', help='Compiler to be used for case', choices=['swcc', 'sw5c'], action='store', default='swcc')
parser.add_argument('case_name')
#parser.print_help()

args = parser.parse_args(sys.argv[1:])
#print(args)
cmdparts = ['./create_newcase -mach sunway -user_mods_dir user_dir -case']
caseroot = '../cases/' + args.case_name
if args.src_branch is not None:
    branch_def = '-b ' + args.src_branch
cmdparts.append(caseroot)

if args.resolution == 's11' and args.compset in ['B', 'F']:
    print('s11 grid only supports G/D compset')
    exit(1)
if args.compset == 'B':
    cmdparts.append('-compset B1850C5')
    if args.resolution == 'ne30':
        cmdparts.append('-res ne30_g16')
    elif args.resolution == 'ne120':
        cmdparts.append('-res ne120_t12')
if args.compset == 'D':
    cmdparts.append('-compset D')
    if args.resolution in ['t12', 'ne120']:
        cmdparts.append('-res T62_t12')
    elif args.resolution in ['g16', 'ne30']:
        cmdparts.append('-res T62_g16')
    elif args.resolution == 's11':
        cmdparts.append('-res T62_s11')
if args.compset == 'G':
    cmdparts.append('-compset G')
    if args.resolution in ['t12', 'ne120']:
        cmdparts.append('-res T62_t12')
    elif args.resolution in ['g16', 'ne30']:
        cmdparts.append('-res T62_g16')
    elif args.resolution == 's11':
        cmdparts.append('-res T62_s11')
if args.compset == 'F':
    cmdparts.append('-compset F1850C5')
    if args.resolution == 'ne30':
        cmdparts.append('-res ne30_ne30')
    elif args.resolution == 'ne120':
        cmdparts.append('-res ne120_ne120')
if args.pes_file is not None:
    cmdparts.append('-pes_file ' + args.pes_file)
if args.pecount is not None:
    cmdparts.append('-pecount ' + args.pecount)
cmdparts.append('-compiler ' + args.compiler)
cmd = ' '.join(cmdparts)
print('Resolution is: ' + args.resolution)
print('Creating case using: ' + cmd)

os.system(cmd)
if args.phys_loadbalance is not None:
    os.system('echo phys_loadbalance=%d >> %s/user_nl_cam' % (args.phys_loadbalance, caseroot))

if not args.output:
    os.system('cd %s && ./xmlchange REST_OPTION=never' % caseroot)
os.system('cd %s && ./xmlchange PIO_TYPENAME=pnetcdf' % caseroot)
os.system('cd %s && rm -rfv SourceMods' % caseroot)
branch_def = ''
if args.src_branch is not None:
    branch_def = '-b ' + args.src_branch
os.system(('cd %s && git clone ' + branch_def + ' file://$HOME/online1/git/SourceMods.git') % caseroot)

os.system('echo profile_depth_limit=99999 >> %s/user_nl_cpl' % caseroot)
os.system('echo profile_detail_limit=99999 >> %s/user_nl_cpl' % caseroot)
os.system('echo profile_outpe_num=512 >> %s/user_nl_cpl' % caseroot)
os.system('echo profile_outpe_stride=4 >> %s/user_nl_cpl' % caseroot)
#os.system("echo \"tavg_freq_opt = 'nmonth' 'never'\" >> %s/user_nl_pop2" % caseroot)
# clinic_distribution_type = 'spacecurve'
# tropic_distribution_type = 'spacecurve'
