#!/usr/bin/python
import os, time, sys

def get_jobid() :
    bjobs = os.popen("bjobs").readlines()
    for i in bjobs:
        s = i.split()
        if len(s) == 14 :
            if int(s[12]) > 10000 : 
                return int(s[0])
    return -1

def check_down(jobid, errnodes) :
    cmd = '/home/export/online1/swmore/release/bin/filternode.py -f ' + '"jobid==' + str(jobid) + " and status!='busy'" + '"'
    # print cmd
    nodes = os.popen(cmd).readlines()
    if len(nodes) > 0 :
        if len(nodes) > 1 :
            exit('Return of filternode should not be great than 1 line!')
        errnodes.extend(nodes[0][:-1].split(','))
        # print errnodes
        return True
    return False

def get_outfile(jobid) :
    jobinfo = os.popen('bjobs -l ' + str(jobid)).readlines()
    return jobinfo[5].replace(' ', '')[9:-2]

def check_abort(outfile, errnodes) :
    cmd = 'grep -o "mpi_rank_[0-9]*" ' + outfile
    # print cmd
    lines = os.popen(cmd).readlines()
    if len(lines) > 0 :
        mpi_ranks = [x[9:-1] for x in lines]
        # print mpi_ranks
        for i in mpi_ranks :
            cmd = 'grep "vn[0-9][0-9][0-9][0-9][0-9][0-9]" '  + outfile + ' | grep "[^1-9]' + i + '[^0-9]" '
            # print cmd
            lines = os.popen(cmd).readlines()
            if len(lines) > 0 :
                hostname = os.popen('echo "' + lines[-1] + \
                        '" | grep -o "vn[0-9][0-9][0-9][0-9][0-9][0-9]"').readlines()[0]
                errnodes.append(int(hostname[2:8]))
                # print errnodes
        return True
    return False

def killjob(jobid) :
    os.popen("bkill " + str(jobid))

def deletenode(errnodes, queue) :
    if len(errnodes) == 0 :
        return
    f = open("/home/export/online1/cesm06/jflfy/errnodes.txt", 'a')
    print 'new error nodes list is', ','.join([str(i) for i in errnodes])
    for i in errnodes :
        f.write(str(i) + '\n')
    f.close()
    cmd = 'qnode -q ' + queue + ' -a -c ' + ','.join([str(i) for i in errnodes])
    # print cmd
    while True:
        ret = os.popen(cmd).readlines()
        if len(ret) != 1 :
            exit("Return of qnode should be only 1 line!")
        ret = ret[0].split()
        if ret[-1] == 'successful' :
            break

def subjob() :
    case_name = os.path.split(os.path.realpath(__file__))[0].split('/')[-1]
    os.system('qname=q_sw_cesm ./' + case_name + '.sub')

def rotate_bar() :
    s = ["\\","|","/","-","|","-"]
    for j in range(10) :
        for i in s : 
            sys.stdout.write("\r")
            sys.stdout.write(i)
            sys.stdout.flush()
            time.sleep(0.5)


jobid = get_jobid()
last_jobid = jobid
if jobid != -1 :
    outfile = get_outfile(jobid)
    print 'jobid is', jobid
    print 'outfile is ', outfile

while True : 
    rotate_bar()
    downnodes = []
    abortnodes = []
    jobid = get_jobid()
    if jobid == -1 :
        subjob()
        continue
    if (last_jobid != jobid) :
        last_jobid = jobid
        outfile = get_outfile(jobid)
        sys.stdout.write("\r")
        sys.stdout.flush()
        print 'jobid is', jobid
        print 'outfile is ', outfile

    if check_down(jobid, downnodes) is True or check_abort(outfile, abortnodes) is True :
        killjob(jobid)
        deletenode(downnodes, 'q_sw_cesm_maintaince')
        deletenode(abortnodes, 'q_sw_cesm_err')
        subjob()

