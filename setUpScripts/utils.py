#!/usr/bin/env python

"purcaro@gmail.com"

import urllib, os, shutil, tarfile, multiprocessing, subprocess, sys
from color_text import ColorText as CT 


class Utils:
    @staticmethod
    def run_in_dir(cmd, d):
        print CT.boldBlack("here")
        cmd = "cd " + Utils.shellquote(d) + " && " + cmd + " && cd -"
        print CT.boldBlack("newcmd")
        print CT.boldGreen(cmd)
        Utils.run(cmd)

    @staticmethod
    def run(cmd):
        # print CT.boldRed("before process")
        # from http://stackoverflow.com/a/4418193
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        #print CT.boldRed("after process")
        # Poll process for new output until finished
        while True:
            nextline = process.stdout.readline()
            if nextline == '' and process.poll() != None:
                break
            sys.stdout.write(nextline)
            sys.stdout.flush()

        output = process.communicate()[0]
        exitCode = process.returncode
        #print "exit code "  + CT.boldRed(str(exitCode))
        if (exitCode == 0):
            return output
        raise Exception(cmd, exitCode, output)

    @staticmethod
    def shellquote(s):
        " from http://stackoverflow.com/a/35857"
        return "'" + s.replace("'", "'\\''") + "'"

    @staticmethod
    def num_cores():
        return multiprocessing.cpu_count()

    @staticmethod
    def mkdir(d):
        '''mkdir if it doesn't already exist '''
        if not os.path.exists(d):
            print "mkdir", d
            os.makedirs(d)

    @staticmethod
    def get_file(url, d):
        '''get file from url and put it into directory d, return new name  '''
        fn = url.split('/')[-1]
        out_fnp = os.path.join(d, fn)
        urllib.urlretrieve(url, out_fnp)
        return out_fnp

    @staticmethod
    def get_file_if_size_diff(url, d):
        '''only download the file if it's needed, not completely fail proof since it is 
        just a size check but fairly likely not to be the same for a difference '''
        fn = url.split('/')[-1]
        out_fnp = os.path.join(d, fn)
        net_file_size = int(urllib.urlopen(url).info()['Content-Length'])
        if os.path.exists(out_fnp):
            fn_size = os.path.getsize(out_fnp)
            if fn_size == net_file_size:
                print "skipping download of", fn
                return out_fnp
            else:
                print "files sizes differed:", "on disk:", fn_size, "from net:", net_file_size
        print "retrieving", fn
        urllib.urlretrieve(url, out_fnp)
        return out_fnp

    @staticmethod
    def rm_rf(d):
        '''remove directory forcibly'''
        if os.path.exists(d):
            print "rm -rf", d
            shutil.rmtree(d)

    @staticmethod
    def untar(fnp, d):
        ''' un pack compressed file, guessing format based on extention '''
        if fnp.endswith(".tar.gz"):
            tar = tarfile.open(fnp, "r:gz")
        elif fnp.endswith(".tar.bz2"):
            tar = tarfile.open(fnp, "r:bz2")
        elif fnp.endswith(".tar"):
            tar = tarfile.open(fnp, "r")
        else:
            raise Exception("invalid file? " + fnp)
        print "untarring", fnp, "to", d
        tar.extractall(d)
        tar.close()

    @staticmethod
    def clear_dir(d):
        ''' forcibly delete directory and then re-make it''' 
        Utils.rm_rf(d)
        Utils.mkdir(d)
