#!/usr/bin/env python2

"purcaro@gmail.com"

import urllib, os, shutil, tarfile, multiprocessing, subprocess, sys, socket
from color_text import ColorText as CT 


class Utils:
    @staticmethod
    def isMac():
        return sys.platform == "darwin"
    
    @staticmethod
    def connectedInternet():
        #from http://stackoverflow.com/questions/20913411/test-if-an-internet-connection-is-present-in-python
        try:
            # see if we can resolve the host name -- tells us if there is
            # a DNS listening
            host = socket.gethostbyname("www.google.com")
            # connect to the host -- tells us if the host is actually
            # reachable
            s = socket.create_connection((host, 80), 2)
            return True
        except:
            pass
        return False
    
    @staticmethod
    def which(program):
        #from http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
        def is_exe(fpath):
            return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    
        fpath, fname = os.path.split(program)
        if fpath:
            if is_exe(program):
                return program
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                path = path.strip('"')
                exe_file = os.path.join(path, program)
                if is_exe(exe_file):
                    return exe_file
        return None
    
    @staticmethod
    def hasProgram(program):
        whichOutput = Utils.which(program);
        return None != whichOutput;
    
    @staticmethod
    def run_in_dir(cmd, d):
        #print CT.boldBlack("here")
        cmd = "cd " + Utils.shellquote(d) + " && " + cmd + " && cd -"
        #print CT.boldBlack("newcmd")
        print CT.boldGreen(cmd)
        Utils.run(cmd)

    @staticmethod
    def run(cmd):
        # print CT.boldRed("before process")
        # from http://stackoverflow.com/a/4418193
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        #print CT.boldRed("after process")
        # Poll process for new output until finished
        actualOutput = ""
        while True:
            nextline = process.stdout.readline()
            if nextline == '' and process.poll() != None:
                break
            sys.stdout.write(nextline)
            actualOutput = actualOutput + nextline
            sys.stdout.flush()

        output = process.communicate()[0]
        exitCode = process.returncode
        #print "exit code "  + CT.boldRed(str(exitCode))
        if (exitCode == 0):
            return actualOutput
        raise Exception(cmd, exitCode, actualOutput)

    @staticmethod
    def runAndCapture(cmd):
        # from http://stackoverflow.com/a/4418193
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        # Poll process for new output until finished
        actualOutput = ""
        while True:
            nextline = process.stdout.readline()
            if nextline == '' and process.poll() != None:
                break
            actualOutput = actualOutput + nextline
        #this is suppose to capture the output but it isn't for some reason so capturing it with the above
        output = process.communicate()[0]
        exitCode = process.returncode
        if (exitCode == 0):
            return actualOutput
            #return output
        raise Exception(cmd, exitCode, actualOutput)

    @staticmethod
    def shellquote(s):
        #from http://stackoverflow.com/a/35857
        return "'" + s.replace("'", "'\\''") + "'"

    @staticmethod
    def num_cores():
        return multiprocessing.cpu_count()

    @staticmethod
    def mkdir(d):
        '''mkdir if it doesn't already exist '''
        if not os.path.exists(d):
            print CT.boldText("mkdir"), CT.boldGreen(d)
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
                print "skipping download of", CT.boldGreen(fn)
                return out_fnp
            else:
                print "files sizes differed:", "on disk:", fn_size, "from net:", net_file_size
        print "retrieving", CT.boldGreen(fn), "from", CT.boldBlue(url)
        urllib.urlretrieve(url, out_fnp)
        return out_fnp

    @staticmethod
    def rm_rf(d):
        '''remove directory forcibly'''
        if os.path.exists(d):
            print CT.boldText("rm -rf"), CT.boldRed(d)
            shutil.rmtree(d)

    @staticmethod
    def untar(fnp, d):
        ''' un pack compressed file, guessing format based on extention '''
        if fnp.endswith(".tar.gz"):
            tar = tarfile.open(fnp, "r:gz")
        elif fnp.endswith(".tgz"):
            tar = tarfile.open(fnp, "r:gz")
        elif fnp.endswith(".tar.bz2"):
            tar = tarfile.open(fnp, "r:bz2")
        elif fnp.endswith(".tar"):
            tar = tarfile.open(fnp, "r")
        else:
            raise Exception("invalid file? " + fnp)
        print "untarring", CT.boldGreen(fnp), "to", CT.boldBlue(d)
        tar.extractall(d)
        tar.close()
        
    @staticmethod
    def getStrFromStrOrList(inputArg):
        if type(inputArg) is list:
            return str(inputArg[0])
        elif type(inputArg) is not str:
            return str(inputArg)
        else:
            return inputArg
          
        

    @staticmethod
    def clear_dir(d):
        ''' forcibly delete directory and then re-make it''' 
        Utils.rm_rf(d)
        Utils.mkdir(d)
        
    @staticmethod
    def fixDyLibOnMac(libDir):
        """
            If a dynamic library's id isn't it's full path name and it isn't in the
            dylib search path it won't be linked in properly, so will modify the id
            of the libraries to be it's full name 
        """
        files = os.listdir(libDir)
        for file in files:
            fullFile = os.path.join(libDir, file)
            if os.path.isfile(fullFile) and str(fullFile).endswith(".dylib"):
                try:
                    cmd = "install_name_tool -id {full_libpath} {full_libpath}".format(full_libpath = os.path.abspath(fullFile))
                    Utils.run(cmd)
                except Exception,e:
                    print (e)
                    print ("Failed to fix dylib for {path}".format(path = os.path.abspath(fullFile)))
            elif os.path.isdir(fullFile):
                Utils.fixDyLibOnMac(fullFile)

    
            
    
    
