#!/usr/bin/env python

import os, sys, time
import json, subprocess

def getFingerOnHomeFolders():
    fnp = '/home/hathawan/users.json'

    def runProcess(exe):    
        p = subprocess.Popen(exe, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        ret = ""
        while(True):
          retcode = p.poll() # returns None while subprocess is running
          line = p.stdout.readline()
          ret += line
          if(retcode is not None):
            break
        return ret

    if os.path.exists(fnp):
        with open(fnp) as f:
            return json.load(f)

    fingers = {}

    for p in os.listdir("/home/"):
        fingers[p] = runProcess("/usr/bin/finger " + p)

    with open(fnp, 'wb') as f:
        json.dump(fingers, f)

    return fingers

def getName(u):
    first = ""
    last = ""
    for line in u.split("\n"):
        toks = line.split("\t")
        for t in toks:
            if t.startswith("Name:"):
                last = t[5:].strip()
            if t.startswith("Office:"):
                first = t[7:].strip()
    return (first, last)

def convertTime(lastLogInTime, addDay):
    '''
    time is the string returned from getLastLoginDatetime
    addDay should be a bool indication whether to append the day
    '''
    # month abbreviations "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", and "Dec"
    monthConverter = {"Jan": "01","Feb": "02", "Mar": "03", "Apr": "04", "May": "05", "Jun": "06", "Jul": "07", "Aug": "08", "Sep": "09", "Oct": "10", "Nov": "11", "Dec": "12"}
    if lastLogInTime == "":
        return "Never_logged_in"
    
    t = lastLogInTime.strip()    
    tSplit = t.split()
    year = ""
    month = monthConverter[str(tSplit[1])]
    day = str(tSplit[2])
    logInTime = str(tSplit[3])
    dayOfTheWeek = str(tSplit[0])
    currentMonth = int((time.strftime("%m")))
    currentYear = int((time.strftime("%Y")))
    if(len(tSplit) == 4):
        monthInt = int(month)
        if(currentMonth <=6 and monthInt > 6):    
            year = str(currentYear - 1)
        else:
            year = str(currentYear)
    elif (len(tSplit) == 5):
        year = str(tSplit[4])
    if(len(day) ==1):
        day = "0" + day
    if addDay:
        return "_".join([year, month, day, logInTime, dayOfTheWeek])
    else:
        return "_".join([year, month, day, logInTime])
    
    
def getLastLoginDatetime(u):
    for line in u.split("\n"):
        toks = line.split("\t")
        for t in toks:
            if t.startswith("Last login"):
                t = t[10:].strip()
                if '(EST)' in t:
                    return t.split('(EST)')[0].strip()
                if '(EDT)' in t:
                    return t.split('(EDT)')[0].strip()
                return t
            elif t.startswith("On since"):
                t = t[8:].strip()
                if '(EST)' in t:
                    return t.split('(EST)')[0].strip()
                if '(EDT)' in t:
                    return t.split('(EDT)')[0].strip()
                return t
    return ""

def get_size(start_path):
    # http://stackoverflow.com/a/1392549
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(start_path):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            try:
                total_size += os.path.getsize(fp)
            except:
                pass
    return total_size

def sizeof_fmt(num, suffix='B'):
    # http://stackoverflow.com/a/1094933
    for unit in ['','Ki','Mi','Gi','Ti','Pi','Ei','Zi']:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, 'Yi', suffix)

def getOldUserIDs():
    ret = []
    with open("/home/hathawan/usersOld.txt") as f:
        for line in f:
            line = line.strip()
            if "" == line:
                continue
            toks = line.split(",")
            ret.append(toks[2].strip())
    return ret

fingers = getFingerOnHomeFolders()
#oldUsers = getOldUserIDs()
print ",".join(["LastName", "FirstName", "userName", "DateOfLastLogin"])
for homeFolder, finger in fingers.iteritems():
    #if homeFolder not in oldUsers:
    #    continue
    if not os.path.exists(os.path.join("/home/", homeFolder)):
        continue
    first, last = getName(finger)
    lastDateTime = getLastLoginDatetime(finger)
    print ','.join([last, first, homeFolder, convertTime(lastDateTime, False)])
    #bytes = get_size(("/home/" + homeFolder).encode("utf-8"))
    #if bytes > 10000:
    #    continue
    #print ', '.join([last, first, homeFolder, lastDateTime, sizeof_fmt(bytes)])
