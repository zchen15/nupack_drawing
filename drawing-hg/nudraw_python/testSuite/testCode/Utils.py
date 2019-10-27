# Some common funtions for tests
import os
import string

#-------------------------------------------------------------------------------
def join_path(plist):
  if len(plist)>0:
    return reduce(lambda x,y:x+os.sep+y, plist)
  else: return plist

#-------------------------------------------------------------------------------
def execTest(exe_name, args):
  pid=os.fork()
  if pid==0: #Child
    os.execv(exe_name,args)
    os.exit(0)
  pid,status=os.waitpid(pid,0)
  #~ print "status=%d, exited=%s, exitstatus=%d"%(status, os.WIFEXITED(status), os.WEXITSTATUS(status))
  if os.WIFEXITED(status) and os.WEXITSTATUS(status)==0:
    return True
  else:
    return False

#-------------------------------------------------------------------------------
def clearData(clearPath, debug=False):
  if debug:
    print "Clearing",clearPath
  for root, dirs, files in os.walk(clearPath, topdown=False):
    for name in files:
        os.remove(os.path.join(root, name))
        if debug:
          print "file=%s"%name
    for name in dirs:
      if debug:
        print "dir=%s"%name
      os.rmdir(os.path.join(root, name))

letters="ACGT"
domains=["A","C","E", "B", "D", "F"]

#-------------------------------------------------------------------------------
def generateSequence(struc):
  global domains
  ll=len(letters)
  lc=len(domains)
  seq=""
  cols=[]
  j=0
  for i in xrange(len(struc)):
    if struc[i]=='+':
      seq+="+"
      cols+="+"
    else:
      seq+=letters[j%ll]
      cols.append("%d"%(j%lc))
      j+=1

  return seq, string.join(cols,","), string.join(domains,",")
