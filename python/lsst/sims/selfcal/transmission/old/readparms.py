# read a catalogue of MODTRAN parameters with their
#  name, format, card sequential number, range in card
# return a list with parameters attributes
#
####


def rdp(parmfile):
        
        parmf = open(parmfile,'r')
        i1 = -1

        parmlis =[]
        # read complete catalogue of parameters that can be changed
        for data in parmf:
                i1 = i1+1
                parmlis = parmlis + [data]#one string per parameter
                print i1 , parmlis[i1]
        parmf.close()
        plist =[]
        for iparm in range(0,len(parmlis)):
                psp=[]
                par=parmlis[iparm]     # isolate one parameter
                psp = par.split('$')   #split parameter attributes
                # format parameter attributes
                psp[0]=psp[0].strip()   #strip name string
                psp[1]=int(psp[1])      # switch card number to numeric
                psp[2]=psp[2].strip()   # strip format string
                psp[3]=int(psp[3])      # switch start position to numeric
                psp[4]=int(psp[4])      # switch end position to numeric
                psp = psp[0:-1]         # remove useless endline
                # store in global parameter list 
                plist = plist + [psp]   
        return plist


# call sequence :
#   import readparms
#   parmfile = 'file name'
#   cardlis = readparms.rdp(parmfile)

