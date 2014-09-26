# read and create a list of of items 
def rdc(myfile):
        
        listfile = open(myfile,'r')
        i1 = -1

        datalis =[]

        for data in listfile:
                i1 = i1+1
                datalis = datalis + [data]
#               print i1 , cardlis[i1]
        listfile.close()
	return datalis


# call sequence :
#   import readcard
#   cardfile = 'file name'
#   cardlis = readcard.rdc(cardfile)

