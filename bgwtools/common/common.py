
numpy_flavors=['','double','complex']
flavors=['NONE','REAL','COMPLEX']
class flavor:
	NONE=0
	REAL=1
	CPLX=2

def get_flavor(str_, die=True):
	try:
		return flavors.index(str_.upper().strip())
	except:
		if not die:
			return flavor.NONE
		else:
			raise ValueError('Invalid flavor "%s"'%(str_))

def get_numpy_flavor(iflavor, die=True):
	try:
		return numpy_flavors[iflavor]
	except:
		if not die:
			return numpy_flavors[flavor.NONE]
		else:
			raise ValueError('Invalid flavor "%s"'%(str_))

ftypes=['UNK','WFN','RHO']
def get_ftype(str_, die=True):
	try:
		return ftypes.index(str_.upper())
	except:
		if not die:
			return 0
		else:
			raise ValueError('Invalid ftype "%s"'%(str_))

