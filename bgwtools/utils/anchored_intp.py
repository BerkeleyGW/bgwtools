
from numpy import *
from scipy.interpolate import interp1d

class anchored_intp:
	def __init__(self, num_anchors):
		self.num_anchors = num_anchors
		self.anchors = []
		self.ts = []
		self.xs = []
		self.ys = []
		self._ready = False

	def add_data(self, t, x, y, anchors):
		'''
		Add a new spectrum to the collection of spectra.
			t: the (unique) contiuous parameter the identifies this data set
			X: array containing the x data
			Y: array containing the y data
			anchors: array of anchors
		'''
		if self._ready:
			raise Exception('Can\'t add more data after interpolation was performed.')
		anchors = array(anchors, dtype=int)
		assert(len(anchors)==self.num_anchors)
		self.xs.append(array(x))
		self.ys.append(array(y))
		self.ts.append(t)
		self.anchors.append(anchors)
		self._ready = False

	def prepare(self):
		'''Order the data so that __call__ can be executed more efficiently'''
		if len(self.ts)<2:
			raise Exception('You must provide at least two spectra before interpolating.')
		order = argsort(self.ts)
		self.ts = array(self.ts, dtype=float)[order]
		self.xs = array(self.xs)[order]
		self.ys = array(self.ys)[order]
		self.anchors = array(self.anchors)[order]
		self._ready = True

	def __call__(self, t, x=None):
		'''Get interpolated spectrum for the particular value (t) of the continuous parameter.'''
		if not self._ready: self.prepare()
		
		rng = arange(len(self.ts))
		#get values from which we interpolate
		if t<self.ts[0]:
			idx1 = 0
			idx2 = 1
		elif t>self.ts[-1]:
			idx1 = len(self.ts)-2
			idx2 = len(self.ts)-1
		else:
			idx1 = rng[self.ts<t][-1]
			idx2 = rng[self.ts>t][0]
		t1 = self.ts[idx1]; t2 = self.ts[idx2]
		x1 = self.xs[idx1]; x2 = self.xs[idx2]
		y1 = self.ys[idx1]; y2 = self.ys[idx2]
		w1 = (t2-t)/(t2-t1)
		w2 = (t-t1)/(t2-t1)

		if x is None:
			if (len(x1)!=len(x2)):
				raise Exception('Provided X for the anchored data were not all the same.'
						'You must therefore specify a X when interpolating.')
			x=x1

		#indices of the anchors
		anc1 = self.anchors[idx1]
		anc2 = self.anchors[idx2]
		#x position of the anchors, and interpolated x position
		x_anc1 = x1[anc1]
		x_anc2 = x2[anc2]
		x_anc = w1*x_anc1 + w2*x_anc2
		#index of interpolated anchor in the output X axis
		anc = array([argmin(fabs(x_-x)) for x_ in x_anc], dtype=int)

		def get_slice(slice_idx, anc_):
			if slice_idx==0:
				return slice(None,anc_[0])
			elif slice_idx==len(anc):
				return slice(anc_[-1],None)
			else:
				return slice(anc_[slice_idx-1], anc_[slice_idx])
			
		#each "f[i]" is a linear interpolation of the data separated by anchors (i-1,i),
		#and the x values are rescalled into the [0,1) interval
		def get_intp_func(slc_, x_, y_):
			x_ = x_[slc_]
			y_ = y_[slc_]
			#print x_, y
			if len(x_)>1:
				x_ = (x_-x_[0])/(x_[-1]-x_[0])
				f = interp1d(x_, y_)
			elif len(x_)==1:
				f = lambda(x): y_[0]
			else:
				raise Exception('Empty slice!')
			return f

		rs = range(len(anc)+1)
		y = empty_like(x)
		for slice_idx in rs:
			slc1 = get_slice(slice_idx, anc1)
			slc2 = get_slice(slice_idx, anc2)
			slc = get_slice(slice_idx, anc)
			f1 = get_intp_func(slc1, x1, y1)
			f2 = get_intp_func(slc2, x2, y2)

			x_ = x[slc]
			if len(x_)>1:
				#map to [0,1) interval
				x_ = (x_-x_[0])/(x_[-1]-x_[0])
				y[slc] = (w1*f1(x_) + w2*f2(x_))
			elif len(x_)==1:
				#take average between anchors (x=0.5)
				y[slc] = (w1*f1(0.5) + w2*f2(0.5))

		return y

