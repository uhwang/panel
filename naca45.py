# naca45.py

import math
import numpy as np

class InvalidNacaDigitException(Exception):
	pass
	
class CheckMemoryException(Exception):
	pass
	
class InvalidSpacingException(Exception):
	pass
	
#-------------------------------------------------
# NACA 4/5
#  
#  4 DIGIT : X1 X2 X3 X 4
#  x1 : Maximum camber (m), one hundredth
#  x2 : Position of maximum camber (p), one theth
#  x3, x4 : Maximum thickness, one hundredth 
#-------------------------------------------------

# Input: const char* naca, 
# Output: int* naca1, DType* fcmax, DType* fcpos, DType* ftmax

c2d = lambda x: ord(x)-ord('0')
X3  = lambda x: (x*x*x)
	
def naca_geom(naca):

	lend = len(naca)
	
	if lend < 2 or lend == 3:
		#print("Error => Insufficient number of digits(%s)"%naca)
		#raise InvalidNacaDigitException("Error => Insufficient number of digits(%s)"%naca)
		return 0,0,0,0

	
	elif lend == 2:

		fcmax, fcpos = 0.0,0.0
		ftmax = 10.0*c2d(naca[0]) + c2d(naca[1])
		naca1 = int(ftmax)
		ftmax *= 0.01

	elif lend == 4:

		fcmax = c2d(naca[0])*0.01
		fcpos = c2d(naca[1])*0.1
		ftmax = 10.0*c2d(naca[2]) + c2d(naca[3])
		naca1 = int(ftmax)
		naca1 += 100*(10*c2d(naca[0]) + c2d(naca[1]))
		ftmax *= 0.01

	return naca1, fcmax, fcpos, ftmax

#------------------------------------------------
# naca : string
# x    : x-coordinate of foil surface
# y    : y-coordinate of foil surface
# np   : number of points (even)
# spc  : 0 is equal spacing
#        1 is cosine spacing
#        2 is half-cosine spacing at L.E.
#        3 is half-cosine spacing at T.E.
#------------------------------------------------


def equal_spacing(x, np1):

	np=np1*2
	
	#--- Normal Spacing ---
	
	dth = 1.0 / (np1-1)
	
	for i in range(np1):
	
		x[i] = 1-dth * i
		x[np-i-1] = x[i]
	

def cosine_spacing(x, np1):

	np=np1*2
	
	#--- Cosine Spacing ---
	
	dth = math.pi / (np1-1)
		
	for i in range(np1):
	
		x[i] = 1-0.5 * ( 1 - math.cos( dth * i ) )
		x[np-i-1] = x[i]
	


def half_cosine_spacing_le(x, np1):

	np=np1*2
	
	#*--- Half-cosine Spacing at L.E. ---*
	
	dth = 0.5*math.pi / (np1-1)
	
	for i in range(np1):
	
		x[i] = 1-math.cos(PI_2-dth*i)
		x[np-i-1] = x[i]
	

def half_cosine_spacing_te(x, np1):

	np=np1*2
	
	# *--- Half-cosine Spacing at T.E. ---*
	
	dth = 0.5 * math.pi / (np1-1)
		
	for i in range(np1):
	
		x[i] = math.cos(dth*i)
		x[np-i-1] = x[i]
	

#
# Input: naca, tmax, m, p, xc
# Output: t, yc, beta
#
#
#   DType* t  , DType* yc, DType* beta)

def naca45_kernel(naca, tmax, m , p, xc):

	r=0.0
	
	if xc < 1.e-10: 
	
		t = 0.0
		
	else:
		t = tmax * 5  * ( 0.2969 * math.sqrt(xc)
		          - xc * ( 0.1260
		          + xc * ( 0.3537
		          - xc * ( 0.2843
		          - xc *   0.1015 ))))
		           
	if m == 0:
	
		yc = 0.0
		beta = 0.0
	
	else:
	
		if naca < 9999:
		
			r = xc/p
			if xc < p:
			
				yc = 2*m*r - m*(r**2)
				beta = math.atan( 2*m*(1-r)/p )
			
			else:
			
				yc = m*(p**2)/(1-p)**2 * ((1-2*p)/(p**2) + 2*r - (r**2))
				beta = math.atan( 2*p*m/(1-p)**2 * (1-r) )
			
		
		else:
		
			r = xc / m
			if xc < m:
			
				yc = p * ( X3(r) - 3*(r**2) + (3-m)*r )
				beta = math.atan( p/m * ( 3*(r**2) - 6*r + 3 - m) )
			
			else:
			
				yc = p * ( 1 - xc )
				beta = math.atan( -p )
			
	return t, yc, beta

# Input: const char* naca, int hnp, int spc
# Output: DType* x, DType* y

def naca45(naca, x, y, hnp, spc): #, cam):

	dth, mc, pc, tc, t, yc, beta, xx=0., 0., 0., 0., 0., 0., 0., 0.
	
	np1, np, naca1 = 0, 0, 0
	
	np1 = hnp
	np = hnp*2
	
	naca1, mc, pc, tc = naca_geom(naca)
	
	if naca1==0 and mc == 0 and pc == 0 and tc == 0:
		raise InvalidNacaDigitException("Error: naca45 => invalid NACA digit(%s)"%naca)
	
	if len(x) ==0 or len(y) == 0:
		raise CheckMemoryException("Error: naca45 => x, y ")
	
	if spc < 0 or spc > 3:
		raise InvalidSpacingException("Error: naca45 => spacing(%d)"%spc)
			
	if   spc == 0: equal_spacing(x, np1)
	elif spc == 1: cosine_spacing(x, np1)
	elif spc == 2: half_cosine_spacingLE(x, np1)
	elif spc == 3: half_cosine_spacingTE(x, np1)
	
	for i in range(np1):
	
		t, yc, beta = naca45_kernel(naca1, tc, mc, pc, x[i])
		#cam[i] = yc
		
		# *--- Lower surface ---*
		xx = x[i]
		x[i] = xx + t*math.sin(beta)
		y[i] = yc - t*math.cos(beta)
		
		#*--- upper surface ---*
		j = np-i-1
		xx = x[j]
		x[j] = xx - t*math.sin(beta)
		y[j] = yc + t*math.cos(beta)
		
		#*--- Initialize ---*
		t = yc = beta = 0
	
def adjpan(x, y, hnp):

	np1=hnp
	np2=hnp+hnp-1
	
	x[0], x[np2] = 1.0, 1.0
	y[0], y[np2] = 0.0, 0.0
	
	x[np1-1], x[np1] = 0.0, 0.0
	y[np1-1], y[np1] = 0.0, 0.0
	
	for i in range(np1-1):
	
		x[np1+i] = x[np1+i+1]
		y[np1+i] = y[np1+i+1]
	
	return hnp+hnp-1

def test(dev):		
	from vgl import symbol, color, symbol
	from vgl import drawaxis, drawtick, drawlabel
	naca_digit="4412"
	npan = 100
	npnt = npan+2
	np1 = int(npan/2+1)
	spc = 1
	x=np.zeros((npnt),dtype=np.float32)
	y=np.zeros((npnt),dtype=np.float32)
	cam=np.zeros((npan),dtype=np.float32)
	
	naca45(naca_digit, x, y, np1, spc) #, cam)
	adjpan(x, y, np1)
	
	w = open("%s.tec"%naca_digit, "wt")
	w.write("variables = x,y\nzone t=\"NACA%s\", i = %d\n"%(naca_digit,x.size))
	for i in range(npnt):
		w.write("%f %f\n"%(x[i],y[i]))
	w.close()
	
	drawaxis.draw_axis(dev)
	drawtick.draw_tick_2d(dev)	
	drawlabel.draw_label_2d(dev)

	dev.polygon(x,y,color.BLUE, color.YELLOW, 0.001*dev.frm.hgt())
	dev.polyline(x,y,color.BLACK,0.003*dev.frm.hgt())
	#dev.polyline(x,cam,color.BLUE,0.003*dev.frm.hgt())
	
	sym = symbol.Circle(0.005, dev.frm.hgt(), 0.001)
	dev.begin_symbol(sym)
	for i in range(0,x.size): dev.symbol(x[i],y[i],sym)
	dev.end_symbol()
	
def main():
	from vgl import color, BBox, Frame, FrameManager, Data
	from vgl import DeviceWindowsMetafile, DeviceCairo
	#from vgl import drawfrm, symbol, drawtick, drawaxis, drawlabel
	
	data = Data(0,1.0,-0.5,0.5)
	fmm = FrameManager()
	frm_x2 = fmm.create(0.0,0.0,4,4, data)
	
	gbbox = fmm.get_gbbox()
	dev_img = DeviceCairo("naca.png", gbbox, 100)
	dev_img.fill_white()
	dev_img.set_plot(frm_x2)
	test(dev_img)
	dev_img.close()
	
	dev_wmf = DeviceWindowsMetafile("naca.wmf", gbbox)
	dev_wmf.set_device(frm_x2)
	test(dev_wmf)
	dev_wmf.close()
	
			
if __name__ == '__main__':
	main()

