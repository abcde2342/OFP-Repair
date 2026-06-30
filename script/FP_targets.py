import math
import mpmath

def fun1(x):
	try:
		return (1-math.cos(x))/(x*x)
	except:
		return 0
def fun2(x):
	return (math.exp(x)-2)+math.exp(-x)
def fun3(x):
	try:
		return math.sqrt((math.exp(2*x)-1)/(math.exp(x)-1))
	except:
		return 0

def fun_high1(x):
	x = mpmath.mpf(x)
	if (abs(x) < 1e-6):
		y = 0.5 + x*x*(-1.0/24.0 + x*x*(1.0/720.0))
	else:
		y = (1-mpmath.cos(x))/(x*x)
	return float(y)
def fun_high2(x):
	x = mpmath.mpf(x)
	if (abs(x) < 1e-6):
		y = x*x*(1.0 + x*x*(1.0/12.0 + x*x*(1.0/360.0)))
	else:
		y = (mpmath.exp(x)-2)+mpmath.exp(-x)
	return float(y)
def fun_high3(x):
	x = mpmath.mpf(x)
	if (abs(x) < 1e-6):
		sqrt2 = mpmath.sqrt(2)
		y = sqrt2 + x*(sqrt2/4.0 + x*(3.0*sqrt2/32.0 + x*(7.0*sqrt2/384.0)))
	else:
		y = mpmath.sqrt((mpmath.exp(2*x)-1) / (mpmath.exp(x)-1))
	# y = mpmath.sqrt((mpmath.exp(2*x)-1) / (mpmath.exp(x)-1))
	return float(y)

def getfp(index):
	fplist = [
		fun1,
		fun2,
		fun3,
	]
	if 0 <= index < len(fplist):
		return fplist[index]
	else:
		print("Invalid function index, please re-check.")

def getoracle(index):
	oraclelist = [
		fun_high1,
		fun_high2,
		fun_high3,
	]
	if 0 <= index < len(oraclelist):
		if mpmath.mp.prec < 256:
			mpmath.mp.prec=256
		return oraclelist[index]
	else:
		print("Invalid oracle function index, please re-check.")
