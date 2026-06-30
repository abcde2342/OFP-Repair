import math
import mpmath
import struct
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import FP_targets

mpmath.mp.prec=256
eps = 1e-30
plt.style.use("seaborn-darkgrid")


def doubleToInt(f):
	return struct.unpack('<Q', struct.pack('<d', f))[0]
def intToDouble(i):
	return struct.unpack('<d', struct.pack('<Q', i))[0]

def acc_sum(xlist):
	s = mpmath.mpf(0.0)
	for x in xlist:
		s += x
	return float(s)

class radiusInfo:
	def __init__(self, centerPoint, radius, func_id):
		self.func=FP_targets.getfp(index=func_id)
		self.orac=FP_targets.getoracle(index=func_id)
		self.centerPoint = centerPoint
		self.radius = radius
		self.insideSampleNum = 3000
		self.pointList=[]
		self.avgPoint=()
		self.xRange=()
		self.yRange=()

	def sampling(self):
		fplb = self.centerPoint - abs(self.centerPoint * self.radius)
		fpub = self.centerPoint + abs(self.centerPoint * self.radius)
		i64lb = doubleToInt(fplb)
		i64ub = doubleToInt(fpub)
		if i64lb > i64ub:
			i64lb, i64ub = i64ub, i64lb
		fpnum = i64ub - i64lb
		for j in range(self.insideSampleNum):
			i64x = random.randint(i64lb, i64ub)
			fpx = intToDouble(i64x)
			fpy = self.func(fpx)
			self.pointList.append((fpx, fpy))
		self.pointList.sort(key=lambda point: point[0])

		# maintain average point
		self._highprec_avg_point()
		# maintain x range and y range
		self._range_maintain()

		self.score_adjacent()
		self.score_adjacent_amend()
		self.score_corrcoef()

		self.scores = {
			'adj': self.score_adj,
			'adj_amend': self.score_adj_amend,
			'corr': self.score_corr,
		}

		return self.scores

	def _highprec_avg_point(self):
		sumx = mpmath.mpf(0.0)
		sumy = mpmath.mpf(0.0)
		for point in self.pointList:
			sumx += point[0]
			sumy += point[1]
		avgx = float(sumx / len(self.pointList))
		avgy = float(sumy / len(self.pointList))
		self.avgPoint = (avgx, avgy)

	def _range_maintain(self):
		ymax = max(self.pointList, key=lambda point: point[1])[1]
		ymin = min(self.pointList, key=lambda point: point[1])[1]
		self.yRange = (ymin, ymax)

		xmax = max(self.pointList, key=lambda point: point[0])[0]
		xmin = min(self.pointList, key=lambda point: point[0])[0]
		self.xRange = (xmin, xmax)

	def score_adjacent(self):
		num = len(self.pointList)
		score = 0.0
		for i in range(num-1):
			ydiff = abs(self.pointList[i][1]-self.pointList[i+1][1])
			score += ydiff
		ymin, ymax = self.yRange
		self.score_adj = score / (ymax-ymin+eps)
		return self.score_adj

	def score_adjacent_amend(self):
		# minus the largest 1 ydiff
		num = len(self.pointList)
		ydifflist = []
		score = 0.0
		for i in range(num-1):
			ydiff = abs(self.pointList[i][1]-self.pointList[i+1][1])
			score += ydiff
			ydifflist.append(ydiff)
		ymin, ymax = self.yRange
		if max(ydifflist) / (ymax-ymin+eps) > 0.4:
			score -= max(ydifflist)
		self.score_adj_amend = score / (ymax-ymin+eps)
		return self.score_adj_amend

	def score_corrcoef(self):
		xlist=[point[0] for point in self.pointList]
		ylist=[point[1] for point in self.pointList]
		self.score_corr = abs(np.corrcoef(xlist, ylist)[0, 1])
		return self.score_corr

	def print(self):
		print("{:+.2e}".format(self.radius), end='\t')
		for key, val in self.scores.items():
			print(key, "{:.4e}".format(val), end=' ')
		print()

class localSampling:
	def __init__(self, func_id):
		self.func_id = func_id
		self.func=FP_targets.getfp(index=func_id)
		self.orac=FP_targets.getoracle(index=func_id)
		self.centerPoint=1.23456789e-03
		self.minRadius = 1e-16
		self.maxRadius = 1e-6
		self.radiusNum = 100
		self.radiusList = []
	def genRadius(self, rList=None):
		if rList is None:
			radiusLnOffset = (math.log(self.maxRadius)-math.log(self.minRadius)) / (self.radiusNum - 1)
			for i in range(self.radiusNum):
				curRadius = math.exp(math.log(self.minRadius) + radiusLnOffset * i)
				self.radiusList.append(radiusInfo(self.centerPoint, curRadius, self.func_id))
		else:
			for curRadius in rList:
				self.radiusList.append(radiusInfo(self.centerPoint, curRadius, self.func_id))

	def sampleRadius(self):
		for neighbour in self.radiusList:
			neighbour.sampling()
			# neighbour.print()
		# minCorr=min(self.radiusList, key=lambda neighbour: neighbour.score_corr)
		# maxAdj =max(self.radiusList, key=lambda neighbour: neighbour.score_adj)
		# maxAdjAmend =max(self.radiusList, key=lambda neighbour: neighbour.score_adj_amend)
		# print("Selected by corr:     ", end=' ')
		# minCorr.print()
		# print("Selected by adj:      ", end=' ')
		# maxAdj.print()
		# print("Selected by adjamend: ", end=' ')
		# maxAdjAmend.print()

		# draw_list = [minCorr, maxAdj, maxAdjAmend]
		# self.drawNeighbour(draw_list)

	def drawNeighbour(self, draw_list, avg=True, oracle=True, index_list=None, savename='micro-structure.png'):
		num = len(draw_list)
		rows = 1
		cols = num
		fig, axs = plt.subplots(rows, cols, figsize=(cols*5, rows*4))
		for j, neighbour in enumerate(draw_list):
			xlist = [point[0] for point in neighbour.pointList]
			ylist = [point[1] for point in neighbour.pointList]
			olist = [self.orac(x) for x in xlist]
			if rows == 1:
				curaxs = axs[j]
			else:
				curaxs = axs[j//cols, j%cols]
			# FP results
			curaxs.scatter(xlist, ylist, s=10, color='tab:blue', label='FP results')
			# oracle results
			if oracle is True:
				curaxs.scatter(xlist, olist, s=10, color='tab:orange', label='Oracle results')
			# average point
			if avg is True:
				curaxs.scatter([neighbour.avgPoint[0]], [neighbour.avgPoint[1]], s=20, marker='x')
			# Set lim
			xmin, xmax = neighbour.xRange
			ymin, ymax = neighbour.yRange
			xr = (xmax - xmin)/2
			xc = xmin + xr

			oracle_y_avg = acc_sum(olist) / len(olist)
			oracle_y_max = max(olist)
			oracle_y_min = min(olist)

			ymin = min(ymin, oracle_y_min)
			ymax = max(ymax, oracle_y_max)
			yr = max(ymax - oracle_y_avg, oracle_y_avg - ymin)
			yc = oracle_y_avg

			# yr = (ymax - ymin)/2
			# yc = ymin + yr

			curaxs.set_xlim([xc-1.2*xr, xc+1.2*xr])
			curaxs.set_ylim([yc-1.2*yr, yc+1.2*yr])
			# Title
			title = ''
			if index_list is not None:
				title = '('+index_list[j]+')  '
			title += 'radius: '+"{:.1e}".format(neighbour.radius)+'\n'
			curaxs.set_title(title, y=1, pad=-20, fontsize=30)
			# ttl = curaxs.title
			# ttl.set_position([.5, 0.05])
			# Remove ticks
			curaxs.set_xticks([])
			curaxs.set_yticks([])

		if rows == 1:
			lastaxs = axs[-1]
		else:
			lastaxs = axs[rows-1, cols-1]
		handles, labels = lastaxs.get_legend_handles_labels()
		lgnd = fig.legend(handles, labels,
			loc='lower left',
			shadow=True,
			frameon=True,
			fancybox=True,
			fontsize=24
			)
		lgnd.legendHandles[0]._sizes = [90]
		lgnd.legendHandles[1]._sizes = [90]

		fig.tight_layout()
		plt.savefig(savename)
		# plt.show()

figurepath = 'figures'
if not os.path.exists(figurepath):
	os.mkdir(figurepath)

def draw_0():
	ls = localSampling(func_id=0)
	rList = [1e-2, 1e-3, 1e-6, 1e-9, 1e-10, 1e-11]
	iList = ['1a','1b','1c','1d','1e','1f']
	ls.genRadius(rList)
	ls.sampleRadius()
	outputpath = os.path.join(figurepath, 'tight-micro-structure-0.png')
	ls.drawNeighbour(ls.radiusList, avg=False, oracle=True, index_list=iList, savename=outputpath)
	print('Micro-Structure Figure Saved to', outputpath)

def draw_1():
	ls = localSampling(func_id=1)
	rList = [1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 4e-14]
	iList = ['2a','2b','2c','2d','2e', '2f']
	ls.genRadius(rList)
	ls.sampleRadius()
	outputpath = os.path.join(figurepath, 'tight-micro-structure-1.png')
	ls.drawNeighbour(ls.radiusList, avg=False, oracle=True, index_list=iList, savename=outputpath)
	print('Micro-Structure Figure Saved to', outputpath)

def draw_11():
	ls = localSampling(func_id=2)
	rList = [1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 3e-14]
	iList = ['3a','3b','3c','3d','3e', '3f']
	ls.genRadius(rList)
	ls.sampleRadius()
	outputpath = os.path.join(figurepath, 'tight-micro-structure-2.png')
	ls.drawNeighbour(ls.radiusList, avg=False, oracle=True, index_list=iList, savename=outputpath)
	print('Micro-Structure Figure Saved to', outputpath)


def main():
	draw_0()
	draw_1()
	draw_11()

if __name__ == '__main__':
	main()