from sys import prefix
from numpy.core.numeric import Inf
import pandas as pd
import numpy as np
import math
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.ticker import NullFormatter

import os, io

sns.set_style('darkgrid')

tableformat = pd.DataFrame(columns=[
	'id',
	'valid'
	'abs-origin',
	'abs-repair',
	'rel-origin',
	'rel-repair',
	])

stable_data  = tableformat.copy()
decayed_data = tableformat.copy()

eval_path   = 'data/evaluation.dat'
valid_path  = 'data/validation.dat'
funcid_path = 'data/funcID.dat'
infolevel_path = 'data/infolevel.dat'
with io.open(infolevel_path) as f:
	line = f.readline()
	if int(line) == 4:
		prefix = 'full-info-'
	elif int(line) == 0:
		prefix = 'no-info-'
	else:
		prefix = ''


def read_csv_data(line):
	items = []

	segs = line.split(',')
	for seg in segs:
		err = float(seg)
		if math.isnan(err) or math.isinf(err):
			err = 10 # Assume a large one
		# if err == 0.0:
		# 	err = 0
		# else:
		# 	err = math.log10(err)
		err = math.log10(err)
		items.append(err)

	line_decayed_data = {}
	line_decayed_data['rel-origin'] = items[0]
	line_decayed_data['rel-repair'] = items[1]
	line_decayed_data['abs-origin'] = items[2]
	line_decayed_data['abs-repair'] = items[3]

	line_stable_data = {}
	line_stable_data['rel-origin'] = items[4]
	line_stable_data['rel-repair'] = items[5]
	line_stable_data['abs-origin'] = items[6]
	line_stable_data['abs-repair'] = items[7]

	return line_decayed_data, line_stable_data

countid = 1
with io.open(eval_path) as f:
	evaldata = f.readlines()
with io.open(valid_path) as f:
	validdata = f.readlines()
with io.open(funcid_path) as f:
	funciddata = f.readlines()

assert len(evaldata) == len(validdata)
assert len(funciddata) == len(evaldata)

nums = len(evaldata)

for index in range(nums):
	funcid = funciddata[index].strip()
	validflag = validdata[index].strip()

	line_decayed_data, line_stable_data = read_csv_data(evaldata[index])
	
	line_stable_data['id'] = funcid
	line_decayed_data['id'] = funcid

	line_stable_data['valid'] = int(validflag)
	line_decayed_data['valid'] = int(validflag)

	stable_data = stable_data.append(line_stable_data, ignore_index=True)
	decayed_data = decayed_data.append(line_decayed_data, ignore_index=True)


# for index in range(30):
# 	chebpath = prefix+str(index)+cheb_suffix
# 	nonepath = prefix+str(index)+none_suffix
# 	if os.path.exists(chebpath):
# 		with io.open(chebpath) as f:
# 			data = f.readlines()
# 		fmt_max_data, fmt_mean_data = read_csv_data(data)
# 		fmt_max_data['id'] = str(countid)
# 		fmt_mean_data['id'] = str(countid)
# 		best_max_data=best_max_data.append(fmt_max_data, ignore_index=True)
# 		best_mean_data=best_mean_data.append(fmt_mean_data, ignore_index=True)

# 	if os.path.exists(nonepath):
# 		with io.open(nonepath) as f:
# 			data = f.readlines()
# 		fmt_max_data, fmt_mean_data = read_csv_data(data)
# 		fmt_max_data['id'] = str(countid)
# 		fmt_mean_data['id'] = str(countid)
# 		none_max_data=none_max_data.append(fmt_max_data, ignore_index=True)
# 		none_mean_data=none_mean_data.append(fmt_mean_data, ignore_index=True)

# 	if os.path.exists(chebpath) and os.path.exists(nonepath):
# 		countid += 1

# arrows
def draw_data(ax, pddata, prevkey, afterkey, fil=True, valid=True):
	#initialize a plot
	# ax = plt.figure(figsize=(5,7))

	#add start points
	# ax = sns.stripplot(data=pddata, x=prevkey, y='id', jitter=0)
	sns.stripplot(data=pddata, 
					x=prevkey, 
					y='id',
					orient='h', 
					order=pddata['id'], 
					size=5,
					jitter=0,
					color='black',
					ax=ax)

	#define arrows
	arrow_starts = pddata[prevkey].values
	arrow_lengths = pddata[afterkey].values - arrow_starts

	#invert x-axis
	xmin = min(np.concatenate((arrow_starts,pddata[afterkey].values)))
	xmax = max(np.concatenate((arrow_starts,pddata[afterkey].values)))
	# don't make xmax too large
	xmax = min(xmax, 5)

	xdiff = xmax-xmin
	# print(arrow_starts)
	# print(xmax+0.1*xdiff, xmin-0.1*xdiff)
	ax.set_xlim([xmax+1, xmin-0.1*xdiff])

	filtered_list = []
	for i, _id in enumerate(pddata['id']):
		# Do not draw invalid patch
		if pddata['valid'][i] == 0:
			continue
		if fil == True and _id in filtered_list:
			continue
		if arrow_lengths[i] == 0:
			continue

		if arrow_lengths[i] > 0:
			arrow_color='#662233'
		elif arrow_lengths[i] < 0:
			arrow_color='#337788'
		else:
			arrow_color='black'

		ax.arrow(arrow_starts[i],		#x start point
		i,								#y start point
		arrow_lengths[i],				#change in x 
		0,								#change in y
		head_width=0.4,					#arrow head width
		head_length=0.03*xdiff,				#arrow head length
		width=0.1,						#arrow stem width
		fc=arrow_color,					#arrow fill color
		ec=arrow_color)					#arrow edge color


	# sns.despine(left=True, bottom=True)                   #remove axes
	ax.grid(axis='y', color='0.8')                        #add a light grid

	# plt.show()

def main():
	figurepath = './figures'
	if not os.path.exists(figurepath):
		os.mkdir(figurepath)

	######### Enable LaTeX Text
	# matplotlib.rc('text', usetex=True)
	######### All Area ########
	fig, axes = plt.subplots(1,4, figsize=(12,6.25))

	# Max error
	draw_data(axes[0], decayed_data, 'abs-origin', 'abs-repair')
	draw_data(axes[1], decayed_data, 'rel-origin', 'rel-repair')
	# Set ticks
	axes[0].set_xlim([1, -21])
	axes[0].set_xticks([0, -5, -10, -15, -20])
	axes[0].set_xticklabels(['1', '1e-5', '1e-10', '1e-15', '1e-20'])
	axes[0].set_xlabel('Maximum Absolute Error\non Decayed Area', fontsize=14)
	axes[0].set_ylabel(None)
	axes[0].xaxis.set_tick_params(labelsize=12)
	axes[0].yaxis.set_tick_params(labelsize=10)
	# Align the ticks / function names
	for tick in axes[0].yaxis.get_majorticklabels():
		tick.set_x(-0.32)
		tick.set_ha("left")

	# only for noinfo
	# axes[1].set_xlim([5, -17])
	axes[1].set_xticks([0, -5, -10, -16])
	axes[1].set_xticklabels(['1', '1e-5', '1e-10', '1e-16'])
	axes[1].set_xlabel('Maximum Relative Error\non Decayed Area', fontsize=14)
	axes[1].set_ylabel(None)
	axes[1].xaxis.set_tick_params(labelsize=12)
	axes[1].yaxis.set_tick_params(labelsize=12)
	axes[1].yaxis.set_major_formatter(NullFormatter())
	# fig.suptitle("Maximum Error")
	# fig.tight_layout()
	# outputpath = os.path.join(figurepath, prefix+'max-err-all.pdf')
	# plt.savefig(outputpath)
	# print('Arrow Figure Saved to', outputpath)

	######## Stable Area ###########
	# fig, axes = plt.subplots(1,2, figsize=(7,4.8))
	# Max error
	draw_data(axes[2], stable_data, 'abs-origin', 'abs-repair')
	draw_data(axes[3], stable_data, 'rel-origin', 'rel-repair')
	# Set ticks
	# axes[0].set_xlim(-3, -24)
	axes[2].set_xticks([0, -5, -10, -15, -20])
	axes[2].set_xticklabels(['1', '1e-5', '1e-10', '1e-15', '1e-20'])
	axes[2].set_xlabel('Maximum Absolute Error\non Stable Area', fontsize=14)
	axes[2].set_ylabel(None)
	axes[2].xaxis.set_tick_params(labelsize=12)
	axes[2].yaxis.set_tick_params(labelsize=12)
	axes[2].yaxis.set_major_formatter(NullFormatter())

	axes[3].set_xticks([0, -5, -10, -16])
	axes[3].set_xticklabels(['1', '1e-5', '1e-10', '1e-16'])
	axes[3].set_xlabel('Maximum Relative Error\non Stable Area', fontsize=14)
	axes[3].set_ylabel(None)
	axes[3].xaxis.set_tick_params(labelsize=12)
	axes[3].yaxis.set_tick_params(labelsize=12)
	axes[3].yaxis.set_major_formatter(NullFormatter())
	# fig.suptitle("Maximum Error")
	fig.tight_layout()
	# outputpath = os.path.join(figurepath, prefix+'max-err-stable.pdf')
	outputpath = os.path.join(figurepath, prefix+'max-err.pdf')
	plt.savefig(outputpath)
	plt.show()
	print('Arrow Figure Saved to', outputpath)

	######## Max for none ########
	# fig, axes = plt.subplots(1,4, figsize=(12,5))
	# draw_data(axes[0], none_max_data, 'abs-origin', 'abs-repair', fil=False)
	# draw_data(axes[1], none_max_data, 'rel-origin', 'rel-repair', fil=False)
	# axes[0].set_xticks([0, -5, -10, -15])
	# axes[0].set_xticklabels(['1', '1e-5', '1e-10', '1e-15'])
	# axes[0].set_xlabel('Maximum Absolute Error', fontsize=11)
	# axes[0].set_ylabel('Function ID')

	# axes[1].set_xticks([10, 0, -10])
	# axes[1].set_xticklabels(['1e10', '1', '1e-10'])
	# axes[1].set_xlabel('Maximum Relative Error', fontsize=11)
	# axes[1].set_ylabel(None)
	# axes[1].yaxis.set_major_formatter(NullFormatter())
	# # fig.tight_layout()
	# # plt.show()

	# # fig, axes = plt.subplots(1,2, figsize=(6,5))
	# draw_data(axes[2], none_mean_data, 'abs-origin', 'abs-repair', fil=False)
	# draw_data(axes[3], none_mean_data, 'rel-origin', 'rel-repair', fil=False)
	# axes[2].set_xticks([0, -5, -10, -15, -20])
	# axes[2].set_xticklabels(['1', '1e-5', '1e-10', '1e-15', '1e-20'])
	# axes[2].set_xlabel('Average Absolute Error', fontsize=11)
	# # axes[2].set_ylabel('Function ID')
	# axes[2].set_ylabel(None)
	# axes[2].yaxis.set_major_formatter(NullFormatter())

	# axes[3].set_xticks([10, 0, -10])
	# axes[3].set_xticklabels(['1e10', '1', '1e-10'])
	# axes[3].set_xlabel('Average Relative Error', fontsize=11)
	# axes[3].set_ylabel(None)
	# axes[3].yaxis.set_major_formatter(NullFormatter())
	# fig.tight_layout()
	# outputpath = os.path.join(figurepath, 'none-info-err-eval.pdf')
	# plt.savefig(outputpath)
	# print('Arrow Figure Saved to', outputpath)
	# # plt.show()

main()
# draw_data(best_mean_data, 'abs-origin', 'abs-repair')
# draw_data(best_mean_data, 'bit-origin', 'bit-repair')
# draw_data(best_max_data, 'abs-origin', 'abs-repair')
# draw_data(best_max_data, 'bit-origin', 'bit-repair')

# draw_data(none_max_data, 'abs-origin', 'abs-repair')

# draw_data(best_mean_data, 'abs-origin', 'abs-repair')
# draw_data(none_mean_data, 'abs-origin', 'abs-repair')