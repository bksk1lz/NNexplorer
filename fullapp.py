import os
import pandas as pd
import numpy as np
from scipy.spatial import KDTree

from bokeh.io import curdoc
from bokeh.layouts import layout, widgetbox
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import Slider, TextInput, RadioButtonGroup, Dropdown
from bokeh.plotting import figure

#a comment!
#and a 2nd comment!!!!
filelist = ['data/' + f for f in os.listdir('data') if f.endswith('.xls')]
#now i added a comment here

def treequerycalc(df, k):
    
    points = list(zip(df.X, df.Y)) # make a list of X-Y pairs (lists) from df.X and df.Y
    tree = KDTree(points)
    [nndist, indices] = tree.query(points, k)

    return nndist, indices

def GetNNs(SampleList, dfList, k):
	dbar = []
	dexp = []
	nni = []

	for sample in SampleList:
		nndlist = []
		npts = 0
		for entry in dfList:
			if entry['name'].startswith(sample['name']):
				nndlist.extend(list(entry['df'].loc[:, 'dist1':].iloc[:, :k].values.flat))
				npts += len(entry['df'].index)
		NNdbar = np.mean(nndlist)
		d_expected = 0.5 * (5 * (3072 * 2304) / npts) ** 0.5
		
		dbar.append(NNdbar)
		dexp.append(d_expected)
		nni.append(NNdbar / d_expected)
	return dbar, nni

def mad(data, axis=None):
    return np.mean(np.absolute(data - np.mean(data, axis)), axis)
	
def angledev(ang1, ang2):
    diff = abs(ang1 - ang2)
    anglediff = min(abs(180 - diff), diff)
    return anglediff	
		
def getMeanQDiffs(SampleList, dfList, qName, k):
	qbar = []
	qi = []
	for sample in SampleList:
		nnqlist = []
		for entry in dfList:
			if entry['name'].startswith(sample['name']):
				quantity = np.array(entry['df'][qName])
				ind0 = np.array(entry['df'].ind0)
				
				for col in np.arange(k):
					indices = np.array(entry['df'].loc[:, 'ind1':].iloc[:, col])
					if qName == 'Orientation':
						qDiffList = [angledev(quantity[i], quantity[indices[i]]) for i in ind0]
					else:
						qDiffList = [abs(quantity[i] - quantity[indices[i]]) for i in ind0]
					nnqlist.extend(qDiffList)
         
       
		NNQbar = np.mean(nnqlist)
		NNQexp = mad(nnqlist)
    
		qbar.append(NNQbar)
		qi.append(NNQbar / NNQexp)
	
	return qbar, qi
	
#Load Data from files
SampleList = [{'name': 'r1c1', 'time': 130, 'color': '#38AB80'},
              {'name': 'r1c2', 'time': 1020, 'color': '#38AB80'},
              {'name': 'r1c3', 'time': 2100, 'color': '#38AB80'},
              {'name': 'r1c4', 'time': 85, 'color': '#E083F9'},
              {'name': 'r2c1', 'time': 380, 'color': '#E083F9'},
              {'name': 'r2c2', 'time': 270, 'color': '#E083F9'},
              {'name': 'r2c3', 'time': 600, 'color': '#E083F9'},
              {'name': 'r2c4', 'time': 1015, 'color': '#E083F9'},
              {'name': 'r2c5', 'time': 1755, 'color': '#E083F9'}]

k = 20
dfList = []
for file in filelist:
    name = file.split('.')[0].split('/')[1]
    df = pd.read_csv(file, delim_whitespace = True)
    [nndist, indices] = treequerycalc(df, k)
    distdf = pd.DataFrame(nndist, columns = ['dist' + str(i) for i in np.arange(k)])
    inddf = pd.DataFrame(indices, columns = ['ind' + str(i) for i in np.arange(k)])
    dfAll = pd.concat([df, distdf, inddf], axis = 1)
    
    dfList.append({'name': name, 'df': dfAll})

#Set up Data
x = [sample['time'] for sample in SampleList]
color = [sample['color'] for sample in SampleList]
dbar, nni = GetNNs(SampleList, dfList, 1)
qbar, qnni = getMeanQDiffs(SampleList, dfList, 'Area', 1)

s1 = ColumnDataSource(data = dict(x = x, y = dbar, color = color))
s2 = ColumnDataSource(data = dict(x = x, y = qbar, color = color))
s3 = ColumnDataSource(data = dict(x = dbar, y = qbar, color = color))
	
#Set up plot
p1 = figure(width = 600, height = 400, title = 'NND',
			x_axis_label = 'Formation time, s')
p1.circle('x', 'y', color = 'color', source = s1, size = 12)

p2 = figure(width = 600, height = 400, title = 'NNQ',
			x_axis_label = 'Formation time, s')
p2.circle('x', 'y', color = 'color', source = s2, size = 12)

p3 = figure(width = 600, height = 400, title = 'NNQ',
			x_axis_label = 'Distance',
			y_axis_label = 'Area')
p3.circle('x', 'y', color = 'color', source = s3, size = 12)

# Set up widgets
kneighb = Slider(title = 'Number of nearest neighbors to average', value = 1, start = 1, end = (k - 1), step = 1)
radio_button_group = RadioButtonGroup(labels=["NN Distance", "NN Index"], active=0)
menu = [('Area', 'Area'), ('Coherency', 'Coherency'), ('Orientaiton', 'Orientation')]
dropdown = Dropdown(label = 'Parameter', menu = menu)
menu2 = [('Area', 'Area'), ('Coherency', 'Coherency'), ('Orientaiton', 'Orientation'), ('Distance', 'Distance')]
x3drop = Dropdown(label = 'X-axis', menu = menu2)
y3drop = Dropdown(label = 'Y-axis', menu = menu2)
w = widgetbox(kneighb, radio_button_group, dropdown, x3drop, y3drop)

#Set up callbacks
def update_data(attrname, old, new):
	ks = kneighb.value
	rbg = radio_button_group.active
	qName = dropdown.value
	x3name = x3drop.value
	y3name = y3drop.value
	
	x = [sample['time'] for sample in SampleList]
	color = [sample['color'] for sample in SampleList]
	dresult = GetNNs(SampleList, dfList, ks)
	s1y = dresult[rbg]
	
	qresult = getMeanQDiffs(SampleList, dfList, qName, ks)
	s2y = qresult[rbg]
	
	s3x = getMeanQDiffs(SampleList, dfList, x3name, ks)[rbg]
	s3y = getMeanQDiffs(SampleList, dfList, y3name, ks)[rbg]
	
	s1.data = dict(x = x, y = s1y, color = color)
	p1.title.text = 'NNI, Average of nearest ' + str(ks) + ' neighbors'
		
	s2.data = dict(x = x, y = s2y, color = color)
	p2.title.text = 'NN ' + qName + ' difference, average of nearest ' + str(ks) + ' neighbors'
	
	s3.data = dict(x = s3x, y = s3y, color = color)
	p3.title.text = 'Plot of ' + y3name + ' vs. ' + x3name
	p3.xaxis.axis_label = x3name
	p3.yaxis.axis_label = y3name
	
kneighb.on_change('value', update_data)
radio_button_group.on_change('active', update_data)
dropdown.on_change('value', update_data)
x3drop.on_change('value', update_data)
y3drop.on_change('value', update_data)

#add to document
curdoc().add_root(layout([[w, p1], [p3, p2]], sizing_mode = 'stretch_both'))
curdoc().title = "test"
