import os
import pandas as pd
import numpy as np
from scipy.spatial import KDTree

from bokeh.io import curdoc
from bokeh.layouts import layout, widgetbox
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import Slider, TextInput, RadioButtonGroup, Dropdown
from bokeh.plotting import figure


filelist = ['data/' + f for f in os.listdir('data') if f.endswith('.xls')]


def treequerycalc(df, k):
    
    points = list(zip(df.X, df.Y)) # make a list of X-Y pairs (lists) from df.X and df.Y
    tree = KDTree(points)
    [nndist, indices] = tree.query(points, k)

    return nndist, indices

def GetNNs(SampleList, dfList, k):
	dbar = []
	nni = []
	dsd = []
	nnisd = []

	for sample in SampleList:
		nndlist = []
		npts = 0
		for entry in dfList:
			if entry['name'].startswith(sample['name']):
				nndlist.extend(list(entry['df'].loc[:, 'dist1':].iloc[:, :k].values.flat))
				npts += len(entry['df'].index)
		
		NNdbar = np.mean(nndlist)
		NNdsd = np.std(nndlist)
		d_expected = 0.5 * (5 * (3072 * 2304) / npts) ** 0.5
		
		dbar.append(NNdbar)
		nni.append(NNdbar / d_expected)
		dsd.append(NNdsd)
		nnisd.append(NNdsd / d_expected)
		
	return dbar, nni, dsd, nnisd

def mad(data, axis=None):
    return np.mean(np.absolute(data - np.mean(data, axis)), axis)
	
def angledev(ang1, ang2):
    diff = abs(ang1 - ang2)
    anglediff = min(abs(180 - diff), diff)
    return anglediff	
		
def getMeanQDiffs(SampleList, dfList, qName, k):
	qbar = []
	qi = []
	qsd = []
	qisd = []
	
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
		NNQsd = np.std(nnqlist)
		NNQexp = mad(nnqlist)
    
		qbar.append(NNQbar)
		qi.append(NNQbar / NNQexp)
		qsd.append(NNQsd)
		qisd.append(NNQsd / NNQexp)
	
	return qbar, qi, qsd, qisd
	
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
dbar, nni, dsd, nnisd = GetNNs(SampleList, dfList, 1)
qbar, qnni, qsd, qisd = getMeanQDiffs(SampleList, dfList, 'Area', 1)

sd_factor = 4

dsdplus = [bar + sd / (2 * sd_factor) for bar, sd in zip(dbar, dsd)]
dsdminus = [bar - sd / (2 * sd_factor) for bar, sd in zip(dbar, dsd)]

qsdplus = [bar + sd / (2 * sd_factor) for bar, sd in zip(qbar, qsd)]
qsdminus = [bar - sd / (2 * sd_factor) for bar, sd in zip(qbar, qsd)]
 

s1 = ColumnDataSource(data = dict(x = x, y = dbar, 
	                              sdplus = dsdplus, sdminus = dsdminus, 
					              color = color))
s2 = ColumnDataSource(data = dict(x = x, y = qbar, 
                                  sdplus = qsdplus, sdminus = qsdminus, 
								  color = color))
s3 = ColumnDataSource(data = dict(x = qbar, xsdplus = qsdplus, xsdminus = qsdminus,
	                              y = qbar, ysdplus = qsdplus, ysdminus = qsdminus,
								  color = color))
	
#Set up plot
p1 = figure(width = 600, height = 400, title = 'NND',
			x_axis_label = 'Formation time, s')
p1.segment('x', 'sdplus', 'x', 'sdminus', source = s1, line_width = 2)
p1.circle('x', 'y', color = 'color', source = s1, size = 12, alpha = 0.65)


p2 = figure(width = 600, height = 400, title = 'NNQ',
			x_axis_label = 'Formation time, s')

p2.segment('x', 'sdplus', 'x', 'sdminus', source = s2, line_width = 2)
p2.circle('x', 'y', color = 'color', source = s2, size = 12, alpha = 0.65)

p3 = figure(width = 600, height = 400, title = 'NNQ',
			x_axis_label = 'Distance',
			y_axis_label = 'Area')
p3.segment('x', 'ysdplus', 'x', 'ysdminus', source = s3, line_width = 2)
p3.segment('xsdplus', 'y', 'xsdminus', 'y', source = s3, line_width = 2)
p3.circle('x', 'y', color = 'color', source = s3, size = 12, alpha = 0.65)

# Set up widgets
kneighb = Slider(title = 'Number of nearest neighbors to average', value = 1, start = 1, end = (k - 1), step = 1)
radio_button_group = RadioButtonGroup(labels=["NN Distance", "NN Index"], active=0)
menu = [('Area', 'Area'), ('Coherency', 'Coherency'), ('Orientaiton', 'Orientation')]
dropdown = Dropdown(label = 'Parameter', menu = menu, default_value = 'Area')
menu2 = [('Area', 'Area'), ('Coherency', 'Coherency'), ('Orientaiton', 'Orientation'), ('Distance', 'Distance')]
x3drop = Dropdown(label = 'X-axis', menu = menu2, default_value = 'Area')
y3drop = Dropdown(label = 'Y-axis', menu = menu2, default_value = 'Area')
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
	s1sd = dresult[(rbg + 2)]
	s1sdplus = [bar + sd / (2 * sd_factor) for bar, sd in zip(s1y, s1sd)]
	s1sdminus = [bar - sd / (2 * sd_factor) for bar, sd in zip(s1y, s1sd)]
	
	qresult = getMeanQDiffs(SampleList, dfList, qName, ks)
	s2y = qresult[rbg]
	s2sd = qresult[(rbg + 2)]
	s2sdplus = [bar + sd / (2 * sd_factor) for bar, sd in zip(s2y, s2sd)]
	s2sdminus = [bar - sd / (2 * sd_factor) for bar, sd in zip(s2y, s2sd)]
	
	if x3name == 'Distance':
		s3x = dresult[rbg]
		s3xsd = dresult[(rbg + 2)]
	else:
		s3xqresult = getMeanQDiffs(SampleList, dfList, x3name, ks)
		s3x = s3xqresult[rbg]
		s3xsd = s3xqresult[(rbg + 2)]
	
	s3xsdplus = [bar + sd / (2 * sd_factor) for bar, sd in zip(s3x, s3xsd)]
	s3xsdminus = [bar - sd / (2 * sd_factor) for bar, sd in zip(s3x, s3xsd)]	
	
	if y3name == 'Distance':
		s3y = dresult[rbg]
		s3ysd = dresult[(rbg + 2)]
	else:
		s3yqresult = getMeanQDiffs(SampleList, dfList, y3name, ks)
		s3y = s3yqresult[rbg]
		s3ysd = s3yqresult[(rbg + 2)]
		
	s3ysdplus = [bar + sd / (2 * sd_factor) for bar, sd in zip(s3y, s3ysd)]
	s3ysdminus = [bar - sd / (2 * sd_factor) for bar, sd in zip(s3y, s3ysd)]	
	
	s1.data = dict(x = x, y = s1y, sdplus = s1sdplus, sdminus = s1sdminus, 
	               color = color)
	p1.title.text = 'NNI, Average of nearest ' + str(ks) + ' neighbors'
		
	s2.data = dict(x = x, y = s2y, sdplus = s2sdplus, sdminus = s2sdminus, 
	               color = color)
	p2.title.text = 'NN ' + qName + ' difference, average of nearest ' + str(ks) + ' neighbors'
	
	s3.data = dict(x = s3x, xsdplus = s3xsdplus, xsdminus = s3xsdminus,
                   y = s3y, ysdplus = s3ysdplus, ysdminus = s3ysdminus, 
				   color = color)
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
