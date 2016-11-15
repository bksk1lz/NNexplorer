import os
import pandas as pd
import numpy as np
from scipy.spatial import KDTree

from bokeh.io import curdoc
from bokeh.layouts import layout, widgetbox
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import Slider, TextInput, RadioButtonGroup, Dropdown, Paragraph
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
	result_names = ['mean', 'devbar', 'index', 'sd', 'devsd', 'indexsd']
	resultDF = pd.DataFrame(columns = result_names)
	
	for sample in SampleList:
		nnqlist = []
		qlist = []
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
					qlist.extend([quantity[i] for i in indices])
         
       
		NNQbar = np.mean(nnqlist)
		NNQsd = np.std(nnqlist)
		NNQexp = mad(nnqlist)
		
		Qbar = np.mean(qlist)
		Qsd = np.std(qlist)

		loopdata = np.array([[Qbar, NNQbar, (NNQbar / NNQexp), Qsd, NNQsd, (NNQsd / NNQexp)]])
		loopdf = pd.DataFrame(data = loopdata, columns = result_names)
		resultDF = resultDF.append(loopdf)
	
	return resultDF
	
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

k = 2
dfList = []
umperpx = 0.03236
for file in filelist:
    name = file.split('.')[0].split('/')[1]
    df = pd.read_csv(file, delim_whitespace = True)
    df.X = df.X * umperpx
    df.Y = df.Y * umperpx
    df.Area = df.Area * (umperpx ** 2)
    [nndist, indices] = treequerycalc(df, k)
    distdf = pd.DataFrame(nndist, columns = ['dist' + str(i) for i in np.arange(k)])
    inddf = pd.DataFrame(indices, columns = ['ind' + str(i) for i in np.arange(k)])
    dfAll = pd.concat([df, distdf, inddf], axis = 1)
    
    dfList.append({'name': name, 'df': dfAll})

#Set up Data
x = [sample['time'] for sample in SampleList]
color = [sample['color'] for sample in SampleList]
dbar, nni, dsd, nnisd = GetNNs(SampleList, dfList, 1)
qdf = getMeanQDiffs(SampleList, dfList, 'Area', 1)
qbar = qdf['mean']
qsd = qdf['sd']

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
p1 = figure(width = 600, height = 400, title = 'Nearest Neighbor Distance',
			x_axis_label = 'Formation time, s')
p1.segment('x', 'sdplus', 'x', 'sdminus', source = s1, line_width = 2)
p1.circle('x', 'y', color = 'color', source = s1, size = 12, alpha = 0.65)


p2 = figure(width = 600, height = 400, title = 'Area difference of Nearest Neighbors',
			x_axis_label = 'Formation time, s')

p2.segment('x', 'sdplus', 'x', 'sdminus', source = s2, line_width = 2)
p2.circle('x', 'y', color = 'color', source = s2, size = 12, alpha = 0.65)

p3 = figure(width = 600, height = 400, title = 'Plot of Area vs. Area',
			x_axis_label = 'Area (um^2)',
			y_axis_label = 'Area (um^2)')
p3.segment('x', 'ysdplus', 'x', 'ysdminus', source = s3, line_width = 2)
p3.segment('xsdplus', 'y', 'xsdminus', 'y', source = s3, line_width = 2)
p3.circle('x', 'y', color = 'color', source = s3, size = 12, alpha = 0.65)

# Set up widgets
kneighb = Slider(title = 'Number of nearest neighbors to average', value = 1, start = 1, end = (k - 1), step = 1)
radio_button_group1 = RadioButtonGroup(labels=["NN Distance", "NN Index"], active = 0)
radio_button_group2 = RadioButtonGroup(labels=["Average", "NN Difference", "NN Difference Index"], active=0)
menu = [('Area', 'Area'), ('Coherency', 'Coherency'), ('Orientaiton', 'Orientation')]
dropdown = Dropdown(label = 'Lower Right Parameter', menu = menu, value = 'Area')
menu2 = [('Area', 'Area'), ('Coherency', 'Coherency'), ('Orientaiton', 'Orientation'), ('Distance', 'Distance')]
x3drop = Dropdown(label = 'Lower Left X-axis', menu = menu2, value = 'Area')
y3drop = Dropdown(label = 'Lower Left Y-axis', menu = menu2, value = 'Area')

pTXT = Paragraph(text = """Welcome to the interactive Nearest Neighbor data explorer. The first toggle controls
 the top right plot, selecting average NN distance in microns, or NN Index (normalized by the number of domains).
 The next menus control the lower two plots: select which parameter to be plotted with the dropdown menus, and
 use the toggle to pick average value, average NN difference, or NN difference index. Green points are constant
 flow series, Pink are constant volume series. Error bars are 1/4 standard deviation.""")

w = widgetbox(pTXT, radio_button_group1, dropdown, x3drop, y3drop, radio_button_group2)

#Set up callbacks
def update_data(attrname, old, new):
	ks = kneighb.value
	rbg1 = radio_button_group1.active
	rbg1labels = ["NN Distance", "NN Index"]
	rbg2 = radio_button_group2.active
	rbg2labels = ["Average", "NN Difference", "NN Difference Index"]
	qName = dropdown.value
	x3name = x3drop.value
	y3name = y3drop.value
	
	x = [sample['time'] for sample in SampleList]
	color = [sample['color'] for sample in SampleList]
	dresult = GetNNs(SampleList, dfList, ks)
	s1y = dresult[rbg1]
	s1sd = dresult[(rbg1 + 2)]
	s1sdplus = [bar + sd / (2 * sd_factor) for bar, sd in zip(s1y, s1sd)]
	s1sdminus = [bar - sd / (2 * sd_factor) for bar, sd in zip(s1y, s1sd)]
	
	qdf = getMeanQDiffs(SampleList, dfList, qName, ks)
	s2y = qdf.iloc[:, rbg2]
	s2sd = qdf.iloc[:, (rbg2 + 3)]
	s2sdplus = [bar + sd / (2 * sd_factor) for bar, sd in zip(s2y, s2sd)]
	s2sdminus = [bar - sd / (2 * sd_factor) for bar, sd in zip(s2y, s2sd)]
	
	if x3name == 'Distance':
		s3x = dresult[rbg1]
		s3xsd = dresult[(rbg1 + 2)]
	else:
		s3xqresult = getMeanQDiffs(SampleList, dfList, x3name, ks)
		s3x = s3xqresult.iloc[:, rbg2]
		s3xsd = s3xqresult.iloc[:, (rbg2 + 3)]
	
	s3xsdplus = [bar + sd / (2 * sd_factor) for bar, sd in zip(s3x, s3xsd)]
	s3xsdminus = [bar - sd / (2 * sd_factor) for bar, sd in zip(s3x, s3xsd)]	
	
	if y3name == 'Distance':
		s3y = dresult[rbg1]
		s3ysd = dresult[(rbg1 + 2)]
	else:
		s3yqresult = getMeanQDiffs(SampleList, dfList, y3name, ks)
		s3y = s3yqresult.iloc[:, rbg2]
		s3ysd = s3yqresult.iloc[:, (rbg2 + 3)]
		
	s3ysdplus = [bar + sd / (2 * sd_factor) for bar, sd in zip(s3y, s3ysd)]
	s3ysdminus = [bar - sd / (2 * sd_factor) for bar, sd in zip(s3y, s3ysd)]	
	
	s1.data = dict(x = x, y = s1y, sdplus = s1sdplus, sdminus = s1sdminus, 
	               color = color)
	p1.title.text = rbg1labels[rbg1]
		
	s2.data = dict(x = x, y = s2y, sdplus = s2sdplus, sdminus = s2sdminus, 
	               color = color)
	p2.title.text = qName + ' ' + rbg2labels[rbg2]
	
	s3.data = dict(x = s3x, xsdplus = s3xsdplus, xsdminus = s3xsdminus,
                   y = s3y, ysdplus = s3ysdplus, ysdminus = s3ysdminus, 
				   color = color)
	p3.title.text = 'Plot of ' + y3name + ' vs. ' + x3name + ' ' + rbg2labels[rbg2]
	p3.xaxis.axis_label = x3name
	p3.yaxis.axis_label = y3name
	
kneighb.on_change('value', update_data)
radio_button_group1.on_change('active', update_data)
radio_button_group2.on_change('active', update_data)
dropdown.on_change('value', update_data)
x3drop.on_change('value', update_data)
y3drop.on_change('value', update_data)

#add to document
curdoc().add_root(layout([[w, p1], [p3, p2]], sizing_mode = 'scale_width'))
curdoc().title = "NN Explorer"
