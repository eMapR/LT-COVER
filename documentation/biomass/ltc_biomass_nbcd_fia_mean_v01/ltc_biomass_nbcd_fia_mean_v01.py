#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  2 15:05:01 2018

@author: braatenj
"""

from osgeo import gdal, ogr, osr
import os
import numpy as np
import numpy.ma as ma
import pandas as pd
import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import geopandas
from PIL import Image

#import seaborn as sns; sns.set()

import seaborn as sns
sns.set(style="ticks")

def get_dims(fileName):
  src = gdal.Open(fileName)
  ulx, xres, xskew, uly, yskew, yres  = src.GetGeoTransform()
  sizeX = src.RasterXSize
  sizeY = src.RasterYSize
  lrx = ulx + (sizeX * xres)
  lry = uly + (sizeY * yres)
  return [ulx,uly,lrx,lry,xres,-yres,sizeX,sizeY]

def get_intersec(files):
  ulxAll=[]
  ulyAll=[]
  lrxAll=[]
  lryAll=[]
  for fn in files:
    dim = get_dims(fn)
    ulxAll.append(dim[0])
    ulyAll.append(dim[1])
    lrxAll.append(dim[2])
    lryAll.append(dim[3])
  return([max(ulxAll),min(ulyAll),min(lrxAll),max(lryAll)])

def get_zone_pixels(feat, input_zone_polygon, input_value_raster, band, coords=[]): #, raster_band
  """
  feat =feature
  input_zone_polygon = shpF
  input_value_raster = trainF
  band = 1
  coords=[commonBox[0], commonBox[2], commonBox[3], commonBox[1]]
  """
  
  # Open data
  raster = gdal.Open(input_value_raster)
  shp = ogr.Open(input_zone_polygon)
  lyr = shp.GetLayer()
  
  # Get raster georeference info
  transform = raster.GetGeoTransform()
  xOrigin = transform[0]
  yOrigin = transform[3]
  pixelWidth = transform[1]
  pixelHeight = transform[5]
  
  sizeX = raster.RasterXSize
  sizeY = raster.RasterYSize
  lrx = xOrigin + (sizeX * pixelWidth)
  lry = yOrigin + (sizeY * pixelHeight)
  
  # Reproject vector geometry to same projection as raster
  #sourceSR = lyr.GetSpatialRef()
  #targetSR = osr.SpatialReference()
  #targetSR.ImportFromWkt(raster.GetProjectionRef())
  #coordTrans = osr.CoordinateTransformation(sourceSR,targetSR)
  #feat = lyr.GetNextFeature()
  #geom = feat.GetGeometryRef()
  #geom.Transform(coordTrans)
  
  # Get extent of feat
  geom = feat.GetGeometryRef()
  if (geom.GetGeometryName() == 'MULTIPOLYGON'):
    count = 0
    pointsX = []; pointsY = []
    for polygon in geom:
      geomInner = geom.GetGeometryRef(count)
      ring = geomInner.GetGeometryRef(0)
      numpoints = ring.GetPointCount()
      for p in range(numpoints):
        lon, lat, z = ring.GetPoint(p)
        pointsX.append(lon)
        pointsY.append(lat)
      count += 1
  elif (geom.GetGeometryName() == 'POLYGON'):
    ring = geom.GetGeometryRef(0)
    numpoints = ring.GetPointCount()
    pointsX = []; pointsY = []
    for p in range(numpoints):
      lon, lat, z = ring.GetPoint(p)
      pointsX.append(lon)
      pointsY.append(lat)

  else:
    sys.exit("ERROR: Geometry needs to be either Polygon or Multipolygon")

  #xmin = min(pointsX)  
  #xmax = max(pointsX)
  #ymin = min(pointsY)
  #ymax = max(pointsY)
  
  
  if len(coords) == 0: 
    xmin = xOrigin if (min(pointsX) < xOrigin) else min(pointsX)
    xmax = lrx if (max(pointsX) > lrx) else max(pointsX)
    ymin = lry if (min(pointsY) < lry) else min(pointsY)
    ymax = yOrigin if (max(pointsY) > yOrigin) else max(pointsY)
  else:
    xmin = coords[0] if (min(pointsX) < coords[0]) else min(pointsX)
    xmax = coords[1] if (max(pointsX) > coords[1]) else max(pointsX)
    ymin = coords[2] if (min(pointsY) < coords[2]) else min(pointsY)
    ymax = coords[3] if (max(pointsY) > coords[3]) else max(pointsY)
    
  # Specify offset and rows and columns to read
  xoff = int((xmin - xOrigin)/pixelWidth)
  yoff = int((yOrigin - ymax)/pixelWidth)
  xcount = int((xmax - xmin)/pixelWidth) #+1 !!!!!!!!!!!!!!!!!!!!! This adds a pixel to the right side
  ycount = int((ymax - ymin)/pixelWidth) #+1 !!!!!!!!!!!!!!!!!!!!! This adds a pixel to the bottom side
  
  #print(xoff, yoff, xcount, ycount)
              
  # Create memory target raster
  target_ds = gdal.GetDriverByName('MEM').Create('', xcount, ycount, 1, gdal.GDT_Byte)
  target_ds.SetGeoTransform((
    xmin, pixelWidth, 0,
    ymax, 0, pixelHeight,
  ))

  # Create for target raster the same projection as for the value raster
  raster_srs = osr.SpatialReference()
  raster_srs.ImportFromWkt(raster.GetProjectionRef())
  target_ds.SetProjection(raster_srs.ExportToWkt())

  # Rasterize zone polygon to raster
  gdal.RasterizeLayer(target_ds, [1], lyr, burn_values=[1])

  # Read raster as arrays
  dataBandRaster = raster.GetRasterBand(band)
  data = dataBandRaster.ReadAsArray(xoff, yoff, xcount, ycount).astype(np.float)
  bandmask = target_ds.GetRasterBand(1)
  datamask = bandmask.ReadAsArray(0, 0, xcount, ycount).astype(np.float)

  # data zone of raster
  dataZone = np.ma.masked_array(data,  np.logical_not(datamask))

  raster_srs = None
  raster = None
  shp = None
  lyr = None
  return [dataZone, [xmin,xmax,ymin,ymax]]

def f2hex(f2rgb, f):
    rgb = f2rgb.to_rgba(f)[:3]
    uint8 = [int(255*fc) for fc in rgb]
    return '#%02x%02x%02x' % (uint8[0],uint8[1],uint8[2])

def cdf(df, trgtCol):
  tot = 0
  cff = []
  cdf = []
  totPixs = sum(df[trgtCol])+0.0

  for i, row in df.iterrows():
    n = row[trgtCol]
    cff.append(n+tot)
    tot += n
    cdf.append(tot/totPixs)
  
  cdfDf = pd.DataFrame({'CDF':cdf, 'CFF':cff})
  return df.join(cdfDf)

def addColorToDf(df, trgField, colorField):
  col = df[trgField]
  minv = min(col)
  maxv = max(col)
  norm = colors.Normalize(vmin=minv, vmax=maxv)
  f2rgb = cm.ScalarMappable(norm=norm, cmap=cm.get_cmap('YlGnBu_r'))  
  color = [f2hex(f2rgb, c) for c in col]
  df[colorField] = color
  return df

def make2dHist(trainP, predP, name, xyLim, xLab, yLab, r, maeM, maeU, annoXY, figSize, hist2dFn):
  # make 2d hist plot 
  #hist2dP = sns.jointplot(trainP, predP, kind="hex", color='blue', xlim=(0,xyLim[0]), ylim=(0,xyLim[1]), size=figSize, ratio=5, space=0)
  hist2dP = sns.jointplot(trainP, predP, kind="hex", color='blue', xlim=(0,xyLim[0]), ylim=(0,xyLim[1]), size=figSize, ratio=5, space=0, marginal_kws=dict(bins=15))
  hist2dP.ax_joint.set_xlabel(xLab)
  hist2dP.ax_joint.set_ylabel(yLab)
  hist2dP.ax_joint.annotate('r: '+str(r)+'\nmae: '+str(int(round(maeM)))+'\n'+u'\u03BCae: '+str(int(round(maeU))), annoXY, fontsize=9)
  #plt.tight_layout()
  hist2dP.savefig(hist2dFn)

def make1dHist(histDf, units, figSize, hist1dFn):
  fig, ax = plt.subplots(figsize=(figSize, figSize))
  sns.lineplot(ax=ax, data=histDf, x="Units", y="Frequency", hue="Layer")
  ax.set(xlabel=units, ylabel='Frequency')
  #ax.legend(loc=1)
  handles, labels = ax.get_legend_handles_labels()
  ax.legend(loc='upper right', handles=handles[1:], labels=labels[1:])
  plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
  plt.tight_layout()
  plt.savefig(hist1dFn)

def makeCDF(histDf, units, figSize, cdfFn):
  fig, ax = plt.subplots(figsize=(figSize,figSize))
  #totPixs = sum(histDf['Frequency'])/2.0
  sns.lineplot(ax=ax, data=histDf, x="Units", y="CDF", hue="Layer")#.set_title(title)
  plt.ylim(0, 1)
  ax.set(xlabel=units, ylabel='Cumulative Density')
  #ax.legend(loc=4)
  handles, labels = ax.get_legend_handles_labels()
  ax.legend(loc='lower right', handles=handles[1:], labels=labels[1:])
  #ax.legend(handles=handles[1:], labels=labels[1:])
  plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
  plt.tight_layout()
  plt.savefig(cdfFn)

def makeHist(ref, pred, minv):
  binMin = minv
  binMax = np.max([np.max(ref), np.max(pred)])+1

  nBins = (binMax-binMin)+1

  theBins = np.linspace(binMin, binMax, nBins)
  tHist, jnk = np.histogram(ref, theBins)
  pHist, jnk = np.histogram(pred, theBins)
  
  tHistDf = pd.DataFrame({'Units':theBins[0:-1], 'Frequency':tHist, 'Layer':['Reference']*len(tHist)})
  pHistDf = pd.DataFrame({'Units':theBins[0:-1], 'Frequency':pHist, 'Layer':['Prediction']*len(pHist)})

  tHistDf = cdf(tHistDf, 'Frequency')
  pHistDf = cdf(pHistDf, 'Frequency')
  
  histDf = tHistDf.append(pHistDf)
  histDf = histDf.reset_index(drop=True)
  
  return histDf
  
def makeStatMaps(df, field, minMax, outDir, outName):
  outFile = os.path.join(outDir, outName)
  df.plot(column=field, cmap='YlGnBu', edgecolor='grey', linewidth=0.25, figsize=(6.5, 4.3), legend=True, vmin=minMax[0], vmax = minMax[1])
  plt.cbar.ax.tick_params(labelsize=20) 
  plt.axis('off')
  plt.tight_layout()
  plt.savefig(outFile, dpi=150)
    
def makeStatHist(df, field, xlab, outDir, outName, figSize):
  vals = df[field]
  binMin = np.min([np.max(vals), np.min(vals)])+1
  binMax = np.max([np.max(vals), np.max(vals)])+1
  nBins = (binMax-binMin)+1
  theBins = np.linspace(binMin, binMax, nBins)
  hist, jnk = np.histogram(vals, theBins)  
  #histDf = pd.DataFrame({'Units':theBins[0:-1], 'Frequency':hist}) #, 'Layer':['Reference']*len(hist)

  outFile = os.path.join(outDir, outName)
  fig, ax = plt.subplots(figsize=(figSize, figSize))
  #sns.lineplot(ax=ax, data=histDf, x="Units", y="Frequency")
  sns.distplot(vals, kde=False)
  ax.set(xlabel=xlab, ylabel='Frequency')
  #handles, labels = ax.get_legend_handles_labels()
  #ax.legend(handles=handles[1:], labels=labels[1:])
  #plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
  plt.tight_layout()
  plt.savefig(outFile)
  return np.median(vals)

def getMinMax(d1, d2):
  vals = d1.append(d2)
  return [min(vals), max(vals)]


# see https://docs.scipy.org/doc/numpy-1.15.1/reference/generated/numpy.histogram.html notes to inform about bin openings
def makeHistDf(fn, band, noData, year, name, binMin, binMax, nBins, outCsv):
  theBins = np.linspace(binMin, binMax, nBins)
  theHist = np.zeros(len(theBins)-1, dtype='uint64')
  bname = os.path.basename(fn)
  print(bname)
  if not os.path.exists(fn):
    sys.exit('Error: file '+fn+' does not exist')
  
  src = gdal.Open(fn)
  srcBand = src.GetRasterBand(band)
  xSize = src.RasterXSize
  ySize = src.RasterYSize
  blockSize = 512
    
  for y in range(0, ySize, blockSize):
    if y + blockSize < ySize:
      rows = blockSize
    else:
      rows = ySize - y
    for x in range(0, xSize, blockSize):
      if x + blockSize < xSize:
        cols = blockSize
      else:
        cols = xSize - x
      
      arr = srcBand.ReadAsArray(x, y, cols, rows)
      arrMa = np.ma.masked_equal(arr, noData)
      htemp, jnk = np.histogram(arrMa, theBins)
      theHist += np.uint64(htemp)
    
  src = None
  
  df = pd.DataFrame({'Units':theBins[0:-1], 'Frequency':theHist, 'Year':[year]*len(theHist), 'Layer':[name]*len(theHist)})
  df = cdf(df, 'Frequency')
  df.to_csv(outCsv, index=False)

def largeAreaMinMax(fn, noData, band):
  srcVrt = gdal.Open(fn)
  xSize = srcVrt.RasterXSize
  ySize = srcVrt.RasterYSize
  blockSize = 512
  srcVrt = None
    
  minV = 50
  maxV = 50 
  for y in range(0, ySize, blockSize):
    if y + blockSize < ySize:
      rows = blockSize
    else:
      rows = ySize - y
    for x in range(0, xSize, blockSize):
      if x + blockSize < xSize:
        cols = blockSize
      else:
        cols = xSize - x
      
      src = gdal.Open(fn)
      data = src.GetRasterBand(band)
      arr = data.ReadAsArray(x, y, cols, rows)
      arrMa = np.ma.masked_equal(arr, noData)
      arrMaMin = np.ma.min(arrMa)
      arrMaMax = np.ma.max(arrMa)
      minV = arrMaMin if arrMaMin < minV else minV
      maxV = arrMaMax if arrMaMax > maxV else maxV
  return [minV, maxV]

def makeLargeAreaHist(refHist, predHist, units, figSize, outDir, outType):
  #outType = 'zero'
  #refHist = csvRefHist0
  #predHist = csvPredHist0
  #outDir = outDirOther
  
  fnReg = os.path.join(outDir, 'hist1d_all_reg_'+outType+'.pdf')
  fnCDF = os.path.join(outDir, 'hist1d_all_cdf_'+outType+'.pdf')
  
  dataHistDf = pd.read_csv(refHist)
  dataHistDf = dataHistDf.append(pd.read_csv(predHist))
  dataHistDf = dataHistDf.reset_index(drop=True)
    
  # plot a normal histogram
  fig, ax = plt.subplots(figsize=(figSize, figSize))
  sns.lineplot(ax=ax, data=dataHistDf, x='Units', y="Frequency", hue="Layer")
  ax.set(xlabel=units, ylabel='Frequency')
  handles, labels = ax.get_legend_handles_labels()
  ax.legend(loc='upper right', handles=handles[1:], labels=labels[1:])
  plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
  plt.tight_layout()
  plt.savefig(fnReg)
  
  # plot the cumulative density function (CDF)
  fig, ax = plt.subplots(figsize=(figSize,figSize))
  sns.lineplot(ax=ax, data=dataHistDf, x="Units", y="CDF", hue="Layer")
  plt.ylim(0, 1)
  ax.set(xlabel=units, ylabel='Cumulative Density')
  handles, labels = ax.get_legend_handles_labels()
  ax.legend(loc='lower right', handles=handles[1:], labels=labels[1:])
  plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
  plt.tight_layout()
  plt.savefig(fnCDF)


##########################################################################
# biomass hexagon
predF = '/vol/v3/lt_stem_v3.1/models/biomassfiaald_20180708_0859/2000/biomassfiaald_20180708_0859_2000_mean.tif'
refF = '/vol/v2/datasets/biomass/nbcd/fia_ald/nbcd_fia_ald_biomass_clipped_to_conus.tif'
shpF = '/vol/v3/lt_stem_v3.1/evaluation/vector/conus_mapzones_1to5million_epsg5070.shp'  #_tester
refND = -32768
predND = -9999
trgField = 'ZONE_NUM'
outDir = '/vol/v3/lt_stem_v3.1/evaluation/ltc_biomass_nbcd/ltc_biomass_nbcd_v01'
xyLim = (500, 500)
xLab = 'Reference (Mg/ha)'
yLab = 'Prediction (Mg/ha)'
units = 'Mg/ha'
annoXY = (15,375)


##########################################################################


# open the shapefile	
vDriver = ogr.GetDriverByName("ESRI Shapefile")
vSrc = vDriver.Open(shpF, 0)
vLayer = vSrc.GetLayer()

commonBox = get_intersec([predF, refF])


# make output dirs for figures
#/vol/v3/lt_stem_v3.1/evaluation/ltc_biomass_nbcd/ltc_biomass_nbcd_v01/figs/zones
outDirOther = os.path.join(outDir, 'other')
outDirFigs = os.path.join(outDir, 'figs')
outDirCdf = os.path.join(outDirFigs, 'cdf')
outDirHist1d = os.path.join(outDirFigs, 'hist1d')
outDirHist2d = os.path.join(outDirFigs, 'hist2d')
outDirMap = os.path.join(outDirFigs, 'map')
for d in [outDirOther, outDirCdf, outDirHist1d, outDirHist2d, outDirMap]:
  if not os.path.isdir(d):
    os.makedirs(d)

# change plt params
params = {'legend.fontsize': 9,
       'axes.labelsize': 9,
       'axes.titlesize':9,
       'xtick.labelsize':9,
       'ytick.labelsize':9}
plt.rcParams.update(params)


for f in range(vLayer.GetFeatureCount()):
  #f = 37
  #print(f)
  # what feature are we working with
  feature = vLayer[f]
  name = feature.GetField(trgField)
  print('f: '+str(name))
  
  # get the pixels for the zone
  predP0, coords = get_zone_pixels(feature, shpF, predF, 1, [commonBox[0], commonBox[2], commonBox[3], commonBox[1]])#.compressed() [commonBox[0], commonBox[2], commonBox[3], commonBox[1]]
  refP0, coords = get_zone_pixels(feature, shpF, refF, 1, [coords[0], coords[1], coords[2], coords[3]])#.compressed()
 
  predP, coords = get_zone_pixels(feature, shpF, predF, 1, [commonBox[0], commonBox[2], commonBox[3], commonBox[1]])#.compressed() [commonBox[0], commonBox[2], commonBox[3], commonBox[1]]
  refP, coords = get_zone_pixels(feature, shpF, refF, 1, [coords[0], coords[1], coords[2], coords[3]])#.compressed()
  
  
  # mask and prep the pixles
  predP0 = ma.masked_equal(predP0, predND)
  refP0 = ma.masked_equal(refP0, refND)

  predP = ma.masked_equal(predP, predND)
  refP = ma.masked_equal(refP, refND)
  refP = ma.masked_equal(refP, 0) # this would remove 0 from the data - I think we need to keep this
  
  combMaskZero = np.logical_not(np.logical_not(predP0.mask) * np.logical_not(refP0.mask))
  combMaskNoZero = np.logical_not(np.logical_not(predP.mask) * np.logical_not(refP.mask))

  predP0[combMaskZero] = ma.masked
  refP0[combMaskZero] = ma.masked
  
  predP[combMaskNoZero] = ma.masked
  refP[combMaskNoZero] = ma.masked
  
  predP0 = predP0.compressed()
  refP0 = refP0.compressed()
  predP = predP.compressed()
  refP = refP.compressed()
  
  if (predP0.shape[0] == 0) | (refP0.shape[0] == 0) | (predP0==0).all() | (refP0==0).all():
    predP0 = np.array([0,0,1,1], dtype='float64')
    refP0 = np.array([0,0,1,1], dtype='float64')

  if (predP.shape[0] == 0) | (refP.shape[0] == 0) | (predP==0).all() | (refP==0).all():
    predP = np.array([0,0,1,1], dtype='float64')
    refP = np.array([0,0,1,1], dtype='float64')

  """  
  # optionally sample the pixels
  sampFrac = 1
  totPixs = trainP.shape[0]
  sampSize = round(totPixs*sampFrac) # todo - this is hard coded to 100 percent - add a parameter to sample a given fraction of the pixels
  if sampFrac != 1:
    sampIndex = np.random.choice(range(sampSize), size=sampSize)
    trainP = trainP[sampIndex]
    predP[sampIndex]
  """
  
  # get stats
  maeM0 = round(np.median(np.absolute(np.subtract(predP0, refP0))),1)
  maeU0 = round(np.mean(np.absolute(np.subtract(predP0, refP0))),1)
  rmseM0 = round(np.sqrt(np.median((predP0-refP0)**2)),1)
  rmseU0 = round(np.sqrt(np.mean((predP0-refP0)**2)),1)
  trainQuan0 = np.quantile(refP0, [0.25, 0.5, 0.75])
  predQuan0 = np.quantile(predP0, [0.25, 0.5, 0.75])
  r0 = round(np.corrcoef(refP0, predP0)[0][1], 2)
  if (maeM0 == 0) & (r0 == 1):
    r0 = 0.0
    
  maeM = round(np.median(np.absolute(np.subtract(predP, refP))),1)
  maeU = round(np.mean(np.absolute(np.subtract(predP, refP))),1)
  rmseM = round(np.sqrt(np.median((predP-refP)**2)),1)
  rmseU = round(np.sqrt(np.mean((predP-refP)**2)),1)
  trainQuan = np.quantile(refP, [0.25, 0.5, 0.75])
  predQuan = np.quantile(predP, [0.25, 0.5, 0.75])
  r = round(np.corrcoef(refP, predP)[0][1], 2)
  if (maeM == 0) & (r == 1):
    r = 0.0
  
  # fill in a dataframe with the stats
  df = pd.DataFrame({
      trgField:name, 
      'r':r, 
      'rmseM':rmseM, 
      'rmseU':rmseU, 
      'maeM':maeM, 
      'maeU':maeU,
      't1q':int(trainQuan[0]),
      't2q':int(trainQuan[1]),
      't3q':int(trainQuan[2]),
      'p1q':int(predQuan[0]),
      'p2q':int(predQuan[1]),
      'p3q':int(predQuan[2]),      
      'r0':r0, 
      'rmseM0':rmseM0,
      'rmseU0':rmseU0,
      'maeM0':maeM0, 
      'maeU0':maeU0, 
      't1q0':int(trainQuan0[0]),
      't2q0':int(trainQuan0[1]),
      't3q0':int(trainQuan0[2]),
      'p1q0':int(predQuan0[0]),
      'p2q0':int(predQuan0[1]),
      'p3q0':int(predQuan0[2])}, index=[f])
  
  if f==0:
    allDf = df
  else:
    allDf = allDf.append(df)
 
  
  # make a histogram df
  histDf0 = makeHist(refP0, predP0, minv=0)
  histDf = makeHist(refP, predP, minv=1)

  
  figSize = 3.15
  # make fig file names
  cdfFn = os.path.join(outDirCdf, 'cdf_nozero_zone_'+str(name)+'.pdf')
  hist1dFn = os.path.join(outDirHist1d, 'hist1d_nozero_zone_'+str(name)+'.pdf')
  hist2dFn = os.path.join(outDirHist2d, 'hist2d_nozero_zone_'+str(name)+'.pdf')
  cdfFn0 = os.path.join(outDirCdf, 'cdf_zero_zone_'+str(name)+'.pdf')
  hist1dFn0 = os.path.join(outDirHist1d, 'hist1d_zero_zone_'+str(name)+'.pdf')
  hist2dFn0 = os.path.join(outDirHist2d, 'hist2d_zero_zone_'+str(name)+'.pdf')
  mapFn = os.path.join(outDirMap, 'map_zone_'+str(name)+'.png')
  #csvFn = os.path.join(outDirOther, 'stats_df.csv')



  # make 2d hist plot 
  make2dHist(refP0, predP0, name, xyLim, xLab, yLab, r0, maeM0, maeU0, annoXY, figSize, hist2dFn0)
  make2dHist(refP, predP, name, xyLim, xLab, yLab, r, maeM, maeU, annoXY, figSize, hist2dFn)
  
   
  # plot a normal histogram
  make1dHist(histDf0, units, figSize, hist1dFn0)
  make1dHist(histDf, units, figSize, hist1dFn)
  
  
  # plot the cumulative density function (CDF)
  makeCDF(histDf0, units, figSize, cdfFn0)
  makeCDF(histDf, units, figSize, cdfFn)



  # plot the map
  #shpFtest = '/vol/v3/lt_stem_v3.1/evaluation/vector/conus_mapzones_1to5million_epsg5070.shp'
  if f == 0:
    shpMap = geopandas.read_file(shpF) #shpF
    shpMap = shpMap.to_crs({'init': 'epsg:5070'})

  color = (shpMap[trgField] == name)
  shpMap['color'] = color
  
  shpMapP = shpMap.plot(column='color', cmap='tab20c_r', edgecolor='grey', linewidth=0.5, figsize=(6.5, 4.3)); #color=white terrain tab20c
  #shpMapP.axes.get_xaxis().set_visible(False)
  #shpMapP.axes.get_yaxis().set_visible(False)
  plt.axis('off')
  #plt.tight_layout(pad=0)
  mapFn = os.path.join(outDirMap, 'map_zone_'+str(name)+'.png')
  plt.savefig(mapFn, dpi=150)
  mapImg = Image.open(mapFn)
  mapImg.convert('L').save(mapFn)

# write out the table after all features are included
#allDf.to_csv(csvFn, ',', index=False)

# work on full area summaries
csvStatsFn = os.path.join(outDirOther, 'stats_summary_df.csv')
csvRefHist0 = os.path.join(outDirOther, 'all_hist_ref_df0.csv')
csvPredHist0 = os.path.join(outDirOther, 'all_hist_pred_df0.csv')
csvRefHist = os.path.join(outDirOther, 'all_hist_ref_df.csv')
csvPredHist = os.path.join(outDirOther, 'all_hist_pred_df.csv')

# add color
allDf = addColorToDf(allDf, 'r0', 'rClr')
allDf = addColorToDf(allDf, 'r', 'rClr')
allDf = addColorToDf(allDf, 'rmseM0', 'rmseM0Clr')
allDf = addColorToDf(allDf, 'rmseM', 'rmseMClr')
allDf = addColorToDf(allDf, 'rmseU0', 'rmseU0Clr')
allDf = addColorToDf(allDf, 'rmseU', 'rmseUClr')
allDf = addColorToDf(allDf, 'maeM0', 'maeM0Clr')
allDf = addColorToDf(allDf, 'maeM', 'maeMClr')
allDf = addColorToDf(allDf, 'maeU0', 'maeU0Clr')
allDf = addColorToDf(allDf, 'maeU', 'maeUClr')

shpMapStats = shpMap.merge(allDf, on=trgField)

#allDf = pd.read_csv('/vol/v3/lt_stem_v3.1/evaluation/ltc_biomass_nbcd/ltc_biomass_nbcd_v01/other/stats_df.csv')

 
# make stat maps 
fields = ['r0', 'r', 'rmseM0', 'rmseM', 'rmseU0', 'rmseU', 'maeM0', 'maeM', 'maeU0', 'maeU', 'p2q', 't2q', 'p2q0', 't2q0']  
minMaxs = [
  getMinMax(shpMapStats['r0'], shpMapStats['r']), getMinMax(shpMapStats['r0'], shpMapStats['r']),
  getMinMax(shpMapStats['rmseM0'], shpMapStats['rmseM']), getMinMax(shpMapStats['rmseM0'], shpMapStats['rmseM']),
  getMinMax(shpMapStats['rmseU0'], shpMapStats['rmseU']), getMinMax(shpMapStats['rmseU0'], shpMapStats['rmseU']),
  getMinMax(shpMapStats['maeM0'], shpMapStats['maeM']), getMinMax(shpMapStats['maeM0'], shpMapStats['maeM']),
  getMinMax(shpMapStats['maeU0'], shpMapStats['maeU']), getMinMax(shpMapStats['maeU0'], shpMapStats['maeU']),
  getMinMax(shpMapStats['p2q'], shpMapStats['t2q']), getMinMax(shpMapStats['p2q'], shpMapStats['t2q']),
  getMinMax(shpMapStats['p2q0'], shpMapStats['t2q0']), getMinMax(shpMapStats['p2q0'], shpMapStats['t2q0']),
]

for field, minMax in zip(fields, minMaxs):
  makeStatMaps(shpMapStats, field, minMax, outDirOther, field+'_stat_map.png')


# summarize the stats
fields = ['r0', 'r', 'rmseM0', 'rmseM', 'rmseU0', 'rmseU', 'maeM0', 'maeM', 'maeU0', 'maeU']
xlabs = [ 'r' , 'r', 'RMSE'  , 'RMSE' , 'RMSE'  , 'RMSE' , 'MAE'  , 'MAE' , u'\u03BCAE'  , u'\u03BCAE']  
statSummary = []
for field, xlab in zip(fields, xlabs):
  statMedian = makeStatHist(shpMapStats, field, xlab, outDirOther, field+'_stat_hist.pdf', figSize)
  tempDf = pd.DataFrame({'statName': [field], 'statVal': [statMedian]})
  statSummary.append(tempDf)
statSummary = pd.concat(statSummary, axis=0)
statSummary.to_csv(csvStatsFn, ',', index=False)



# get the min and max of the ref and pred data
refMinMax = largeAreaMinMax(refF, refND, 1)
predMinMax = largeAreaMinMax(predF, predND, 1)
dataMin = min([refMinMax[0],predMinMax[0]]) 
dataMax = min([refMinMax[1],predMinMax[1]])+1
dataNbins = (dataMax-dataMin)+1 

# make a histograms
makeHistDf(refF, 1, refND, 'NA', 'Reference', dataMin, dataMax, dataNbins, csvRefHist0)
makeHistDf(predF, 1, predND, 'NA', 'Prediction', dataMin, dataMax, dataNbins, csvPredHist0)
makeHistDf(refF, 1, refND, 'NA', 'Reference', 1, dataMax, dataNbins-1, csvRefHist)
makeHistDf(predF, 1, predND, 'NA', 'Prediction', 1, dataMax, dataNbins-1, csvPredHist)

makeLargeAreaHist(csvRefHist0, csvPredHist0, units, figSize, outDirOther, 'zero')
makeLargeAreaHist(csvRefHist, csvPredHist, units, figSize, outDirOther, 'nozero')






 
