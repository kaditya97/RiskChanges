#!/usr/bin/env python
# coding: utf-8

# In[2]:


from osgeo import gdal
import numpy
from osgeo import ogr
from osgeo import osr
import os
import sys
import pandas as pd
import numpy as np
import copy
#from sqlalchemy import *
import pandas as pd
import re
from sqlalchemy import create_engine


# In[21]:


def zonalPoly(lyr, input_value_raster, Ear_Table_PK, exposure_id, agg_col=None):
    tempDict = {}
    print(agg_col, 'agg_col')
    featlist = range(lyr.GetFeatureCount())
    raster = gdal.Open(input_value_raster)

    projras = osr.SpatialReference(wkt=raster.GetProjection())
    print(projras, lyr.__dict__, 'porjras')

    # projear = lyr.GetSpatialRef()
    # print(projear, 'projear')

    # epsgear = projear.GetAttrValue('AUTHORITY', 1)
    # # print(epsgear,epsgras)
    # if not epsgras == epsgear:
    #     toEPSG = "EPSG:"+str(epsgear)
    #     output_raster = input_value_raster.replace(".tif", "_projected.tif")
    #     gdal.Warp(output_raster, input_value_raster, dstSRS=toEPSG)
    #     raster = None
    #     raster = gdal.Open(output_raster)
    # else:
    #     pass
    # Get raster georeference info
    raster_srs = osr.SpatialReference()
    raster_srs.ImportFromWkt(raster.GetProjectionRef())
    gt = raster.GetGeoTransform()
    xOrigin = gt[0]
    yOrigin = gt[3]
    pixelWidth = gt[1]
    pixelHeight = gt[5]
    rb = raster.GetRasterBand(1)

    df = pd.DataFrame()
    add_1 = True

    try:
        a = lyr.GetFeature(0)

    except:
        add_1 = True

    print(featlist, add_1, '65')

    for FID in featlist:
        if add_1:
            FID += 1

        feat = lyr.GetFeature(FID)
        geom = feat.GetGeometryRef()
        area = geom.GetArea()
        extent = geom.GetEnvelope()
        xmin = extent[0]
        xmax = extent[1]
        ymin = extent[2]
        ymax = extent[3]
        print(FID)

        xoff = int((xmin - xOrigin)/pixelWidth)
        yoff = int((yOrigin - ymax)/pixelWidth)
        xcount = int((xmax - xmin)/pixelWidth)+1
        ycount = int((ymax - ymin)/pixelWidth)+1

        target_ds = gdal.GetDriverByName('MEM').Create(
            '', xcount, ycount, 1, gdal.GDT_Byte)
        target_ds.SetGeoTransform((
            xmin, pixelWidth, 0,
            ymax, 0, pixelHeight,
        ))
        # Create for target raster the same projection as for the value raster
        target_ds.SetProjection(raster_srs.ExportToWkt())

        gdal.RasterizeLayer(target_ds, [1], lyr, burn_values=[1])

        # Read raster as arrays
        banddataraster = raster.GetRasterBand(1)
        dataraster = banddataraster.ReadAsArray(
            xoff, yoff, xcount, ycount).astype(numpy.float)

        bandmask = target_ds.GetRasterBand(1)
        datamask = bandmask.ReadAsArray(
            0, 0, xcount, ycount).astype(numpy.float)
        # Mask zone of raster
        zoneraster = numpy.ma.masked_array(
            dataraster,  numpy.logical_not(datamask))
       # print(zoneraster)
        (unique, counts) = numpy.unique(zoneraster, return_counts=True)
        unique[unique.mask] = 9999
        if 9999 in unique:
            falsedata = numpy.where(unique == 9999)[0][0]
            ids = numpy.delete(unique, falsedata)
            cus = numpy.delete(counts, falsedata)
        else:
            ids = unique
            cus = counts
        # print(ids)
        frequencies = numpy.asarray((ids, cus)).T
        len_ras = zoneraster.count()
        for i in range(len(frequencies)):
            frequencies[i][1] = (frequencies[i][1]/len_ras)*100

        df_temp = pd.DataFrame(frequencies, columns=['class', 'exposed'])
        df_temp['geom_id'] = feat[Ear_Table_PK]
        df_temp['exposure_id'] = exposure_id
        df_temp['areaOrLen'] = area
        if agg_col is not None:
            df_temp['admin_unit'] = feat[agg_col]
            print(feat[agg_col], 'featagg')
        df = df.append(df_temp, ignore_index=True)

    raster = None
    print(df, 'df')
    return df


# In[4]:


# lead co-ordinate system should be EAR
def zonalPoint(lyr, input_value_raster, exposure_id, Ear_Table_PK, agg_col):
    tempDict = {}
    featlist = range(lyr.GetFeatureCount())
    raster = gdal.Open(input_value_raster)

    projras = osr.SpatialReference(wkt=raster.GetProjection())
    epsgras = projras.GetAttrValue('AUTHORITY', 1)

    projear = lyr.GetSpatialRef()
    epsgear = projear.GetAttrValue('AUTHORITY', 1)
    print(epsgear, epsgras)
    if not epsgras == epsgear:
        toEPSG = "EPSG:"+str(epsgear)
        output_raster = input_value_raster.replace(".tif", "_projected.tif")
        gdal.Warp(output_raster, input_value_raster, dstSRS=toEPSG)
        raster = None
        raster = gdal.Open(output_raster)
    else:
        pass
    # Get raster georeference info

    gt = raster.GetGeoTransform()
    xOrigin = gt[0]
    yOrigin = gt[3]
    pixelWidth = gt[1]
    pixelHeight = gt[5]
    rb = raster.GetRasterBand(1)

    df = pd.DataFrame()
    add_1 = False

    try:
        a = lyr.GetFeature(0)

    except:
        add_1 = True

    for FID in featlist:
        if add_1:
            FID += 1
        feat = lyr.GetFeature(FID)
        geom = feat.GetGeometryRef()
        mx, my = geom.GetX(), geom.GetY()
        px = int((mx - gt[0]) / gt[1])
        py = int((my - gt[3]) / gt[5])
        intval = rb.ReadAsArray(px, py, 1, 1)
        df_temp = pd.DataFrame([[intval[0][0]]], columns=['class'])
        df_temp['geom_id'] = FID
        df_temp['exposure_id'] = exposure_id
        #df_temp['class'] = intval[0][0]
        df_temp['exposed'] = 100
        print(df_temp)
        if agg_col is not None:
            df_temp['admin_unit'] = feat[agg_col]
        df = df.append(df_temp, ignore_index=True)
    raster = None
    return df


# In[5]:


def zonalLine(lyr, input_value_raster, Ear_Table_PK, exposure_id, agg_col=None):
    tempDict = {}
    featlist = range(lyr.GetFeatureCount())
    raster = gdal.Open(input_value_raster)

    projras = osr.SpatialReference(wkt=raster.GetProjection())
    epsgras = projras.GetAttrValue('AUTHORITY', 1)

    projear = lyr.GetSpatialRef()
    epsgear = projear.GetAttrValue('AUTHORITY', 1)
    # print(epsgear,epsgras)
    if not epsgras == epsgear:
        toEPSG = "EPSG:"+str(epsgear)
        output_raster = input_value_raster.replace(".tif", "_projected.tif")
        gdal.Warp(output_raster, input_value_raster, dstSRS=toEPSG)
        raster = None
        raster = gdal.Open(output_raster)
    else:
        pass
    # Get raster georeference info
    raster_srs = osr.SpatialReference()
    raster_srs.ImportFromWkt(raster.GetProjectionRef())
    gt = raster.GetGeoTransform()
    xOrigin = gt[0]
    yOrigin = gt[3]
    pixelWidth = gt[1]
    pixelHeight = gt[5]
    rb = raster.GetRasterBand(1)

    df = pd.DataFrame()
    add_1 = False
    try:
        a = lyr.GetFeature(0)

    except:
        add_1 = True

    for FID in featlist:
        if add_1:
            FID += 1
        feat = lyr.GetFeature(FID)
        geom = feat.GetGeometryRef()
        length = geom.Length()
        extent = geom.GetEnvelope()
        xmin = extent[0]
        xmax = extent[1]
        ymin = extent[2]
        ymax = extent[3]

        xoff = int((xmin - xOrigin)/pixelWidth)
        yoff = int((yOrigin - ymax)/pixelWidth)
        xcount = int((xmax - xmin)/pixelWidth)+1
        ycount = int((ymax - ymin)/pixelWidth)+1

        target_ds = gdal.GetDriverByName('MEM').Create(
            '', xcount, ycount, 1, gdal.GDT_Byte)
        target_ds.SetGeoTransform((
            xmin, pixelWidth, 0,
            ymax, 0, pixelHeight,
        ))
        # Create for target raster the same projection as for the value raster
        target_ds.SetProjection(raster_srs.ExportToWkt())

        gdal.RasterizeLayer(target_ds, [1], lyr, burn_values=[1])

        # Read raster as arrays
        banddataraster = raster.GetRasterBand(1)
        dataraster = banddataraster.ReadAsArray(
            xoff, yoff, xcount, ycount).astype(numpy.float)

        bandmask = target_ds.GetRasterBand(1)
        datamask = bandmask.ReadAsArray(
            0, 0, xcount, ycount).astype(numpy.float)
        # Mask zone of raster
        zoneraster = numpy.ma.masked_array(
            dataraster,  numpy.logical_not(datamask))
       # print(zoneraster)
        (unique, counts) = numpy.unique(zoneraster, return_counts=True)
        unique[unique.mask] = 9999
        if 9999 in unique:
            falsedata = numpy.where(unique == 9999)[0][0]
            ids = numpy.delete(unique, falsedata)
            cus = numpy.delete(counts, falsedata)
        else:
            ids = unique
            cus = counts
        # print(ids)
        frequencies = numpy.asarray((ids, cus)).T
        len_ras = zoneraster.count()
        for i in range(len(frequencies)):
            frequencies[i][1] = (frequencies[i][1]/len_ras)*100

        df_temp = pd.DataFrame(frequencies, columns=['class', 'exposed'])
        df_temp['geom_id'] = feat[Ear_Table_PK]
        df_temp['exposure_id'] = exposure_id
        df_temp['areaOrLen'] = length
        if agg_col is not None:
            df_temp['admin_unit'] = feat[agg_col]
        #df_temp['admin_unit'] =feat[agg_col]
        df = df.append(df_temp, ignore_index=True)
    raster = None
    return df


# In[6]:


def ExposurePgAg(input_zone, admin_unit, agg_col, input_value_raster, connString, exposure_id, Ear_Table_PK, schema):
    conn = ogr.Open(connString)
    sql = '''
    SELECT  "{0}".*, "{1}"."{2}" FROM "{3}"."{0}", "{3}"."{1}" WHERE ST_Intersects("{0}".geom, "{1}".geom)
    '''.format(input_zone, admin_unit, agg_col, schema)

    print(sql)
    lyr = conn.ExecuteSQL(sql)

    featList = range(lyr.GetFeatureCount())

    print(featList, 'flist')
    statDict = {}
    feat = lyr.GetNextFeature()

    geom = feat.GetGeometryRef()
    print(geom, 'geom', geom.GetArea())
    geometrytype = geom.GetGeometryName()
    print(geometrytype)
    i = 1
    df = pd.DataFrame()
    if (geometrytype == 'POLYGON' or geometrytype == 'MULTIPOLYGON'):
        df = zonalPoly(lyr, input_value_raster,
                       Ear_Table_PK, exposure_id, agg_col=agg_col)
        return df
    elif(geometrytype == 'POINT' or geometrytype == 'MULTIPOINT'):
        return zonalPoint(lyr, input_value_raster, exposure_id, Ear_Table_PK, agg_col)
    elif(geometrytype == 'LINESTRING' or geometrytype == 'MULTILINESTRING'):
        df = zonalLine(lyr, input_value_raster,
                       Ear_Table_PK, exposure_id, agg_col)
        return df


# In[ ]:


# In[7]:


def ExposurePG(input_zone, input_value_raster, connString, exposure_id, Ear_Table_PK):
    # print(input_zone)
    conn = ogr.Open(connString)
    lyr = conn.GetLayer(input_zone)
    featList = range(lyr.GetFeatureCount())
    # print(featList)
    statDict = {}
    feat = lyr.GetNextFeature()
    # print(lyr.GetLayerDefn())
    geom = feat.GetGeometryRef()
    geometrytype = geom.GetGeometryName()
    # print(geometrytype)
    i = 1
    df = pd.DataFrame()
    if (geometrytype == 'POLYGON' or geometrytype == 'MULTIPOLYGON'):
        df = zonalPoly(lyr, input_value_raster, Ear_Table_PK,
                       exposure_id, agg_col=None)
        return df
    elif(geometrytype == 'POINT' or geometrytype == 'MULTIPOINT'):
        return zonalPoint(lyr, input_value_raster, exposure_id, Ear_Table_PK, agg_col=None)
    elif(geometrytype == 'LINESTRING' or geometrytype == 'MULTILINESTRING'):
        df = zonalLine(lyr, input_value_raster, Ear_Table_PK,
                       exposure_id, agg_col=None)
        return df


# In[8]:


def reclassify(in_image, out_image, out_dir, classification):

    driver = gdal.GetDriverByName('GTiff')
    file = gdal.Open(in_image)
    band = file.GetRasterBand(1)
    lista = band.ReadAsArray()
    classess = classification.keys()
    reclass = copy.copy(lista)
    for key in classess:
        reclass[np.where(lista < classification[key][1])
                ] = classification[key][0]
    # reclassification
    # for j in  range(file.RasterXSize):
     #   for i in  range(file.RasterYSize):
      #      for key in classess:
       #         if lista[i,j] < classification[key][1]:
        #            lista[i,j] = classification[key][0]
        #           break

    export = out_dir+"//"+out_image+".tif"
    # create new file
    file2 = driver.Create(export, file.RasterXSize, file.RasterYSize, 1)
    file2.GetRasterBand(1).WriteArray(lista)

    # spatial ref system
    proj = file.GetProjection()
    georef = file.GetGeoTransform()
    file2.SetProjection(proj)
    file2.SetGeoTransform(georef)
    file2.FlushCache()


# In[9]:


def aggregate(df, agg_col):
    try:
        df['exposed_areaOrLen'] = df['exposed']*df['areaOrLen']/100
        df_aggregated = df.groupby([agg_col, 'class'], as_index=False).agg(
            {'exposed_areaOrLen': 'sum', 'exposed': 'count'})
    except:
        df_aggregated = df.groupby(['admin_unit', 'class'], as_index=False).agg({
            'exposed': 'count', 'exposure_id': 'mean'})
    return df_aggregated


# In[10]:


def todatabase(df, connstr, table_name, schema):
  # Creating SQLAlchemy's engine to use
    engine = create_engine(connstr)

    # ... [do something with the geodataframe]

   # Use 'dtype' to specify column's type
    # For the geom column, we will use GeoAlchemy's type 'Geometry'
    try:
        df.to_sql(table_name, engine, schema, if_exists='append', index=False)
    except:
        return("error, trying to append in non related table,please store in same table as EAR")
    engine.dispose()


# In[11]:


def CalculateExposure(Ear_Table, Ear_Table_PK, haz_dir, connString,
                      connSchema, exposure_id, exposure_table,
                      admin_unit=None, agg_col=None, aggregation=False):
    a = re.split(':|//|/|@', connString)
    print(a, 'a')
    databaseServer = a[4]
    databaseName = a[6]
    databaseUser = a[2]
    databasePW = a[3]
    connStringOGR = "PG: host=%s dbname=%s user=%s password=%s  schemas=%s" % (
        databaseServer, databaseName, databaseUser, databasePW, connSchema)
    input_zone = Ear_Table
    input_value_raster = haz_dir

    if aggregation is True:
        if admin_unit is not None:
            df = ExposurePgAg(input_zone, admin_unit, agg_col,
                              input_value_raster, connStringOGR, exposure_id, Ear_Table_PK, connSchema)
            df_aggregated = aggregate(df, agg_col)
            todatabase(df, connString, exposure_table, connSchema)
            todatabase(df_aggregated, connString,
                       exposure_table+'_agg', connSchema)
            return ('completed')
        elif admin_unit is None:
            return ("Please provide the admin unit and aggregation column for the aggregation")
    else:
        df = ExposurePG(input_zone, input_value_raster,
                        connStringOGR, exposure_id, Ear_Table_PK)
        todatabase(df, connString, exposure_table, connSchema)
        return ('completed')


# '''    databaseServer = "127.0.0.1"
#    databaseName = "postgres"
#    databaseUser = "ashok"
#    databasePW = "gicait123"
#    connString = "PG: host=%s dbname=%s user=%s password=%s" % (databaseServer,databaseName,databaseUser,databasePW)
# '''


# TEST

# In[13]:


# In[16]:
databaseServer = "203.159.29.45"
databaseName = "sdssv2"
databaseUser = "postgres"
databasePW = "gicait123"
connString = "postgresql://postgres:gicait123@203.159.29.45:5432/sdssv2"

print(connString)

input_zone = "bu78EPK"
admin_unit = "panjakent_suFC6QW"
agg_col = 'id'
input_value_raster = r"C:\Users\tek\Downloads\Landslide_map.tif"
exposure_id = 800
CalculateExposure(input_zone, 'id', input_value_raster, connString,
                  'itc', exposure_id, 'demo_exp',
                  admin_unit=admin_unit, agg_col='id', aggregation=True)

# %%
