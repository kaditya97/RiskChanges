import geopandas as gpd
from . import readmeta
import psycopg2
import pandas as pd 


def readear(connstr,earid):
    engine=psycopg2.connect(connstr)
    metatable=readmeta.earmeta(connstr,earid)
    eartablename=metatable.layer_name[0]
    schema=metatable.workspace[0]
    
    sql=f'SELECT * FROM {schema}."{eartablename}";'
    print(sql)
    ear_table=gpd.read_postgis(sql,con=engine)
    engine.close()
    return ear_table


def readexposure(connstr,exposureid,schema):
    engine=psycopg2.connect(connstr)
    metatable=readmeta.exposuremeta(connstr,exposureid)
    exposuretable=metatable.exposure_table[0]
    sql=f'SELECT * FROM {schema}."{exposuretable}" as exposuretable WHERE exposuretable.exposure_id={exposureid};'
    #print(sql)
    exposure_table=pd.read_sql(sql,con=engine)
    engine.close()
    return exposure_table

def prepareExposureForLoss(connstr,exposureid):
    metadict=readmeta.computeloss_meta(connstr,exposureid)
    schema=metadict["Schema"]
    earid=metadict["earID"]
    pk=metadict["earPK"]
    exposuredata=readexposure(connstr,exposureid,schema)
    eardatageom=readear(connstr,earid)
    eardata=pd.DataFrame(eardatageom.drop(columns='geom'))
    exposure_all=pd.merge(left=exposuredata, right=eardata, how='left', left_on=['geom_id'], right_on=[pk])
    assert not exposure_all.empty , f"The exposure data  {exposureid} returned empty from database"
    #exposure=readvulnerability.linkvulnerability(connstr,exposure_all)
    return exposure_all

def readLoss(connstr,lossid):
    #write more on loss reading
    engine=psycopg2.connect(connstr)
    metatable=readmeta.readLossMeta(connstr,lossid)
    losstable=metatable.loss_table[0]
    schema=metatable.workspace[0] ##TEK add workspace in loss index

    sql=f'SELECT * FROM {schema}."{losstable}" as losstable WHERE losstable.loss_id={lossid};'
    #print(sql)
    loss_table=pd.read_sql(sql,con=engine)
    engine.close()
    return loss_table

def readRiskGeometry(connstr,riskid):
    engine=psycopg2.connect(connstr)
    meta=readmeta.getRiskMeta(connstr,riskid)
    earid=meta.earid[0]
    schema=meta.workspace[0]
    risktable=meta.risk_table[0]
    sql=f'SELECT * FROM {schema}."{risktable}" as risktable WHERE risktable.loss_id={riskid};'
    #print(sql)
    risk_table=pd.read_sql(sql,con=engine)
    ear=readear(connstr,earid)
    engine.close()
    risk_table=pd.merge(left=risk_table, right=ear[['id','geom']], how='left', left_on=['Unit_ID'], right_on=['id'],right_index=False)
    risk_table=gpd.GeoDataFrame(risk_table,crs=ear.crs,geometry='geom')
    return risk_table

def readAdmin(connstr,adminunit):
    meta=readmeta.getAdminMeta(connstr,adminunit)
    schema=meta.workspace[0]
    engine=psycopg2.connect(connstr)
    admintablename=meta.layer_name[0]    
    sql=f'SELECT * FROM {schema}."{admintablename}";'
    #print(sql)
    ear_table=gpd.read_postgis(sql,con=engine)
    engine.close()
    ear_table=ear_table.rename(columns={'id':'ADMIN_ID'})
    return ear_table
    



