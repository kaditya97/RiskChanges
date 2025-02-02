from datetime import MINYEAR
from types import MappingProxyType
import RiskChangesOps.readmeta as readmeta
import RiskChangesOps.readvector as readvector
from sklearn.metrics import auc
import RiskChangesOps.writevector as writevector
import RiskChangesOps.AggregateData as aggregator
import geopandas as gpd
import pandas as pd
def dutch_method(xx,yy):
        #compute risk based on dutch method where xx is value axis and yy is probability axis
        AAL=auc(yy,xx)+(xx[0]*yy[0])
        return AAL

def checkUniqueHazard(con, lossids):
    blank_list=[]
    for lossid in lossids:
        haztype=readmeta.getHazardType(connstr,lossid)
        blank_list.append(haztype)
    assert pd.Series(blank_list).nunique()==1, "Only multiple return periods of single hazard is supported"

def predict_loss(prepared_loss,rps,probs,extensions,hazard):
    Lmin_x=extensions['left'][0][0]
    Lmin_y=extensions['left'][0][1]
    Rmin_x=extensions['right'][0][0]
    Rmin_y=extensions['right'][0][1]
    Lmax_x=extensions['left'][1][0]
    Lmax_y=extensions['left'][1][1]
    Rmax_x=extensions['right'][1][0]
    Rmax_y=extensions['right'][1][1]
    m_left=abs(Lmin_y-Lmax_y)/abs(Lmin_x-Lmax_x)
    m_right=abs(Rmin_y-Rmax_y)/abs(Rmin_x-Rmax_x)
    min_rp=int(1/Lmin_x)
    max_rp=int(1/Rmax_x)
    return_periods=np.linspace(min_rp, max_rp, 6).astype(int)
    return_periods=np.around(return_periods, decimals=-1).tolist() 
    return_periods[0]=min_rp
    return_periods[-1]=max_rp
    low_rp='loss_rp_'+str(rps.min())
    high_rp='loss_rp_'+str(rps.max())
    colname_min='loss_rp_'+str(min_rp)
    colname_max='loss_rp_'+str(max_rp)
    left_dx=abs(probs.min()-min_rp)
    right_dx=abs(probs.max()-min_rp)
    prepared_loss[colname_max]=prepared_loss.apply(lambda row: row[high_rp]+left_dx*m_left, axis=1) #negative because extending in left //negative slope
    prepared_loss[colname_min]=prepared_loss.apply(lambda row: row[low_rp]-right_dx*m_right, axis=1) #negative because extending in right
    # perform interpolation for loss values using following methods
    rps=rps.sort(reverse=True)
    probs=probs.sort()
    probs.insert(0,Lmin_x)
    probs.insert(-1,Rmax_x)
    interp_cols=['loss_rp_'+str(x) for x in rps]
    interp_cols.insert(-1,colname_min)
    interp_cols.insert(0,colname_max)
    interp_probs=probs

    predict_probs=1/np.array(return_periods,dtype=int)
    new_cols=[hazard+'_lrp_'+str(x) for x in return_periods]
    new_cols.append('Unit_ID')
    new_losses=[]
    for index, row in prepared_loss.iterrows():
        y=row[interp_cols].values 
        x=interp_probs
        loss=np.interp(predict_probs, x, y, left=y.max(), right=y.min()).tolist()
        loss.append(row.Unit_ID)
        new_losses.append(loss)
    prepared_loss=pd.DataFrame(new_losses, columns =new_cols)
    return prepared_loss,new_cols[:-1],predict_probs

def getMaxcost(con,lossid):
    meta=readmeta.readLossMeta(con,lossid)
    ear_id=meta['ear_index_id'][0]
    cost_column=meta['ear_cost_column'][0]
    cost_data=readvector.readear(con,ear_id)
    return cost_data,cost_column

def PrepareLossForRisk(con, lossids,extensions,hazard):
    i=True
    rps=[]
    probs=[]
    for id in lossids:
        lossdata=readvector.readLoss(con,id)
        return_period=float(readmeta.getReturnPeriod(con,id))
        colname='loss_rp_'+str(return_period)
        rps.append(return_period)
        probs.append(1.0/return_period)
        lossdata=lossdata.rename(columns={'loss': colname,'geom_id':'Unit_ID'})
        if i:
            prepared_loss=lossdata
            i=False
        else:
            prepared_loss=prepared_loss.merge(lossdata, on='Unit_ID')
    modified_loss,cols, probs=predict_loss(prepared_loss,rps,probs,extensions,hazard)
    cost_data,cost_col=getMaxcost(con,lossids[0])
    prepared_loss=pd.merge(left=modified_loss, right=cost_data[['id',cost_col]], how='outer', left_on=['Unit_ID'], right_on=['id'],right_index=False).rename(columns={cost_col: 'MAX_COST'})
    return modified_loss,cols,probs

def combineLosses(lossess,interacting_hazards,interaction='independent',haz_prob=(1,1)):
    if interaction== 'cascading':
        interaction_function=mul_cascading
    elif interaction== 'compounding':
        interaction_function=mul_compounding
    elif interaction== 'conditional':
        interaction_function=mul_conditional
    elif interaction== 'coupled':  
        interaction_function=mul_coupled
    elif interaction=='independent':  
        interaction_function=mul_independent
    else: 
        raise ValueError("The interaction names do not match the possible hazard interaction. It can be only independent, compounding, cascading, conditional or coupled")
    haz1=interacting_hazards[0]
    haz2=interacting_hazards[1]
    
    loss1=lossess[haz1]['normalized_loss']
    loss2=lossess[haz2]['normalized_loss']
    first_data=True
    #max_cost obtain from database
    cols=[w.replace(haz1, '') for w in lossess[haz1]['cols']]
    new_columns=[]
    for col in cols:
        col1=haz1+col
        col2=haz2+col
        df1=loss1[['UNIT_ID','MAX_COST',col1]].rename(columns={col1: 'loss1'})
        df2=loss2[['UNIT_ID',col2]].rename(columns={col2: 'loss2'})
        df=pd.merge(left=df1, right=df2, how='outer', left_on=['Unit_ID'], right_on=['Unit_ID'],right_index=False) #also merge maximum value
        df['loss1']=haz_prob[0]*df['loss1']
        df['loss2']=haz_prob[1]*df['loss2']
        combined_loss_step=interaction_function(df,haz_prob)
        new_name='combined'+col
        combined_loss_step=combined_loss_step.rename(columns={'combined':new_name})
        new_columns.append(new_name)
        if first_data:
            combined_loss=combined_loss_step
            first_data=False
        else:
            combined_loss=pd.merge(left=combined_loss,right=combined_loss_step,how='outer',left_on=['Unit_ID'], right_on=['Unit_ID'],right_index=False)
    return combined_loss

def mul_independent(loss_data):
    loss_data['combined']=loss_data[['loss1','loss2']].sum(axis=1)
    independent_loss=loss_data[['Unit_ID','combined']]
    return independent_loss

    #return combination of risks

def mul_compounding(loss_data):
    loss_data['combined_nf']=loss_data[['loss1','loss2']].sum(axis=1)
    loss_data['combined']=loss_data[['combined_nf','maxvalue']].min(axis=1)
    compounding_loss=loss_data[['Unit_ID','combined']]
    return compounding_loss
    #return combination of risks

def mul_coupled(loss_data):
    loss_data['combined']=loss_data[['loss1','loss2']].max(axis=1)
    coupled_loss=loss_data[['Unit_ID','combined']]
    return coupled_loss

    #return combination of risks

def mul_cascading(loss_data):
    loss_data['combined_nf']=loss_data[['loss1','loss2']].sum(axis=1)
    loss_data['combined']=loss_data[['combined_nf','maxvalue']].min(axis=1)
    cascading_loss=loss_data[['Unit_ID','combined']]
    return cascading_loss
    #return combination of risks

def mul_conditional(individual_losses,value,triggering,probability_conditional=1):
    loss_data['combined']=loss_data[['loss1','loss2']].sum(axis=1)
    conditional_loss=loss_data[['Unit_ID','combined']]
    return conditional_loss
    #return combination of risks

def dutch_method(xx,yy):
        #compute risk based on dutch method where xx is value axis and yy is probability axis
        AAL=auc(yy,xx)+(xx[0]*yy[0])
        return AAL

def calculateRisk(lossdf,columns, probs):
    risktable=pd.DataFrame(columns=['Unit_ID','AAL'])
    for index, row in lossdf.iterrows():
        xx=row[columns].values.tolist()
        yy=probs
        aal=dutch_method(xx,yy)
        #print('ear',aal)
        ear_id=row['Unit_ID']
        new_row = {'Unit_ID':ear_id, 'AAL':aal}
        #append row to the dataframe
        risktable =risktable.append(new_row, ignore_index=True)
    assert not risktable.empty , f"The Risk calculation failed"
    return risktable


def computeMulRisk(connstr,losscombinations,hazardinteractions,extensions,riskid,**kwargs): #pass dictionaries of loss and hazard combinations. losscombiations can be like this {'flood':[5,6,1],'earthquake':[8,4,9]} interactions can be {'independent':[('eartqhuake','flood'),('earthquake','debrisflow')],'cascading':[('landslide','flood')],'cascading_weights':[(1,0.5)]}. extensions provide information on how to extend the datasets it should be like this {'earthquake':{'left':[(x1,y1),(x2,y2)], 'right':[(a1,b1),(a2,b2)]}, 'flood':{'left':[(x1,y1),(x2,y2)], 'right':[(a1,b1),(a2,b2)]}}
    try:
        is_aggregated=kwargs['is_aggregated']
        onlyaggregated=kwargs['only_aggregated']
        adminid=kwargs['adminunit_id']
    except:
        is_aggregated= False
        onlyaggregated= False
    
    loss_normalization={}
    first=True
    
    for hazard in losscombinations:
        Name_hazard=hazard
        lossids=losscombinations[hazard]
        extention=extensions[hazard]
        normalized_loss,cols,probs=PrepareLossForRisk(connstr, lossids,extention,hazard)
        loss_haz={'normalized_loss':normalized_loss,'cols':cols,'probs':probs}
        loss_normalization[hazard]=loss_haz
    for interaction in hazardinteractions:
        interacting_hazards=hazardinteractions[interaction]
        if (interaction=='cascading') | (interaction=='conditional'):
            try:
                haz_prob=hazardinteractions[interaction+'_weights']
            except:
                raise ValueError('the hazard weights for cascading and conditional interaction in all the hazards are not provided use 1 if it is not available')
        else:
            haz_prob = [(1,1)]
            haz_prob = haz_prob*len(interacting_hazards)

        if len(interacting_hazards)==0:
            continue
        elif len(interacting_hazards)==1:
            combined_loss=combineLosses(loss_normalization,interacting_hazards[0],interaction,haz_prob[0])
        else :
            second=True
            for single_interaction in interacting_hazards:
                index_nm=interacting_hazards.index(single_interaction)
                combined_loss_step=combineLosses(loss_normalization,single_interaction,interaction,haz_prob[index_nm])
                if second:
                    combined_loss=combined_loss_step
                    second=False
                else:
                    cols=combined_loss_step.columns.tolist()
                    cols_selected=[word for word in cols if ('combined' in word)]
                    merged=pd.merge(left=combined_loss,right=combined_loss_step,how='outer',left_on=['Unit_ID'], right_on=['Unit_ID'],right_index=False,suffixes=('_left', '_right'))
                    del combined_loss,combined_loss_step
                    merged_cols=merged.columns.tolist()
                    for col in cols_selected:
                        merged_selected=[word for word in merged_cols if (col in word)]
                        merged[col]=merged[merged_selected].sum(axis=1)
                        merged=merged.drop(columns=[merged_selected])
                    combined_loss=merged
                # do same as in above steps
        if first:
            final_loss=combined_loss
            first=False
        else:
            cols=combined_loss.columns.values.tolist()
            cols_selected=[word for word in cols if ('combined' in word)]
            merged=pd.merge(left=final_loss,right=combined_loss,how='outer',left_on=['Unit_ID'], right_on=['Unit_ID'],right_index=False,suffixes=('_left', '_right'))
            del final_loss,combined_loss
            merged_cols=merged.columns.values.tolist()
            for col in cols_selected:
                merged_selected=[word for word in merged_cols if (col in word)]
                merged[col]=merged[merged_selected].sum(axis=1)
                merged=merged.drop(columns=[merged_selected])
            final_loss=merged
    
    # prepare information for the calculation of risk
    probabilities=loss_haz[probs]
    columns=loss_haz[cols]
    columns=[w.replace(Name_hazard, 'combined') for w in columns]
    risk=calculateRisk(final_loss,columns,probabilities)
    metatable=readmeta.readLossMeta(connstr,lossids[0])
    schema= metatable.workspace[0]
    risk['risk_id']=riskid
    
    if not onlyaggregated:
        writevector.writeRisk(risk,con,schema)
    if is_aggregated:
        admin_unit=readvector.readAdmin(con,adminid)
        ear_id=metatable['ear_index_id'][0]
        ear= readvector.readear(con,ear_id)
        risk=pd.merge(left=risk, right=ear['id','geom'], left_on='Unit_ID',right_on='id',right_index=False)
        risk= gpd.GeoDataFrame(risk,geometry='geom')
        risk=aggregator.aggregaterisk(risk,admin_unit)
        writevector.writeRiskAgg(risk,con,schema)
    #what to do with the over sampling of the computation. other than that already done. 