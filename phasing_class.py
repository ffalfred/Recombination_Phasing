import math
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import operator
pd.options.mode.chained_assignment = None

class phasing_class_tool:

    def __init__(self,BEDfile,min_SNPS):
        """Instants of the files"""
        self.inputfile = pd.read_csv(BEDfile,sep='\t')   #Input
        self.min_snps = min_SNPS

    def limit_snps(self,array_,pos_array):
        pos_pd=pd.DataFrame(pos_array).set_index('Position')
        array_=np.append(array_,np.zeros((1,np.shape(array_)[1])),axis=0)
        for i in  range(np.shape(array_)[1]):
            num_snps=pd.notna(pos_pd.loc[array_[1,i]:array_[2,i]][pos_pd.columns[0]]).sum()
            array_[3,i]=num_snps
        haplo_count=pd.DataFrame(array_.T,columns=['Haplotype','Start','End','NumSNPS'])
        del_haplo_count=haplo_count[(haplo_count['Haplotype']!=0.5)&(haplo_count['NumSNPS']>=self.min_snps)].reset_index()
        deleting_rows=[]
        for index, row in del_haplo_count.iterrows():
            if index!=0:
                if del_haplo_count.loc[index]['Haplotype']==del_haplo_count.loc[index-1]['Haplotype']:
                    del_haplo_count.at[index,'Start']=del_haplo_count.loc[index-1]['Start']
                    del_haplo_count.at[index,'NumSNPS']=del_haplo_count.loc[index-1]['NumSNPS']+del_haplo_count.loc[index]['NumSNPS']
                    deleting_rows.append(index-1)
        del_haplo_count.drop(deleting_rows, inplace=True)
        del_haplo_count=del_haplo_count.reset_index()
        del_haplo_count.drop('level_0',axis=1,inplace=True)
        del_haplo_count.drop('index',axis=1,inplace=True)
        return(del_haplo_count)


    def treat_PB1(self,dict_cells,mark_pb1,i):
        regions=[]
        pos_pd=mark_pb1[['Position',str(i)+'.PB1.Phase']].set_index('Position')
        for d in dict_cells:
            regions.extend(dict_cells[d]['Start'].values)
            regions.extend(dict_cells[d]['End'].values)
        regions=list(set(regions))
        regions.sort()
        new_coord=np.zeros((len(regions)-1,2))
        list_=[]
        for n in range(len(regions)):
            if n+1<len(regions):
                if regions[n]!=regions[n+1]:
                    num_snps=pos_pd.loc[regions[n]:regions[n+1],str(i)+'.PB1.Phase'].count()
                    try:
                        value_0=pos_pd.loc[regions[n]:regions[n+1],str(i)+'.PB1.Phase'].value_counts(normalize=True).loc[0]
                    except KeyError:
                        value_0=0
                    try:
                        value_05=pos_pd.loc[regions[n]:regions[n+1],str(i)+'.PB1.Phase'].value_counts(normalize=True).loc[0.5]
                    except KeyError:
                        value_05=0
                    try:
                        value_1=pos_pd.loc[regions[n]:regions[n+1],str(i)+'.PB1.Phase'].value_counts(normalize=True).loc[1]
                    except KeyError:
                        value_1=0
                    if value_0>0.66:
                        haplo=0
                    elif value_1>0.66:
                        haplo=1
                    else:
                        haplo=0.5
                    list_.append([regions[n],regions[n+1],value_0,value_05,value_1,num_snps,haplo])
        first_haplo=pd.DataFrame(list_,columns=['Start','End','0','0.5','1','SNPS','Haplotype'])
        last_list=[]
        done_rows=[]
        for index,row in first_haplo.iterrows():
            if (row['SNPS']<self.min_snps) & (index not in done_rows):
                if index==0 or row['Haplotype']==first_haplo['Haplotype'][index+1]:
                    value_0=((row['0']*row['SNPS'])+(first_haplo['0'][index+1]*first_haplo['SNPS'][index+1]))/(row['SNPS']+first_haplo['SNPS'][index+1])
                    value_05=((row['0.5']*row['SNPS'])+(first_haplo['0.5'][index+1]*first_haplo['SNPS'][index+1]))/(row['SNPS']+first_haplo['SNPS'][index+1])
                    value_1=((row['1']*row['SNPS'])+(first_haplo['1'][index+1]*first_haplo['SNPS'][index+1]))/(row['SNPS']+first_haplo['SNPS'][index+1])
                    if value_0>0.66:
                        haplo=0
                    elif value_1>0.66:
                        haplo=1
                    else:
                        haplo=0.5
                    last_list.append([row['Start'],first_haplo['End'][index+1],value_0,value_05,value_1,row['SNPS']+first_haplo['SNPS'][index+1],haplo])
                    done_rows.append(index+1)
                elif index==len(first_haplo.index) or row['Haplotype']==first_haplo['Haplotype'][index-1] or first_haplo['Haplotype'][index+1]==first_haplo['Haplotype'][index-1]:
                    value_0=((row['0']*row['SNPS'])+(first_haplo['0'][index-1]*first_haplo['SNPS'][index-1]))/(row['SNPS']+first_haplo['SNPS'][index-1])
                    value_05=((row['0.5']*row['SNPS'])+(first_haplo['0.5'][index+1]*first_haplo['SNPS'][index-1]))/(row['SNPS']+first_haplo['SNPS'][index-1])
                    value_1=((row['1']*row['SNPS'])+(first_haplo['1'][index-1]*first_haplo['SNPS'][index-1]))/(row['SNPS']+first_haplo['SNPS'][index-1])
                    if value_0>0.66:
                        haplo=0
                    elif value_1>0.66:
                        haplo=1
                    else:
                        haplo=0.5
                    last_list.append([first_haplo['Start'][index-1],row['End'],value_0,value_05,value_1,row['SNPS']+first_haplo['SNPS'][index-1],haplo])
                    done_rows.append(index-1)
                else:
                    continue
            elif index not in done_rows:
                last_list.append(row.tolist())
                done_rows.append(index)
        del_haplo_count=pd.DataFrame(last_list,columns=['Start','End','0','0.5','1','SNPS','Haplotype'])
        deleting_rows=[]
        for index, row in del_haplo_count.iterrows():
            if index!=0:
                if del_haplo_count.loc[index]['Haplotype']==del_haplo_count.loc[index-1]['Haplotype']:
                    del_haplo_count.at[index,'Start']=del_haplo_count.loc[index-1]['Start']
                    del_haplo_count.at[index,'SNPS']=del_haplo_count.loc[index-1]['SNPS']+del_haplo_count.loc[index]['SNPS']
                    deleting_rows.append(index-1)
        del_haplo_count.drop(deleting_rows, inplace=True)
        del_haplo_count=del_haplo_count.reset_index() 
        del_haplo_count.drop('index',axis=1,inplace=True)
        return(del_haplo_count[['Haplotype', 'Start', 'End', 'SNPS']])


    def mark_snps(self,HetMatSign,i,cell):
        list_phase=[]
        list_init_phase=[]
        list_end_phase=[]
        chromosome=[]
        previous_phase=float('nan')
        for index, row in HetMatSign[['Chr','Position',str(i)+'.'+str(cell)+'.Phase']].iterrows():
            if math.isnan(row[str(i)+'.'+str(cell)+'.Phase'])==False:
                if math.isnan(previous_phase): 
                    previous_phase=row[str(i)+'.'+str(cell)+'.Phase']
                    list_phase.append(row[str(i)+'.'+str(cell)+'.Phase'])
                    list_init_phase.append(row['Position'])
                    chromosome.append(row['Chr'])
                if row[str(i)+'.'+str(cell)+'.Phase']!=previous_phase:
                    previous_phase=row[str(i)+'.'+str(cell)+'.Phase']
                    list_phase.append(row[str(i)+'.'+str(cell)+'.Phase'])
                    list_init_phase.append(row['Position'])
                    meh=True
                    previous_step=-1
                    while meh:
                        if math.isnan(HetMatSign.iloc[index+previous_step][str(i)+'.'+str(cell)+'.Phase']):
                            previous_step=previous_step-1
                        else:
                            list_end_phase.append(HetMatSign.iloc[index+previous_step]['Position'])
                            meh=False
            if index==HetMatSign[['Chr','Position',str(i)+'.'+str(cell)+'.Phase']].index.values[-1]: 
                meh=True
                previous_step=0
                while meh:
                    if math.isnan(HetMatSign.iloc[index+previous_step][str(i)+'.'+str(cell)+'.Phase']):
                        previous_step=previous_step-1
                    else:
                        list_end_phase.append(HetMatSign.iloc[index+previous_step]['Position'])
                        meh=False
        array_=np.array((list_phase,list_init_phase,list_end_phase))
        return(array_)

    def save_bed(self, dict, i):
        for d in dict:
            dict[d]['Chr']= 'chr' + dict[d]['Chr'].astype(str)
            saving_file=pd.concat([dict[d]['Chr'], dict[d]['Start'],dict[d]['End'],dict[d]['Haplotype']], axis=1, keys=['Chr', 'Start','End','Haplotype'])
            saving_file.to_csv('./Trio'+str(i)+'Cell'+str(d)+'.bed',index=False,header=False,sep='\t')



    def phasing(self):
        """Calculates the phasing of the cells in the input file"""

        inputfile=self.inputfile
        HetMatIndex=inputfile.index[(inputfile['gDNA.GType']=='AB') | (inputfile['gDNA.GType']=='BA')]
        HetMatIndexPost=HetMatIndex+1
        HetMat=inputfile.iloc[HetMatIndex.tolist(),:]
        list_strings_egg=("\n".join(s for s in list(HetMat) if '.egg.GType'.lower() in s.lower()))
        list_num_trios=[int(s) for s in list(list_strings_egg) if s.isdigit()]
        HetMatSign=HetMat[(HetMat['1.egg.GType']!='NC')&(HetMat['1.egg.GType']!='AB')]
        HetMatSign=HetMatSign.reset_index()
        HetMatSign=HetMatSign.rename(columns = {'index':'SNP.num'})
        self.reference=HetMatSign['1.egg.GType']
        for i in list_num_trios:
            HetMatSign[str(i)+'.PB1.Phase']=pd.Series(float('nan'))
            HetMatSign.loc[:,str(i)+'.PB1.Phase'][(HetMatSign[str(i)+'.PB1.GType']!=self.reference)&(HetMatSign[str(i)+'.PB1.GType']!='NC')]=0.
            HetMatSign[str(i)+'.PB1.Phase'][(HetMatSign[str(i)+'.PB1.GType']==self.reference)&(HetMatSign[str(i)+'.PB1.GType']!='NC')]=1.
            HetMatSign[str(i)+'.PB1.Phase'][(HetMatSign[str(i)+'.PB1.GType']=='AB')&(HetMatSign[str(i)+'.PB1.GType']!='NC')]=0.5
            HetMatSign[str(i)+'.PB2.Phase']=pd.Series(float('nan'))
            HetMatSign[str(i)+'.PB2.Phase'][(HetMatSign[str(i)+'.PB2.GType']!=self.reference)&(HetMatSign[str(i)+'.PB2.GType']!='NC')]=0.
            HetMatSign[str(i)+'.PB2.Phase'][(HetMatSign[str(i)+'.PB2.GType']==self.reference)&(HetMatSign[str(i)+'.PB2.GType']!='NC')]=1.
            HetMatSign[str(i)+'.PB2.Phase'][(HetMatSign[str(i)+'.PB2.GType']=='AB')&(HetMatSign[str(i)+'.PB2.GType']!='NC')]=float('nan')   
            HetMatSign[str(i)+'.egg.Phase']=pd.Series(float('nan'))
            HetMatSign[str(i)+'.egg.Phase'][(HetMatSign[str(i)+'.egg.GType']!=self.reference)&(HetMatSign[str(i)+'.egg.GType']!='NC')]=0.
            HetMatSign[str(i)+'.egg.Phase'][(HetMatSign[str(i)+'.egg.GType']==self.reference)&(HetMatSign[str(i)+'.egg.GType']!='NC')]=1.
            HetMatSign[str(i)+'.egg.Phase'][(HetMatSign[str(i)+'.egg.GType']=='AB')&(HetMatSign[str(i)+'.egg.GType']!='NC')]=float('nan')   
            cells_results={}
            for cell in ['PB2','egg']:
                chr_results=[]
                for chr in HetMatSign['Chr'].drop_duplicates().values:
                    chr_HetMatSign=HetMatSign[HetMatSign['Chr']==chr]
                    array_=self.mark_snps(chr_HetMatSign,i,cell)
                    min_results=self.limit_snps(array_,chr_HetMatSign[[str(i)+'.'+str(cell)+'.Phase','Position']])
                    final_results=min_results.sort_values(by=['Start'])
                    final_results['Chr']=chr
                    chr_results.append(final_results)
                cells_results[cell]=pd.concat(chr_results)
            chr_results=[]
            for chr in HetMatSign['Chr'].drop_duplicates().values:
                chr_HetMatSign=HetMatSign[HetMatSign['Chr']==1]
                min_results=self.treat_PB1(cells_results,chr_HetMatSign,i)
                final_results=min_results.sort_values(by=['Start'])
                final_results['Chr']=chr
                chr_results.append(final_results)
            cells_results['PB1']=pd.concat(chr_results)
            self.save_bed(cells_results,i)




#end=phasing_class_tool('./input_data_trios.txt',20)

#end.phasing()
