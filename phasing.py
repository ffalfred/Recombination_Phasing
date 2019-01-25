import math
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

class phasing_class_tool:

    def __init__(self,BEDfile):
        """Instants of the files"""
        self.inputfile = pd.read_csv(BEDfile,sep='\t')   #Input

    def limit_snps(self,min_snps,array_,pos_array):
        pos_pd=pd.DataFrame(pos_array).set_index('Position')
        array_=np.append(array_,np.zeros((1,np.shape(array_)[1])),axis=0)
        for i in  range(np.shape(array_)[1]):
            num_snps=pd.notna(pos_pd.loc[array_[1,i]:array_[2,i]][pos_pd.columns[0]]).sum()
            array_[3,i]=num_snps
        haplo_count=pd.DataFrame(array_.T,columns=['Haplotype','Start','End','NumSNPS'])
        del_haplo_count=haplo_count[(haplo_count['Haplotype']!=0.5)&(haplo_count['NumSNPS']>=min_snps)].reset_index()
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

    def detect_centromer(self, data):
        data.dropna(inplace=True)
        data['dA'] = data['Position']-data['Position'].shift(-1)
        end_cent=data.loc[data['dA'].idxmin()+1]['Position']
        start_cent=data.loc[data['dA'].idxmin()]['Position']
        return start_cent,end_cent

    def mark_snps(self,HetMatSign,i,cell):
        list_phase=[]
        list_init_phase=[]
        list_end_phase=[]
        chromosome=[]
        previous_phase=float('nan')
        for index, row in HetMatSign[['Chr','Position',str(i)+'.'+str(cell)+'.Phase']].iterrows():
            if math.isnan(row[str(i)+'.'+str(cell)+'.Phase'])==False:
                if math.isnan(previous_phase): #initial
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
            if index==HetMatSign[['Chr','Position',str(i)+'.'+str(cell)+'.Phase']].index.values[-1]: #end
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


    def phasing(self, min_snps):
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
        for i in [2,3,4,1]:
            HetMatSign[str(i)+'.PB1.Phase']=pd.Series(float('nan'))
            HetMatSign[str(i)+'.PB1.Phase'][(HetMatSign[str(i)+'.PB1.GType']!=self.reference)&(HetMatSign[str(i)+'.PB1.GType']!='NC')]=0.
            HetMatSign[str(i)+'.PB1.Phase'][(HetMatSign[str(i)+'.PB1.GType']==self.reference)&(HetMatSign[str(i)+'.PB1.GType']!='NC')]=1.
            HetMatSign[str(i)+'.PB1.Phase'][(HetMatSign[str(i)+'.PB1.GType']=='AB')&(HetMatSign[str(i)+'.PB1.GType']!='NC')]=0.5
            HetMatSign[str(i)+'.PB2.Phase']=pd.Series(float('nan'))
            HetMatSign[str(i)+'.PB2.Phase'][(HetMatSign[str(i)+'.PB2.GType']!=self.reference)&(HetMatSign[str(i)+'.PB2.GType']!='NC')]=0.
            HetMatSign[str(i)+'.PB2.Phase'][(HetMatSign[str(i)+'.PB2.GType']==self.reference)&(HetMatSign[str(i)+'.PB2.GType']!='NC')]=1.
            HetMatSign[str(i)+'.PB2.Phase'][(HetMatSign[str(i)+'.PB2.GType']=='AB')&(HetMatSign[str(i)+'.PB2.GType']!='NC')]=float('nan')   #heterozygosis
            HetMatSign[str(i)+'.egg.Phase']=pd.Series(float('nan'))
            HetMatSign[str(i)+'.egg.Phase'][(HetMatSign[str(i)+'.egg.GType']!=self.reference)&(HetMatSign[str(i)+'.egg.GType']!='NC')]=0.
            HetMatSign[str(i)+'.egg.Phase'][(HetMatSign[str(i)+'.egg.GType']==self.reference)&(HetMatSign[str(i)+'.egg.GType']!='NC')]=1.
            HetMatSign[str(i)+'.egg.Phase'][(HetMatSign[str(i)+'.egg.GType']=='AB')&(HetMatSign[str(i)+'.egg.GType']!='NC')]=float('nan')   #heterozygosis
            cells_results={}
            for cell in ['PB2','egg']:
                array_=self.mark_snps(HetMatSign,i,cell)
                start_cent,end_cent=self.detect_centromer(HetMatSign[['Position',str(i)+'.'+str(cell)+'.Phase']])
                min_results=self.limit_snps(min_snps,array_,HetMatSign[[str(i)+'.'+str(cell)+'.Phase','Position']])
                new_add=pd.Series([min_results[(min_results['Start']<start_cent)&(min_results['End']>end_cent)]['Haplotype'].values[0],end_cent,min_results[(min_results['Start']<start_cent)&(min_results['End']>end_cent)]['End'].values[0],float('nan')],index=list(min_results))
                min_results=min_results.append(new_add,ignore_index=True)
                min_results.loc[(min_results['Start']<start_cent)&(min_results['End']>end_cent),'End']=start_cent
                min_results.loc[(min_results['Start']<start_cent)&(min_results['End']>end_cent),'NumSNPS']=float('nan')
                final_results=min_results.sort_values(by=['Start'])
                print(final_results)
                cells_results[cell]=final_results
            print(cells_results)
            exit()




end=phasing_class_tool('./input_data_trios.txt')

end.phasing(20)
