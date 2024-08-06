#!/usr/bin/env python
# coding: utf-8

# In[37]:


import numpy
import scipy
import pandas
import sqlalchemy
import patsy
import statsmodels
import h5py

get_ipython().system('git clone https://github.com/hakyimlab/MetaXcan')
get_ipython().system('cd MetaXcan/software')


# In[9]:


#harmonize data set for s-predixcan
get_ipython().system('python MetaXcan/software/M03_betas.py  --snp_map_file /Users/halk/Downloads/data/coordinate_map/map_snp150_hg38.txt.gz  --gwas_file /Users/halk/liftoveredgwas.txt  --snp_column rsID  --non_effect_allele_column NEA  --effect_allele_column EA  --beta_column BETA  --pvalue_column P  --keep_non_rsid  --throw  --output output/newharmonizedgwas.txt.gz')


# In[12]:


#TWAS on whole blood
get_ipython().system('python MetaXcan/software/SPrediXcan.py  --gwas_file /Users/halk/output/newharmonizedgwas.txt.gz  --snp_column snp  --effect_allele_column effect_allele  --non_effect_allele_column non_effect_allele  --zscore_column zscore  --model_db_path /Users/halk/Downloads/eqtl/mashr/mashr_Whole_Blood.db  --covariance /Users/halk/Downloads/eqtl/mashr/mashr_Whole_Blood.txt.gz  --keep_non_rsid  --additional_output  --model_db_snp_key varID  --throw  --output_file output/Wholebloodresults.csv')


# In[13]:


#TWAS on skin no sun
get_ipython().system('python MetaXcan/software/SPrediXcan.py  --gwas_file /Users/halk/output/newharmonizedgwas.txt.gz  --snp_column snp  --effect_allele_column effect_allele  --non_effect_allele_column non_effect_allele  --zscore_column zscore  --model_db_path /Users/halk/Downloads/eqtl/mashr/mashr_Skin_Not_Sun_Exposed_Suprapubic.db  --covariance /Users/halk/Downloads/eqtl/mashr/mashr_Skin_Not_Sun_Exposed_Suprapubic.txt.gz  --keep_non_rsid  --additional_output  --model_db_snp_key varID  --throw  --output_file output/Skinnosunresults.csv')


# In[14]:


#TWAS on Skin Sun exposed
get_ipython().system('python MetaXcan/software/SPrediXcan.py  --gwas_file /Users/halk/output/newharmonizedgwas.txt.gz  --snp_column snp  --effect_allele_column effect_allele  --non_effect_allele_column non_effect_allele  --zscore_column zscore  --model_db_path /Users/halk/Downloads/eqtl/mashr/mashr_Skin_Sun_Exposed_Lower_leg.db  --covariance /Users/halk/Downloads/eqtl/mashr/mashr_Skin_Sun_Exposed_Lower_leg.txt.gz  --keep_non_rsid  --additional_output  --model_db_snp_key varID  --throw  --output_file output/skinsunexposed.csv')


# In[15]:


#TWAS on Fibroblasts
get_ipython().system('python MetaXcan/software/SPrediXcan.py  --gwas_file /Users/halk/output/newharmonizedgwas.txt.gz  --snp_column snp  --effect_allele_column effect_allele  --non_effect_allele_column non_effect_allele  --zscore_column zscore  --model_db_path /Users/halk/Downloads/data/models/eqtl/mashr/mashr_Cells_Cultured_fibroblasts.db --covariance /Users/halk/Downloads/eqtl/mashr/mashr_Cells_Cultured_fibroblasts.txt.gz --keep_non_rsid  --additional_output  --model_db_snp_key varID  --throw  --output_file output/fibroblastresults.csv')


# In[18]:


#TWAS on Lymphocytes
get_ipython().system('python MetaXcan/software/SPrediXcan.py  --gwas_file /Users/halk/output/newharmonizedgwas.txt.gz  --snp_column snp  --effect_allele_column effect_allele  --non_effect_allele_column non_effect_allele  --zscore_column zscore  --model_db_path /Users/halk/Downloads/data/models/sqtl/mashr/mashr_Cells_EBV-transformed_lymphocytes.db  --covariance /Users/halk/Downloads/eqtl/mashr/mashr_Cells_EBV-transformed_lymphocytes.txt.gz  --keep_non_rsid  --additional_output  --model_db_snp_key varID  --throw  --output_file output/ebvlymphocytesresults.csv')


# In[19]:


get_ipython().system('pip install qmplot')


# In[13]:


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns 
import numpy as np
#Upload whole blood results
WB= pd.read_csv('/Users/halk/output/Wholebloodresults.csv')


# In[59]:


#create bar plot for lead genes from WB twas
def lead_genes_plot(WB, lead_n=10, gene_col='gene_name', p_value_col= 'pvalue'):
    lead_genes= WB.nsmallest(lead_n, p_value_col)
    lead_genes['-log10(pvalue)']= -np.log10(lead_genes[p_value_col])
    
    plt.figure(figsize=(12,10))
    sns.barplot(x='-log10(pvalue)', y=gene_col, data= lead_genes)
    plt.xlabel('-log10(p-value)')
    plt.ylabel('Gene')
    plt.title('Top 10 Significant Genes for Whole Blood')
    plt.savefig('WBtop10.png')
    plt.show()
    
lead_genes_plot(WB)


# In[89]:


#Create QQ plot for all 4 tissue
from scipy import stats
def qq_plot(dfs, pvalue_cols,labels):
    plt.figure(figsize=(10,10))
    
    max_expected=0
    max_observed =0
    
    
    for i, df in enumerate(dfs):
    
        pvalues= df[pvalue_cols[i]]
   
    #Calculate expected qs
        sorted_pvalues = np.sort(pvalues)
        expected_qs = -np.log10(np.linspace(1/ len(sorted_pvalues),1, len(sorted_pvalues)))
    
    #calculate observed -log10ps
        observed_qs= -np.log10(sorted_pvalues)
    
       #create values for reference line
        max_expected= max(max_expected, max(expected_qs))
        max_observed= max(max_observed, max(observed_qs))
    #plot
        plt.plot(expected_qs,observed_qs, 'o', markersize=2,label=labels[i])
    
    
    
    #Create line
    max_val = max(max_expected, max_observed)
    plt.plot([0, max_val],[0, max_val], 'r--')
    
    #label graph
    plt.xlabel('Expected -log10(p-values)')
    plt.ylabel('Observed -log10(p-values)')
    plt.title('QQ Plot for All Tissues')
    plt.legend()
    plt.xlim([0,5])
    plt.savefig('QQplot.png')
    plt.show()
    
dfs= [WB,CF,sun,nosun]
pvalue_cols= ['pvalue','pvalue','pvalue','pvalue']
labels=['Whole Blood','Fibroblasts','Sun Exposed Skin','Sun-free Skin']
qq_plot(dfs, pvalue_cols,labels)


# In[11]:


#Upload fibroblast results
CF=pd.read_csv('/Users/halk/output/fibroblastresults.csv')


# In[62]:


#create bar plot for lead genes from CF twas
def lead_genes_plot(CF, lead_n=10, gene_col='gene_name', p_value_col= 'pvalue'):
    lead_genes= CF.nsmallest(lead_n, p_value_col)
    lead_genes['-log10(pvalue)']= -np.log10(lead_genes[p_value_col])
    
    plt.figure(figsize=(12,10))
    sns.barplot(x='-log10(pvalue)', y=gene_col, data= lead_genes)
    plt.xlabel('-log10(p-value)')
    plt.ylabel('Gene')
    plt.title('Top 10 Significant Genes for Cultured Fibroblasts')
    plt.savefig('CFtop10.png')
    plt.show()
    
lead_genes_plot(CF)


# In[3]:


#save sun exposed skin results as sun
sun=pd.read_csv('/Users/halk/output/skinsunexposed.csv')


# In[65]:


#create bar plot for lead genes from Sun exposed skin twas
def lead_genes_plot(sun, lead_n=10, gene_col='gene_name', p_value_col= 'pvalue'):
    lead_genes= sun.nsmallest(lead_n, p_value_col)
    lead_genes['-log10(pvalue)']= -np.log10(lead_genes[p_value_col])
    
    plt.figure(figsize=(12,10))
    sns.barplot(x='-log10(pvalue)', y=gene_col, data= lead_genes)
    plt.xlabel('-log10(p-value)')
    plt.ylabel('Gene')
    plt.title('Top 10 Significant Genes for Sun Exposed Skin')
    plt.savefig('suntop10.png')
    plt.show()
    
lead_genes_plot(sun)


# In[7]:


#save non-sun exposed skin results as nosun
nosun=pd.read_csv('/Users/halk/output/Skinnosunresults.csv')


# In[71]:


#create bar plot for lead genes from Sun exposed skin twas
def lead_genes_plot(nosun, lead_n=10, gene_col='gene_name', p_value_col= 'pvalue'):
    lead_genes= nosun.nsmallest(lead_n, p_value_col)
    lead_genes['-log10(pvalue)']= -np.log10(lead_genes[p_value_col])
    
    plt.figure(figsize=(12,10))
    sns.barplot(x='-log10(pvalue)', y=gene_col, data= lead_genes)
    plt.xlabel('-log10(p-value)')
    plt.ylabel('Gene')
    plt.title('Top 10 Significant Genes for Sun-free Skin')
    plt.savefig('2nosuntop10.png')
    plt.show()
    
lead_genes_plot(nosun)


# In[6]:


#create bar plot of top 10 zscores for sun exposed skin
data_sorted=sun.sort_values('zscore',ascending= False)
top_n=10
top_data= data_sorted.head(top_n)

#plot zscores
plt.figure(figsize=(12,10))
plt.barh(top_data['gene_name'],top_data['zscore'])
plt.xlabel('Z-score')
plt.ylabel('Gene')
plt.title('Top 10 Significant Z-scores for Sun Exposed Skin')
plt.savefig('zsuntop10.png')
plt.show()


# In[10]:


#create bar plot of top 10 zscores for nonsun exposed skin
data_sorted=nosun.sort_values('zscore',ascending= False)
top_n=10
top_data= data_sorted.head(top_n)

#plot zscores
plt.figure(figsize=(12,10))
plt.barh(top_data['gene_name'],top_data['zscore'])
plt.xlabel('Z-score')
plt.ylabel('Gene')
plt.title('Top 10 Significant Z-scores for Sun-Free Skin')
plt.savefig('znosuntop10.png')
plt.show()


# In[12]:


#create bar plot of top 10 zscores for cultured fibroblasts
data_sorted=CF.sort_values('zscore',ascending= False)
top_n=10
top_data= data_sorted.head(top_n)

#plot zscores
plt.figure(figsize=(12,10))
plt.barh(top_data['gene_name'],top_data['zscore'])
plt.xlabel('Z-score')
plt.ylabel('Gene')
plt.title('Top 10 Significant Z-scores for Cultured Fibroblasts')
plt.savefig('zCFtop10.png')
plt.show()


# In[17]:


#create bar plot of top 10 zscores for whole blood
data_sorted=WB.sort_values('zscore',ascending= False)
top_n=10
top_data= data_sorted.head(top_n)

#plot zscores
plt.figure(figsize=(12,10))
plt.barh(top_data['gene_name'],top_data['zscore'])
plt.xlabel('Z-score')
plt.ylabel('Gene')
plt.title('Top 10 Significant Z-scores for Whole Blood')
plt.savefig('zWBtop10.png')
plt.show()


# In[ ]:




