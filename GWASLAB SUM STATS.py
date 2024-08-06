#!/usr/bin/env python
# coding: utf-8

# In[1]:


#GWASLab Reference:
# He, Y., Koido, M., Shimmori, Y., Kamatani, Y. (2023). GWASLab: a Python package for processing and visualizing GWAS summary statistics. Preprint at Jxiv, 2023-5. https://doi.org/10.51094/jxiv.370


# In[2]:


#import gwaslab and upload gwas sum stats
import gwaslab as gl

mysumstats = gl.Sumstats("EAGLE_AD_no23andme_results_29072015.txt",
             rsid="rsID",
             chrom="chromosome",
             pos="position",
             ea="reference_allele",
             nea="other_allele",
             eaf="eaf",
             n="AllEthnicities_N",
             beta="beta",
             se="se",
             p="p.value",
             i2= "i2",
             build= "19")
mysumstats.data


# In[3]:


#Fill missing columns
mysumstats.fill_data(to_fill=["Z","MAF","MLOG10P"])
mysumstats.data


# In[4]:


#remove variants in HLA region and run a check
mysumstats.exclude_hla()
mysumstats.basic_check(remove=True)
mysumstats.data


# In[4]:


mysumstats.data.to_csv("cleandatabuild37.txt",sep='\t', index=False)


# In[6]:


mysumstats.data


# In[5]:


#Fix ids
mysumstats.fix_id(fixsep=True)


# In[6]:


#normalize alleles
mysumstats.normalize_allele()


# In[8]:


mysumstats.liftover(n_cores=3, from_build="19", to_build="38")
mysumstats.data.to_csv("liftoveredgwas.txt",sep='\t', index=False)


# In[10]:


#Get summary
mysumstats.summary()


# In[11]:


#Download reference data
gl.download_ref("ensembl_hg19_gtf", overwrite=True)
gl.get_path("ensembl_hg19_gtf")


# In[12]:


output_df = mysumstats.get_lead(windowsizekb=500,
                                sig_level=5e-08,
                                anno=True,
                                build="19",
                                source="ensembl",
                                verbose=True)


# In[13]:


get_ipython().system('pip install kaleido')


# In[27]:


#store table as image
import matplotlib.pyplot as plt
import plotly.graph_objects as go


fig= go.Figure(data=[go.Table(
    columnwidth=[400]*len(output_df.columns),
    header= dict(values= list (output_df.columns),
                                           align='left', font=dict(size=20)),
    cells=dict(values=[output_df[col] for col in output_df.columns],
                                        align='left', font=dict(size=15)))
])
fig.update_layout(width=2000, height=900)
fig.write_image("leadsnptable.png")


# In[5]:


#Calculate Novel SNPs using EFO
mysumstats.get_novel(efo="EFO_0000274")


# In[ ]:


#Create MMQ PLot
mysumstats.plot_mqq(skip=3,cut=20,anno=True,mode= "mqq",check=False,verbose=False,save="GWAS_MMQ_plot.png", 
                    save_args={"dpi":400,"facecolor":"white"})


# In[ ]:


#Create FLG regional plot
mysumstats.plot_mqq(mode="r",
                    region=(1,151302165,153325239),                    
                    gtf_path=gl.get_path("ensembl_hg19_gtf"),
                    anno="GENENAME",
                    save="flg_MMQ_plot.png", 
                    save_args={"dpi":400,"facecolor":"white"})
                   


# In[ ]:


#Create Power plot
a = gl.plot_power_x(mafs=[0.5,0.4,0.3,0.2,0.9,0.1], 
                    betas= 0.1,  
                    xscale="nonlog",
                    save="power_plot.png", 
                    save_args={"dpi":400,"facecolor":"white"})


# In[ ]:




