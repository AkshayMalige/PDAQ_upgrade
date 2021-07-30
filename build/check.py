import binascii
import os,sys

from datetime import datetime
import numpy
from array import array
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

import numpy as np
import pandas as pd


soft_df = pd.read_csv("lab_marking_track_44584_soft_dump.txt",sep='\t')
hard_df = pd.read_csv("lab_marking_track_20of256bin_44584_hard_dump.txt",sep='\t')
soft_df.columns = ['trigger', 's_slope', 's_const']
hard_df.columns = ['trigger', 'h_slope', 'h_const']


soft_df=soft_df.drop_duplicates(subset =['trigger','s_slope','s_const'],keep = 'first').reset_index(drop = True)
hard_df=hard_df.drop_duplicates(subset =['trigger','h_slope','h_const'],keep = 'first').reset_index(drop = True)

soft_multiplicity = list(soft_df.pivot_table(index=['trigger'], aggfunc='size'))
hard_multiplicity = list(hard_df.pivot_table(index=['trigger'], aggfunc='size'))

soft_df.sort_values("trigger", inplace = True)
hard_df.sort_values("trigger", inplace = True) 

soft_trig = list(soft_df.trigger.unique())
hard_trig = list(hard_df.trigger.unique())


soft_trigger = ['0x' + st for st in soft_trig]
hard_trigger = ['0x' + ht for ht in hard_trig]

softint_trigger = [int(x, 16) for x in soft_trigger]
hardint_trigger = [int(x, 16) for x in hard_trigger]

merge_df = pd.merge(soft_df,hard_df,how='inner',on='trigger')

merge_df['slope_diff']= round(abs(merge_df['s_slope'] - merge_df['h_slope']),3)
merge_df['const_diff']= round(abs(merge_df['s_const'] - merge_df['h_const']),3)

merge_df = merge_df.sort_values(by='const_diff', ascending=False)
merge_df = merge_df.drop_duplicates(subset='trigger', keep="last")

merge_df['slope_diff_r']= round((merge_df['s_slope'] - merge_df['h_slope']),3)
merge_df['const_diff_r']= round((merge_df['s_const'] - merge_df['h_const']),3)

merge_df.info()
merge_df.isna().sum()
print(merge_df)

merge_df.to_csv(r'pandas.txt', header=None, index=None, sep='\t', mode='a')

slope_diff = list(merge_df.slope_diff)
const_diff = list(merge_df.const_diff)

'''
plt.hist(const_diff, bins=np.arange(-5, 5 + 0.5, 0.5),density = True) 
#plt.hist(const_diff, bins=np.arange(min(const_diff),max(const_diff)+1,1))
plt.title('const_diff')
plt.xlabel('Const_diff [Cm]')
plt.ylim(0,10)
plt.show()

plt.hist(slope_diff, bins=np.arange(-5, 5+0.5, 0.5),density=True)
#plt.hist(slope_diff, bins=np.arange(min(slope_diff),max(slope_diff)+1,1))
plt.title('slope_diff')
plt.xlabel('Slope_diff [Cm]')
plt.ylim(0,10)
plt.show()

plt.hist(hard_multiplicity, bins=np.arange(min(hard_multiplicity), max(hard_multiplicity) + 1, 1))
plt.title('hard_track_mult')
plt.xlabel('hard_track_mult')
plt.show()

plt.hist(soft_multiplicity,2, range=[1, 3],density=True)
#plt.hist(soft_multiplicity, bins=np.arange(min(soft_multiplicity), max(soft_multiplicity) + 1, 1))
plt.title('soft_track_mult')
plt.xlabel('soft_track_mult')
plt.yscale('log')
plt.show() '''

'''
######## plot trigger no match ###########
unique_soft = 0
unique_hard = 0
match0 =0
match1 =0

for i in softint_trigger:
	#print("soft  :",i)
	if i not in hardint_trigger:
		unique_soft = unique_soft+1
		print(hex(i))
	else:
		match0 = match0+1
		
for i in hardint_trigger:
	if i not in softint_trigger:
		unique_hard = unique_hard+1
	else:
		match1 = match1+1		
print("#Software_Events :",unique_soft)
print("#Hardware_Events :",unique_hard)
print("#Match :",match0)
#print(match1)
#venn2(subsets = (unique_hard, unique_soft,match0 ), set_labels = ('Hardware', 'Software'))
plt.title('Enents id of tracks')
plt.show() '''


