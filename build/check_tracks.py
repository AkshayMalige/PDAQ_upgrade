import binascii
import os,sys

from datetime import datetime
import numpy
from array import array
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

import numpy as np
import pandas as pd


with open("44854soft_dump.txt") as soft:
    soft_con = soft.read().splitlines()
#print(len(soft_con))


with open("44854hard_dump.txt") as hard:
    hard_con = hard.read().splitlines()
#print(len(hard_con))

soft_list = list(dict.fromkeys(soft_con))
hard_list = list(dict.fromkeys(hard_con))

print(len(soft_list))
print(len(hard_list))

unique_soft = 0
unique_hard = 0
match0 =0
match1 =0

for i in soft_list:
	#print("soft  :",i)
	if i not in hard_list:
		unique_soft = unique_soft+1
		print(i)
	else:
		match0 = match0+1
		
for i in hard_list:
	if i not in soft_list:
		unique_hard = unique_hard+1
	else:
		match1 = match1+1
		
print("#Software_Events :",unique_soft)
print("#Hardware_Events :",unique_hard)
print("#Match :",match0)
#print(match1)
venn2(subsets = (unique_hard, unique_soft,match0 ), set_labels = ('Hardware', 'Software'))
plt.title('Enents id of tracks')
plt.show()


