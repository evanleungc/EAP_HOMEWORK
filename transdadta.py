import pandas as pd
import numpy as np
import os
import datetime
os.chdir('C:/Users/Kanghua/Desktop/Im')
data = pd.read_csv('Q6.csv')
date = data.date
date = list(map(lambda x : datetime.datetime.strptime(x, '%Y/%m/%d').strftime('%Y%m%d'), date))
year = list(map(lambda x : x[0:4], date))
month = list(map(lambda x : x[4:6], date))
day = list(map(lambda x : x[6:8], date))
retlist = data.RET
ret2list = []
for i in retlist:
    try:
        i + 1
        ret2list.append(i)
    except:
        ret2list.append(np.nan)
data['RET'] = ret2list
data['date'] = date
data['year'] = year
data['month'] = month
data['day'] = day
data.to_csv('Q6_2.csv', index = False, header = None)