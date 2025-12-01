import pandas as pd
import numpy as np
import sys
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt

ff_info=sys.argv[1]
autoff_chryff=pd.read_csv(ff_info,sep='\t',index_col=None)

autoff_chryff=autoff_chryff[autoff_chryff['chry_ff']<6]

autoff_chryff['idx']=range(len(autoff_chryff))

poly_features=PolynomialFeatures(degree=3,include_bias=True)

X_poly=poly_features.fit_transform(autoff_chryff['auto_ff'].values.reshape(len(autoff_chryff),1))

lin_reg=LinearRegression()
lin_reg.fit(X_poly,autoff_chryff['chry_ff'].values)

y_predict=lin_reg.predict(X_poly)

plt.scatter(autoff_chryff['auto_ff'].values,autoff_chryff['chry_ff'].values,color='blue')
plt.scatter(autoff_chryff['auto_ff'].values,y_predict,color='red')
plt.xlabel('auto_ff')
plt.ylabel('chry_ff')
plt.savefig('lowff_polynomial_predict.png')
