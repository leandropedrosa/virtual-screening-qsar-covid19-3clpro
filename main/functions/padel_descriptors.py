import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from math import floor
#Rdkit: coleção de quiminformática e software de aprendizado de máquina escrito em C++ e Python de Código Aberto.
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import PandasTools
from rdkit.Chem.Draw import IPythonConsole
from collections import Counter
from sklearn.feature_selection import VarianceThreshold

## Carregar dados
def padel_descriptors(moldf):
    desc = pd.read_csv('descriptors/descritores_padel.csv', sep=',')
    desc.drop(desc.columns[0:1], axis=1,inplace=True)
    descriptors = desc.columns.difference(moldf.columns).tolist()
    desc.head()
    moldf_desc = pd.concat([moldf,desc], axis=1)
    balance_data = 'no'

    if balance_data == 'yes':
        # Equilibre os dados usando 1/2 similaridade e 1/2 aleatória
        moldf_desc = BalanceBySim(moldf_desc, 'Outcome', 2)
        # Forma de impressão
        print('Forma do conjunto de treinamento: %s' % Counter(moldf_desc['Outcome'].loc[moldf_desc['Set'] == 'train']))
        print('Forma externa definida: %s' % Counter(moldf_desc['Outcome'].loc[moldf_desc['Set'] == 'ext']))

    else:
        moldf_desc['Set'] = 'train'
        # Forma de impressão
        print('Forma do conjunto de treinamento: %s' % Counter(moldf_desc['Outcome']))
        print('Forma externa definida: %s' % Counter(moldf_desc['Outcome'].loc[moldf_desc['Set'] == 'ext']))
        
    moldf_train = moldf_desc[(moldf_desc['Set'] == 'train')]

    y_train = moldf_train['Outcome'].to_numpy()
    X_train = moldf_train[descriptors]
    X_train.shape        
    
    ##### Remover variáveis constantes e quase constantes
    X_train = X_train.select_dtypes(exclude=['object'])
    X_train = X_train.dropna(axis=1, how='any')
    X_train = X_train.fillna(0)

    # Definir filtro de baixa variação (limite de 10%)
    def variance_filter(data, threshold=0.1):
        selector = VarianceThreshold(threshold)
        selector.fit(data)
        return data[data.columns[selector.get_support(indices=True)]]

    # Aplicar filtro
    X = variance_filter(X_train)
    
    ##### Remover variáveis correlacionadas
    correlated_features = set()  
    correlation_matrix = X_train.corr()

    for i in range(len(correlation_matrix.columns)):  
        for j in range(i):
            if abs(correlation_matrix.iloc[i, j]) > 0.9:
                colname = correlation_matrix.columns[i]
                correlated_features.add(colname)

    X_train.drop(labels=correlated_features, axis=1, inplace=True)

    X_train.shape
    X_train.to_csv('descriptors/padel-chembl-sars-cov-3C-like-proteinase-processed.txt', sep='\t', index=False)
    
    data_train = {'moldf_desc': moldf_desc, 'moldf_train': moldf_train, 'Y_train': moldf_train['Outcome'].to_numpy(), 'X_train': moldf_train[descriptors]}
    return data_train
    