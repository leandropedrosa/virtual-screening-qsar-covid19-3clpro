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

## Carregar dados
def carregar_dados():
    # Definir caminho do arquivo
    file = '../dataset/formats/sdf/COVID-Sulkowsk.sdf'

    # Novo dicionário inicializado a partir de um objeto de mapeamento
    sdfInfo = dict(smilesName='SMILES', molColName='ROMol')

    # Carregando o arquivo SDF com os dicionarios mapeados
    moldf = PandasTools.LoadSDF(file, **sdfInfo)
    print('Original data: ', moldf.shape)

    # Renomear ROMol
    moldf = moldf.rename(columns={'ROMol': 'Mol'})

    # Remover moléculas RDKit ausentes
    moldf = moldf[pd.notnull(moldf['Mol'])]
    if 'StandardizerResult' in moldf.columns:
        moldf = moldf.drop(columns='StandardizerResult')

    # Colunas
    print('Dados mantidos: ', moldf.shape)
    
    ## Forma dos dados
    # (27 ativos e 64 inativos) 91 compostos utilizando o software ChemAxon Standardizer 
    # (13 ativos e 09 inativos) 22 compostos obtidos de empresas encontradas do PDB
    # 2 Classes criadas Classe 1: 40 (Ativos) e Classe 0: 73 (Inativos)

    moldf['Outcome'] = moldf['Outcome'].replace('Active', 1)
    moldf['Outcome'] = moldf['Outcome'].replace('Inactive', 0)

    classes = Counter(moldf['Outcome'])
    print('\033[1m' + 'Forma do conjunto de treinamento:' + '\n' + '\033[0m')
    for key, value in classes.items():
        print('\t\t Classe %d: %d' % (key, value))
    print('\t\t Número total de compostos: %d' % (len(moldf['Outcome'])))

    print('Class labels:', np.unique(classes))
  
    return moldf
    