# virtual-screening-covid19
Modelos para protease principal 3C-like protease (Mpro) de SARS-CoV
Os modelos nesse fluxo de trabalho foram criados utilizando os decritores:

## Artigos e tutoriais
- [De novo design and bioactivity prediction of SARS‑CoV‑2 main protease inhibitors using recurrent neural network‑based transfer learning](https://www.arca.fiocruz.br/bitstream/icict/46077/2/Santana_Marcos_etal_IOC_2021_COVID-19.pdf)
- [Virtual Screening on Indonesian Herbal Compounds as COVID-19 Supportive Therapy: Machine Learning and Pharmacophore Modelling Approaches](https://assets.researchsquare.com/files/rs-29119/v2_stamped.pdf?c=1603225813)
- [Computational Models Identify Several FDA Approved or Experimental Drugs as Putative Agents Against SARS-CoV-2](https://chemrxiv.org/engage/chemrxiv/article-details/60c74a200f50dba0b23969fa)
- [QSAR tutorial. Dr. Pavel Polishchuk] (http://www.qsar4u.com/files/qsar_rdkit_tutorial/qsar-rdkit.html)

## Base de dados de treinamento e testes (potenciais inibidores de SARS-CoV-2 3CLpro.)

- [950 compostos com atividade ativa] (https://plpro-inhibitors.cent.uw.edu.pl/)
- [1662 compostos parte dos dados com atividade Inativo] (https://github.com/marcossantanaioc/De_novo_design_SARSCOV2/blob/master/data/ChEMBL_v1.zip)

## Descritores
- [PaDEL-Descriptor](https://github.com/ecrl/padelpy)
- [Impressões digitais de Morgan](https://www.rdkit.org/docs/GettingStartedInPython.html#morgan-fingerprints-circular-fingerprints)
- [SiRMS (Simplex Representation of Molecular Structure)](http://www.qsar4u.com/pages/sirms.php)
- [Descritos gerados a partir do RDKIT](https://www.rdkit.org/docs/GettingStartedInPython.html#list-of-available-descriptors)

## Algoritmos de classificação
E os algorítimos utilizando as bibliotecas da scikit-learn:
- [Random Forest classifier](https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html)
- [Support Vector Machine classifier](https://scikit-learn.org/stable/modules/svm.html)
- [Multi-layer Perceptron classifier](https://scikit-learn.org/stable/modules/neural_networks_supervised.html)

### Modelos QSAR criados 
No total de 12 moledos criados com a combinação dos descritores e algorítimos:
1. Random Forest
 - RF_Morgan
 - RF_SiRMS
 - RF_Rdkit
 - RF_PaDEL
2. Support Vector Machine
 - SVM_Morgan
 - SVM_SiRMS
 - SVM_Rdkit
 - SVM_PaDEL

3. Multi-layer Perceptron
 - MLP_Morgan
 - MLP_SiRMS
 - MLP_Rdkit
 - MLP_PaDEL
