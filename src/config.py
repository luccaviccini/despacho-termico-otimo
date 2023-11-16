import pandas as pd

# BARRA  CUSTO($/MWh)  MAX(MW) MIN(MW) INICIAL(MW) RAMPA(MW)
barras = [[1,10,30,5,0,10],
               [2,30,40,15,0,5],
               [3,100,40,0,0,3]]

DADOS_BARRAS = pd.DataFrame(barras, columns=['NUM_BARRA', "CUSTO($/MWh)", 'MAX(MW)',  "MIN(MW)","INICIAL (MW)", "RAMPA(MW)"])

# DE   PARA   SUSCEPTÂNCIA(OHMS) CONDUTÂNCIA(OHMS) LIMITES (MW)
linhas = [[1,2,1,33,25,20],
          [1,3,1,50,20,25],
          [1,3,2,50,20,25],
          [2,3,1,50,20,30]]

DADOS_LINHA = pd.DataFrame(linhas, columns=['DE', 'PARA','NLINHA','SUSCEPTÂNCIA(OHMS)', 'CONDUTÂNCIA(OHMS)', 'LIMITES (MW)'])

DADOS_DEMANDA = [[1,0,40,30],
               [2,0,43,25],
               [3,0,25,25]]



