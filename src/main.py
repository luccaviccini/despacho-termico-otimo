import numpy as np
from cvxopt.modeling import variable, op
from config import DADOS_BARRAS, DADOS_LINHA, DADOS_DEMANDA
import time

def otimizacao_despacho(DADOS_BARRAS, DADOS_LINHA, DADOS_DEMANDA):
  
  n_horas = np.size(DADOS_DEMANDA,0)
  n_barras = np.size(DADOS_BARRAS,0) 

  v_gt = variable(size = n_horas*n_barras,name='Geração Térmica (MW)')
  v_theta = variable(size = n_horas*n_barras,name='Ângulo das barras (rad)')
  v_deficit = variable(size = n_horas*n_barras, name='Déficit (MW)')

  r_balanco_carga = []  
  for hora in range(n_horas):
    for ute in range(len(DADOS_BARRAS)):
      barra = DADOS_BARRAS['NUM_BARRA'][ute]

      balanco_potencia = 0
      balanco_potencia = balanco_potencia + v_gt[ute + (hora)*(n_barras)] + v_deficit[ute + (hora)*(n_barras)]

      for lin in range(len(DADOS_LINHA)):
        de = int(DADOS_LINHA['DE'][lin])
        para = int(DADOS_LINHA['PARA'][lin])
        susceptancia = float(DADOS_LINHA['SUSCEPTÂNCIA(OHMS)'][lin])

        if de==int(barra):
          balanco_potencia = balanco_potencia - ((v_theta[(de-1) + (hora)*(n_barras)] - v_theta[(para-1) + (hora)*(n_barras)])*susceptancia)
        elif para==int(barra):
          balanco_potencia = balanco_potencia - ((v_theta[(para-1) + (hora)*(n_barras)] - v_theta[(de-1) + (hora)*(n_barras)])*susceptancia)


      r_balanco_carga.append(balanco_potencia==float((DADOS_DEMANDA[hora][barra])))


  r_gt_max_min, r_rampa, r_deficit = list(), list(), list()
  for hora in range(n_horas):
    for ute in range(len(DADOS_BARRAS)):
      barra = DADOS_BARRAS['NUM_BARRA'][ute]

      r_gt_max_min.append(v_gt[ute + (hora)*(n_barras)]<=float(DADOS_BARRAS['MAX(MW)'][ute]))
      r_gt_max_min.append(v_gt[ute + (hora)*(n_barras)]>=float(DADOS_BARRAS['MIN(MW)'][ute]))

      r_deficit.append(v_deficit[ute + (hora)*(n_barras)]>=0)
      r_deficit.append(v_deficit[ute + (hora)*(n_barras)]<=float(DADOS_DEMANDA[hora][barra]))

      if (hora==0):
        r_rampa.append(v_gt[ute + (hora)*(n_barras)] <=100000)
        r_rampa.append(v_gt[ute + (hora)*(n_barras)] >=-100000)
      else:
        r_rampa.append(v_gt[ute + (hora)*(n_barras)] - v_gt[ute + (hora-1)*(n_barras)]<=float(DADOS_BARRAS['RAMPA(MW)'][ute]))
        r_rampa.append(v_gt[ute + (hora)*(n_barras)] - v_gt[ute + (hora-1)*(n_barras)]>=-float(DADOS_BARRAS['RAMPA(MW)'][ute]))

  r_linha = list()
  for hora in range(n_horas):
    for lin in range(len(DADOS_LINHA)):
      de = int(DADOS_LINHA['DE'][lin])
      para = int(DADOS_LINHA['PARA'][lin])
      susceptancia = float(DADOS_LINHA['SUSCEPTÂNCIA(OHMS)'][lin])

      r_linha.append((v_theta[(de-1) + (hora)*(n_barras)] - v_theta[(para-1) + (hora)*(n_barras)])<=float((1/susceptancia)*DADOS_LINHA['LIMITES (MW)'][lin]))
      r_linha.append((v_theta[(de-1) + (hora)*(n_barras)] - v_theta[(para-1) + (hora)*(n_barras)])>=-float((1/susceptancia)*DADOS_LINHA['LIMITES (MW)'][lin]))

  r_angular = list()
  for hora in range(n_horas):
    for bus in range(len(DADOS_BARRAS)):
      if bus==0: 
        r_angular.append(v_theta[bus + (hora)*n_barras]<=0)
        r_angular.append(v_theta[bus + (hora)*n_barras]>=0)

      else:
        r_angular.append(v_theta[bus + (hora)*n_barras]<=float(np.pi))
        r_angular.append(v_theta[bus + (hora)*n_barras]>=-float(np.pi))

  restricoes = r_balanco_carga + r_gt_max_min + r_rampa + r_linha + r_angular + r_deficit

  fob = 0
  for hora in range(n_horas):
    for ute in range(len(DADOS_BARRAS)):
      fob = fob + float(DADOS_BARRAS['CUSTO($/MWh)'][ute])*v_gt[ute + (hora)*n_barras]
      fob = fob + 50000000*v_deficit[ute + (hora)*(n_barras)]


  # Solução da função objetivo
  problema = op(fob,restricoes)
  problema.solve('dense','glpk') 
   
  # Resultados das barras
  geracao, corte, angulo = dict(), dict(), dict()
  for bus in DADOS_BARRAS.index:

    geracao[str(DADOS_BARRAS['NUM_BARRA'][bus])] = [round(float(v_gt[bus + (h-1)*n_barras].value()[0]), 2) for h in range(1,n_barras+1)]
    corte[str(DADOS_BARRAS['NUM_BARRA'][bus])] = [round(float(v_deficit[bus + (h-1)*n_barras].value()[0]), 2) for h in range(1,n_barras+1)]
    angulo[str(DADOS_BARRAS['NUM_BARRA'][bus])] = [round(float(v_theta[bus + (h-1)*n_barras].value()[0]), 2) for h in range(1,n_barras+1)]


  fluxo = dict()

  for lin in DADOS_LINHA.index:

    de = int(DADOS_LINHA['DE'][lin])
    para = int(DADOS_LINHA['PARA'][lin])
    susceptancia = float(DADOS_LINHA['SUSCEPTÂNCIA(OHMS)'][lin])

    de_onde_para_onde = str(de)+'-'+str(para)+'('+str(DADOS_LINHA['NLINHA'][lin])+')'

    fluxo[de_onde_para_onde] = [((float(v_theta[(de-1) + (h-1)*(n_barras)].value()[0]) - float(v_theta[(para-1) + (h-1)*(n_barras)].value()[0]))*susceptancia) for h in range(1,n_barras+1)]

  return geracao,corte,angulo,fluxo

start_time = time.time()
sol = otimizacao_despacho( DADOS_BARRAS, DADOS_LINHA, DADOS_DEMANDA)
end_time = time.time()
print("Tempo de execução: ", end_time - start_time)

print(sol[0])
print(sol[1])   
print(sol[2])                    
print(sol[3])

def verificar_fluxos(geracao, demanda, fluxos):
    # Assumindo que geracao e fluxos são dicionários com as chaves sendo as barras ou linhas e os valores sendo listas com os valores por hora.
    # Demanda é uma lista de listas no formato [[hora, demanda_barra1, demanda_barra2, demanda_barra3], ...]

    for h in range(len(demanda)):  # itera através de cada hora
        hora, *demanda_barras = demanda[h]  # separa a hora e a demanda de cada barra
        
        for i, barra in enumerate(geracao.keys(), start=1):  # itera pelas barras
            # Calcula o balanço de energia para a barra na hora h
            balance = geracao[barra][h] - demanda_barras[i-1]

            # Soma os fluxos entrantes e subtrai os fluxos saindo
            for linha, fluxo in fluxos.items():
                if str(barra) == linha.split('-')[0]:  # se a barra atual é a origem na linha
                    balance -= fluxo[h]  # subtrai o fluxo saindo
                elif str(barra) == linha.split('-')[1].split('(')[0]:  # se a barra atual é o destino na linha
                    balance += fluxo[h]  # soma o fluxo entrante

            # Checa se o balanço é zero (considerando uma pequena margem de erro)
            if abs(balance) > 1e-6:  # 1e-6 é a margem de erro que você está disposto a aceitar
                print(f"Desbalanceamento detectado na barra {barra} na hora {hora}: {balance}")
                return False

    print("Todos os fluxos estão corretos.")
    return True
  
verificar_fluxos(sol[0], DADOS_DEMANDA, sol[3])

import matplotlib.pyplot as plt

def plot_geracao_deficit(geracao, deficit):
    # Supondo que 'geracao' e 'deficit' sejam dicionários com chaves representando as barras (1, 2, 3, ...)
    # e os valores sejam listas com os valores de geração e déficit por hora.

    barras = sorted(geracao.keys())  # Ordena as chaves para plotar na ordem
    n_horas = range(1, len(next(iter(geracao.values()))) + 1)  # Começa de 1 em vez de 0

    # Cria uma figura com três subplots lado a lado
    fig, axes = plt.subplots(1, len(barras), figsize=(15, 5))
    
    for i, barra in enumerate(barras):
        # Plota os valores de geração e déficit para cada barra com marcadores em cada ponto
        axes[i].plot(n_horas, geracao[barra], '-o', label='Geração')
        axes[i].plot(n_horas, deficit[barra], '-o', label='Déficit', linestyle='--')
        
        axes[i].set_title(f'Barra {barra}')
        axes[i].set_xlabel('Hora')
        axes[i].set_ylabel('MW')
        axes[i].legend()
        axes[i].grid(True)  # Adiciona o grid
        axes[i].set_ylim(-2, 42)  # Define o limite inferior e superior do eixo y

    plt.tight_layout()
    #plt.show()
    
plot_geracao_deficit(sol[0], sol[1])


def plot_angulos(angulos):
    n_horas = list(range(1, len(next(iter(angulos.values()))) + 1))  # assume todas as barras tem o mesmo nº de horas
    fig, ax = plt.subplots()
    for barra, angulos_barra in angulos.items():
        ax.plot(n_horas, angulos_barra, '-o', label=f'Barra {barra}')
    ax.set_xlabel('Hora')
    ax.set_ylabel('Ângulo (rad)')
    ax.set_title('Ângulos por Barra')
    ax.legend()
    ax.grid(True)
    #plt.show()
plot_angulos(sol[2])
def plot_fluxos(fluxos):
    n_horas = list(range(1, len(next(iter(fluxos.values()))) + 1))  # assume todas as linhas tem o mesmo nº de horas
    fig, ax = plt.subplots()
    for linha, fluxos_linha in fluxos.items():
        ax.plot(n_horas, fluxos_linha, '-o', label=f'Linha {linha}')
    ax.set_xlabel('Hora')
    ax.set_ylabel('Fluxo (em MW)')
    ax.set_title('Fluxo por Linha de Transmissão')
    ax.legend()
    ax.grid(True)
    #plt.show()
    
plot_fluxos(sol[3])

def plot_custo_total(DADOS_BARRAS, sol):
    # Extrair custos e geração por barra
    custos = DADOS_BARRAS['CUSTO($/MWh)'].values
    horas = len(next(iter(sol.values())))  # Número de horas igual para todas as barras
    custo_total_por_hora = np.zeros(horas)
    
    # Calcular custo total por hora
    for barra in sol:
        geracao_barra = np.array(sol[barra])
        custo_barra = custos[int(barra) - 1]  # Ajustando o índice da barra para o índice do Python
        custo_total_por_hora += geracao_barra * custo_barra

    # Plotar o gráfico de custo total
    plt.figure(figsize=(10, 5))
    plt.plot(range(1, horas+1), custo_total_por_hora, '-o')
    plt.title('Custo Total de Geração por Hora')
    plt.xlabel('Hora')
    plt.ylabel('Custo Total ($)')
    plt.grid(True)
    
plot_custo_total(DADOS_BARRAS, sol[0])
plt.show()