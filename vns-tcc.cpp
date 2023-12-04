/*
 * Descrição: Código sequencial para solução do problema da mochila multidimensional
 *            utilizando o algoritmo VNS (Variable Neighborhood Search).
 * Programador: Klelber Dias Januário
 * Orientadora: Bianca de Almeida Dantas
 * Data de Modificação: 15/11/2023
 * Algoritmo final para o TCC

	ALGORITMO:
    01: Encontra uma solução inicial S
    02: para k = 1 até maxIterações faça
    04:      Gere um vizinho qualquer S’ ∈ N(k)(S)
    05:      s’’ ← Faz uma busca local em N(k)(S’)
    06:
    07:      Se f (S”) > f (S) então
    08:          S ← S”;  
    09:          k ← 1;
    10:      senão 
    11:          k ← k + 1;
 */


#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <limits.h>
#include <queue>
#include <utility> 

using namespace std; 


#define EXECS 30


// DEFINIÇÃO DA ESTRUTURA DA MOCHILA
struct Mochila
{
	int qtdObjs;	 						// Quantidade de objetos
	int numComparts; 						// Numero de compartimentos da mochila
	float* valores;	 						// Beneficio associado a escolha de um objeto
	
	// Onde sera armazenado os pesos referente o quanto cada item gasta (matriz em forma de vetor)
	int* pesosItens;

	// Capacidade limite de cada compartimento da mochila
	int* capacidadeCompartimentoMochila;
};

struct Solucao
{
	float valor;
	int* itensLevados;
	int* rc; // Quantidade que os itens levados estão gastando da mochila
};


// PROTÓTIPO DE FUNÇÕES E VARIÁVEIS
void geraPopulacaoInicialGrasp(int*, int*, float*, float);
void calculoPseudoUtilidades(int*, float*, int*);
bool consegueAdicionar(int, int*, int*);
int pegarIndexMaxPsUtilidades(float*);
int pegarIndexMinPsUtilidades(float*);

/* Gera dois vizinhos baseado nas trocas */
void geraVizinho(int);
void n1();
void n2();
void n3();
void n4();
void n5();

/* Busca e funções que utilizamos dentro dela */
void buscaLocal();
void atualizarRC(Solucao*);
bool ehPossivel(Solucao*);
void removeItem(Solucao*);
void adicionarItem(Solucao*, int);
void adicionaAleatorio(Solucao*);
void calculaFitness(Solucao*);
bool allTested(int*);

/* Funções para imprimir no arquivo e gerar tabelas */
void escreveMelhorResultado(FILE *, char*, double, int);
void escritaFinal(FILE*, char*, double, double, double, double, int, int);
void tabelaDados(char*);
void imprimeTabela(FILE*, char*, int, int);

Mochila m;
Solucao* s;
Solucao* s1;
Solucao* s2;
int maxIteracoes = 5;
float alpha;

int optSol;

int* vetSols;
float* vetGaps;
float* vetTimes;
double avgSol;
double mediumGap;
double mediumTime;
double dpSols;
double dpTimes;
double dpGaps;
int bestOverallSolValue;
int bestSolValue;
double bestOverallTime;
double bestOverallGap;
int numOpt;
int numBetter;
int bestIndex;

int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		printf("Uso:\n./vns.exe arq_teste alpha\n");

		return 0;
	}
	
	// Abertura do arquivo de entrada
	FILE* fileIn = NULL;
	FILE* fileOut = NULL;

	timeval start, finish;
	double totalTime;
	char nomeArq[40];

	// Lê nome do arquivo e a quantidade de iterações
	strcpy(nomeArq, argv[1]);
	alpha = atof(argv[2]);

	// Definindo a saída
	char fileOutName[50];
	sprintf(fileOutName, "%s_cb%s", nomeArq, ".out");

	// Abre arquivo de entrada
	fileIn = fopen(nomeArq, "r");

	// Inicia a leitura do arquivo de entrada
	fscanf(fileIn, "%d", &m.qtdObjs);
	fscanf(fileIn, "%d", &m.numComparts);

	fscanf(fileIn, "%d", &optSol);

	// Alocação do vetor de valores
	m.valores = (float*) malloc(m.qtdObjs * sizeof(float));

	// Leitura do vetor de valores
	for (int j = 0; j < m.qtdObjs; j++)
		fscanf(fileIn, "%f", &(m.valores[j]));

	// Alocação da matriz de pesos como um vetor
	m.pesosItens = (int*) malloc(m.numComparts * m.qtdObjs * sizeof(int));

	// Leitura da matriz de pesos
	for (int j = 0; j < m.numComparts; j++)
		for (int k = 0; k < m.qtdObjs; k++)
			fscanf(fileIn, "%d", &(m.pesosItens[j * m.qtdObjs + k]));

	// Alocação do vetor de restrições
	m.capacidadeCompartimentoMochila = (int*) malloc(m.numComparts * sizeof(int));

	// Leitura do vetor de restrições
	for (int j = 0; j < m.numComparts; j++)
		fscanf(fileIn, "%d", &(m.capacidadeCompartimentoMochila[j]));

	// Alocando estruturas e atribuindo valores
	vetSols = (int*) malloc(EXECS * sizeof(int));
	vetGaps = (float*) malloc(EXECS * sizeof(float));
	vetTimes = (float*) malloc(EXECS * sizeof(float));

	avgSol = 0.0;
	mediumTime = 0.0;
	mediumGap = 0.0;
	bestOverallSolValue = 0;
	numOpt = 0;


	srand(time(NULL));


	s = (Solucao*) malloc(sizeof(Solucao));
	s1 = (Solucao*) malloc(sizeof(Solucao));
	s2 = (Solucao*) malloc(sizeof(Solucao));


	s->itensLevados = (int*) malloc(m.qtdObjs * sizeof(int));
	s1->itensLevados = (int*) malloc(m.qtdObjs * sizeof(int));
	s2->itensLevados = (int*) malloc(m.qtdObjs * sizeof(int));


	s->rc = (int*) malloc(m.numComparts * sizeof(int));
	s1->rc = (int*) malloc(m.numComparts * sizeof(int));
	s2->rc = (int*) malloc(m.numComparts * sizeof(int));
	

	/*******************************************************************
	 * Laço de execuções                                               *
	 ******************************************************************/
	for (int exec = 0; exec < EXECS; exec++)
	{
		// Inicio da tomada de tempo
		gettimeofday(&start, NULL);

	  	// Zerando as posições dos vetores
	  	for (int p = 0; p < m.qtdObjs; p++) {
			s->itensLevados[p] = 0;
			s1->itensLevados[p] = 0;
			s2->itensLevados[p] = 0;
		}

		for (int p = 0; p < m.numComparts; p++) {
			s->rc[p] = 0;
			s1->rc[p] = 0;
			s2->rc[p] = 0;
		}

		s->valor = 0;
		s1->valor = 0;
		s2->valor = 0;	


		geraPopulacaoInicialGrasp(s->itensLevados, s->rc, &s->valor, alpha);


		// Laço que controla o número de iterações
		for (int k = 1; k <= maxIteracoes;)
		{
			geraVizinho(k);
			
			buscaLocal();

			if (s1->valor > s->valor || s2->valor > s->valor)
			{
				k = 1;

				if(s1->valor > s2->valor) {
					// s = sDuasLinhas1;
					for (int i = 0; i < m.qtdObjs; i++)
						s->itensLevados[i] = s1->itensLevados[i];

					for (int j = 0; j < m.numComparts; j++)
						s->rc[j] = s1->rc[j];

					s->valor = s1->valor;
				} else {
					// s = sDuasLinhas2;
					for (int i = 0; i < m.qtdObjs; i++)
						s->itensLevados[i] = s2->itensLevados[i];

					for (int j = 0; j < m.numComparts; j++)
						s->rc[j] = s2->rc[j];

					s->valor = s2->valor;
				}
			}
			else
				k = k + 1;
		}
		
		// Fim da tomada de tempo
		gettimeofday(&finish, NULL);

		long start_usecs, finish_usecs, diff_usecs;

		finish_usecs = finish.tv_usec + (1000000 * finish.tv_sec);
		start_usecs = start.tv_usec + (1000000 * start.tv_sec);
		diff_usecs = finish_usecs - start_usecs;

		totalTime = (double)(diff_usecs) / 1000000;

		// Escreve o resultado final da execução
		bestSolValue = s->valor;

		vetGaps[exec] = (float)(optSol - bestSolValue) / optSol * 100;
		vetTimes[exec] = totalTime;
		vetSols[exec] = bestSolValue;

		escreveMelhorResultado(fileOut, fileOutName, totalTime, exec);

		avgSol += bestSolValue;
		mediumTime += totalTime;
		mediumGap += vetGaps[exec];

		if (bestSolValue >= bestOverallSolValue)
		{
			bestOverallSolValue = bestSolValue;
			bestOverallTime = totalTime;
		}

		if (vetGaps[exec] == 0.0)
			numOpt++;
		else if (vetGaps[exec] < 0.0)
			numBetter++;

	} // Fim do laço de execuções

	// Calcular média e desvio padrão das execuções
	avgSol /= EXECS;
	mediumTime /= EXECS;
	mediumGap /= EXECS;

	dpSols = 0.0;
	dpTimes = 0.0;
	dpGaps = 0.0;

	for (int i = 0; i < EXECS; i++)
	{
		dpSols += pow((vetSols[i] - avgSol), 2);
		dpTimes += pow((vetTimes[i] - mediumTime), 2);
		dpGaps += pow((vetGaps[i] - mediumGap), 2);
	}

	dpSols /= EXECS;
	dpTimes /= EXECS;
	dpGaps /= EXECS;

	dpSols = sqrt(dpSols);
	dpTimes = sqrt(dpTimes);
	dpGaps = sqrt(dpGaps);

	escritaFinal(fileOut, fileOutName, mediumTime, mediumGap, dpTimes,
				 dpGaps, bestOverallSolValue, numOpt);

	// Insere os dados na tabela que foi criada no arquivo testes.cpp
	tabelaDados(nomeArq);

	// Libera a memória utilizada
	free(m.valores);
	free(m.pesosItens);
	free(m.capacidadeCompartimentoMochila);

	free(vetSols);
	free(vetGaps);
	free(vetTimes);

	free(s->itensLevados);
	free(s1->itensLevados);
	free(s2->itensLevados);

	free(s->rc);
	free(s1->rc);
	free(s2->rc);

	free(s);
	free(s1);
	free(s2);

	// Fecha o arquivo de entrada
	fclose(fileIn);

	return 0;
}


// Gera o vetor com a população inicial utilizando o mesmo algoritmo do GRASP
void geraPopulacaoInicialGrasp(int* curSol, int* curRc, float* valor, float alpha) {
	int i, count;
	float rclFloor;

	// Construir a lista RCL explicitamente
	// Nesta implementação, ela contém os índices dos elementos
	int* rclIndexes = (int*) malloc(m.qtdObjs * sizeof(int));
	float* curPsUt = (float*) malloc(m.qtdObjs * sizeof(float));

	bool full = false;
	

	while (!full)
	{
		//Calcula as pseudoutilidades dos itens
		calculoPseudoUtilidades(curSol, curPsUt, curRc);

		// Baseado no valor de alpha considera a RCL composta pelos 
		// alpha * m.numItems primeiros elementos do vetor auxiliar
		int max = pegarIndexMaxPsUtilidades(curPsUt);
		int min = pegarIndexMinPsUtilidades(curPsUt);	

		// QUANTO MENOR O ALPHA, MAIS RIGOSO FICA O TETO
		rclFloor = curPsUt[max] - alpha * (curPsUt[max] - curPsUt[min]);


		// Preencher a RCL com os índices dos candidatos
		i = 0;
		count = 0;
		while (i < m.qtdObjs)
		{
			//Não está na solução
			if (curSol[i] == 0 && curPsUt[i] >= rclFloor)
				rclIndexes[count++] = i;

			i++;
		}

		int actualRclSize = count;


		//Gera um número aleatório para escolher um dos elementos contidos na RCL
		int index = rand() % actualRclSize;
		index = rclIndexes[index];

		// Tenta adicionar o item escolhido, verificando se ele não torna a mochila inviável
		// Nessa implementação, se o primeiro não puder ser colocado,
		// considera-se a solução como completada
		if (consegueAdicionar(index, curSol, curRc))
		{
			curSol[index] = 1;

			// Atualiza o vetor de recursos disponíveis
			for (int k = 0; k < m.numComparts; k++)
				curRc[k] += m.pesosItens[k * m.qtdObjs + index];

			// Atualiza o valor
			*valor += m.valores[index];
		}
		else
			full = true;
	}
	
	free(rclIndexes);
	free(curPsUt);
}


// Faz o cálculo de Pseudo-Utilidades
void calculoPseudoUtilidades(int* curSol, float* curPsUt, int* curRc)
{
	for (int i = 0; i < m.qtdObjs; i++)
	{
		if (curSol[i] == 0)
		{
			float sum = 0.0;
			// curRc[j] = o quanto já foi usado (como é início, zero)
			for (int j = 0; j < m.numComparts; j++) {
				sum += (float) m.pesosItens[j * m.qtdObjs + i] / (m.capacidadeCompartimentoMochila[j] - curRc[j]);

				printf("CUR <- %d\n", curRc[j]);
			}

			curPsUt[i] = (float) m.valores[i] / sum;
		}
		else
			curPsUt[i] = 0.0;
	}
}


// Pega o índice do vetor que possui a MAIOR pseudoutilidade
int pegarIndexMaxPsUtilidades(float* psUt)
{
	int i;

	int iMax = 0;

	for (i = 1; i < m.qtdObjs; i++)
	{
		if (psUt[i] != 0 && psUt[i] > psUt[iMax])
			iMax = i;
	}

	return iMax;
}


// Pega o índice do vetor que possui a MENOR pseudoutilidade
int pegarIndexMinPsUtilidades(float* psUt)
{
	int i;

	int iMin = 0;

	for (i = 1; i < m.qtdObjs; i++)
	{
		if (psUt[i] != 0 && psUt[i] < psUt[iMin])
			iMin = i;
	}

	return iMin;
}


// Função de adicionar no vetor, UTILIZADO APENAS NA INICIALIZAÇÃO DO GRASP
bool consegueAdicionar(int index, int* curSol, int* curRc)
{
	if (curSol[index] == 1)
		return false;
	
	for (int i = 0; i < m.numComparts; i++)
	{
		if (m.capacidadeCompartimentoMochila[i] - (curRc[i] + m.pesosItens[i * m.qtdObjs + index]) < 0)
			return false;
	}

	return true;
}


// Encontra um vizinho qualquer de S
void geraVizinho(int vizinhanca) {
	if(vizinhanca == 1)
		n1();
	else if(vizinhanca == 2)
		n2();
	else if(vizinhanca == 3)
		n3();
	else if(vizinhanca == 4)
		n4();
	else if(vizinhanca == 5)
		n5();
}


// Troca Simples
void n1() {
	// s1 = s;
	for (int i = 0; i < m.qtdObjs; i++) {
		s1->itensLevados[i] = s->itensLevados[i];
	}
	for (int j = 0; j < m.numComparts; j++) {
		s1->rc[j] = s->rc[j];
	}
	s1->valor = s->valor;

	// Realiza a troca e se não puder, tenta mais uma vez
	int* valores0 = (int*) malloc(m.qtdObjs * sizeof(int));
	int* valores1 = (int*) malloc(m.qtdObjs * sizeof(int));

	for(int x = 0; x < 2; x++) {
		// Verificando os valores e guardando suas posições
		int qtdZe = 0;
		int qtdUm = 0;

		for(int p = 0; p < m.qtdObjs; p++) {
			if(s1->itensLevados[p] == 1)
				valores1[qtdUm++] = p;
			else
				valores0[qtdZe++] = p;
		}

		// Sorteando i e j
		int i = rand() % qtdZe;
		int j = rand() % qtdUm;

		// Retirando o a posição que o índice j guarda
		s1->itensLevados[valores1[j]] = 0;
		s1->valor -= m.valores[valores1[j]];

		for (int l = 0; l < m.numComparts; l++)
			s1->rc[l] -= m.pesosItens[l * m.qtdObjs + valores1[j]];

		// Adicionando o item (obs.: função já adiciona RC e incrementa s1->valor)
		adicionarItem(s1, valores0[i]);
		
		// Conseguiu adicionar
		if(s->itensLevados[valores0[i]] == 1)
			break;
	}

	free(valores0);
	free(valores1);
}


// Troca Excessiva
void n2() {
	// s1 = s;
	for (int i = 0; i < m.qtdObjs; i++) {
		s1->itensLevados[i] = s->itensLevados[i];
	}
	for (int j = 0; j < m.numComparts; j++) {
		s1->rc[j] = s->rc[j];
	}
	s1->valor = s->valor;

	// Avalia
	int melhorFitness = 0;
	int melhorIndice1 = -1;
	int melhorIndice2 = -1;
	int possui1;
	int possui0;

	for(int i = 0; i < m.qtdObjs; i++) 
	{
		for(int j = i + 1; j < m.qtdObjs; j++) 
		{
			// Trocar essas duas posicoes e verificar se é possível
			if(s1->itensLevados[i] != s1->itensLevados[j]) {
				// Preciso verificar qual dos dois possui 1
				if(s1->itensLevados[i] == 1) {
					possui1 = i;
					possui0 = j;
				} else {
					possui1 = j;
					possui0 = i;
				}

				// Removendo o item
				s1->itensLevados[possui1] = 0;
				s1->valor -= m.valores[possui1];

				for (int l = 0; l < m.numComparts; l++)
					s1->rc[l] -= m.pesosItens[l * m.qtdObjs + possui1];

				// Adicionando o item (obs.: função já adiciona RC e incrementa s1->valor)
				adicionarItem(s1, possui0);

				if(s1->valor > melhorFitness) {
					melhorFitness = s1->valor;
					melhorIndice1 = i;
					melhorIndice2 = j;
				}

				// SE EU INSERI, desfaço as modificações
				if(s1->itensLevados[possui0] == 1)
				{
					s1->itensLevados[possui0] = 0;
				}

				s1->itensLevados[possui1] = 1;
				s1->valor = s->valor;

				for (int l = 0; l < m.numComparts; l++)
					s1->rc[l] = s->rc[l];
			}
		}
	}

	if(melhorIndice1 != -1 && melhorIndice2 != -1) {
		adicionarItem(s1, melhorIndice1);
		adicionarItem(s1, melhorIndice2);
	}
}


// Path Relinking 
void n3() {
	Solucao* newSol = (Solucao*) malloc(sizeof(Solucao));
	newSol->itensLevados = (int*) malloc(m.qtdObjs * sizeof(int));
	newSol->rc = (int*) malloc(m.numComparts * sizeof(int));

	// Criando uma solução inicial e melhorando um pouco (START)
	for (int i = 0; i < m.qtdObjs; i++)
		newSol->itensLevados[i] = rand() % 2;

	atualizarRC(newSol);

	while (!ehPossivel(newSol))
		removeItem(newSol);

	adicionaAleatorio(newSol);
	calculaFitness(newSol);

	// Criando o mapeameento dos itens que faltam adicionar em START para virar s
	int* itensEm1 = (int*) malloc(m.qtdObjs * sizeof(int));
	int qtd = 0;

	for(int p = 0; p < m.qtdObjs; p++) {
		if(s->itensLevados[p] == 1 && newSol->itensLevados[p] == 0) {
			itensEm1[qtd++] = p;
		}
	}

	// primeiro valor = valor do item
	// segundo valor = índice do item
	priority_queue<pair<int, int>> mapRelation;

	for(int k = 0; k < qtd; k++) {
		// Adiciona na lista de prioridade o valor do item e seu índice
		mapRelation.push(make_pair(m.valores[itensEm1[k]], itensEm1[k]));
	}

	// Enquanto estiver elemento
	while (!mapRelation.empty()) {
        adicionarItem(newSol, mapRelation.top().second);

        mapRelation.pop();
    }

	// Transfere de newSol para s1
	for (int k = 0; k < m.qtdObjs; k++)
		s1->itensLevados[k] = newSol->itensLevados[k];
	
	for (int k = 0; k < m.numComparts; k++)
		s1->rc[k] = newSol->rc[k];

	s1->valor = newSol->valor;

	free(newSol->itensLevados);
	free(newSol->rc);
	free(newSol);
	free(itensEm1);
}


// Rotação Circular
void n4() {
	int k = 4;

	for (int i = 0; i < m.qtdObjs; i++)
		s1->itensLevados[(i + k) % m.qtdObjs] = s->itensLevados[i];

	// Adicionando os valores dos itens no rc
	atualizarRC(s1);

	// Com a mudança no vetor, não podemos garantir que é uma solução viavél
	while (!ehPossivel(s1))
		removeItem(s1);

	adicionaAleatorio(s1);
	calculaFitness(s1);
}


// Empty Bucket
void n5() {
	// s1 = s;
	for (int i = 0; i < m.qtdObjs; i++) {
		s1->itensLevados[i] = s->itensLevados[i];
	}
	for (int j = 0; j < m.numComparts; j++) {
		s1->rc[j] = s->rc[j];
	}
	s1->valor = s->valor;

	int tamBucket = 2;

	// Definindo o tamanho do balde
	for (int i = 1; i <= 4; i++) {
		if (m.qtdObjs % i == 0 && i > tamBucket)
			tamBucket = i;
    }

	// Analisando a quantidade de itens adicionados em cada balde
	priority_queue<pair<int, int>> qtdItensBucket;

	for (int bucket = 0; bucket < m.qtdObjs / tamBucket; bucket++) {
		int qtd = 0;

		// Percorrendo os itens do balde
		for(int j = bucket * tamBucket; j < m.qtdObjs && j < ((bucket + 1) * tamBucket); j++) {
			if(s1->itensLevados[j] == 1)
				qtd++;
		}	

		// Adicionando a quantidade (primeiro valor = qtd de itens & segundo valor = índice do balde)
		qtdItensBucket.push(make_pair(qtd, bucket));
    }
	
	// printf("Tamanho de cada balde: %d\n", tamBucket);
	// printf("Total Baldes: %d\n", (int) qtdItensBucket.size());
	// printf("Itens no Balde Maximo: %d\n", qtdItensBucket.top().first);
	// printf("Indice do Balde Maximo: %d\n", qtdItensBucket.top().second);

	if (!qtdItensBucket.empty()) {
		// Escolher um dos máximos
		int* maximalBuckets = (int*) malloc(sizeof(int) * qtdItensBucket.size());
		int qtdMax = 0;

		int maior = qtdItensBucket.top().first;

		// Extrair os valores iguais ao primeiro
		while (!qtdItensBucket.empty() && qtdItensBucket.top().first == maior) {
			maximalBuckets[qtdMax++] = qtdItensBucket.top().second;

			qtdItensBucket.pop();
		}
		
		// Se existir baldes máximos
		if(qtdMax > 0) {
			int bucketSortedIndex = rand() % qtdMax;
			
			// printf("BALDE SORTEADO: %d\n", bucketSorted);
			// printf("=================================\n\n");

			// Percorrendo os itens do balde para retirá-los
			int bucketIndex = maximalBuckets[bucketSortedIndex];

			for(int j = bucketIndex * tamBucket; j < m.qtdObjs && j < (bucketIndex + 1) * tamBucket; j++) {
				// Necessário verificar, pois nem todos os itens do balde precisam estar na mochila
				if(s1->itensLevados[j] == 1) {

					s1->valor -= m.valores[j];
					s1->itensLevados[j] = 0;
					
					// Retirando itens do RC
					for (int l = 0; l < m.numComparts; l++)
						s1->rc[l] -= m.pesosItens[l * m.qtdObjs + j];
				}
			}

		}

		free(maximalBuckets);
	}
}


// Faz uma busca local para tentar encontrar um máximo "local"
void buscaLocal() {
	double temp = 500;
	double finalTemp = 0.00001;
	int maxChainLength = m.qtdObjs * 10;
	int factor = 0.85;

	// solucao = s1
	// melhorSolucao = s2
	// novaSolucao = newSol
	
	Solucao* newSol = (Solucao*) malloc(sizeof(Solucao));
	newSol->itensLevados = (int*) malloc(m.qtdObjs * sizeof(int));
	newSol->rc = (int*) malloc(m.numComparts * sizeof(int));


	for (int k = 0; k < m.qtdObjs; k++)
		s2->itensLevados[k] = s1->itensLevados[k];
	
	for (int k = 0; k < m.numComparts; k++)
		s2->rc[k] = s1->rc[k];

	s2->valor = s1->valor;
	
	
	/*******************************************
	 * Laço principal do SA.
	 * Executado até alcançar a temperatura final
	 *******************************************/
	while (temp >= finalTemp)
	{	
		for (int chainLength = 0; chainLength < maxChainLength; chainLength++)
		{
			for (int k = 0; k < m.qtdObjs; k++)
				newSol->itensLevados[k] = s1->itensLevados[k];
				
			for (int k = 0; k < m.numComparts; k++)
				newSol->rc[k] = s1->rc[k];

			newSol->valor = s1->valor;
			
			//Escolhe um item aleatoriamente
			int index = rand() % m.qtdObjs;
			
			//Se o item ainda não estiver na solução, é adicionado
			if (newSol->itensLevados[index] == 0)
			{
				newSol->itensLevados[index] = 1;
				
				//Atualizar rc para refletir a adição do novo item
				for (int l = 0; l < m.numComparts; l++)
					newSol->rc[l] += m.pesosItens[l * m.qtdObjs + index];
				
				newSol->valor += m.valores[index];
				

				int drop;
				
				// Deixa a solução viável
				while (!ehPossivel(newSol))
				{						
					// Laço para garantir que o item retirado é diferente do que acabamos de adicionar e que ele realmente está na mochila	
					// E SE NÃO TIVER NINGUÉM NA MOCHILA???	NÃO TEM COMO, SE NÃO A MOCHILA SERIA VIÁVEL				
					do {
						drop = rand() % m.qtdObjs;						
					} while(drop == index || newSol->itensLevados[drop] == 0);
					
					newSol->itensLevados[drop] = 0;
					
					// Atualiza rc para refletir a retirada do item
					for (int l = 0; l < m.numComparts; l++)
						newSol->rc[l] -= m.pesosItens[l * m.qtdObjs + drop];

					newSol->valor -= m.valores[drop];
				}
				

				float delta = newSol->valor - s1->valor;
				

				if (delta >= 0 || ((double) rand() / (RAND_MAX)) < ((double) exp((double) delta / temp)))
				{
					// Atualiza solução
					for (int k = 0; k < m.qtdObjs; k++)
						s1->itensLevados[k] = newSol->itensLevados[k];
						
					for (int k = 0; k < m.numComparts; k++)
						s1->rc[k] = newSol->rc[k];
						
					s1->valor = newSol->valor;
				}
			}
			else // Se o item já estava na solução, ele é retirado
			{
				newSol->itensLevados[index] = 0;
				
				// Atualizar rc para refletir na remoção do item
				for (int l = 0; l < m.numComparts; l++)
					newSol->rc[l] -= m.pesosItens[l * m.qtdObjs + index];
				
				newSol->valor -= m.valores[index];
				
				// Outro item é adicionado aleatoriamente
				int* tested = (int*) malloc(m.qtdObjs * sizeof(int));
				int add;
				bool found = false;
				
				// Código potencialmente problemático!
				// E se ninguém puder ser adicionado?
				// Pode ser necessário estabelecer um limite para a busca
				bool goon = true;					
				
				for (int k = 0; k < m.qtdObjs; k++)
					tested[k] = 0;
					
				tested[index] = 1;
				
				while (!found && goon)
				{
					// E SE TODOS OS ITENS JÁ ESTIVEREM NA MOCHILA?
					do
					{
						add = rand() % m.qtdObjs;
						tested[add] = 1;

					} while(add == index || newSol->itensLevados[add] == 1);

					// Tenta adicionar
					adicionarItem(newSol, add);

					// Realmente adicionou
					if(newSol->itensLevados[add] == 1)
						found = true;										
					else
					{
						// Testa se todos os elementos do vetor já estão marcados com 1
						if (allTested(tested))
							goon = false;
					}
				}

				// Libera vetor de teste
				free(tested);


				float delta = newSol->valor - s1->valor;
				

				if (delta >=0 || ((double) rand() / (RAND_MAX)) < ((double) exp((double) delta / temp)))
				{
					// Atualiza solução
					for (int k = 0; k < m.qtdObjs; k++)
						s1->itensLevados[k] = newSol->itensLevados[k];
						
					for (int k = 0; k < m.numComparts; k++)
						s1->rc[k] = newSol->rc[k];
						
					s1->valor = newSol->valor;
				}				
			}
			// Fim if-else
			
			if (s1->valor > s2->valor)
			{
				for (int k = 0; k < m.qtdObjs; k++)
					s2->itensLevados[k] = s1->itensLevados[k];
					
				for (int k = 0; k < m.numComparts; k++)
					s2->rc[k] = s1->rc[k];
					
				s2->valor = s1->valor;
			}				
		}
		// Fim do for da cadeia de Markov
		
		// Atualiza a temperatura de acordo com o fator de controle de temperatura
		temp = factor * temp;			
	}
	// Fim do laço principal do SA

	free(newSol->itensLevados);
	free(newSol->rc);
	free(newSol);
}


// Reseta e adiciona os valores dos itens da mochila
void atualizarRC(Solucao* p) {
	// Zerando o vetor p->rc
	for (int i = 0; i < m.numComparts; i++)
		p->rc[i] = 0;

	// Para cada item levado, somar
	for (int i = 0; i < m.qtdObjs; i++)
	{
		if (p->itensLevados[i] == 1)
		{
			/****************************************************************/
			// Vai atribuindo para cada compartimento o quanto cada item ocupa
			// da restrição
			/****************************************************************/
			for (int j = 0; j < m.numComparts; j++)
				p->rc[j] += m.pesosItens[j * m.qtdObjs + i];
		}
	}
}


// Verifica se a solução é possível
bool ehPossivel(Solucao* p) {
	//Para cada restrição
	for (int i = 0; i < m.numComparts; i++)
		if (m.capacidadeCompartimentoMochila[i] - p->rc[i] < 0)
			return false;

	return true;
}


// Remove um item da solucao
void removeItem(Solucao* k) {
	// Vetor que contém as posições preenchidas
	int posicoes[m.qtdObjs];
	int tamanho = 0;

	// Percorre os indices, verificando as posicoes que estão ocupadas
	for(int i = 0; i < m.qtdObjs; i++) {
		if(k->itensLevados[i] == 1) {
			posicoes[tamanho] = i;

			tamanho++;
		}
	}

	
	if (tamanho != 0)
	{
		int i = rand() % tamanho;
		i = posicoes[i];
		
		k->itensLevados[i] = 0;

		for (int j = 0; j < m.numComparts; j++)
			k->rc[j] -= m.pesosItens[j * m.qtdObjs + i];

		k->valor -= m.valores[i];
	}
}


// Adiciona item i na solução
void adicionarItem(Solucao* k, int i) {
	if (k->itensLevados[i] != 1)
	{
		k->itensLevados[i] = 1;

		// ADICIONANDO item no vetor RC
		for (int j = 0; j < m.numComparts; j++)
			k->rc[j] += m.pesosItens[j * m.qtdObjs + i];

		k->itensLevados[i] = (ehPossivel(k)) ? 1 : 0;

		// ITEM NÃO FOI ADICIONADO
		if (k->itensLevados[i] == 0) {
			// RETIRANDO item no vetor RC
			for (int j = 0; j < m.numComparts; j++)
				k->rc[j] -= m.pesosItens[j * m.qtdObjs + i];
		} else
			k->valor += m.valores[i];
	}
}


// "Tenta" adicionar dois itens para tentar melhorar a solucao
void adicionaAleatorio(Solucao* k) {
	// Vetor que contém as posições disponíveis
	int posicoes[m.qtdObjs];
	int tamanho = 0;

	// Percorre os indices, verificando as posicoes que estão livre
	for(int i = 0; i < m.qtdObjs; i++) {
		if(k->itensLevados[i] == 0) {
			posicoes[tamanho] = i;

			tamanho++;
		}
	}

	// Sorteia dois números e vai tentando adicionar um por vez
	int indice1;
	int indice2;

	indice1 = rand() % tamanho;

	do {
		indice2 = rand() % tamanho;
	} while(indice1 == indice2);

	adicionarItem(k, posicoes[indice1]);
	adicionarItem(k, posicoes[indice2]);
}


// Calcula Fitness = soma dos valores
void calculaFitness(Solucao* k) {
	k->valor = 0;

	for (int i = 0; i < m.qtdObjs; i++)
	{
		k->valor += m.valores[i] * k->itensLevados[i];
	}
}


// Testa se todos os elementos do vetor já estão marcados com 1
bool allTested(int* tested)
{
	for (int i = 0; i < m.qtdObjs; i++)
	{
		if (tested[i] == 0)
			return false;
	}
			
	return true;
}


// Função que escreve o melhor resultado de todas as iterações de uma execução
void escreveMelhorResultado(FILE *f, char *nome, double time, int exec) {
	FILE* fGaps;
	FILE* fSols;
	FILE* fTimes;

	f = fopen(nome, "a+");
	fGaps = fopen("seq_allgaps_nopr.out", "a+");
	fSols = fopen("seq_allsols.out", "a+");
	fTimes = fopen("seq_alltimes_nopr.out", "a+");

	float bestGap = (float)(optSol - bestSolValue) / optSol * 100;

	fprintf(f, "Execução: %d\nTempo Total: %.4f\n Valor da Melhor Solução: %d",
			exec, time, bestSolValue);
	fprintf(f, "\nNúmero de gerações = %d\nSolução ótima = %d\nMelhor solução = %d\nMelhor índice = %d\nGap: %f%%\n",
			maxIteracoes, optSol, bestSolValue, bestIndex, bestGap);

	fprintf(fGaps, "%.4f\n", bestGap);
	fprintf(fSols, "%d\n", bestSolValue);
	fprintf(fTimes, "%.4f\n", time);

	fprintf(f, "Solução:\n");

	for (int i = 0; i < m.qtdObjs; i++)
		fprintf(f, "%d", s->itensLevados[i]);

	fprintf(f, "\n\n");

	fclose(f);
	fclose(fGaps);
	fclose(fSols);
	fclose(fTimes);
}


// Função para escrever no arquivo de saída
void escritaFinal(FILE* f, char* nome, double mediumTime, double mediumGap,
				  double dpTimes, double dpGaps,
				  int bestOverallSolValue, int numOpt)
{
	FILE* f2;
	double bestGap = (float)(optSol - bestOverallSolValue) / optSol * 100;

	f = fopen(nome, "a+");
	f2 = fopen("seqdados_ag.out", "a+");

	fprintf(f, "Execuções: %d \nTempo Medio: %.4f \nSolução Média: %.4lf \nGap Médio: %.4lf \nDesvio-Padrão (Soluções): %.4lf\nDesvio-Padrão (Gaps): %.4lf\nDesvio-Padrão (Tempos): %.4lf\nMelhor Solução: %d, Melhor Gap: %.4lf, Melhor Tempo: %.4lf\nÓtimas: %d\nMelhores: %d\n",
			EXECS, mediumTime, avgSol, mediumGap, dpSols, dpGaps, dpTimes,
			bestOverallSolValue, bestGap, bestOverallTime, numOpt, numBetter);

	fprintf(f2, "%d %d %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n", numOpt,
			bestOverallSolValue, avgSol, dpSols, bestGap, mediumGap,
			dpGaps, mediumTime, dpTimes);

	fclose(f);
	fclose(f2);
}


// Imprime os dados
void tabelaDados(char* arquivoIn) {
	FILE* arquivo = fopen("tabela.in", "a+");

	char nomeArquivoIn[20] = "";
	
	// ARQUIVO
	if(arquivoIn[0] == 'g') {
		strcat(nomeArquivoIn, "| GK");
	}
	else if(arquivoIn[0] == 'O') {
		strcat(nomeArquivoIn, "| OR");
	}
	else if(arquivoIn[0] == 'p') {
		if(arquivoIn[1] == 'b')
			strcat(nomeArquivoIn, "| PB");
		else
			strcat(nomeArquivoIn, "| PET");
	}
	else if(arquivoIn[0] == 's') {
		strcat(nomeArquivoIn, "| SENTO");
	}
	else if(arquivoIn[0] == 'w') {
		if(arquivoIn[3] == 'n')
			strcat(nomeArquivoIn, "| WEING");
		else
			strcat(nomeArquivoIn, "| WEISH");
	}
	else {
		strcat(nomeArquivoIn, "| HP");
	}

	imprimeTabela(arquivo, nomeArquivoIn, 11, 0);

	// M * N
	char dimensoes[20] = "";
	sprintf(dimensoes, "| %3d x %4d", m.numComparts, m.qtdObjs);

	imprimeTabela(arquivo, dimensoes, 13, 0);

	// ALPHA
	char alpha[20] = "";

	if(arquivoIn[0] == 'O') {
		int i = 0;

		while(arquivoIn[i] != '-')
			i++;
		
		sprintf(alpha, "| %c.%c%c", arquivoIn[i+1], arquivoIn[i+3], arquivoIn[i+4]);
	}
	else
		strcat(alpha, "| X");

	imprimeTabela(arquivo, alpha, 8, 0);

	// INSTANCIA
	char instancia[20] = "";
	int j = 0;

	while(arquivoIn[j] != '\0')
		j++;
		
	if(	arquivoIn[j-6] == '1' || arquivoIn[j-6] == '2' || arquivoIn[j-6] == '3')
		sprintf(instancia, "| %c%c", arquivoIn[j-6], arquivoIn[j-5]);
	else
		sprintf(instancia, "| 0%c", arquivoIn[j-5]);

	imprimeTabela(arquivo, instancia, 8, 0);

	// ÓTIMOS
	char otimos[20] = "";
	sprintf(otimos, "| %d", numOpt);

	imprimeTabela(arquivo, otimos, 9, 0);

	// GAP - MELHOR
	double bestGap = (float)(optSol - bestOverallSolValue) / optSol * 100;

	char gapMelhor[20] = "";
	sprintf(gapMelhor, "| %.4lf", bestGap);

	imprimeTabela(arquivo, gapMelhor, 13, 0);

	// GAP MÉDIO
	char gapMedio[20] = "";
	sprintf(gapMedio, "| %.4lf", mediumGap);

	imprimeTabela(arquivo, gapMedio, 12, 0);

	// GAP DESVIO
	char gapDesvio[20] = "";
	sprintf(gapDesvio, "| %.4lf", dpGaps);

	imprimeTabela(arquivo, gapDesvio, 13, 0);

	// TEMPO MEDIO
	char tempoMedio[20] = "";
	sprintf(tempoMedio, "| %.4lf", mediumTime);

	imprimeTabela(arquivo, tempoMedio, 14, 0);
	
	// TEMPO DESVIO
	char tempoDesvio[20] = "";
	sprintf(tempoDesvio, "| %.4lf", dpTimes);
	
	imprimeTabela(arquivo, tempoDesvio, 15, 1);

	
	fclose(arquivo);
}


// Função auxiliar para imprimir
void imprimeTabela(FILE* arquivo, char* p, int k, int terOuNao) {
    int cont = 0;
    
    while(p[cont] != '\0')
        cont++;

    fprintf(arquivo, "%s", p);

    for(int i = cont; i < k; i++) {
        fprintf(arquivo, " ");
    }
    
    if(terOuNao == 1)
		fprintf(arquivo, "|\n");
}
