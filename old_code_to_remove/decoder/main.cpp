/*
 * BRKGA_SNW.cpp
 *
 *  Created on: 28 de mai de 2016
 *      Author: guest-2FMpQr
 */

#include "BRKGA.h"
#include "MTRand.h"
#include "SampleDecoder.h"
#include <float.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits.h>
#include <new> //bad alloc do comando new
#include <sstream>
#include <string>
#include <vector>

#include "Utility.h"

using namespace std;

int main(int argc, char *argv[]) {

    const unsigned p = 100;   // size of population //***Antes era 1000 o valor desse parâmetro
    const double pe = 0.25;   // fraction of population to be the elite-set
    const double pm = 0.05;   // fraction of population to be replaced by mutants
    const double rhoe = 0.70; // probability that offspring inherit an allele from elite parent
    const unsigned K =
        1; // number of independent populations         //***Antes era 3 o valor desse parâmetro
    const unsigned MAXT =
        1; // number of threads for parallel decoding   //***Antes era 3 o valor desse parâmetro

    const unsigned X_INTVL = 100;   // exchange best individuals at every 100 generations
    const unsigned X_NUMBER = 2;    // exchange top 2 best
    const unsigned MAX_GENS = 1000; // run for 1000 gens

    ifstream listaConjuntos, filaNomeArquivos;
    ofstream dataFileWrite, dataFileSolutions, dataFileBestSolutions;
    string Buffer, fileName, fileNameSolutions;
    clock_t TempoFO_StarInic; // TempoInicial
    vector<double> auxConnections;
    vector<int> dataFeasibleSolutionsPoplt;
    vector<pair<double, int>> dataFeasibleBestSolutions;
    vector<string> listSetInstances, listFileInstances;
    double TempoExecTotal = 0.0, TempoFO_Star = 0.0, FO_Star = DBL_MAX, FO_Min = DBL_MIN;
    int bestGeneration = 0, minGeneration = 0;
    stringstream converterIntToString;

    vector<double> chromosomeSolution;

    int iterSemMelhora, iterMax = 10, quantIteracoes = 0, bestIteration = 0;
    double seed;

    Utility utilities;

    // Abertura do arquivo que contém a lista de cada conjunto de instâncias (cada conjunto possui
    // sua própria lista de instâncias) listaConjuntos.open("ListaConjuntos.txt");
    listaConjuntos.open("ListaConjuntos.txt");
    if (listaConjuntos.is_open()) {
        // cout<<"Operation succesfully performed: listaConjuntos!"<<endl;
    } else {
        cout << "Not was possible open the file: listaConjuntos!" << endl;
        listaConjuntos.clear();
        exit(1);
    }

    // Loop em que se é carregada a lista de conjuntos de instâncias
    listaConjuntos >> Buffer;
    do {
        listSetInstances.push_back(Buffer);
        listaConjuntos >> Buffer;
    } while (listaConjuntos.good());

    listaConjuntos.clear();
    listaConjuntos.close();
    Buffer.clear();

    // Execução dos teste por conjunto de instâncias   //listSetInstances.size()
    for (unsigned i = 0; i < listSetInstances.size(); i++) {

        fileName = fileName + listSetInstances[i] + ".txt";

        // Abertura da lista de instâncias do conjunto i
        filaNomeArquivos.open(fileName.c_str());
        if (filaNomeArquivos.is_open()) {
            // cout<<"Operation successfully performed: filaNomeArquivos"<<endl;
        } else {
            cout << "Não foi possivel abrir o arquivo!: " << fileName << endl;
            filaNomeArquivos.clear();
            exit(1);
        }

        fileName.clear();

        // Carrega uma lista com os nomes das instâncias
        filaNomeArquivos >> Buffer;
        do {
            listFileInstances.push_back(Buffer);
            filaNomeArquivos >> Buffer;
        } while (filaNomeArquivos.good());

        filaNomeArquivos.clear();
        filaNomeArquivos.close();

        // Abertura do arquivo que armazena os dados da solução de cada instância de um determinado
        // conjunto
        fileName = "Dados_Solucoes/" + listSetInstances[i] + ".txt";

        dataFileWrite.open(fileName.c_str());

        if (dataFileWrite.is_open()) {
            // cout<<"Operation successfully performed: dataFileWrite"<<endl;
        } else {
            cout << "Não foi possivel abrir o arquivo dataFileWrite!: " << fileName << endl;
            dataFileWrite.clear();
            exit(1);
        }

        fileName.clear();

        // Execução dos Testes para as instâncias do conjunto i atual   //listFileInstances.size()
        for (unsigned j = 10; j < listFileInstances.size(); j++) {

            fileName = "Instancias/" + listSetInstances[i] + "/";
            fileName = fileName + listFileInstances[j].c_str();

            //====== Inicialização do Decoder e do Gerador de Números Aleatórios
            SampleDecoder decoder(fileName.c_str()); // initialize the decoder

            seed = time(NULL);

            MTRand rng(seed);

            fileName.clear();

            fileNameSolutions =
                "Solucoes/" + listSetInstances[i] + "/Solucao_" + listFileInstances[j].c_str();

            dataFileSolutions.open(fileNameSolutions.c_str());

            if (dataFileSolutions.is_open()) {
                // cout<<"Operation successfully performed: dataFileWrite"<<endl;
            } else {
                cout << "Não foi possivel abrir o arquivo dataFileSolutions!: " << fileNameSolutions
                     << endl;
                dataFileSolutions.clear();
                exit(1);
            }

            fileNameSolutions.clear();
            fileNameSolutions = "Melhores_Solucoes/" + listSetInstances[i] + "/Dados_Tempo_" +
                                listFileInstances[j].c_str();

            dataFileBestSolutions.open(fileNameSolutions.c_str());

            if (dataFileBestSolutions.is_open()) {
                // cout<<"Operation successfully performed: dataFileWrite"<<endl;
            } else {
                cout << "Não foi possivel abrir o arquivo dataFileBestSolutions!: "
                     << fileNameSolutions << endl;
                dataFileBestSolutions.clear();
                exit(1);
            }

            const unsigned n = decoder.getQuantConnections() * 2; // size of chromosomes

            // initialize the BRKGA-based heuristic
            BRKGA<SampleDecoder, MTRand> algorithm(n, p, pe, pm, rhoe, decoder, rng, K, MAXT);

            TempoFO_StarInic = clock();
            decoder.setInitialTime();

            FO_Star = DBL_MAX;
            FO_Min = DBL_MAX * -1;
            bestGeneration = 0;
            minGeneration = 0;

            iterSemMelhora = 0;

            iterMax = 1;

            quantIteracoes = 0;
            bestIteration = 0;

            unsigned generation = 0; // current generation
            do {
                algorithm.evolve(); // evolve the population for one generation

                if ((++generation) % X_INTVL == 0) {
                    algorithm.exchangeElite(X_NUMBER); // exchange top individuals
                }

                if (algorithm.getBestFitness() < FO_Star) {
                    TempoFO_Star = (((double)(clock() - TempoFO_StarInic)) / CLOCKS_PER_SEC);
                    FO_Star = algorithm.getBestFitness();
                    bestGeneration = generation;

                    dataFileBestSolutions
                        << " " << setprecision(12) << (algorithm.getBestFitness() / 10.0) * -1
                        << "   " << setprecision(12) << TempoFO_Star << "  "
                        << algorithm.getBestFitness() << "  " << quantIteracoes << "  "
                        << "\n";

                    bestIteration = quantIteracoes;
                }
                quantIteracoes++;
            } while ((((double)(clock() - TempoFO_StarInic)) / CLOCKS_PER_SEC) < 600);

            TempoExecTotal = (((double)(clock() - TempoFO_StarInic)) / CLOCKS_PER_SEC);

            dataFileWrite << " " << listFileInstances[j].c_str() << "   " << setprecision(12)
                          << (FO_Star / 10.0) * -1 << "   " << setprecision(12) << TempoExecTotal
                          << "  " << algorithm.getBestFitness() << "  " << quantIteracoes << "  "
                          << bestIteration << "  " << TempoFO_Star << "  " << seed << "\n";

            decoder.solutionDecode(algorithm.getBestChromosome(), dataFileSolutions);

            fileNameSolutions.clear();

            dataFileSolutions.clear();
            dataFileSolutions.close();

            dataFileBestSolutions.clear();
            dataFileBestSolutions.close();

            decoder.cleanMemory();
        }

        listFileInstances.clear();
        dataFileWrite.clear();
        dataFileWrite.close();
    }

    listSetInstances.clear();
    return 0;
}
