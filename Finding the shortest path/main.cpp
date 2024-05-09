#include <bits/stdc++.h>

using namespace std;

#define isz(x) (int)x.size()
#define all(x) x.begin(),x.end()
#define rall(x) x.rbegin(),x.rend()
#define fr first
#define sc second

vector < vector <int> > pathsLists;
void precalc(int vertexsCount){
    int fact = 1;
    vector <int> curPath;
    for (int v = 1;v <= vertexsCount;v++){
        curPath.push_back(v);
        fact *= v;
    }
    while(fact--){
        pathsLists.push_back(curPath);
        next_permutation(all(curPath));
    }
}

const int INF = 1000;
const vector < vector <int> > distsMatrix = {
        {0,12,INF,INF,1,INF},
        {12,0,INF,1,INF,8},
        {INF,INF,0,20,INF,2},
        {INF,1,20,0,INF,INF},
        {1,INF,INF,INF,0,8},
        {INF,8,2,INF,8,0}
};

int populationMalloc = 0;
int individMalloc = 0;
static std::mt19937 genn(std::chrono::high_resolution_clock::now().time_since_epoch().count());

int get_inumber(int from,int to){
    std::uniform_int_distribution<int> dist(from, to);
    return dist(genn);
}

double get_dnumber(double from,double to){
    std::uniform_real_distribution<double> dist(from,to);
    return dist(genn);
}

class CIndivid{
public:
    explicit CIndivid(int vertexsCount){
        chromosomeLen = vertexsCount;
        chromosomePaths.resize(chromosomeLen);
        fact = 1;
        for (int i = 2;i <= vertexsCount;i++) fact *= i;
        individMalloc++;
    }
    ~CIndivid(){
        chromosomePaths.clear();
        individMalloc--;
    }
    void generate(){
        for (int i = 0;i < chromosomeLen;i++) {
            int pathNumber = get_inumber(1, fact);
            chromosomePaths[i] = pathsLists[pathNumber - 1];
        }
    }
    void print(int vertexFrom) const{
        for (int v = 1,path = 0;path < chromosomeLen;path++,v++){
            std::cout << "shortest path from " << vertexFrom << " to " << v << ": " << vertexFrom << " -> ";
            for (auto& p : chromosomePaths[path]){
                std::cout << p << ' ';
                if (p == v) break;
                std::cout << "-> ";
            }
            std::cout << "\n";
        }
    }
    int get_fitness_value(int vertexFrom) const{
        int fitnessValue = 0;
        for (int i = 0,u = 1;u <= chromosomeLen;u++,i++){
            int pred = vertexFrom;
            for (auto& v : chromosomePaths[i]){
                fitnessValue -= distsMatrix[pred - 1][v - 1];
                if (u == v) break;
                pred = v;
            }
        }
        return fitnessValue;
    }
    int get_chromosome_len() const{
        return chromosomeLen;
    }
    int get_gen(int path,int gen) const{
        return chromosomePaths[path][gen];
    }
    void set_gen(int path,int gen,int newValue){
        chromosomePaths[path][gen] = newValue;
    }
    auto *clone() const{
        auto *clone = new CIndivid(this->chromosomeLen);
        clone->fact = this->fact;
        clone->chromosomePaths = this->chromosomePaths;
        return clone;
    }
    vector <int> get_path(int path) const{
        return chromosomePaths[path];
    };
private:
    int fact;
    int chromosomeLen;
    vector < vector <int> > chromosomePaths;
};

class CPopulation{
public:
    explicit CPopulation(int populationSize,
                         int vertexsCount){
        this->populationSize = populationSize;
        this->vertexsCount = vertexsCount;
        populationMalloc++;
    }
    ~CPopulation(){
        for (auto& individ : population){
            if (individ){
                delete individ; individ = nullptr;
            }
        }
        population.clear();
        populationMalloc--;
    }
    void gen_population(){
        for (int individNumber = 0;individNumber < populationSize;individNumber++){
            auto *individ = new CIndivid(vertexsCount);
            individ->generate();
            population.push_back(individ);
        }
    }
    void add_individ(CIndivid *individ){
        population.push_back(individ);
        populationSize++;
    }
    int size() const{
        return populationSize;
    }
    auto *get_individ(int individNumber) const{
        return population[individNumber];
    }
private:
    int populationSize;
    int vertexsCount;
    vector <CIndivid*> population;
};

auto *new_population(CPopulation *population,int vertexsCount,int vertexFrom){
    auto *newPopulation = new CPopulation(0,vertexsCount);
    while(newPopulation->size() < population->size()){
        int individNumber1 = 1,individNumber2 = 1,individNumber3 = 1;
        while(individNumber1 == individNumber2 || individNumber1 == individNumber3 ||
              individNumber2 == individNumber3){
            individNumber1 = get_inumber(0,population->size() - 1);
            individNumber2 = get_inumber(0,population->size() - 1);
            individNumber3 = get_inumber(0,population->size() - 1);
        }
        auto individ1 = population->get_individ(individNumber1);
        auto individ2 = population->get_individ(individNumber2);
        auto individ3 = population->get_individ(individNumber3);
        vector < pair <int,CIndivid*> > individsFitnessValue = {
                std::make_pair(individ1->get_fitness_value(vertexFrom),individ1),
                std::make_pair(individ2->get_fitness_value(vertexFrom),individ2),
                std::make_pair(individ3->get_fitness_value(vertexFrom),individ3)
        };
        stable_sort(rall(individsFitnessValue));
        newPopulation->add_individ(individsFitnessValue.front().second->clone());
    }
    delete population;
    population = nullptr;
    return newPopulation;
}

void update_individ(CIndivid *individ,
                    int path,int to,int vertexsCount,
                    set <int> &individStableGens,
                    const vector <int> &notUpdatedPath){
    int gen = (to + 1) % vertexsCount;
    int ptr = gen;
    while(isz(individStableGens) != vertexsCount){
        if (!individStableGens.count(notUpdatedPath[ptr])){
            individ->set_gen(path,gen,notUpdatedPath[ptr]);
            individStableGens.insert(individ->get_gen(path,gen));
            gen = (gen + 1) % vertexsCount;
        }
        ptr = (ptr + 1) % vertexsCount;
    }
}

void cross(CIndivid *individ1,CIndivid *individ2,int vertexsCount,double crossChance){
    int chromosomeLen = individ1->get_chromosome_len();
    for (int path = 0;path < chromosomeLen;path++){
        double curChance = get_dnumber(0,1);
        if (curChance >= crossChance) continue;
        int from = get_inumber(2,chromosomeLen - 1);
        int to = get_inumber(2,chromosomeLen - 1);
        from--; to--;
        if (from > to) std::swap(from,to);
        set <int> firstIndividStableGens;
        set <int> secondIndividStableGens;
        vector <int> notUpdatedPath1 = individ1->get_path(path);
        vector <int> notUpdatedPath2 = individ2->get_path(path);
        for (int gen = from;gen <= to;gen++){
            int temp = individ1->get_gen(path,gen);
            individ1->set_gen(path,gen,individ2->get_gen(path,gen));
            individ2->set_gen(path,gen,temp);
            firstIndividStableGens.insert(individ1->get_gen(path,gen));
            secondIndividStableGens.insert(individ2->get_gen(path,gen));
        }
        update_individ(individ1,path,to,vertexsCount,firstIndividStableGens,notUpdatedPath1);
        update_individ(individ2,path,to,vertexsCount,secondIndividStableGens,notUpdatedPath2);
    }
}

void mutate(CIndivid *individ,double mutateChance){
    int chromosomeLen = individ->get_chromosome_len();
    for (int path = 0;path < chromosomeLen;path++){
        double curChance = get_dnumber(0,1);
        if (curChance >= mutateChance) continue;
        int firstGen = 1,secondGen = 1;
        while(firstGen == secondGen){
            firstGen = get_inumber(0,chromosomeLen - 1);
            secondGen = get_inumber(0,chromosomeLen - 1);
        }
        int temp = individ->get_gen(path,firstGen);
        individ->set_gen(path,firstGen,individ->get_gen(path,secondGen));
        individ->set_gen(path,secondGen,temp);
    }
}

class CGeneticAlgorithm{
public:
    CGeneticAlgorithm(){
        population = new CPopulation(populationSize,vertexsCount);
    }
    ~CGeneticAlgorithm(){
        delete population;
    }
    void gen_start_population(){
        population->gen_population();
    }
    void start(){
        int generationNumber = 1;
        while(generationNumber <= cntGenerations){
            population = new_population(population,vertexsCount,vertexFrom);
            for (int individNumber = 0;individNumber < population->size() - 1;individNumber++){
                auto individ1 = population->get_individ(individNumber);
                auto individ2 = population->get_individ(individNumber + 1);
                double krossChance = get_dnumber(0,1);
                if (krossChance < krossingoverChance){
                    cross(individ1,individ2,vertexsCount,double(krossingoverChance / vertexsCount));
                }
            }
            for (int individNumber = 0;individNumber < population->size();individNumber++){
                double curChance = get_dnumber(0,1);
                auto individ = population->get_individ(individNumber);
                if (curChance < mutateChance) mutate(individ,double(mutateChance / individ->get_chromosome_len()));
            }
            int maxFitnesValue = -1e9;
            for (int individNumber = 0;individNumber < population->size();individNumber++){
                maxFitnesValue = max(maxFitnesValue,population->get_individ(individNumber)->get_fitness_value(vertexFrom));
            }
            generationNumber++;
        }
    }
    auto *get_best_individ() const{
        vector < std::pair <int,CIndivid*> > individsFitnessValue;
        for (int individNumber = 0;individNumber < population->size();individNumber++){
            auto individ = population->get_individ(individNumber);
            individsFitnessValue.emplace_back(individ->get_fitness_value(vertexFrom),individ);
        }
        stable_sort(rall(individsFitnessValue));
        return individsFitnessValue.front().second;
    }
private:
    CPopulation *population;
    const int vertexsCount = 6;
    const int vertexFrom = 1;
    const int populationSize = 500;
    const int cntGenerations = 150;
    const double krossingoverChance = 0.9;
    const double mutateChance = 0.1;
};

int32_t main(){
    {
        const int vertexsCount = 6;
        const int vertexFrom = 1;
        precalc(vertexsCount);
        CGeneticAlgorithm geneticAlgorithm;
        geneticAlgorithm.gen_start_population();
        geneticAlgorithm.start();
        auto bestIndivid = geneticAlgorithm.get_best_individ();
        std::cout << "shortest path sum - " << -bestIndivid->get_fitness_value(vertexFrom) << "\n";
        bestIndivid->print(vertexFrom);
    }
    if (individMalloc + populationMalloc != 0) std::cout << "memory leak!!!";
}
