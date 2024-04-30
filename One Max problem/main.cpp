#include <bits/stdc++.h>

using namespace std;

#define isz(x) (int)x.size()
#define all(x) x.begin(),x.end()
#define rall(x) x.rbegin(),x.rend()
#define fr first
#define sc second

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
    CIndivid(){
        this->gens.resize(gensLen);
        individMalloc++;
    }
    ~CIndivid(){
        gens.clear();
        individMalloc--;
    }
    void generate(){
        for (int i = 0;i < gensLen;i++) gens[i] = get_inumber(0,1);
    }
    void print() const{
        for (const auto& ch : gens) std::cout << ch;
        std::cout << "\n";
    }
    int get_gen(int pos) const{
        return gens[pos];
    }
    int set_gen(int pos,bool value){
        gens[pos] = value;
    }
    int get_fitness_value() const{
        return std::accumulate(all(gens),0);
    }
    int get_gens_len() const{
        return gensLen;
    }
    vector <bool> get_gens() const{
        return gens;
    }
    void mutate_gen(int gen){
        gens[gen] = !gens[gen];
    }
    CIndivid *clone() const{
        CIndivid *newIndivid = new CIndivid();
        newIndivid->gens = this->gens;
        return newIndivid;
    }
private:
    const int gensLen = 100;
    vector <bool> gens;
};

class CPopulation{
public:
    CPopulation(){
        populationSize = 0;
        populationMalloc++;
    }
    explicit CPopulation(int populationSize){
        this->populationSize = populationSize;
        populationMalloc++;
    }
    ~CPopulation(){
        for (auto &individ : population){
            if (individ) {
                delete individ; individ = nullptr;
            }
        }
        population.clear();
        populationMalloc--;
    }
    void gen_population(){
        for (int i = 0;i < populationSize;i++){
            auto individ = new CIndivid();
            individ->generate();
            population.push_back(individ);
        }
    }
    void add_individ(CIndivid *newIndivid){
        population.push_back(newIndivid);
        populationSize = isz(population);
    }
    int size() const{
        return populationSize;
    }
    CIndivid *get_individ(int pos) const{
        return population[pos];
    }
private:
    int populationSize;
    vector <CIndivid*> population;
};

CPopulation *new_population(CPopulation *population){
    auto newPopulation = new CPopulation();
    while(newPopulation->size() < population->size()){
        int individPos1 = 1,individPos2 = 1,individPos3 = 1;
        while(individPos1 == individPos2 || individPos1 == individPos3
              || individPos2 == individPos3){
            individPos1 = get_inumber(0,population->size() - 1);
            individPos2 = get_inumber(0,population->size() - 1);
            individPos3 = get_inumber(0,population->size() - 1);
        }
        auto individ1 = population->get_individ(individPos1);
        auto individ2 = population->get_individ(individPos2);
        auto individ3 = population->get_individ(individPos3);
        vector < pair <int,CIndivid *> > individsRang = {
                std::make_pair(individ1->get_fitness_value(),individ1),
                std::make_pair(individ2->get_fitness_value(),individ2),
                std::make_pair(individ3->get_fitness_value(),individ3)
        };
        stable_sort(rall(individsRang));
        newPopulation->add_individ(individsRang.front().second->clone());
    }
    delete population;
    population = nullptr;
    return newPopulation;
}

void cross(CIndivid *left,CIndivid *right){
    int crossPos = get_inumber(2,left->get_gens_len() - 1);
    crossPos--;
    for (int gen = crossPos;gen < left->get_gens_len();gen++){
        int leftGen = left->get_gen(gen);
        left->set_gen(gen,right->get_gen(gen));
        right->set_gen(gen,leftGen);
    }
}

void mutate(CIndivid *individ,double mutateChance){
    for (int gen = 0;gen < individ->get_gens_len();gen++){
        double chance = get_dnumber(0,1);
        if (chance < mutateChance) individ->mutate_gen(gen);
    }
}

class CGenericAlgorithm{
public:
    CGenericAlgorithm(){
        population = new CPopulation(populationSize);
    }
    ~CGenericAlgorithm(){
        delete population;
    }
    void gen_start_population(){
        population->gen_population();
        for (int individ = 0; individ < population->size();individ++){
            fitnessValues.push_back(population->get_individ(individ)->get_fitness_value());
        }
    }
    void start(){
        int generationNumber = 0;
        while(*max_element(all(fitnessValues)) < population->get_individ(0)->get_gens_len()
              && generationNumber < cntGeneration){
            population = new_population(population);
            fitnessValues.clear();
            for (int individ = 0;individ < population->size() - 1;individ++){
                auto leftIndivid = population->get_individ(individ);
                auto rightIndivid = population->get_individ(individ + 1);
                double chance = get_dnumber(0,1);
                if (chance < krossingoverChance) cross(leftIndivid,rightIndivid);
            }
            for (int individ = 0;individ < population->size();individ++){
                double chance = get_dnumber(0,1);
                if (chance < this->mutateChance) mutate(population->get_individ(individ),1. / population->get_individ(individ)->get_gens_len());
            }
            int maxFitnessValue = 0;
            for (int individPos = 0;individPos < population->size();individPos++){
                auto individ = population->get_individ(individPos);
                maxFitnessValue = max(maxFitnessValue,individ->get_fitness_value());
                fitnessValues.push_back(individ->get_fitness_value());
            }
            std::cout << "generation - " << generationNumber << "\n";
            std::cout << "max fitness value - " << maxFitnessValue << "\n";
            generationNumber++;
        }
    }
    CIndivid *get_best_individ() const{
        pair <int,CIndivid*> best = std::make_pair(population->get_individ(0)->get_fitness_value(),
                                                   population->get_individ(0));
        for (int individPos = 1;individPos < population->size();individPos++){
            auto individ = population->get_individ(individPos);
            if (individ->get_fitness_value() > best.first){
                best.first = individ->get_fitness_value();
                best.second = individ;
            }
        }
        return best.second;
    }
private:
    const int cntGeneration = 50;
    const int populationSize = 100;
    const double krossingoverChance = 0.6456;
    const double mutateChance = 0.1;
    CPopulation *population;
    vector <int> fitnessValues;
};

int32_t main() {
    {
        CGenericAlgorithm genericAlgorithm;
        genericAlgorithm.gen_start_population();
        genericAlgorithm.start();
        std::cout << "final result - " << genericAlgorithm.get_best_individ()->get_fitness_value() << "\n";
    }
    if (populationMalloc + individMalloc != 0) std::cout << "memory leak!!!";
}
