//Projekt 3 Projektowanie efektwnych algorytmów
//Autor Cyprian Kozubek nr.indeksu 272959
//Problem komiwojazera:  Algorytmu Symulowanego Wyżarzania
// Implementacja i analiza efektywności algorytmu Symulowanego Wyżarzania dla problemu komiwojażera
//Badane pliki: ftv47.atsp (1776), ftv170.atsp (2755) , rgb403.atsp (2465).
#include <iostream>
#include <windows.h>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <sstream>
#include <random>
#include <algorithm>

using namespace std;

//Elementy zajmujące się mierzeniem czasu
double t1 = 0.0;
double PCFreq = 0.0;
__int64 CounterStart = 0;
void StartCounter(){
    LARGE_INTEGER li;
    if( !QueryPerformanceFrequency( & li ) )
        cout << "QueryPerformanceFrequency failed!\n";
    PCFreq = double( li.QuadPart ) ;
    QueryPerformanceCounter( & li );
    CounterStart = li.QuadPart;
}
double GetCounter(){
    LARGE_INTEGER li;
    QueryPerformanceCounter( & li );
    return double( li.QuadPart - CounterStart ) / PCFreq;
}

int ** Matrix = nullptr; //Tablica dwuwymiarowa dla grafu
int matrix_size; //Rozmiar grafu
// Funkcja alokująca pamięć dla macierzy
void allocateMatrix(int rows, int cols) {
    Matrix = new int*[rows];
    for (int i = 0; i < rows; i++) {
        Matrix[i] = new int[cols];
    }
}
//Funkcja zwalniajaca pamiec macierzy
void freeMatrix(int rows) {
    for (int i = 0; i < rows; ++i) {
        delete[] Matrix[i];
    }
    delete[] Matrix;
}
//Wczytywanie danych z pliku xml
void loadGraphFromFile(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Unable to open file: " << filename << endl;
        exit(1);
    }
    string line;
    int current_vertex = -1;
    matrix_size = 0;
    // Zliczanie liczby wierzchołków
    while (getline(file, line)) {
        if (line.find("<vertex>") != string::npos) {
            ++matrix_size;
        }
    }
    allocateMatrix(matrix_size, matrix_size);
    //Wczytywanie danych
    file.clear();
    file.seekg(0, ios::beg);
    while (getline(file, line)) {
        if (line.find("<vertex>") != string::npos) {
            ++current_vertex;
        } else if (line.find("<edge") != string::npos) {
            size_t cost_pos = line.find("cost=\"");
            if (cost_pos != string::npos) {
                cost_pos += 6;
                size_t end_pos = line.find("\"", cost_pos);
                string cost_str = line.substr(cost_pos, end_pos - cost_pos);
                double cost = atof(cost_str.c_str());
                size_t close_tag_pos = line.find(">", end_pos);
                size_t vertex_end_pos = line.find("<", close_tag_pos);
                string vertex_str = line.substr(close_tag_pos + 1, vertex_end_pos - close_tag_pos - 1);
                int target_vertex = atoi(vertex_str.c_str());
                Matrix[current_vertex][target_vertex] = static_cast<int>(cost);
            }
        }
    }
    file.close();
    cout << "Graf prawidlowo wczytany. Rozmiar macierzy: " << matrix_size << "x" << matrix_size << endl;
}
// Generowanie losowego grafu (asymetrycznego)
void generate_graph(int graphSize) {
    matrix_size = graphSize;
    allocateMatrix(graphSize, graphSize);
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> distrib(1, 100);
    for (int i = 0; i < graphSize; i++) {
        for (int j = 0; j < graphSize; j++) {
            if (i == j) {
                Matrix[i][j] = -1;
            } else {
                int liczba = distrib(gen);
                Matrix[i][j] = liczba;
            }
        }
    }
}
// Generowanie losowego grafu (symetrycznego)
void generate_graph_sym(int graphSize) {
    matrix_size = graphSize;
    allocateMatrix(graphSize, graphSize);
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> distrib(1, 100);
    for (int i = 0; i < graphSize; i++) {
        for (int j = 0; j < graphSize; j++) {
            if (i == j) {
                Matrix[i][j] = -1;
            } else if (i < j) {
                int liczba = distrib(gen);
                Matrix[i][j] = liczba;
                Matrix[j][i] = liczba;
            }
        }
    }
}
//Wyswietlanie macierzy
void print_matrix(int size){
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            cout<<Matrix[i][j]<<" ";
        }
        cout<<endl;
    }
}
//******************************************************************************************************

// Funkcja obliczająca koszt cyklu Hamiltona
int calculate_cost(int* solution, int size) {
    int cost = 0;
    for (int i = 0; i < size - 1; ++i) {
        cost += Matrix[solution[i]][solution[i + 1]];
    }
    cost += Matrix[solution[size - 1]][solution[0]];
    return cost;
}
//Funkcja najblizszego sasiada obliczajaca wartosc poczatkowa dla symulowanego wyzarzania
void nearest_neighbor(int* solution, int matrix_size) {
    int* best_solution = new int[matrix_size];
    int best_cost = INT_MAX;
    for (int start = 0; start < matrix_size; ++start) {
        bool* visited = new bool[matrix_size]();
        int* temp_solution = new int[matrix_size];
        temp_solution[0] = start;
        visited[start] = true;
        for (int i = 1; i < matrix_size; ++i) {
            int current = temp_solution[i - 1];
            int next_city = -1;
            int min_cost = INT_MAX;
            for (int j = 0; j < matrix_size; ++j) {
                if (!visited[j] && Matrix[current][j] != -1 && Matrix[current][j] < min_cost) {
                    min_cost = Matrix[current][j];
                    next_city = j;
                }
            }
            if (next_city != -1) {
                temp_solution[i] = next_city;
                visited[next_city] = true;
            } else {
                for (int k = 0; k < matrix_size; ++k) {
                    if (!visited[k]) {
                        temp_solution[i] = k;
                        visited[k] = true;
                        break;
                    }
                }
            }
        }
        int cost = calculate_cost(temp_solution, matrix_size);
        if (cost < best_cost) {
            best_cost = cost;
            copy(temp_solution, temp_solution + matrix_size, best_solution);
        }
        delete[] temp_solution;
        delete[] visited;
    }
    copy(best_solution, best_solution + matrix_size, solution);
    delete[] best_solution;
}
void random_swap(int* solution, int matrix_size, mt19937& gen) {
    uniform_int_distribution<int> dist(0, matrix_size - 1);
    int i = dist(gen);
    int j;
    do {
        j = dist(gen);
    } while (i == j);
    swap(solution[i], solution[j]);
}
//Funkcja wyboru sasiada metoda 2-opt
void two_opt(int* solution, int matrix_size, mt19937& gen) {
    uniform_int_distribution<int> dist(0, matrix_size - 1);
    int i = dist(gen);
    int j;
    do {
        j = dist(gen);
    } while (i == j || abs(i - j) < 2);

    if (i > j) swap(i, j);
    while (i < j) {
        swap(solution[i], solution[j]);
        i++;
        j--;
    }
}
//Funkcja wyboru sasiada metoda 3-opt
void three_opt(int* solution, int matrix_size, mt19937& gen) {
    uniform_int_distribution<int> dist(0, matrix_size - 1);
    int i = dist(gen);
    int j = dist(gen);
    int k = dist(gen);
    while (j == i) j = dist(gen);
    while (k == i || k == j) k = dist(gen);

    if (i > j) swap(i, j);
    if (j > k) swap(j, k);
    if (i > j) swap(i, j);

    int temp[matrix_size];
    for (int x = 0; x < i; ++x) temp[x] = solution[x];
    for (int x = i, y = j; y >= i; ++x, --y) temp[x] = solution[y];
    for (int x = j + 1, y = k; y >= j + 1; ++x, --y) temp[x] = solution[y];
    for (int x = k + 1; x < matrix_size; ++x) temp[x] = solution[x];

    copy(temp, temp + matrix_size, solution);
}
//Funkcja wyboru sasiada metoda insert
void insert(int* solution, int matrix_size, mt19937& gen) {
    uniform_int_distribution<int> dist(0, matrix_size - 1);
    int i = dist(gen);
    int j = dist(gen);
    while (i == j) j = dist(gen);

    int city = solution[i];
    if (i < j) {
        for (int k = i; k < j; ++k) {
            solution[k] = solution[k + 1];
        }
        solution[j] = city;
    } else {
        for (int k = i; k > j; --k) {
            solution[k] = solution[k - 1];
        }
        solution[j] = city;
    }
}
//Funkcja wykonujaca algorytm symulowanego wyzarzania
void simulated_annealing(double temperature, int initial_neighborhood_type, double stop_time) {
    StartCounter();
    int* current_solution = new int[matrix_size];
    int* best_solution = new int[matrix_size];
    nearest_neighbor(current_solution, matrix_size);
    int current_cost = calculate_cost(current_solution, matrix_size);
    copy(current_solution, current_solution + matrix_size, best_solution);
    int best_cost = current_cost;
    double best_time = 0.0;
    mt19937 gen(random_device{}());
    uniform_real_distribution<> real_distrib(0, 1);
    uniform_int_distribution<int> neighborhood_distrib(0, 3);
    int* temp_solution = new int[matrix_size];
    double cooling_rate = 0.995;
    double initial_temperature = temperature;
    int no_improvement_count = 0;

    while (GetCounter() < stop_time) {
        copy(current_solution, current_solution + matrix_size, temp_solution);
        // Dynamiczny wybór typu sąsiedztwa
        int neighborhood_type = neighborhood_distrib(gen);
        if (initial_neighborhood_type != -1) {
            neighborhood_type = initial_neighborhood_type;
        }
        if (neighborhood_type == 0) {
            random_swap(temp_solution, matrix_size, gen);
        } else if (neighborhood_type == 1) {
            two_opt(temp_solution, matrix_size, gen);
        } else if (neighborhood_type == 2) {
            insert(temp_solution, matrix_size, gen);
        } else if (neighborhood_type == 3) {
            three_opt(temp_solution, matrix_size, gen);
        }
        int new_cost = calculate_cost(temp_solution, matrix_size);
        int cost_diff = new_cost - current_cost;
        // Kryterium Metropolisa z adaptacyjną akceptacją
        if (cost_diff < 0 || real_distrib(gen) < exp(-cost_diff / temperature)) {
            copy(temp_solution, temp_solution + matrix_size, current_solution);
            current_cost = new_cost;
            if (current_cost < best_cost) {
                copy(current_solution, current_solution + matrix_size, best_solution);
                best_cost = current_cost;
                best_time = GetCounter();
                // Wypisanie nowego najlepszego rozwiązania
                cout << "Nowe najlepsze rozwiazanie znalezione!" << endl;
                cout << "Koszt: " << best_cost << ", Czas: " << best_time << " sekund" << endl;
                no_improvement_count = 0;
            }
        } else {
            no_improvement_count++;
        }
        temperature *= cooling_rate;
        if (temperature < 1e-4) {
            temperature = initial_temperature;
        }
    }
    // Wyświetlenie najlepszego wyniku
    cout << "Najlepsze rozwiazanie: ";
    for (int i = 0; i < matrix_size; ++i) {
        cout << best_solution[i] << "-";
    }
    cout << best_solution[0] << endl;
    cout << "Najlepszy koszt: " << best_cost << endl;
    cout << "Czas odnalezienia najlepszego rozwiazania: " << best_time << " sekund" << endl;
    cout << "Koncowa temperatura: " << temperature << endl;
    delete[] current_solution;
    delete[] best_solution;
    delete[] temp_solution;
}

// Funkcja do wczytania danych z pliku konfiguracyjnego
bool loadConfig(const string& configFile, int& graphType, int& algorithmType,
 int& load_graph, string& graphFile, int& repeatCount, int& graphSize,
 int& displayGraph, double& temp, int& neighbour, double& stop_time)
 {
    ifstream config(configFile);
    if (!config.is_open()) {
        cout << "Blad otwarcia pliku konfiguracyjnego!" << endl;
        return false;
    }
    string line;
    while (getline(config, line)) {
        stringstream ss(line);
        string key, value;
        if (line.empty() || line[0] == '#') continue;
        getline(ss, key, '=');
        getline(ss, value);
        key.erase(0, key.find_first_not_of(" \t"));
        key.erase(key.find_last_not_of(" \t") + 1);
        value.erase(0, value.find_first_not_of(" \t"));
        value.erase(value.find_last_not_of(" \t") + 1);
        if (key == "graph_type") {
            graphType = stoi(value);
        } else if (key == "algorithm") {
            algorithmType = stoi(value);
        } else if (key == "load_graph") {
            load_graph = stoi(value);
        } else if (key == "graph_file") {
            graphFile = value;
        } else if (key == "repeat_count") {
            repeatCount = stoi(value);
        } else if (key == "graph_size") {
            graphSize = stoi(value);
        } else if (key == "display_graph") {
            displayGraph = stoi(value);
        } else if (key == "temp"){
            temp = stod(value);
        } else if (key == "neighbour"){
            neighbour = stoi(value);
        }else if (key == "stop_time") {
            stop_time = stod(value);
        }
    }
    config.close();
    return true;
}
// Menu wyboru typu generowanego grafu
void Menu_Type(int graphType, int graphSize) {
    if (graphType == 1) {
        generate_graph(graphSize);
    } else if (graphType == 2) {
        generate_graph_sym(graphSize);
    }
}
// Menu wyboru algorytmów
void Menu_Alg(int algorithmType, int repeatCount, double temp, int neighbour, double stop_time) {
    if ( algorithmType == 1) {
        for(int i=0; i<repeatCount; i++) {
            StartCounter();
            simulated_annealing(temp, neighbour, stop_time);
            t1 = GetCounter();
        }
        cout<< "Czas dzialania simulated annealing w sekundach: " << t1/repeatCount <<endl;
    }
}
// Menu główne programu
void Menu() {
    int graphType = 0, algorithmType = 0, load_graph = 0, repeatCount = 1, graphSize = 0, displayGraph = 0, neighbour = 0;
    double temp = 0,  stop_time = 0;
    string graphFile;
    //Należy podać ścieżkę bezwzględną pliku
    if (!loadConfig(R"(C:\Users\cypri\OneDrive\Pulpit\I\Aizo_projekt1\PEA3\config.txt)",
    graphType, algorithmType, load_graph, graphFile, repeatCount,
    graphSize,displayGraph, temp, neighbour, stop_time)) {
        return;
    }
    if (load_graph == 1 && !graphFile.empty()) {
        loadGraphFromFile(graphFile);
    } else if (load_graph == 2) {
        Menu_Type(graphType, graphSize);
    } else {
        cout << "Nieprawidlowy wybor opcji wczytywania grafu!" << endl;
        return;
    }
    if (displayGraph == 1) {
        cout << "Aktualny graf:" << endl;
        print_matrix(matrix_size);
    }
    Menu_Alg(algorithmType, repeatCount, temp, neighbour, stop_time);
}
int main() {
    Menu();
    freeMatrix(matrix_size);
    return 0;
}
