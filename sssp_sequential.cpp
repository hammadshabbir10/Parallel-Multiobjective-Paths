#include <iostream>
#include <vector>
#include <queue>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <string>
#include <climits>
#include <chrono>
#include <random>  // For random edge generation
#include <algorithm>  // For shuffle function

using namespace std;
using namespace std::chrono;

const int INF = INT_MAX;

struct Edge {
    int to;
    int weight;
};

struct Change {
    int u, v, weight;
    bool inserted;
};

vector<vector<Edge>> graph;
unordered_map<int, int> nodeToIndex;
vector<int> indexToNode;
int currentIndex = 0;
int N = 0;

vector<int> dist, parent;
vector<bool> affected;

// Load graph from unweighted edge list format
void loadUnweightedEdgeList(const string& filename) {
    ifstream infile(filename);
    if (!infile.is_open()) {
        cerr << "Error opening input file.\n";
        exit(1);
    }

    string line;
    int total_edges = 0;
    while (getline(infile, line)) {
        if (line[0] == '#') continue; // Skip comments
        stringstream ss(line);
        int u_orig, v_orig;
        ss >> u_orig >> v_orig;

        // Ensure all nodes are indexed
        if (!nodeToIndex.count(u_orig)) {
            nodeToIndex[u_orig] = currentIndex++;
            indexToNode.push_back(u_orig);
        }
        if (!nodeToIndex.count(v_orig)) {
            nodeToIndex[v_orig] = currentIndex++;
            indexToNode.push_back(v_orig);
        }

        int u = nodeToIndex[u_orig];
        int v = nodeToIndex[v_orig];

        if (graph.size() <= max(u, v))
            graph.resize(max(u, v) + 1);

        graph[u].push_back({v, 1});  // weight = 1
        total_edges++;
    }

    N = currentIndex;
    dist.resize(N, INF);
    parent.resize(N, -1);
    affected.resize(N, false);
    cout << "Loaded " << N << " nodes and " << total_edges << " edges from file.\n";
}

// Initialize for Dijkstra
void initialize(int src) {
    fill(dist.begin(), dist.end(), INF);
    fill(parent.begin(), parent.end(), -1);
    dist[src] = 0;
}

// Run full Dijkstra
void dijkstra(int src) {
    auto start = high_resolution_clock::now();
    
    initialize(src);
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<>> pq;
    pq.push({0, src});

    while (!pq.empty()) {
        auto [d, u] = pq.top(); pq.pop();
        if (d > dist[u]) continue;

        for (auto edge : graph[u]) {
            int v = edge.to;
            if (dist[v] > dist[u] + edge.weight) {
                dist[v] = dist[u] + edge.weight;
                parent[v] = u;
                pq.push({dist[v], v});
            }
        }
    }
    
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start).count();
    cout << "Full Dijkstra's algorithm completed in " << duration << " ms\n";
}

// Check if edge exists
bool edgeExists(int u, int v) {
    for (const auto& edge : graph[u]) {
        if (edge.to == v) {
            return true;
        }
    }
    return false;
}

// Process a single edge change and run Dijkstra after each change
bool processChangeWithFullDijkstra(const Change& c, int src) {
    auto start = high_resolution_clock::now();
    bool success = true;
    
    if (c.inserted) {
        if (edgeExists(c.u, c.v)) {
            cout << "Edge already exists: " << indexToNode[c.u] << " -> " << indexToNode[c.v] << "\n";
            success = false;
        } else {
            cout << "Inserting edge: " << indexToNode[c.u] << " -> " << indexToNode[c.v] << " (weight: " << c.weight << ")\n";
            graph[c.u].push_back({c.v, c.weight});
        }
    } else {
        if (!edgeExists(c.u, c.v)) {
            cout << "Edge doesn't exist: " << indexToNode[c.u] << " -> " << indexToNode[c.v] << "\n";
            success = false;
        } else {
            cout << "Removing edge: " << indexToNode[c.u] << " -> " << indexToNode[c.v] << "\n";
            for (auto it = graph[c.u].begin(); it != graph[c.u].end(); ++it) {
                if (it->to == c.v) {
                    graph[c.u].erase(it);
                    break;
                }
            }
        }
    }
    
    auto modifyEnd = high_resolution_clock::now();
    auto modifyDuration = duration_cast<microseconds>(modifyEnd - start).count();
    
    // Only run Dijkstra if the graph was actually modified
    if (success) {
        // Run full Dijkstra after each change
        dijkstra(src);
        return true;
    }
    
    return false;
}

// Process edge changes
int processChanges(vector<Change>& changes) {
    int affectedCount = 0;
    
    for (auto& c : changes) {
        if (c.inserted) {
            graph[c.u].push_back({c.v, c.weight});
            if (dist[c.u] != INF && dist[c.v] > dist[c.u] + c.weight) {
                dist[c.v] = dist[c.u] + c.weight;
                parent[c.v] = c.u;
                affected[c.v] = true;
                affectedCount++;
            }
        } else {
            for (auto it = graph[c.u].begin(); it != graph[c.u].end(); ++it) {
                if (it->to == c.v) {
                    graph[c.u].erase(it);
                    break;
                }
            }
            if (parent[c.v] == c.u) {
                dist[c.v] = INF;
                parent[c.v] = -1;
                affected[c.v] = true;
                affectedCount++;
            }
        }
    }
    
    cout << "Number of nodes affected by changes: " << affectedCount << " (" 
         << (float)affectedCount / N * 100 << "% of total nodes)\n";
    
    return affectedCount;
}

// Incremental update for affected nodes
void updateAffected() {
    auto start = high_resolution_clock::now();
    
    int iterationCount = 0;
    int totalUpdates = 0;
    
    bool change = true;
    while (change) {
        change = false;
        iterationCount++;
        int updatesThisIteration = 0;
        
        for (int u = 0; u < N; ++u) {
            if (!affected[u]) continue;
            affected[u] = false;

            for (auto edge : graph[u]) {
                int v = edge.to;
                if (dist[v] > dist[u] + edge.weight) {
                    dist[v] = dist[u] + edge.weight;
                    parent[v] = u;
                    affected[v] = true;
                    change = true;
                    updatesThisIteration++;
                }
            }
        }
        
        totalUpdates += updatesThisIteration;
        cout << "  - Iteration " << iterationCount << ": Updated " << updatesThisIteration << " nodes\n";
    }
    
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start).count();
    cout << "Incremental update completed in " << duration << " ms\n";
    cout << "Total iterations: " << iterationCount << ", Total node updates: " << totalUpdates << "\n";
}

// Generate random changes - guarantees the exact requested number of valid changes
vector<Change> generateRandomChanges(int count, bool insertions) {
    vector<Change> changes;
    
    // Random number generator
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> nodeDist(0, N-1);
    
    if (insertions) {
        // For large graphs, don't try to enumerate all possible edges
        // Instead, generate random edges until we find enough valid ones
        int attempts = 0;
        int maxAttempts = count * 100; // Limit attempts to avoid infinite loop
        
        cout << "Generating " << count << " random insertions...\n";
        while (changes.size() < count && attempts < maxAttempts) {
            attempts++;
            
            // Pick two random nodes
            int u = nodeDist(gen);
            int v = nodeDist(gen);
            
            // Don't create self-loops
            if (u == v) continue;
            
            // Make sure edge doesn't exist
            if (!edgeExists(u, v)) {
                changes.push_back({u, v, 1, true});
                if (changes.size() % 10 == 0 || changes.size() == count) {
                    cout << "  Progress: " << changes.size() << "/" << count << " edges found\r";
                    cout.flush();
                }
            }
        }
        cout << endl;
        
        if (changes.size() < count) {
            cout << "Warning: Could only generate " << changes.size() << " valid insertions after " 
                 << attempts << " attempts.\n";
        }
    } else {
        // For deletions, use reservoir sampling to select random existing edges
        vector<pair<int, int>> existingEdges;
        
        cout << "Collecting existing edges for deletion...\n";
        // Collect existing edges
        for (int u = 0; u < N; u++) {
            for (const auto& edge : graph[u]) {
                existingEdges.push_back({u, edge.to});
            }
        }
        
        cout << "Found " << existingEdges.size() << " existing edges.\n";
        
        // Check if we have enough existing edges
        if (existingEdges.size() < count) {
            cout << "Warning: Only " << existingEdges.size() << " edges exist for deletion. ";
            cout << "Using all available edges.\n";
            count = existingEdges.size();
        }
        
        // Shuffle the existing edges
        cout << "Shuffling edges...\n";
        shuffle(existingEdges.begin(), existingEdges.end(), gen);
        
        // Take the first 'count' edges
        for (int i = 0; i < count; i++) {
            changes.push_back({existingEdges[i].first, existingEdges[i].second, 1, false});
        }
    }
    
    cout << "Successfully generated " << changes.size() << " valid " 
         << (insertions ? "insertions" : "deletions") << ".\n";
    
    return changes;
}

// Display shortest path tree
void printSSSP(int limit = 100) {
    cout << "\nShortest paths from source:\n";
    for (int i = 0; i < min(N, limit); ++i) {
        cout << "Node " << indexToNode[i]
             << ": Distance = " << (dist[i] == INF ? -1 : dist[i])
             << ", Parent = " << (parent[i] == -1 ? -1 : indexToNode[parent[i]])
             << endl;
    }
}

// Main program
int main() {
    loadUnweightedEdgeList("graph.txt");

    int srcID;
    cout << "Enter source node ID: ";
    cin >> srcID;

    if (!nodeToIndex.count(srcID)) {
        cerr << "Invalid source node ID.\n";
        return 1;
    }

    int src = nodeToIndex[srcID];
    dijkstra(src);
    printSSSP();

    int choice;
    do {
        cout << "\n1. Random Batch Updates\n2. Manual Updates\n3. Exit\nChoice: ";
        cin >> choice;

        if (choice == 1) {
            // Random batch updates
            int count, updateType;
            cout << "Enter number of updates: ";
            cin >> count;
            
            cout << "Update type (1: Insert, 2: Delete): ";
            cin >> updateType;
            
            bool insertions = (updateType == 1);
            
            vector<Change> changes = generateRandomChanges(count, insertions);
            
            cout << "\nGenerating " << count << " random " 
                 << (insertions ? "insertions" : "deletions") << "...\n";
            
            auto totalStart = high_resolution_clock::now();
            
            // Process each change individually and run Dijkstra after each
            int successfulChanges = 0;
            long long totalDijkstraTime = 0;
            
            for (const auto& change : changes) {
                auto changeStart = high_resolution_clock::now();
                bool success = processChangeWithFullDijkstra(change, src);
                auto changeEnd = high_resolution_clock::now();
                
                if (success) {
                    successfulChanges++;
                    totalDijkstraTime += duration_cast<milliseconds>(changeEnd - changeStart).count();
                }
            }
            
            auto totalEnd = high_resolution_clock::now();
            auto totalDuration = duration_cast<milliseconds>(totalEnd - totalStart).count();
            
            cout << "\nBatch processing summary:\n";
            cout << "- Successful changes: " << successfulChanges << "/" << changes.size() << "\n";
            cout << "- Total execution time: " << totalDuration << " ms\n";
            if (successfulChanges > 0) {
                cout << "- Average time per successful update: " << totalDuration / successfulChanges << " ms\n";
            }
            
            printSSSP();
            
        } else if (choice == 2) {
            // Manual updates
            int updateChoice;
            cout << "\n1. Insert Batch\n2. Delete Batch\nChoice: ";
            cin >> updateChoice;

            if (updateChoice == 1 || updateChoice == 2) {
                int count;
                cout << "Enter number of changes: ";
                cin >> count;
                vector<Change> changes;
                
                for (int i = 0; i < count; ++i) {
                    int uID, vID;
                    cout << "Change " << i + 1 << " (from to): ";
                    cin >> uID >> vID;

                    if (!nodeToIndex.count(uID)) {
                        nodeToIndex[uID] = currentIndex++;
                        indexToNode.push_back(uID);
                        graph.resize(currentIndex);
                        dist.push_back(INF);
                        parent.push_back(-1);
                        affected.push_back(false);
                    }
                    if (!nodeToIndex.count(vID)) {
                        nodeToIndex[vID] = currentIndex++;
                        indexToNode.push_back(vID);
                        graph.resize(currentIndex);
                        dist.push_back(INF);
                        parent.push_back(-1);
                        affected.push_back(false);
                    }

                    changes.push_back({nodeToIndex[uID], nodeToIndex[vID], 1, updateChoice == 1});
                }

                auto totalStart = high_resolution_clock::now();
                
                // Process each change individually and run Dijkstra after each
                int successfulChanges = 0;
                long long totalDijkstraTime = 0;
                
                for (const auto& change : changes) {
                    auto changeStart = high_resolution_clock::now();
                    bool success = processChangeWithFullDijkstra(change, src);
                    auto changeEnd = high_resolution_clock::now();
                    
                    if (success) {
                        successfulChanges++;
                        totalDijkstraTime += duration_cast<milliseconds>(changeEnd - changeStart).count();
                    }
                }
                
                auto totalEnd = high_resolution_clock::now();
                auto totalDuration = duration_cast<milliseconds>(totalEnd - totalStart).count();
                
                cout << "\nBatch processing summary:\n";
                cout << "- Successful changes: " << successfulChanges << "/" << changes.size() << "\n";
                cout << "- Total execution time: " << totalDuration << " ms\n";
                if (successfulChanges > 0) {
                    cout << "- Average time per successful update: " << totalDuration / successfulChanges << " ms\n";
                }

                printSSSP();
            }
        }

    } while (choice != 3);

    return 0;
}
