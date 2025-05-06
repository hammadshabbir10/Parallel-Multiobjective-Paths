#include <iostream>
#include <fstream>
#include <climits>
#include <vector>
#include <map>
#include <set>
#include <queue>
#include <mpi.h>
#include <metis.h>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <chrono>
#include <omp.h>

using namespace std;
using namespace std::chrono;

const int INF = INT_MAX;

struct AdjNode {
    int dest;
    int weight;
    AdjNode* next;
};

struct Graph {
    int V;
    vector<AdjNode*> adj;

    Graph(int vertices) : V(vertices), adj(vertices, nullptr) {}

    void addEdge(int src, int dest, int weight) {
        removeEdge(src, dest);
        adj[src] = new AdjNode{dest, weight, adj[src]};
    }

    void removeEdge(int src, int dest) {
        AdjNode* curr = adj[src];
        AdjNode* prev = nullptr;
        while (curr) {
            if (curr->dest == dest) {
                if (prev) prev->next = curr->next;
                else adj[src] = curr->next;
                delete curr;
                break;
            }
            prev = curr;
            curr = curr->next;
        }
    }

    void clear() {
        for (auto& node : adj) {
            while (node) {
                AdjNode* temp = node;
                node = node->next;
                delete temp;
            }
        }
    }
};

void exchangeBoundaryData(Graph& g, vector<int>& dist, vector<int>& parent, 
                         const vector<int>& part, int rank, int size) {
    vector<int> boundary_nodes;
    #pragma omp parallel
    {
        vector<int> local_boundary_nodes;
        #pragma omp for schedule(dynamic)
        for (int u = 0; u < g.V; u++) {
            if (part[u] != rank) continue;
            
            bool is_boundary = false;
            for (AdjNode* curr = g.adj[u]; curr; curr = curr->next) {
                if (part[curr->dest] != rank) {
                    is_boundary = true;
                    break;
                }
            }
            
            if (is_boundary) {
                local_boundary_nodes.push_back(u);
            }
        }
        #pragma omp critical
        boundary_nodes.insert(boundary_nodes.end(), 
                            local_boundary_nodes.begin(), 
                            local_boundary_nodes.end());
    }

    int local_count = boundary_nodes.size();
    vector<int> boundary_counts(size);
    
    for (int i = 0; i < size; i++) {
        int count_to_share = (i == rank) ? local_count : 0;
        MPI_Bcast(&count_to_share, 1, MPI_INT, i, MPI_COMM_WORLD);
        boundary_counts[i] = count_to_share;
    }
    
    for (int i = 0; i < size; i++) {
        if (boundary_counts[i] > 0) {
            vector<int> proc_data;
            
            if (i == rank) {
                for (int u : boundary_nodes) {
                    proc_data.push_back(u);
                    proc_data.push_back(dist[u]);
                    proc_data.push_back(parent[u]);
                }
            } else {
                proc_data.resize(boundary_counts[i] * 3);
            }
            
            MPI_Bcast(proc_data.data(), proc_data.size(), MPI_INT, i, MPI_COMM_WORLD);
            
            if (i != rank) {
                #pragma omp parallel for schedule(static)
                for (size_t j = 0; j < proc_data.size(); j += 3) {
                    int u = proc_data[j];
                    int new_dist = proc_data[j+1];
                    int new_parent = proc_data[j+2];
                    
                    #pragma omp critical
                    {
                        if (new_dist < dist[u]) {
                            dist[u] = new_dist;
                            parent[u] = new_parent;
                        }
                    }
                }
            }
        }
    }
}

void parallelIncrementalSSSP(Graph& g, vector<int>& dist, vector<int>& parent, 
                           const vector<int>& part, int rank, int size, 
                           long long& total_sssp_time) {
    int updated = 1;
    int iterations = 0;
    const int MAX_ITERATIONS = 1000;
    
    auto start_time = high_resolution_clock::now();
    
    if (rank == 0) {
        cout << "Starting Dijkstra's algorithm..." << endl;
    }
    
    while (updated && iterations++ < MAX_ITERATIONS) {
        updated = 0;
        
        int local_updated = 0;
        #pragma omp parallel
        {
            int thread_updated = 0;
            #pragma omp for schedule(dynamic)
            for (int u = 0; u < g.V; u++) {
                if (part[u] != rank || dist[u] == INF) continue;
                
                for (AdjNode* curr = g.adj[u]; curr; curr = curr->next) {
                    int v = curr->dest;
                    int new_dist = dist[u] + curr->weight;
                    
                    if (new_dist < dist[v]) {
                        #pragma omp critical
                        {
                            if (new_dist < dist[v]) {
                                dist[v] = new_dist;
                                parent[v] = u;
                                thread_updated = 1;
                            }
                        }
                    }
                }
            }
            #pragma omp critical
            local_updated |= thread_updated;
        }
        
        exchangeBoundaryData(g, dist, parent, part, rank, size);
        
        if (rank == 0 && iterations % 10 == 0) {
            cout << "Iteration " << iterations << " completed." << endl;
        }
        
        int global_updated = 0;
        MPI_Allreduce(&local_updated, &global_updated, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
        updated = global_updated;
    }
    
    auto end_time = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end_time - start_time);
    total_sssp_time += duration.count();
    
    if (rank == 0) {
        cout << "Dijkstra's algorithm completed in " << iterations << " iterations." << endl;
        cout << "This execution took: " << duration.count() << " milliseconds" << endl;
    }
    
    if (iterations >= MAX_ITERATIONS && rank == 0) {
        cerr << "Warning: SSSP computation reached maximum iterations" << endl;
    }
}

void printResults(const vector<int>& dist, const vector<int>& parent, 
                 const vector<int>& idx_to_id, int src, int rank) {
    ofstream fout("output_" + to_string(rank) + ".txt");
    fout << "Process " << rank << " results:\n";
    
    int reachable_nodes = 0;
    int total_distance = 0;
    int max_distance = 0;
    
    #pragma omp parallel for schedule(static) reduction(+:reachable_nodes,total_distance) reduction(max:max_distance)
    for (size_t i = 0; i < dist.size(); i++) {
        if (dist[i] != INF) {
            reachable_nodes++;
            total_distance += dist[i];
            max_distance = max(max_distance, dist[i]);
            
            stringstream path_stream;
            path_stream << "Node " << idx_to_id[i] << ": Distance = " << dist[i] << ", Path = ";
            
            vector<int> path;
            for (int v = i; v != -1; v = parent[v]) {
                path.push_back(idx_to_id[v]);
            }
            
            for (auto it = path.rbegin(); it != path.rend(); ++it) {
                path_stream << *it << " ";
            }
            path_stream << "\n";
            
            #pragma omp critical
            fout << path_stream.str();
        }
    }
    
    fout << "\nSummary Statistics:\n";
    fout << "Reachable nodes: " << reachable_nodes << " out of " << dist.size() << "\n";
    if (reachable_nodes > 0) {
        fout << "Average distance: " << fixed << setprecision(2) << (double)total_distance / reachable_nodes << "\n";
        fout << "Maximum distance: " << max_distance << "\n";
    }
    
    fout.close();
    
    int global_reachable = 0;
    int global_total_dist = 0;
    int global_max_dist = 0;
    
    MPI_Reduce(&reachable_nodes, &global_reachable, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_distance, &global_total_dist, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&max_distance, &global_max_dist, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        cout << "\n===== Dijkstra's Algorithm Results =====" << endl;
        cout << "Source node: " << idx_to_id[src] << endl;
        cout << "Reachable nodes: " << global_reachable << " out of " << dist.size() << endl;
        if (global_reachable > 0) {
            cout << "Average distance: " << fixed << setprecision(2) 
                 << (double)global_total_dist / global_reachable << endl;
            cout << "Maximum distance: " << global_max_dist << endl;
        }
        cout << "Detailed results written to output_[rank].txt files" << endl;
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Initialize OpenMP
    omp_set_nested(true);
    omp_set_dynamic(true);
    
    auto overall_start_time = high_resolution_clock::now();
    long long total_sssp_time = 0;
    long long update_time = 0;
    long long graph_read_time = 0;
    long long partition_time = 0;

    map<int, int> id_to_idx;
    vector<int> idx_to_id;
    Graph g(0);

    // Graph reading
    auto graph_read_start = high_resolution_clock::now();
    if (rank == 0) {
        cout << "Reading graph file..." << endl;
        ifstream fin("graph.txt");
        if (!fin) {
            cerr << "Error: Could not open graph.txt" << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        string line;
        int from, to;
        
        while (getline(fin, line)) {
            if (line.empty() || line[0] == '#') continue;
            sscanf(line.c_str(), "%d%d", &from, &to);
            if (id_to_idx.find(from) == id_to_idx.end()) {
                id_to_idx[from] = idx_to_id.size();
                idx_to_id.push_back(from);
            }
            if (id_to_idx.find(to) == id_to_idx.end()) {
                id_to_idx[to] = idx_to_id.size();
                idx_to_id.push_back(to);
            }
        }
        
        g = Graph(idx_to_id.size());
        fin.clear();
        fin.seekg(0);
        
        while (getline(fin, line)) {
            if (line.empty() || line[0] == '#') continue;
            int from, to;
            sscanf(line.c_str(), "%d%d", &from, &to);
            g.addEdge(id_to_idx[from], id_to_idx[to], 1);
        }
        
        auto graph_read_end = high_resolution_clock::now();
        graph_read_time = duration_cast<milliseconds>(graph_read_end - graph_read_start).count();
        cout << "Graph reading completed. Time taken: " << graph_read_time << " ms" << endl;
        cout << "Number of vertices: " << g.V << endl;
    }

    int num_vertices = 0;
    if (rank == 0) {
        num_vertices = g.V;
    }
    MPI_Bcast(&num_vertices, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if (rank != 0) {
        g = Graph(num_vertices);
        idx_to_id.resize(num_vertices);
    }
    
    if (num_vertices > 0) {
        vector<int> temp_id_buffer;
        if (rank == 0) {
            temp_id_buffer = idx_to_id;
        } else {
            temp_id_buffer.resize(num_vertices);
        }
        
        MPI_Bcast(temp_id_buffer.data(), num_vertices, MPI_INT, 0, MPI_COMM_WORLD);
        
        if (rank != 0) {
            idx_to_id = temp_id_buffer;
        }
    }

    // Graph partitioning
    auto partition_start = high_resolution_clock::now();
    vector<idx_t> part(num_vertices);
    if (num_vertices > 0) {
        if (rank == 0) {
            cout << "Building CSR format for graph partitioning..." << endl;
        }
        
        vector<idx_t> xadj(num_vertices + 1);
        vector<idx_t> adjncy;
        
        if (rank == 0) {
            int edge_pos = 0;
            for (int i = 0; i < num_vertices; i++) {
                xadj[i] = edge_pos;
                for (AdjNode* curr = g.adj[i]; curr; curr = curr->next) {
                    adjncy.push_back(curr->dest);
                    edge_pos++;
                }
            }
            xadj[num_vertices] = edge_pos;
        }
        
        int edge_count = 0;
        if (rank == 0) {
            edge_count = adjncy.size();
            cout << "Number of edges: " << edge_count << endl;
        }
        MPI_Bcast(&edge_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        if (rank != 0) {
            adjncy.resize(edge_count);
        }
        
        MPI_Bcast(xadj.data(), num_vertices + 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        if (edge_count > 0) {
            MPI_Bcast(adjncy.data(), edge_count, MPI_INT, 0, MPI_COMM_WORLD);
        }
        
        if (rank != 0 && edge_count > 0) {
            #pragma omp parallel for schedule(static)
            for (int i = 0; i < num_vertices; i++) {
                for (int j = xadj[i]; j < xadj[i+1]; j++) {
                    g.addEdge(i, adjncy[j], 1);
                }
            }
        }
        
        if (rank == 0) {
            cout << "Partitioning graph with METIS..." << endl;
        }
        
        if (num_vertices > 0 && edge_count > 0) {
            idx_t ncon = 1;
            idx_t nparts = size;
            idx_t objval;
            idx_t options[METIS_NOPTIONS];
            METIS_SetDefaultOptions(options);
            options[METIS_OPTION_CONTIG] = 0;
            
            int result = METIS_PartGraphKway(
                &num_vertices, &ncon, xadj.data(), adjncy.data(),
                NULL, NULL, NULL, &nparts, NULL, NULL, options, &objval, part.data()
            );
            
            if (result != METIS_OK) {
                cerr << "METIS partitioning failed on rank " << rank << endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        } else {
            #pragma omp parallel for schedule(static)
            for (int i = 0; i < num_vertices; i++) {
                part[i] = i % size;
            }
        }
    }
    
    auto partition_end = high_resolution_clock::now();
    partition_time = duration_cast<milliseconds>(partition_end - partition_start).count();
    if (rank == 0) {
        cout << "Partitioning time: " << partition_time << " ms" << endl;
    }

    // Source node setup
    int src_id = -1, src = -1;
    if (rank == 0) {
        cout << "Enter source node ID: ";
        cin >> src_id;
        
        if (id_to_idx.find(src_id) != id_to_idx.end()) {
            src = id_to_idx[src_id];
            cout << "Source node " << src_id << " mapped to internal index " << src << endl;
        } else {
            cout << "Source node not found in graph." << endl;
        }
    }
    MPI_Bcast(&src, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if (src == -1) {
        MPI_Finalize();
        return 0;
    }

    // Initialize distances
    vector<int> dist(num_vertices, INF);
    vector<int> parent(num_vertices, -1);
    if (src >= 0 && src < num_vertices && part[src] == rank) {
        dist[src] = 0;
        parent[src] = -1;
    }

    // Initial SSSP computation
    if (num_vertices > 0) {
        if (rank == 0) {
            cout << "\n===== Initial SSSP Computation =====" << endl;
        }
        parallelIncrementalSSSP(g, dist, parent, part, rank, size, total_sssp_time);
    }

    // Dynamic updates
    int num_updates = 0;
    if (rank == 0) {
        cout << "\nEnter number of updates: ";
        cin >> num_updates;
    }
    MPI_Bcast(&num_updates, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if (num_vertices > 0) {
        srand(time(0) + rank * 17);
        
        for (int i = 0; i < num_updates; i++) {
            int update_type = 0, u = -1, v = -1;
            
            if (rank == 0) {
                cout << "\n===== Update " << i+1 << " =====" << endl;
            }
            
            auto update_start = high_resolution_clock::now();
            
            if (rank == 0) {
                update_type = rand() % 2;
                u = rand() % num_vertices;
                do {
                    v = rand() % num_vertices;
                } while (v == u);
                
                cout << (update_type ? "Delete" : "Insert") 
                     << " edge " << idx_to_id[u] << "->" << idx_to_id[v] << endl;
            }
            
            MPI_Bcast(&update_type, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&u, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&v, 1, MPI_INT, 0, MPI_COMM_WORLD);
            
            if (update_type == 0) {
                g.addEdge(u, v, 1);
            } else {
                g.removeEdge(u, v);
            }
            
            auto update_end = high_resolution_clock::now();
            update_time += duration_cast<milliseconds>(update_end - update_start).count();
            
            parallelIncrementalSSSP(g, dist, parent, part, rank, size, total_sssp_time);
        }
    }

    printResults(dist, parent, idx_to_id, src, rank);
    
    auto overall_end_time = high_resolution_clock::now();
    auto overall_duration = duration_cast<milliseconds>(overall_end_time - overall_start_time);
    
    if (rank == 0) {
        cout << "\n===== Performance Summary =====" << endl;
        cout << "Total execution time: " << overall_duration.count() << " milliseconds" << endl;
        cout << "  Graph reading time: " << graph_read_time << " milliseconds" << endl;
        cout << "  Partitioning time: " << partition_time << " milliseconds" << endl;
        cout << "  Dijkstra computation time: " << total_sssp_time << " milliseconds" << endl;
        cout << "  Graph update time: " << update_time << " milliseconds" << endl;
        cout << "  Other overhead: " << overall_duration.count() - graph_read_time - partition_time 
             - total_sssp_time - update_time << " milliseconds" << endl;
        cout << "Number of processes: " << size << endl;
        cout << "Number of threads per process: " << omp_get_max_threads() << endl;
        cout << "Number of vertices: " << num_vertices << endl;
        cout << "Number of dynamic updates: " << num_updates << endl;
    }

    g.clear();
    MPI_Finalize();
    return 0;
}