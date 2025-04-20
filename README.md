# Single-Source Shortest Path  (SSSP)

This project implements a **parallel algorithm** to update **single-source shortest paths (SSSP)** in large-scale **dynamic graphs**, based on the research paper:  
 **"A Parallel Algorithm Template for Updating Single-Source Shortest Paths in Large-Scale Dynamic Networks"**  
 *Authors: Khanda, Srinivasan, Bhowmick, Norris, Das (IEEE TPDS, 2022)*

---

## Objective

To design and implement a **scalable parallel algorithm** using **OpenMP** that efficiently updates SSSP trees in dynamic graphs where edge insertions and deletions occur — without recomputing everything from scratch.

---

## Technologies & Tools

- **C++** – Language used for implementation  
- **OpenMP** – For shared-memory, intra-node parallelism  
- **Dynamic Scheduling** – For workload balancing during updates  
- **CSR Format** – For graph representation  
- **Real-world Datasets** – From Network Repository and synthetic R-MAT graphs  
- **GitHub** – For collaboration and progress tracking

---

## Project Flow

1. **Graph Input**: Read static graph in CSR or adjacency list format  
2. **Edge Updates**: Load batch of inserted/deleted edges  
3. **Subgraph Detection**: Mark affected vertices  
4. **Parallel Updates**: Relax distances using OpenMP threads  
5. **Asynchronous Rounds**: Optionally control level of synchronization  
6. **Final Output**: Updated SSSP tree and performance stats

---

## Key Concepts from the Research

- **Dynamic Graphs**: Graphs where edges are frequently inserted or deleted
- **SSSP Tree Representation**: Uses a rooted tree storing distance and parent per vertex
- **Affected Subgraph Identification**: Only vertices impacted by changes are processed
- **Iterative Updates Without Locks**: Convergence via repeated distance checks instead of synchronization
- **Parallelism**:
  - Intra-node using **OpenMP**
  - Dynamic scheduling for load balancing
- **Avoid Cycle Formation**: During insertions, maintain acyclic SSSP tree
- **Asynchronous Updates**: Reduces synchronization cost by processing neighbors multiple levels deep
- **Batch Processing**: Handles massive edge changes in manageable chunks for performance

---

## 👨‍💻 Contributors

- Hammad Shabbir  
- Haider Zia 
- Muhammad Iqrash Qureshi
