#  Parallel Multi-Objective Shortest Path Updater (PMOSPU)

This project implements a parallel algorithm to update multi-objective shortest paths in large dynamic graphs, based on the research paper:  
**"A Parallel Algorithm for Updating a Multi-objective Shortest Path in Large Dynamic Networks" (Khanda et al.)**

---

##  Objective

To design and implement a scalable, parallel solution using MPI, OpenMP, and METIS that efficiently updates shortest paths in dynamic graphs without recomputing them from scratch.

---

## Technologies & Tools

- **MPI** – For distributed inter-node parallelism  
- **OpenMP** – For shared-memory (intra-node) parallelism  
- **METIS** – For efficient graph partitioning  
- **C++** – Language for implementation  
- **GitHub** – Version control and progress tracking  
- **MPI Performance Analyzer** – For performance evaluation

---

## Key Concepts from the Research

- Incremental shortest path updates using parallel processing  
- Reduces redundant computations by localizing updates  
- Uses multi-threaded processing within and across graph partitions  
- Efficient message passing between partitions for synchronization

---

## Project Flow

1. **Graph Partitioning** with METIS  
2. **Parallel Updates** in each subgraph using OpenMP  
3. **MPI Synchronization** across partitions  
4. **Performance Evaluation** on different datasets and machines

---

## 📂 Repository Structure

